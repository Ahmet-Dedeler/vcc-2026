"""158 — the per-target scale lever: reproduce the oracle, then price what is buildable.

THE LEVER
---------
`expr_mse_unbiased_capped_norm`'s uncorrected core is a RATIO OF SUMS over targets,

    sse(s) = sum_p ||s_p d_p - delta_p||^2 / sum_p ||delta_p||^2

with `d_p` our delivered (a23 mass-fixed, rendered) effect and `delta_p` the truth.
`pds_cosine` unit-normalises each row before the cosine, so ANY s_p > 0 leaves every
rank exactly unchanged. That makes `s_p` a free lever on `mse` at zero `pds` cost --
provided s_p stays strictly positive (s_p = 0 destroys the row, s_p < 0 flips it, and
neither is a pure scaling; §5 of RESULTS.md verifies this numerically).

Minimiser: s_p* = <d_p, delta_p> / ||d_p||^2, which needs the truth. The whole
experiment is whether a BUILD-TIME key predicts it.

Everything here is deterministic: no emitter, no seeds, exact-space pseudobulk truth
from experiment 156 (`work/truth-*.npz`), production's own `mass_consistent_effect`.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
E156 = os.path.join(REPO, "experiments/156-delivered-geometry")
WORK = os.path.join(HERE, "work")
sys.path.insert(0, E156)

from massfix import mass_fix, production_effect  # noqa: E402
from screen import render  # noqa: E402

SURFACES = ("hct116", "hek293t")


# ------------------------------------------------------------------ the object
def delivered(surface: str):
    """The a23 delivered effect and the truth, on the scored gene mask.

    `D` is production's log-space effect (frozen a15 baseline, CPM-10 gated exactly as
    the builder gates it); `mass_fix(...,'minnorm')` is a23; `render` is the
    builder-plus-scorer round trip, which is the identity to 1e-12 after the fix but is
    applied anyway so screen and build cannot drift (156's own discipline).
    """
    D, T, b, b_ref, gate, mask, ptar = production_effect(surface)
    assert np.abs(b - b_ref).max() == 0.0
    Pd = render(mass_fix(D, b, gate, mode="minnorm"), b)
    return Pd[:, mask], T[:, mask], ptar, dict(D=D, b=b, gate=gate, mask=mask)


# -------------------------------------------------------------------- endpoints
def sse_of(scales: np.ndarray, P: np.ndarray, T: np.ndarray) -> float:
    """sum_p ||s_p d_p - delta_p||^2 / sum_p ||delta_p||^2, exactly as written."""
    R = scales[:, None] * P - T
    return float((R * R).sum() / (T * T).sum())


def global_star(P: np.ndarray, T: np.ndarray) -> float:
    """The single scale that minimises sse over all targets jointly."""
    return float((P * T).sum() / (P * P).sum())


def per_target_star(P: np.ndarray, T: np.ndarray) -> np.ndarray:
    """s_p* = <d_p, delta_p> / ||d_p||^2, the oracle."""
    num = (P * T).sum(1)
    den = (P * P).sum(1)
    return num / np.maximum(den, 1e-300)


def cos_agg(sse: float) -> float:
    """The aggregate cosine implied by an sse, i.e. sqrt(1 - sse). At the GLOBAL optimum
    this is exactly the aggregate cosine 156 reports; at a per-target optimum it is the
    norm-weighted quadratic mean of the per-target cosines, which is why the same
    identity is used for both (156's `sse*` column is `1 - cos_agg^2`)."""
    return float(np.sqrt(max(0.0, 1.0 - sse)))


def per_target_cos(P: np.ndarray, T: np.ndarray) -> np.ndarray:
    pn = np.linalg.norm(P, axis=1)
    tn = np.linalg.norm(T, axis=1)
    return (P * T).sum(1) / np.maximum(pn * tn, 1e-300)


def pds_of(P: np.ndarray, T: np.ndarray) -> float:
    """156's exact rank statistic, ties at half."""
    A = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
    B = T / np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-300)
    S = A @ B.T
    n = S.shape[0]
    own = np.diag(S).copy()
    beaten = (S < own[:, None]).sum(1)
    ties = (S == own[:, None]).sum(1) - 1
    return float(np.mean((beaten + 0.5 * ties) / (n - 1)))


# ----------------------------------------------------------------------- report
def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    out = {}
    for surface in SURFACES:
        P, T, ptar, _ = delivered(surface)
        n = P.shape[0]

        g = global_star(P, T)
        sse_g = sse_of(np.full(n, g), P, T)
        s = per_target_star(P, T)
        sse_o = sse_of(s, P, T)

        # oracle MAGNITUDE only: the true row norm is known, the per-target cosine is not,
        # so each row is rescaled to ||delta_p|| and then ONE global scale is re-optimised.
        rn = np.linalg.norm(T, axis=1) / np.maximum(np.linalg.norm(P, axis=1), 1e-300)
        Pm = rn[:, None] * P
        sse_m = sse_of(np.full(n, global_star(Pm, T)), Pm, T)

        # oracle SIGN only: the same global scale, but rows whose s_p* < 0 flipped.
        flip = np.where(s < 0, -1.0, 1.0)
        Pf = flip[:, None] * P
        sse_f = sse_of(np.full(n, global_star(Pf, T)), Pf, T)

        row = dict(
            n=n,
            sse_global=sse_g, cos_global=cos_agg(sse_g), g_star=g,
            sse_oracle=sse_o, cos_oracle=cos_agg(sse_o),
            sse_mag_oracle=sse_m, gain_mag=sse_g - sse_m,
            sse_sign_oracle=sse_f, gain_sign=sse_g - sse_f,
            frac_neg=float(np.mean(s <= 0)),
            corr_s_norm=float(np.corrcoef(s, np.linalg.norm(P, axis=1))[0, 1]),
            s_star=s.tolist(), targets=[str(x) for x in ptar],
            per_target_cos=per_target_cos(P, T).tolist(),
            row_norm=np.linalg.norm(P, axis=1).tolist(),
            truth_norm=np.linalg.norm(T, axis=1).tolist(),
        )
        out[surface] = row

        print(f"\n{'=' * 78}\n{surface}: {n} targets, {P.shape[1]} scored genes")
        print(f"  global-optimal single scale      sse {sse_g:.6f}  "
              f"cos_agg {cos_agg(sse_g):.5f}   (g* = {g:.4f})")
        print(f"  ORACLE per-target scale          sse {sse_o:.6f}  "
              f"cos_agg {cos_agg(sse_o):.5f}   (+{100*(cos_agg(sse_o)/cos_agg(sse_g)-1):.0f}%)")
        print(f"  ORACLE MAGNITUDE ONLY            sse {sse_m:.6f}  "
              f"gain vs global {sse_g - sse_m:+.6f}")
        print(f"  ORACLE SIGN ONLY                 sse {sse_f:.6f}  "
              f"gain vs global {sse_g - sse_f:+.6f}")
        print(f"  fraction of targets s_p* <= 0    {np.mean(s <= 0):.3f}")
        print(f"  corr(s_p*, ||d_p||)              {np.corrcoef(s, np.linalg.norm(P, axis=1))[0,1]:+.3f}")

        # pds invariance, verified rather than asserted
        base = pds_of(P, T)
        rng = np.random.default_rng(0)
        for label, sc in (("random positive", np.exp(rng.normal(0, 1.5, n))),
                          ("oracle clipped >0", np.maximum(s, 1e-6)),
                          ("oracle raw (has <=0)", s)):
            print(f"  pds under {label:22s} {pds_of(sc[:, None] * P, T):.10f} "
                  f"(base {base:.10f}, delta {pds_of(sc[:, None]*P, T) - base:+.2e})")

    with open(os.path.join(WORK, "oracle.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {os.path.join(WORK, 'oracle.json')}")


if __name__ == "__main__":
    main()
