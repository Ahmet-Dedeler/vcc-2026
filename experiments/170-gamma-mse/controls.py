"""170 — controls that would FAIL if the thing being tested were absent.

Two vacuous checks have already been caught in this project in two days (an algebraic identity
sold as a mass-fix verification, and a row permutation that is invariant by construction), so
each control below is paired with what it reads when the mechanism is removed.

1. RENDER IDENTITY. With `MASS_FIX=minnorm` the renderer's normalisation is exactly the
   identity, so `log1p(TS*p) - b_ctrl` must equal the gated, mass-fixed effect.
   WITHOUT the mass fix the same check reads ~2-3% of the row norm (156 measured
   2.83-3.35% on the challenge frame). So it is not an identity by construction -- it is
   the thing the mass fix buys, and the control is run BOTH ways to show the contrast.

2. EMITTED-PROFILE FIDELITY. The quasi emitter repairs each gene's panel total exactly, so
   the pseudobulk the scorer reads back from the emitted cells must match the intended one.
   a21 measured rel L2 0.0014 for a correct build and 0.0155 for the wrong repair order.
   If the emission were mis-wired this reads percent-level.

3. GAMMA IS ACTUALLY DOING SOMETHING. The predicted effect at gamma 0.35 vs 0.80 must differ
   by much more than the arm-to-arm noise floor. Reported as row-wise cosine between the two
   delivered effect matrices; if the operator were being silently skipped this reads 1.000.

    VCC_MASS_FIX=minnorm uv run python controls.py
    uv run python controls.py            # same checks with the mass fix OFF
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from h1kit import (  # noqa: E402
    DISPERSION, KEEP_SPREAD, KEEP_TOP, MIN_CPM, NCELL, OVERDISPERSION, TAIL_SPREAD, TS,
    Surface170, builder, deliver, to_h1_axis)

WORK = HERE / "work"


def main() -> None:
    from cell_eval2.prep import bulk_lognorm_means
    from vcc2026.donor_transfer import emit_counts
    from scipy.sparse import csr_matrix

    import importlib.util
    path = HERE.parents[1] / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("prod_builder_170c", path)
    B = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_170c"] = B
    spec.loader.exec_module(B)
    print(f"[ctrl] MASS_FIX={B.MASS_FIX}")

    S = Surface170()
    pre = np.load(WORK / "pre.npz", allow_pickle=True)
    pos = {str(t): i for i, t in enumerate(pre["h1_targets"])}
    order = np.asarray([pos[str(t)] for t in S.targets])

    delivered = {}
    for gamma in (0.35, 0.65, 0.80):
        eff_prod = np.load(WORK / "eff" / f"h1_g{gamma:.2f}.npz")["effect"]
        eff, _ = to_h1_axis(eff_prod, pre["genes"], pre["h1_genes"])
        eff = eff[order]

        # --- 1. render identity
        gated = np.where(S.gate[None, :], eff, 0.0)
        if B.MASS_FIX is not None:
            fixed = np.vstack([B.mass_consistent_effect(gated[i], S.b_ctrl, S.gate,
                                                        mode=B.MASS_FIX)
                               for i in range(gated.shape[0])])
        else:
            fixed = gated
        P, dh = deliver(B, eff, S.b_ctrl)
        rel = np.median(np.linalg.norm(dh - fixed, axis=1)
                        / np.linalg.norm(fixed, axis=1))
        print(f"[ctrl] gamma {gamma:.2f}  render residual vs the intended effect: "
              f"{100 * rel:.4f}% of the row norm "
              f"({'mass fix ON -- expect ~0' if B.MASS_FIX else 'mass fix OFF -- expect ~2-3%'})")
        delivered[gamma] = dh

        # --- 2. emitted-profile fidelity, on 12 targets
        sel = np.linspace(0, P.shape[0] - 1, 12).astype(int)
        errs = []
        for i in sel:
            spread = B.emission_spread_vector(
                eff[i], S.b_ctrl, min_cpm=MIN_CPM, tail_spread=TAIL_SPREAD,
                keep_top=KEEP_TOP, keep_spread=KEEP_SPREAD, spread_map="block")
            v, cl, ip = emit_counts(P[i], S.libs, DISPERSION,
                                    np.random.default_rng(7919 * 0 + int(i)), mode="quasi",
                                    overdispersion=OVERDISPERSION, spread=spread)
            X = csr_matrix((v.astype(np.float64), cl.astype(np.int64), ip.astype(np.int64)),
                           shape=(NCELL, P.shape[1]))
            got = bulk_lognorm_means(np.asarray(X.sum(0)).ravel()[None, :], TS)[0]
            want = S.b_ctrl + dh[i]
            errs.append((np.linalg.norm(got - want) / np.linalg.norm(want),
                         np.linalg.norm(got - want) / np.linalg.norm(want - S.b_ctrl)))
        e = np.asarray(errs)
        print(f"[ctrl] gamma {gamma:.2f}  emitted-pseudobulk rel L2 vs intended: "
              f"{np.median(e[:, 0]):.6f} of the PROFILE norm "
              f"(a21's convention: 0.0014 correct, 0.0155 mis-ordered repair)"
              f" | {np.median(e[:, 1]):.5f} of the EFFECT norm -- that second number is "
              f"the floor() truncation of `n*L*profile`, identical across arms")

    # --- 3. gamma is not a no-op
    a, b = delivered[0.35], delivered[0.80]
    ca = a / np.maximum(np.linalg.norm(a, axis=1, keepdims=True), 1e-300)
    cb = b / np.maximum(np.linalg.norm(b, axis=1, keepdims=True), 1e-300)
    cos = np.einsum("ij,ij->i", ca, cb)
    print(f"[ctrl] row-wise cosine between the delivered effects at gamma 0.35 and 0.80: "
          f"median {np.median(cos):.5f}, min {cos.min():.5f}  (1.000 would mean the "
          f"operator is not being applied)")


if __name__ == "__main__":
    main()
