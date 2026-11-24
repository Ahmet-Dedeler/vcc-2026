#!/usr/bin/env python3
"""625 — the readout DECODE, committed with synthetic known answers BEFORE any read (coordinator's condition b).

Inputs (later, from captured board records): the base B's members and the 18 arms' members (raw nmae, fid, reach),
arm order exactly as modules625.npz's `flips` (16 Hadamard rows) then `check_flips` (modules 0, 15).
decode(base, arms) -> per member: x_M (16), the additivity check (§5), and the §6 transfer mask.

Algebra (PREREGISTRATION §4): Δ_k = m(arm k) - m(B). S = Δ_1 (all 16 flipped). For k >= 2: sum_M H[k-1,M] x_M =
S - 2 Δ_k (arm k flips the 8 modules with H = -1). With y = (S, S - 2Δ_2, ..., S - 2Δ_16): x = H^T y / 16 (H H^T = 16 I).
Additivity (§5): |Δ_check_j - x_{M_j}| <= 0.10 |x_{M_j}| + 1e-5 for the two single-module arms, on nmae and fid.
Transfer (§6): flip M iff x_M^nmae < 0 AND x_M^fid > 0 AND |x_M^nmae| >= 0.25 * median_M |x_M^nmae|.
If §5 fails: the fallback set = modules whose nmae and fid signs agree as in §6 AND whose |x^nmae| ranks in the top half.

    uv run python experiments/625-module-readout/decode625.py --selftest        # synthetic KAs; must print SELFTEST PASS
    uv run python experiments/625-module-readout/decode625.py --reads <json>   # the real decode (after the arms are scored)
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DESIGNS = {"v16": "work/modules625.npz", "v8": "work/modules625_v8.npz"}   # Amendment 2: v8 (8 CLEAN4 octiles) is current
H = FLIPS = CHECKS = None
CHECK_MODULES: list[int] = []


def use_design(name: str = "v8") -> None:
    """Load a design's Hadamard, flip rows and check rows (module index == Hadamard column)."""
    global H, FLIPS, CHECKS, CHECK_MODULES
    MOD = np.load(HERE / DESIGNS[name], allow_pickle=True)
    H = MOD["hadamard"].astype(float)
    FLIPS = MOD["flips"].astype(bool)
    CHECKS = MOD["check_flips"].astype(bool)
    CHECK_MODULES = [int(np.flatnonzero(c)[0]) for c in CHECKS]


use_design("v8")


def decode_member(base: float, arms: list[float]) -> np.ndarray:
    n = H.shape[0]
    d = np.asarray(arms[:n], float) - float(base)
    S = d[0]
    y = np.empty(n)
    y[0] = S
    y[1:] = S - 2.0 * d[1:]
    return H.T @ y / float(n)


def decode(reads: dict) -> dict:
    """reads = {"base": {member: value}, "arms": [{member: value}, ... 18 in FLIPS then CHECKS order]}."""
    n = H.shape[0]
    assert len(reads["arms"]) == n + 2, f"need {n} Hadamard arms + 2 check arms"
    out = {}
    for m in ("nmae", "fid", "reach"):
        x = decode_member(reads["base"][m], [a[m] for a in reads["arms"]])
        chk = []
        for j, M in enumerate(CHECK_MODULES):
            dj = reads["arms"][n + j][m] - reads["base"][m]
            chk.append(dict(module=M, delta=float(dj), x=float(x[M]), ok=bool(abs(dj - x[M]) <= 0.10 * abs(x[M]) + 1e-5)))
        out[m] = dict(x=x.tolist(), additivity=chk)
    xn, xf = np.asarray(out["nmae"]["x"]), np.asarray(out["fid"]["x"])
    additive = all(c["ok"] for m in ("nmae", "fid") for c in out[m]["additivity"])
    bar = 0.25 * float(np.median(np.abs(xn)))
    mask = (xn < 0) & (xf > 0) & (np.abs(xn) >= bar)
    if not additive:
        top = np.abs(xn) >= np.median(np.abs(xn))
        mask = (xn < 0) & (xf > 0) & top
    out["additive"] = additive
    out["transfer_mask"] = mask.astype(int).tolist()
    out["n_flip"] = int(mask.sum())
    out["reading"] = ("ABORT-additivity (fallback mask)" if not additive else ("GO" if mask.any() else "NULL"))
    return out


def simulate(x_n, x_f, x_r, base=(0.9, 0.51, 0.28), nonadd=0.0, pair=None, noise=0.0, rng=None):
    """Board reads for given per-module effects; nonadd adds an interaction on arms flipping both modules of `pair`;
    noise adds iid N(0, noise) to every arm's members (the fused design's per-entry leak)."""
    rows = list(FLIPS) + list(CHECKS)
    p0, p1 = pair if pair is not None else (0, len(x_n) - 1)
    arms = []
    for r in rows:
        inter = nonadd * (r[p0] and r[p1])
        e = (rng.normal(0, noise, 3) if (noise and rng is not None) else np.zeros(3))
        arms.append({"nmae": base[0] + float(np.dot(r, x_n)) + inter + e[0], "fid": base[1] + float(np.dot(r, x_f)) - inter + e[1],
                     "reach": base[2] + float(np.dot(r, x_r)) + e[2]})
    return {"base": {"nmae": base[0], "fid": base[1], "reach": base[2]}, "arms": arms}


def selftest() -> None:
    for design in ("v16", "v8"):
        use_design(design)
        n = H.shape[0]
        c0, c1 = CHECK_MODULES
        rng = np.random.default_rng(6253)
        # KA-1: random additive effects are recovered exactly
        for t in range(200):
            x_n, x_f, x_r = rng.normal(0, 1e-3, n), rng.normal(0, 1e-3, n), rng.normal(0, 1e-3, n)
            r = decode(simulate(x_n, x_f, x_r))
            assert np.allclose(r["nmae"]["x"], x_n, atol=1e-12) and np.allclose(r["fid"]["x"], x_f, atol=1e-12)
            assert np.allclose(r["reach"]["x"], x_r, atol=1e-12) and r["additive"]
        # KA-2: the transfer rule on a constructed case
        a, b, c, d = 2, n - 3, n - 4, n - 1
        x_n = np.full(n, 2e-4); x_f = np.full(n, -1e-4)        # CLEAN4 agrees everywhere -> flips hurt: NULL
        x_n[[a, b]] = -6e-4; x_f[[a, b]] = 3e-4                # modules a, b disagree on both members -> flip
        x_n[c] = -6e-4; x_f[c] = -3e-4                         # nmae says flip, fid disagrees -> no flip
        x_n[d] = -1e-6; x_f[d] = 1e-4                          # below the 0.25 * median bar -> no flip
        r = decode(simulate(x_n, x_f, np.zeros(n)))
        assert r["transfer_mask"] == [1 if k in (a, b) else 0 for k in range(n)], (design, r["transfer_mask"])
        assert r["reading"] == "GO"
        r0 = decode(simulate(np.full(n, 2e-4), np.full(n, -1e-4), np.zeros(n)))
        assert r0["reading"] == "NULL" and r0["n_flip"] == 0
        # KA-3: a non-additive interaction between the two CHECK modules trips the additivity check
        r = decode(simulate(rng.normal(0, 1e-3, n), rng.normal(0, 1e-3, n), np.zeros(n), nonadd=5e-3, pair=(c0, c1)))
        assert not r["additive"] and r["reading"].startswith("ABORT"), (design, r["reading"])
        # KA-4 (coordinator's arithmetic): decoded sd under iid per-entry noise sigma: column 0 -> sqrt(n^2+4(n-1))/n,
        # other columns -> 2 sqrt(n-1)/n
        sims = np.array([decode(simulate(np.zeros(n), np.zeros(n), np.zeros(n), noise=1.0, rng=rng))["nmae"]["x"]
                         for _ in range(4000)])
        sd = sims.std(0)
        assert abs(sd[0] - np.sqrt(n * n + 4 * (n - 1)) / n) < 0.06 and np.all(np.abs(sd[1:] - 2 * np.sqrt(n - 1) / n) < 0.06), sd
        print(f"SELFTEST PASS [{design}]: 200 additive cases to 1e-12; transfer rule exact; check-module interaction trips "
              f"additivity; decoded sd col0 {sd[0]:.3f} (theory {np.sqrt(n * n + 4 * (n - 1)) / n:.3f}), "
              f"others {sd[1:].mean():.3f} (theory {2 * np.sqrt(n - 1) / n:.3f})")
    use_design("v8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--reads", type=Path)
    ap.add_argument("--design", choices=("v8", "v16"), default="v8")
    a = ap.parse_args()
    use_design(a.design)
    if a.selftest:
        selftest()
        return
    r = decode(json.loads(a.reads.read_text()))
    print(json.dumps(r, indent=1))


if __name__ == "__main__":
    main()
