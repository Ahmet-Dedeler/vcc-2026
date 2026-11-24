"""165 — error bars on the PROFILE members, which is where the placement family actually
separates.

`pds` and the `mse` proxy are deterministic in the delivered profile (no emitter, no seed), but
they are still means over 85 perturbations, so a placement difference needs a standard error
before it is a finding.  Both are computed here as PAIRED per-perturbation contributions
against `a15`, with the decoy set held at all 85 rows, and clustered on perturbation.

`pds` per perturbation is its own rank contribution `(beaten + ties/2)/(n-1)`; the `mse` proxy
per perturbation is `1 - cos^2` at the panel's own best global amplitude, decomposed as the
per-row squared-error share.  The known-answer row is `cst`: 156 measured that the uniform
log shift LOSES pds on both held Orion lines, so it must lose here too.

Reads only `work/profiles_all.npz`.  No emitter, no cell-eval2, no build, no submission.
"""
from __future__ import annotations

import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
TS = 5.0e4


def per_pert(P, T):
    """(pds contribution, own cosine, squared-error share) per perturbation."""
    A = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
    B = T / np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-300)
    S = A @ B.T
    n = S.shape[0]
    own = np.diag(S).copy()
    beaten = (S < own[:, None]).sum(1)
    ties = (S == own[:, None]).sum(1) - 1
    pds = (beaten + 0.5 * ties) / (n - 1)
    # sse* is 1 - (sum P.T)^2 / (sum P.P * sum T.T): a ratio of panel sums, so the honest
    # per-row object is the row's contribution to each of the three sums.
    return pds, own, (P * P).sum(1), (P * T).sum(1)


def main() -> None:
    z = np.load(f"{WORK}/profiles_all.npz", allow_pickle=True)
    b = z["b_ctrl"]
    perts = np.asarray([str(p) for p in z["perts"]])
    genes = np.asarray([str(g) for g in z["genes"]])
    rt = z["real_tot"]
    T = np.log1p(TS * rt / rt.sum(1, keepdims=True)) - b[None, :]
    gpos = {g: i for i, g in enumerate(genes)}
    drop = {gpos[p] for p in perts if p in gpos}
    mask = np.asarray([i for i in range(genes.size) if i not in drop])
    T = T[:, mask]
    dd = float((T * T).sum())
    names = [k[2:] for k in z.files if k.startswith("p_")]
    names = [n for n in ["a15", "mn", "mn_half", "cst", "hp_p500", "hp_p2000", "hp_p6000",
                         "hp_p9000", "hp_w2000", "eps_keep", "eps_tail", "inveff", "eff",
                         "mn_c10"]
             if n in names]

    cache = {}
    for a in names:
        P = (np.log1p(TS * z[f"p_{a}"].astype(np.float64)) - b[None, :])[:, mask]
        cache[a] = per_pert(P, T)
    p0, c0, a0, x0 = cache["a15"]
    sse0 = 1.0 - x0.sum() ** 2 / (a0.sum() * dd)
    n = len(perts)

    print("=" * 100)
    print("profile members, paired per-perturbation, clustered on perturbation "
          f"(n = {n}, {mask.size} scored genes)")
    print("=" * 100)
    print(f"{'arm':<11}{'d pds':>10}{'se':>9}{'t':>7}{'d cos':>11}{'se':>10}{'t':>7}"
          f"{'d sse*':>11}{'d cos %':>10}")
    for a in names:
        if a == "a15":
            continue
        p1, c1, a1, x1 = cache[a]
        dp, dc = p1 - p0, c1 - c0
        sep = dp.std(ddof=1) / np.sqrt(n)
        sec = dc.std(ddof=1) / np.sqrt(n)
        sse1 = 1.0 - x1.sum() ** 2 / (a1.sum() * dd)
        print(f"{a:<11}{dp.mean():>+10.5f}{sep:>9.5f}{dp.mean() / sep:>7.2f}"
              f"{dc.mean():>+11.6f}{sec:>10.6f}{dc.mean() / sec:>7.2f}"
              f"{sse1 - sse0:>+11.6f}{100 * (c1.mean() / c0.mean() - 1):>+9.2f}%")
    print("\n  KNOWN-ANSWER: `cst` (the uniform log shift) must LOSE pds — 156 measured "
          "-0.00555 / -0.00312 on the two held Orion lines.")


if __name__ == "__main__":
    main()
