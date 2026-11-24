"""190 Phase 2/3 — price the regulatory prior as an ORDERING KEY on h1-test.

THE SCORE, DEFINED BEFORE ANYTHING WAS MEASURED
-----------------------------------------------
A network is reduced to one signed weighted matrix

    N[t, g]     sign = the direction the NETWORK expects gene g to move when target t is
                       KNOCKED DOWN   (CollecTRI: -mor;  co-expression: -sign(r))
                |.|  = the network's evidence weight for that pair
                0    = the network says nothing about this pair

The per-(target, gene) REGULON-AGREEMENT SCORE is

    R(t, g) = w(t, g) * sign(N[t, g]) * sign(H[t, g])          if N != 0 and H != 0
            = 0                                                 otherwise
    w(t, g) = |N[t, g]| / median(|N| over its nonzeros)          (scale normalisation)

R > 0  the network CONFIRMS the sign we already predict  -> we should be MORE confident
R < 0  the network CONTRADICTS it                        -> we should be LESS confident
R = 0  no information

and the key is

    key_lambda(t, g) = |H[t, g]| * exp(lambda * R(t, g))

which at lambda = 0 is exactly |H|, the known-answer row (reach 0.23533).  A rank-additive
variant is reported alongside because the multiplicative one is sensitive to |H|'s scale.

`reach` reads only the order of the reference's DE pool; it never reads the profile.  So a key
change is FREE on pds / mse / nmae by construction -- that is the whole point of the mandate,
and it is why this is priced as a key and not as a direction blend (Phase 4 does that).
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.join(ROOT, "experiments", "190-regulon-ordering")
WORK = os.path.join(HERE, "work")
sys.path.insert(0, os.path.join(ROOT, "experiments", "172-sign-confidence-key"))
sys.path.insert(0, HERE)

BIG = 1e9
LAMBDAS = (-4.0, -2.0, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0)


# --------------------------------------------------------------------------------------
def agreement(N: np.ndarray, H: np.ndarray) -> np.ndarray:
    """R(t, g), as defined in the module docstring."""
    nz = N != 0
    if not nz.any():
        return np.zeros_like(N, dtype=np.float64)
    w = np.abs(N).astype(np.float64)
    med = float(np.median(w[nz]))
    w = w / (med if med > 0 else 1.0)
    return np.where(nz, w * np.sign(N) * np.sign(H), 0.0)


def key_mult(H: np.ndarray, R: np.ndarray, lam: float) -> np.ndarray:
    return np.abs(H) * np.exp(np.clip(lam * R, -50, 50))


def midrank_within(A: np.ndarray) -> np.ndarray:
    """Within-row midrank of A, scaled to (0, 1].  Monotone in A, so rank-0 == |H| order."""
    out = np.empty(A.shape, dtype=np.float64)
    n = A.shape[1]
    for i in range(A.shape[0]):
        order = np.argsort(A[i], kind="stable")
        r = np.empty(n, dtype=np.float64)
        r[order] = np.arange(1, n + 1, dtype=np.float64)
        # midrank for ties
        v = A[i][order]
        s = 0
        while s < n:
            e = s
            while e + 1 < n and v[e + 1] == v[s]:
                e += 1
            if e > s:
                r[order[s:e + 1]] = 0.5 * (s + 1 + e + 1)
            s = e + 1
        out[i] = r / n
    return out


def key_rank(rH: np.ndarray, R: np.ndarray, lam: float) -> np.ndarray:
    return rH + lam * R


# --------------------------------------------------------------------------------------
def main() -> None:
    import harness172
    import prior190

    T = harness172.H1Test()
    print("KNOWN-ANSWER ROWS (rule 10) -- abort on any failure")
    harness172.check(T)
    H = T.H.astype(np.float64)
    pool = T.pool
    sT = np.sign(T.S.real_lfc)
    rH = midrank_within(np.abs(H))
    base = T.reach(np.abs(H))
    assert abs(base - 0.23533) < 5e-6
    assert abs(T.reach(rH) - base) < 1e-9, "rank transform must be order-preserving"
    print(f"  OK  reach(rank(|H|)) == reach(|H|)        {T.reach(rH):.5f}\n")

    out: dict = {"base_absH": base, "base_production": T.reach()}

    # ---- the networks ----------------------------------------------------------------
    nets: dict[str, np.ndarray] = {}
    N_ct, s_ct = prior190.collectri_matrix(T.perts, T.genes, weighted=True)
    nets["collectri"] = N_ct.astype(np.float64)
    N_ct1, _ = prior190.collectri_matrix(T.perts, T.genes, weighted=False)
    nets["collectri_unweighted"] = N_ct1.astype(np.float64)
    Z = np.load(os.path.join(WORK, "coex190.npz"), allow_pickle=True)
    assert np.array_equal(np.asarray([str(g) for g in Z["genes"]]), T.genes)
    assert np.array_equal(np.asarray([str(p) for p in Z["perts"]]), T.perts)
    for k in (0, 2, 5):
        nets[f"coexpr_pc{k}"] = Z[f"C_pc{k}"].astype(np.float64)
    # co-expression restricted to each target's strongest module (top-K by |r|) -- the
    # rest of the row is noise (174) and would swamp the informative part
    for K in (100, 500):
        C = nets["coexpr_pc0"].copy()
        M = np.zeros(C.shape, bool)
        for p in range(C.shape[0]):
            M[p, np.argsort(-np.abs(C[p]))[:K]] = True
        nets[f"coexpr_pc0_top{K}"] = np.where(M, C, 0.0)

    # ---- A. does the network's SIGN beat ours where it speaks? -----------------------
    print("A. NETWORK SIGN ACCURACY ON THE REFERENCE DE POOL, where the network speaks")
    print(f"{'network':24} {'pool pairs':>11} {'% of pool':>10} {'net sign':>9} "
          f"{'H sign there':>13} {'H sign all':>11}")
    accs = {}
    for name, N in nets.items():
        m = pool & (N != 0) & (H != 0)
        n = int(m.sum())
        if n == 0:
            continue
        a_net = float(np.mean(np.sign(N[m]) == sT[m]))
        a_H = float(np.mean(np.sign(H[m]) == sT[m]))
        accs[name] = dict(n=n, frac_pool=n / int(pool.sum()), net=a_net, H_there=a_H)
        print(f"{name:24} {n:11,} {100*n/pool.sum():9.2f}% {a_net:9.4f} {a_H:13.4f} "
              f"{float(np.mean(np.sign(H[pool]) == sT[pool])):11.4f}")
    out["sign_accuracy"] = accs

    # ---- B. the CEILING: perfect use of the covered pairs ----------------------------
    print("\nB. CEILING -- an ORACLE restricted to the pairs the network covers")
    print("   (promote covered pairs that are truly right to the front, demote the wrong")
    print("    ones to the back; leave every uncovered pair at |H|.  No prior built on this")
    print("    network can beat this.)")
    correct = (np.sign(T.S.pred_lfc) == np.sign(T.S.real_lfc))
    ceil = {}
    for name, N in nets.items():
        cov = N != 0
        k = np.abs(H).copy()
        k = np.where(cov & correct, BIG + np.abs(H), k)
        k = np.where(cov & ~correct, -BIG + np.abs(H), k)
        r = T.reach(k)
        ceil[name] = r
        print(f"   {name:24} ceiling reach {r:.5f}   (+{r - base:+.5f} over |H|)")
    out["ceiling"] = ceil

    # ---- C. the key sweep -------------------------------------------------------------
    print("\nC. KEY SWEEP -- reach(|H| * exp(lambda R))   [lambda=0 IS the known-answer row]")
    hdr = "   " + f"{'network':24}" + "".join(f"{l:>9.2f}" for l in LAMBDAS)
    print(hdr)
    sweeps = {}
    for name, N in nets.items():
        R = agreement(N, H)
        row = [T.reach(key_mult(H, R, l)) for l in LAMBDAS]
        sweeps[name] = {"mult": row}
        print(f"   {name:24}" + "".join(f"{v:9.5f}" for v in row))
    print(f"\n   rank-additive variant: reach(midrank(|H|) + lambda R)")
    print(hdr)
    for name, N in nets.items():
        R = agreement(N, H)
        row = [T.reach(key_rank(rH, R, l)) for l in LAMBDAS]
        sweeps[name]["rank"] = row
        print(f"   {name:24}" + "".join(f"{v:9.5f}" for v in row))
    out["sweep"] = {"lambdas": list(LAMBDAS), "reach": sweeps}

    # ---- D. hard promotion / demotion -------------------------------------------------
    print("\nD. HARD variants (the most aggressive thing a binary prior can do)")
    hard = {}
    for name, N in nets.items():
        R = agreement(N, H)
        rows = {
            "promote_agree": T.reach(np.where(R > 0, BIG + np.abs(H), np.abs(H))),
            "demote_disagree": T.reach(np.where(R < 0, -BIG + np.abs(H), np.abs(H))),
            "both": T.reach(np.where(R > 0, BIG + np.abs(H),
                                     np.where(R < 0, -BIG + np.abs(H), np.abs(H)))),
        }
        hard[name] = rows
        print(f"   {name:24} " + "  ".join(f"{k} {v:.5f}" for k, v in rows.items()))
    out["hard"] = hard

    # ---- E. top-X% precision, the objective 172 says reach actually reads -------------
    print("\nE. TOP-X% SIGN ACCURACY over the pool (pooled), best lambda per network")
    print(f"   {'key':32}" + "".join(f"{f'{int(f*100)}%':>8}" for f in harness172.TOPX))
    tx = {"absH": T.topx(np.abs(H))}
    print(f"   {'|H| (benchmark)':32}" + "".join(f"{v:8.3f}" for v in tx["absH"]))
    txm = {"absH": T.topx_macro(np.abs(H))}
    best = {}
    for name, N in nets.items():
        R = agreement(N, H)
        row = sweeps[name]["mult"]
        j = int(np.argmax(row))
        best[name] = {"lambda": LAMBDAS[j], "reach": row[j]}
        k = key_mult(H, R, LAMBDAS[j])
        tx[name] = T.topx(k)
        txm[name] = T.topx_macro(k)
        print(f"   {name + f' (lam={LAMBDAS[j]:g})':32}" + "".join(f"{v:8.3f}" for v in tx[name]))
    out["topx_pooled"] = tx
    out["topx_macro"] = txm
    out["best_lambda"] = best

    with open(os.path.join(WORK, "key190.json"), "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\nwrote {WORK}/key190.json")


if __name__ == "__main__":
    main()
