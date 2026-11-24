"""155 — blends, and the two nulls the raw numbers have to clear.

1. BLEND.  Does adding a fraction of the control-covariance signal to the incumbent buy sign
   accuracy or head purity?  The candidate is L2-matched to the incumbent row by row first,
   so `f` is a genuine mixing fraction and not a hidden magnitude knob.

2. STRATUM NULL.  `probe.py` §B shows the k=0 candidate reaching 0.5487 sign accuracy on the
   33 best-measured targets.  Before that is called signal it has to beat (a) its own
   permutation and (b) the marginal, i.e. the fraction of that stratum's reference DE pool
   that is simply down-regulated.

3. THE MARGINAL ITSELF.  A constant all-down vector scores sign accuracy 0.5974 against the
   incumbent's 0.5366.  Unlike whitening or reweighting, an additive offset DOES flip signs,
   so this is the one cheap "different signs" source on the table and it is priced here.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bench import HDR, ROW, head_purity, load_frame, score  # noqa: E402

s, inc, gdst, pok = load_frame()


def frame(M):
    P = np.zeros_like(s.pred_lfc)
    P[:, gdst] = M[:, gdst]
    return P


def dbl(M):
    M = M - M.mean(0, keepdims=True)
    return M - M.mean(1, keepdims=True)


def match(P):
    """Rescale each row of P to the incumbent row's L2 norm on the scored columns."""
    a = np.linalg.norm(inc[:, gdst], axis=1)
    b = np.maximum(np.linalg.norm(P[:, gdst], axis=1), 1e-30)
    return P * (a / b)[:, None]


def main() -> None:
    base = score(s, inc)
    print(HDR)
    print(ROW.format(name="incumbent (104 baseline, g=0.65)", **base))
    print()

    cands = {}
    for k, mode in ((0, "raw"), (0, "double"), (20, "double")):
        B = np.load(os.path.join(HERE, "work", f"coef_k{k}.npz"),
                    allow_pickle=True)["beta"].astype(np.float64)
        cands[f"k={k} {mode}"] = match(frame(-(dbl(B) if mode == "double" else B)))

    print("1. BLEND  (1-f)*incumbent + f*candidate, candidate L2-matched per row")
    for name, C in cands.items():
        print(f"   candidate: control-cov -beta {name}")
        for f in (0.0, 0.02, 0.05, 0.10, 0.20, 0.40, 1.0):
            r = score(s, (1 - f) * inc + f * C)
            print(ROW.format(name=f"     f = {f:.2f}", **r)
                  + f"   dacc {r['sign_acc']-base['sign_acc']:+.4f}"
                    f"  dhead {r['head_purity']-base['head_purity']:+.4f}"
                    f"  dreach {r['reach']-base['reach']:+.5f}")
        print()

    print("2. STRATUM NULL -- the 33 best-measured targets, where k=0 looked positive")
    cpm = np.load(os.path.join(HERE, "work", "coef_k20.npz"), allow_pickle=True)["target_cpm"]
    top = np.argsort(-cpm)[:33]
    keep = np.zeros_like(s.real_sig)
    keep[top] = True
    rs, ad_ = s.real_sig.copy(), s.adj.copy()
    s.real_sig, s.adj = rs & keep, ad_ & keep
    pool = np.concatenate([s.real_lfc[p][s.real_sig[p] & s.adj[p]] for p in top])
    print(f"   stratum reference pool: {pool.size} gene-calls, {100*(pool<0).mean():.1f}% down")

    def q(P):
        s.pred_lfc = P
        return s.sign_accuracy(), head_purity(s, P)[0]

    B0 = np.load(os.path.join(HERE, "work", "coef_k0.npz"),
                 allow_pickle=True)["beta"].astype(np.float64)
    P0 = frame(-dbl(B0))
    a, h = q(P0)
    print(f"   candidate  (k=0, double-centered)   sign acc {a:.4f}   head purity {h:.4f}")
    accs, hps = [], []
    for seed in range(20):
        perm = np.random.default_rng(seed).permutation(len(s.perts))
        aa, hh = q(P0[perm])
        accs.append(aa)
        hps.append(hh)
    print(f"   permutation null (20 seeds)         sign acc {np.mean(accs):.4f} "
          f"+- {np.std(accs):.4f}   head purity {np.mean(hps):.4f} +- {np.std(hps):.4f}")
    z = (a - np.mean(accs)) / max(np.std(accs), 1e-9)
    print(f"   -> z = {z:+.2f} against its own permutation null")
    s.real_sig, s.adj = rs, ad_
    print()

    print("3. GLOBAL DOWN-SHIFT  incumbent - c * (row sd), the one additive lever that "
          "flips signs")
    sd = inc[:, gdst].std(1)
    shift = np.zeros_like(inc)
    shift[:, gdst] = -1.0
    shift *= sd[:, None]
    for c in (0.0, 0.1, 0.25, 0.5, 1.0, 2.0):
        r = score(s, inc + c * shift)
        print(ROW.format(name=f"   c = {c:.2f}", **r)
              + f"   dacc {r['sign_acc']-base['sign_acc']:+.4f}"
                f"  dhead {r['head_purity']-base['head_purity']:+.4f}"
                f"  dreach {r['reach']-base['reach']:+.5f}"
                f"  dpds {r['pds']-base['pds']:+.4f}")


if __name__ == "__main__":
    main()
