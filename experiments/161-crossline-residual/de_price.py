"""161 SS6 — price the DE-side members for the arms this experiment considered.

AGENTS.md §11.1 rule 12: price every arm on all six members, or state which are unmeasured
AND why they cannot move.  a19 was selected on a three-member screen and lost -0.01352 on the
member the screen never computed.

WHAT CAN AND CANNOT BE READ HERE
--------------------------------
This surface is a donor-table pseudobulk surface with no cells, so there is no Wilcoxon and no
call set.  Rule 9 is the harder constraint: `fid`, `jac` and `reach`'s normalisation all carry
`n_real` in a denominator or a union, and `n_real` is a property of the reference's cell count,
which no surface we hold reproduces.  So the LEVELS are unreadable by construction.

What IS readable is the ORDERING and SIGN content of the delivered effect, which is the part
of those members an arm can actually move, and it is read exactly as 156 §4 read it:

    sign accuracy   over the reference's own top-|truth| pool, per target, then averaged
    head purity     sign accuracy inside each target's top decile by |predicted|
    reach proxy     deepest prefix of the |predicted| ordering that stays >=90% sign-pure,
                    as a fraction of the pool -- per target, NEVER pooled (rule 8)
    nmae proxy      mean |pred - truth| over the pool, normalised by mean |truth|

Everything is per target and then averaged, because every DE member is scored per
perturbation and then averaged; pooling all (target, gene) pairs is a different and much
easier statistic (rule 8, exp 154).
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np

from residual import Surface, deliver, restore_norms
from common import WORK

POOL = 200          # per-target reference DE pool: the truth's strongest genes


def price(P, T, pool=POOL):
    """Per-target ordering/sign statistics of a delivered effect P against truth T."""
    nT = P.shape[0]
    sacc, head, reach, nmae = [], [], [], []
    for t in range(nT):
        tr, pr = T[t], P[t]
        idx = np.argsort(-np.abs(tr))[:pool]          # the reference's own DE pool
        s_t, s_p = np.sign(tr[idx]), np.sign(pr[idx])
        ok = s_t == s_p
        sacc.append(ok.mean())
        nmae.append(np.abs(pr[idx] - tr[idx]).sum() / max(np.abs(tr[idx]).sum(), 1e-12))
        order = np.argsort(-np.abs(pr[idx]))          # OUR ordering of that pool
        oks = ok[order]
        head.append(oks[:max(len(oks) // 10, 1)].mean())
        cum = np.cumsum(oks) / np.arange(1, len(oks) + 1)
        deep = np.where(cum >= 0.90)[0]
        reach.append((deep.max() + 1) / len(oks) if len(deep) else 0.0)
    return dict(sign_acc=float(np.mean(sacc)), head_purity=float(np.mean(head)),
                reach_proxy=float(np.mean(reach)), nmae_proxy=float(np.mean(nmae)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--surface", default="core5")
    ap.add_argument("--resid", default=os.path.join(WORK, "residual-core5-R1.json"))
    a = ap.parse_args()

    S = Surface(a.surface, ctrl=False)
    out = {}
    print(f"surface {a.surface}: {S.nT} targets x {S.nG} genes, pool = top {POOL} |truth| genes")
    for L in S.lines:
        pool_l = [x for x in S.lines if x != L]
        C, _, _ = S.consensus(pool_l)
        b, gate, T = S.bL[L], S.gateL[L], S.EL[L]
        dn = np.median([np.median(np.linalg.norm(S.R[d], axis=1)) for d in pool_l])
        alpha = dn / np.median(np.linalg.norm(C, axis=1))
        arms = {
            "consensus (baseline)": C,
            "power p=2 (context-free)": restore_norms(np.sign(C) * np.abs(C) ** 2, C),
            "ORACLE truth": T,
        }
        print(f"\n{'=' * 100}\nHELD LINE {L}")
        print(f"  {'arm':<28} {'sign acc':>9} {'head pur':>9} {'reach*':>9} {'nmae*':>9}")
        row = {}
        for nm, D in arms.items():
            P = deliver(alpha * D, b, gate) if nm != "ORACLE truth" else deliver(D, b, gate)
            m = price(P, T)
            row[nm] = m
            print(f"  {nm:<28} {m['sign_acc']:9.4f} {m['head_purity']:9.4f} "
                  f"{m['reach_proxy']:9.4f} {m['nmae_proxy']:9.4f}")
        out[L] = row
    with open(os.path.join(WORK, "de_price.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print(f"\n{'=' * 100}\nMEAN OVER HELD LINES")
    names = list(next(iter(out.values())))
    print(f"  {'arm':<28} {'sign acc':>9} {'head pur':>9} {'reach*':>9} {'nmae*':>9}")
    for nm in names:
        v = {k: np.mean([out[L][nm][k] for L in out]) for k in
             ("sign_acc", "head_purity", "reach_proxy", "nmae_proxy")}
        print(f"  {nm:<28} {v['sign_acc']:9.4f} {v['head_purity']:9.4f} "
              f"{v['reach_proxy']:9.4f} {v['nmae_proxy']:9.4f}")


if __name__ == "__main__":
    main()
