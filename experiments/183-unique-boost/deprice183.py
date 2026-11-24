"""183 — what the boost does to `q` and `k`, the two things the DE-side members read.

`nmae`, `reach` and `fid` are not locally priceable as levels (rule 9: `n_real`; rule 13:
2/6 and 3/6 sign hits).  But the 2026-09-03 02:15 landscape note pins the MECHANISM:

    nmae(k, q) = q|k - 1| + (1 - q)(k + 1)

and `nmae`, `reach` and `fid` all move only through the sign accuracy `q` and the DE-pool
magnitude ratio `k`.  Both are measurable on h1-test's own reference DE pool without any
`n_real` normalisation, so this script measures them directly instead of declaring the three
members unmeasurable.  A transform that leaves `q` and the pool magnitude alone cannot move
them, which is the rule-12 "mechanism that pins it" exemption rather than a hand-wave.

Deterministic.  No emitter, no cell-eval2, no build, no submission.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "experiments/168-sharpening"))
sys.path.insert(0, os.path.join(REPO, "experiments/156-delivered-geometry"))
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
sys.path.insert(0, HERE)
from h1_tiebreak import group_profiles, H1, CTRL, TS                   # noqa: E402
from massfix import mass_fix                                           # noqa: E402
from screen import render                                              # noqa: E402
from debench import Surface                                            # noqa: E402
from uboost import unique_boost                                        # noqa: E402

ALPHAS = (0.00, 0.25, 0.50, 0.75, 1.00)

tg, genes_t, acc_t = group_profiles(H1)
ci = int(np.flatnonzero(tg == CTRL)[0])
prof_t = np.log1p(TS * acc_t / acc_t.sum(1, keepdims=True))
b_ctrl = prof_t[ci]
targets_all = tg[np.asarray([i for i in range(len(tg)) if i != ci])]

E = np.load(os.path.join(REPO, "experiments/169-sharpening-six-member/work/h1test-effect.npz"),
            allow_pickle=True)
tpos = {str(t): i for i, t in enumerate(E["targets"])}
targets = np.asarray([t for t in targets_all if str(t) in tpos])
H = E["H"][np.asarray([tpos[str(t)] for t in targets])].astype(np.float64)
gate = E["gate"].astype(bool)
gpos = {str(g): i for i, g in enumerate(genes_t)}

S = Surface(os.path.join(REPO, "experiments/164-callset-artefact/de164.a15.npz"))
de_gi = np.asarray([gpos[g] for g in S.genes])
de_ti = np.asarray([list(S.perts).index(str(t)) for t in targets])
pool = S.real_sig & S.adj                       # the reference's own DE pool
sT = np.sign(S.real_lfc)
tmag = np.abs(S.real_lfc)

Hg = np.where(gate[None, :], H, 0.0)
D0 = render(mass_fix(Hg, b_ctrl, gate), b_ctrl)

print(f"h1-test reference DE pool: {int(pool.sum())} (perturbation, gene) pairs "
      f"over {len(targets)} targets")
print(f"\n{'alpha':>6} {'q (sign acc)':>13} {'d q':>9} {'mean |p| on pool':>18} "
      f"{'r |p|':>9} {'k = med(|p|/|t|)':>18} {'pool energy share':>18}")
base = None
for alpha in ALPHAS:
    D = unique_boost(D0, alpha)
    P = np.zeros(S.real_lfc.shape)
    P[de_ti] = D[:, de_gi]
    q = float(np.mean(np.sign(P)[pool] == sT[pool]))
    mp = float(np.mean(np.abs(P[pool])))
    k = float(np.median(np.abs(P[pool]) / np.maximum(tmag[pool], 1e-12)))
    share = float((P[pool] ** 2).sum() / max((D * D).sum(), 1e-300))
    if base is None:
        base = (q, mp)
    print(f"{alpha:>6.2f} {q:13.5f} {q - base[0]:+9.5f} {mp:18.6f} "
          f"{mp / base[1]:9.6f} {k:18.6f} {share:18.6f}")
