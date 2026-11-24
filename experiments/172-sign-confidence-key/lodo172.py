"""172 — leave-one-donor-out training data for the sign-confidence key.

FOLDS.  Production's donors are `k562_gwps`, `orion_hct116`, `orion_hek293t`.  Each fold holds
one out as the RECIPIENT and transfers from the other two exactly as production would
(`MultiDonorModel(..., gene_norm=False, eb=False)`, uniform donor weights, `beta = 0`):

    held out          sources                              == deployment?
    orion_hct116      k562_gwps + orion_hek293t            YES - this is 095's h1-test pair
    orion_hek293t     k562_gwps + orion_hct116             one Orion line swapped
    k562_gwps         orion_hct116 + orion_hek293t         both sources from the same atlas

The third fold is the odd one out and is reported separately: two X-Atlas/Orion lines are far
less independent than K562+Orion, so "cross-source agreement" means something weaker there.

`data/donors/h1_2025` IS NOT USED ANYWHERE.  It covers 85/85 of h1-test's targets
(`donors172.py --overlap`) and is pure oracle leakage (AGENTS.md trap #1, and 157 §7 records
this exact leak firing at reach 0.404).  The audit is printed by `donors172.py` and reproduced
in RESULTS.md.

`orion_hct116` and `orion_hek293t` each cover 84/85 of h1-test's targets and ARE used as
training donors.  That is not leakage: the label in a fold is whether the transferred sign
matched the HELD-OUT LINE's own response, the model carries no gene identity and no target
identity, and it therefore has no channel through which "target X, gene Y goes up" could be
memorised.  What a fold can teach is only the shape of the map from measurable scalars to
sign reliability.

THE TRAINING POPULATION, AND ITS ONE APPROXIMATION
--------------------------------------------------
`reach` walks the REFERENCE's DE pool, so the rows the model should be fitted on are the rows
it will be asked to order: the recipient's own differentially-expressed genes.  A donor screen
has no p-values, only `delta`, so the pool is proxied by a Wald z built from the same Poisson
noise model `eb_shrink` already uses in production:

    z(t, g)  ~  |delta_L[t, g]| * sqrt(n_cells * depth * p_g)      ~  |delta| * sqrt(cpm)

within a target (`n_cells`, `depth` are per-target constants and cancel in a within-target
ranking).  `pool_check()` validates this proxy ON h1-test, where the true BH-significant set
IS known, and reports the recall of the real pool by the proxy's top-`n_real`.  A proxy that
recovered the pool at chance would make every fold label a different problem from the one
h1-test poses, so this check is load-bearing and its number is in RESULTS.md.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from donors172 import DonorHandle                     # noqa: E402
from feats172 import FEATURES, Recipient, build, stack  # noqa: E402

WORK = os.path.join(HERE, "work")
PROD = ("k562_gwps", "orion_hct116", "orion_hek293t")
FOLDS = {h: tuple(n for n in PROD if n != h) for h in PROD}
MIN_CPM = 5.0                 # the gate the 095 h1-test artifact shipped (169 probe_effect.py)
POOL_FRAC = 0.00826           # h1-test's median n_real (89) / its gated axis (10,771)
N_TARGETS = 3000              # per fold
BOTH_FRAC = 0.95              # h1-test has 81/85 targets carried by BOTH sources
CHUNK = 375


def pool_check(T) -> None:
    """Does `|lfc| * sqrt(cpm)` recover the reference's BH-significant set on h1-test?"""
    print("pool proxy validation on h1-test (where the true DE pool is known):")
    hits = {"|lfc|*sqrt(cpm)": [], "|lfc| alone": [], "cpm alone": [], "random": []}
    rng = np.random.default_rng(0)
    live = T.gate
    for p in range(len(T.perts)):
        n = int(T.pool[p].sum())
        if n < 10:
            continue
        idx = np.flatnonzero(live & T.S.adj[p])
        lfc = np.abs(np.nan_to_num(T.S.real_lfc[p, idx]))
        cpm = T.cpm[idx]
        truth = T.pool[p, idx]
        for name, key in (("|lfc|*sqrt(cpm)", lfc * np.sqrt(cpm)),
                          ("|lfc| alone", lfc),
                          ("cpm alone", cpm),
                          ("random", rng.random(idx.size))):
            sel = np.argpartition(-key, n - 1)[:n]
            hits[name].append(float(truth[sel].mean()))
    for name, v in hits.items():
        print(f"  recall of the real DE pool by the proxy's top-n_real, {name:<16} "
              f"median {np.median(v):.3f}  mean {np.mean(v):.3f}")


def build_fold(held: str, handles: dict, n_targets: int = N_TARGETS,
               seed: int = 0, verbose: bool = True) -> dict:
    src = FOLDS[held]
    L = handles[held]
    e = np.expm1(L.b_ctrl)
    cpm = 1e6 * e / e.sum()
    live = (cpm > MIN_CPM) & (L.b_ctrl > 0)
    gidx = np.flatnonzero(live)
    axis = L.genes[gidx]                      # recipient axis, gated -- the deployable columns
    cpm_live = cpm[gidx]
    b_live = L.b_ctrl[gidx]
    if verbose:
        print(f"\n=== fold: held out {held}  <-  {' + '.join(src)}")
        print(f"    recipient axis {len(L.genes)} genes -> {len(axis)} gated at "
              f"{MIN_CPM} CPM")

    have = {n: set(handles[n].pindex) for n in src}
    both = [str(p) for p in L.perts if all(str(p) in have[n] for n in src)]
    one = [str(p) for p in L.perts
           if any(str(p) in have[n] for n in src) and str(p) not in set(both)]
    rng = np.random.default_rng(seed)
    n_both = min(len(both), int(round(n_targets * BOTH_FRAC)))
    n_one = min(len(one), n_targets - n_both)
    targets = ([both[i] for i in rng.choice(len(both), n_both, replace=False)] +
               [one[i] for i in rng.choice(len(one), n_one, replace=False)])
    targets = np.asarray(sorted(targets))
    if verbose:
        print(f"    {len(targets)} targets ({n_both} carried by both sources, {n_one} by one) "
              f"out of {len(both)}+{len(one)} available")

    K = max(10, int(round(POOL_FRAC * len(axis))))
    gpos = {str(g): i for i, g in enumerate(axis)}
    X, y, tg, gg, cons_all = [], [], [], [], []
    centers = None
    for c0 in range(0, len(targets), CHUNK):
        tt = targets[c0:c0 + CHUNK]
        rec = Recipient(held, axis, b_live, tt, np.ones(len(axis), bool))
        B = build(rec, src, handles={n: handles[n] for n in src}, verbose=False)
        if centers is None:
            centers = True
        truth, carried = L.raw_rows(tt, axis)
        cons = B["cons"]
        F = B["F"]
        z = np.abs(truth) * np.sqrt(cpm_live)[None, :]
        for r, t in enumerate(tt):
            zz = z[r].copy()
            own = gpos.get(str(t))
            if own is not None:
                zz[own] = -np.inf                       # the target's own row is never scored
            zz = np.where(np.isfinite(zz) & (truth[r] != 0), zz, -np.inf)
            sel = np.argpartition(-zz, K - 1)[:K]
            sel = sel[np.isfinite(zz[sel])]
            if sel.size == 0:
                continue
            rows_t = np.full(sel.size, r)
            X.append(stack(F, rows_t, sel))
            y.append((np.sign(cons[r, sel]) == np.sign(truth[r, sel])).astype(np.uint8))
            tg.append(np.full(sel.size, c0 + r, dtype=np.int32))
            gg.append(sel.astype(np.int32))
            cons_all.append(cons[r, sel])
        del B, F, truth, cons, z
    out = dict(X=np.concatenate(X), y=np.concatenate(y), t=np.concatenate(tg),
               g=np.concatenate(gg), cons=np.concatenate(cons_all),
               targets=targets, axis=axis, held=held, sources=np.asarray(src))
    if verbose:
        print(f"    {out['X'].shape[0]} rows x {out['X'].shape[1]} features, "
              f"pool size {K}/target, base sign accuracy {out['y'].mean():.4f}")
    return out


def main() -> None:
    from harness172 import H1Test, check
    T = H1Test()
    check(T, verbose=False)
    print("known-answer rows OK (harness172.check)\n")
    pool_check(T)
    os.makedirs(WORK, exist_ok=True)
    handles = {n: DonorHandle(n) for n in PROD}
    for held in PROD:
        d = build_fold(held, handles)
        np.savez_compressed(os.path.join(WORK, f"fold.{held}.npz"),
                            X=d["X"], y=d["y"], t=d["t"], g=d["g"], cons=d["cons"],
                            targets=d["targets"], axis=d["axis"],
                            features=np.asarray(FEATURES), sources=d["sources"])
        print(f"    wrote work/fold.{held}.npz")


if __name__ == "__main__":
    main()
