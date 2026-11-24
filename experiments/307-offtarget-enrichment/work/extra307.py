"""307 step 6 -- POST-HOC extensions of declared secondary S3, clearly labelled as post-hoc.

S3 asked for `GRB10_P2 -> SNX8` with its partner construct.  This generalises that specificity
check to all 34 constructs: every candidate gene of construct c is re-scored in the cells of c's
PARTNER construct (the other construct of the same intended target, different protospacers), using
c's own matched pools.  If the association is spacer/construct-specific the partner sits at the
wrong-recipient level; if it is intended-target downstream biology the partner is repressed too.

Also: the construct-level sign count, and the bias-corrected bounded estimate.

    PYTHONDONTWRITEBYTECODE=1 uv run python extra307.py
"""
from __future__ import annotations

import collections
import json
import pickle

import numpy as np
from scipy import stats

import lib307 as X
import enrich_pools307 as EP

SP = {s: X.Split(s) for s in X.SPLITS}
GIX = SP["train"].gix
cons_rows = X.load_constructs()
rowmap = {r["construct"]: r for r in cons_rows}
TARG = [r["construct"] for r in cons_rows if r["unit"] == "targeting"]
gene_meta = X.load_genes()
lookup = pickle.load(open(X.HERE / "seq307_lookup.pkl", "rb"))
pi_off = pickle.load(open(X.HERE / "poolhits307.pkl", "rb"))["pi_off"]
EP.setup(SP, cons_rows, gene_meta, lookup, pi_off)
SPLIT_OF = {r["construct"]: (r["split"] if r["unit"] == "targeting" else "train") for r in cons_rows}
EFF = {c: SP[SPLIT_OF[c]].effect(c, loo_control=(rowmap[c]["unit"] != "targeting"))
       for c in [r["construct"] for r in cons_rows]}
POOLS = {c: EP.build_pools(c, rowmap[c]) for c in [r["construct"] for r in cons_rows]}
ZOBS = {c: EP.z_of(EFF[c], POOLS[c]) for c in POOLS}
log = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    log.append(s)


# ---------------------------------------------------------------- partners
by_target = collections.defaultdict(list)
for c in TARG:
    by_target[rowmap[c]["intended_target"]].append(c)
partner = {}
for t, cs in by_target.items():
    if len(cs) == 2:
        partner[cs[0]], partner[cs[1]] = cs[1], cs[0]
P(f"paired intended targets: {sum(1 for v in by_target.values() if len(v) == 2)} of "
  f"{len(by_target)}; constructs with a partner: {len(partner)}")


def excl_ix(c):
    bad = {rowmap[c]["intended_target"]} if rowmap[c]["intended_target"] else set()
    for L in (12, 11, 10, 9):
        bad |= set(lookup["cand_L"][c][L])
    return np.asarray(sorted(GIX[g] for g in bad if g in GIX))


def z_cross(scorer: str, setowner: str) -> float:
    """z of setowner's candidate set, evaluated in scorer's cells, with scorer's exclusions."""
    bad = excl_ix(scorer)
    e, qs = EFF[scorer], []
    for p in POOLS[setowner]:
        if p["gi"] in bad:
            continue
        mem = p["members"][~np.isin(p["members"], bad)]
        if mem.size < 20:
            continue
        v, pv = e[p["gi"]], e[mem]
        qs.append((np.sum(pv < v) + 0.5 * np.sum(pv == v)) / pv.size)
    return float(np.sqrt(12 * len(qs)) * (0.5 - np.mean(qs))) if qs else 0.0


zp = np.asarray([z_cross(partner[c], c) for c in TARG if c in partner])
zo = np.asarray([ZOBS[c] for c in TARG if c in partner])
rng = np.random.default_rng(X.SEED307 + 21)
bp = zp[rng.integers(0, zp.size, (20000, zp.size))].mean(axis=1)
diff = np.asarray([zo[rng.integers(0, zo.size, zo.size)].mean()
                   - zp[rng.integers(0, zp.size, zp.size)].mean() for _ in range(20000)])
P(f"\nPOST-HOC partner specificity (n = {zp.size} constructs with a partner, unit = construct):")
P(f"   own construct   z = {zo.mean():+.5f}")
P(f"   PARTNER construct, same intended target, different protospacers: "
  f"z = {zp.mean():+.5f}  95% [{np.quantile(bp, 0.025):+.5f}, {np.quantile(bp, 0.975):+.5f}]")
P(f"   own - partner = {zo.mean() - zp.mean():+.5f}  95% "
  f"[{np.quantile(diff, 0.025):+.5f}, {np.quantile(diff, 0.975):+.5f}]  "
  f"(paired Wilcoxon p = {stats.wilcoxon(zo, zp).pvalue:.5f})")

# ---------------------------------------------------------------- sign count
zt = np.asarray([ZOBS[c] for c in TARG])
npos = int(np.sum(zt > 0))
sgn = stats.binomtest(npos, zt.size, 0.5, alternative="greater")
P(f"\nsign count: {npos}/{zt.size} targeting constructs have z > 0 "
  f"(exact binomial p = {sgn.pvalue:.6g}); median z {np.median(zt):+.5f}")

# ---------------------------------------------------------------- bias-corrected bounded estimate
sens = json.load(open(X.HERE / "sens307.json"))
bias, bse = sens["P1"]["mean"], sens["P1"]["se"]
res = json.load(open(X.HERE / "results307.json"))["primary"]
S, lo, hi = res["S"], *res["ci95_construct_bootstrap"]
P(f"\nbias-corrected bounded estimate (subtracting the {sens['P1']['n']:,}-derangement "
  f"wrong-recipient level {bias:+.5f} +- {bse:.5f}):")
P(f"   S - bias = {S - bias:+.5f}, construct-bootstrap 95% "
  f"[{lo - bias:+.5f}, {hi - bias:+.5f}]  (n = {res['n']} targeting constructs)")

out = dict(partner=dict(n=int(zp.size), z_own=float(zo.mean()), z_partner=float(zp.mean()),
                        partner_ci=[float(np.quantile(bp, 0.025)), float(np.quantile(bp, 0.975))],
                        diff=float(zo.mean() - zp.mean()),
                        diff_ci=[float(np.quantile(diff, 0.025)), float(np.quantile(diff, 0.975))],
                        wilcoxon_p=float(stats.wilcoxon(zo, zp).pvalue),
                        per_construct={c: dict(own=ZOBS[c], partner=float(z_cross(partner[c], c)))
                                       for c in TARG if c in partner}),
           sign=dict(n_pos=npos, n=int(zt.size), binom_p=float(sgn.pvalue),
                     median_z=float(np.median(zt))),
           bias_corrected=dict(S_minus_bias=S - bias, ci=[lo - bias, hi - bias],
                               bias=bias, bias_se=bse))
(X.HERE / "extra307.json").write_text(json.dumps(out, indent=1, default=float))
(X.HERE / "extra307.log").write_text("\n".join(log) + "\n")
