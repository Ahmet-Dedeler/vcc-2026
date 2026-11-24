"""307 step 7 -- P6, POST-HOC and labelled: does the aggregate reproduce across 305's two cell
halves, which share no cell AND no control construct?

305 sec.5: within each batch a construct's cells alternate into half 0 / half 1 (deterministic, no
RNG); the 31 non-targeting constructs are split by parity into control half A (16 constructs,
21,581 cells) and control half B (15 constructs, 16,595 cells).  Half 0 is referenced to A, half 1
to B.  This is a CELL-SAMPLING split, not a biological replicate -- 305 and 307's preregistration
both say so, and nothing here is labelled replication.  What it does test is whether the aggregate
is a stable feature of the measurement or an accident of which cells and which controls were used.

    PYTHONDONTWRITEBYTECODE=1 uv run python halves307.py
"""
from __future__ import annotations

import json
import pickle

import numpy as np
from scipy import stats

import lib307 as X
import enrich_pools307 as EP

SP = {s: X.Split(s) for s in X.SPLITS}
cons_rows = X.load_constructs()
rowmap = {r["construct"]: r for r in cons_rows}
TARG = [r["construct"] for r in cons_rows if r["unit"] == "targeting"]
gene_meta = X.load_genes()
lookup = pickle.load(open(X.HERE / "seq307_lookup.pkl", "rb"))
pi_off = pickle.load(open(X.HERE / "poolhits307.pkl", "rb"))["pi_off"]
EP.setup(SP, cons_rows, gene_meta, lookup, pi_off)
SPLIT_OF = {r["construct"]: r["split"] for r in cons_rows if r["unit"] == "targeting"}
POOLS = {c: EP.build_pools(c, rowmap[c]) for c in TARG}
log = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    log.append(s)


rng = np.random.default_rng(X.SEED307 + 31)
z0, z1, zp = [], [], []
for c in TARG:
    sp = SP[SPLIT_OF[c]]
    e0, e1 = sp.half_effects(c)
    z0.append(EP.z_of(e0, POOLS[c]))
    z1.append(EP.z_of(e1, POOLS[c]))
    zp.append(EP.z_of(sp.effect(c), POOLS[c]))
z0, z1, zp = map(np.asarray, (z0, z1, zp))


def boot(v):
    b = v[rng.integers(0, v.size, (20000, v.size))].mean(axis=1)
    return float(np.quantile(b, 0.025)), float(np.quantile(b, 0.975))


P("P6 POST-HOC -- 305's two cell halves (disjoint cells AND disjoint control constructs; a "
  "cell-sampling split, NOT a biological replicate):")
P(f"   pooled  S = {zp.mean():+.5f}  95% [{boot(zp)[0]:+.5f}, {boot(zp)[1]:+.5f}]")
P(f"   half 0  S = {z0.mean():+.5f}  95% [{boot(z0)[0]:+.5f}, {boot(z0)[1]:+.5f}]   "
  f"(vs control half A, 16 NTC constructs / 21,581 cells)")
P(f"   half 1  S = {z1.mean():+.5f}  95% [{boot(z1)[0]:+.5f}, {boot(z1)[1]:+.5f}]   "
  f"(vs control half B, 15 NTC constructs / 16,595 cells)")
r_p = stats.pearsonr(z0, z1)
r_s = stats.spearmanr(z0, z1)
agree = int(np.sum(np.sign(z0) == np.sign(z1)))
P(f"   per-construct agreement (n = {len(TARG)} targeting constructs): Pearson "
  f"{r_p.statistic:+.4f} (p {r_p.pvalue:.4g}), Spearman {r_s.statistic:+.4f}; "
  f"same sign in {agree}/{len(TARG)}, both positive in "
  f"{int(np.sum((z0 > 0) & (z1 > 0)))}/{len(TARG)}")
P(f"   both halves positive AND pooled positive: "
  f"{int(np.sum((z0 > 0) & (z1 > 0) & (zp > 0)))}/{len(TARG)}")

out = dict(S_pooled=float(zp.mean()), S_half0=float(z0.mean()), S_half1=float(z1.mean()),
           ci_pooled=list(boot(zp)), ci_half0=list(boot(z0)), ci_half1=list(boot(z1)),
           pearson=float(r_p.statistic), pearson_p=float(r_p.pvalue),
           spearman=float(r_s.statistic), same_sign=agree, n=len(TARG),
           both_positive=int(np.sum((z0 > 0) & (z1 > 0))),
           per_construct={c: dict(z0=float(a), z1=float(b), z=float(d))
                          for c, a, b, d in zip(TARG, z0, z1, zp)})
(X.HERE / "halves307.json").write_text(json.dumps(out, indent=1))
(X.HERE / "halves307.log").write_text("\n".join(log) + "\n")
