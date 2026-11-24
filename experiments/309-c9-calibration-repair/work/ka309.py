"""309 step 1 -- the known-answer rows that do not touch the repaired statistic.

KA-0  split composition of the 34 targeting constructs (reported, not a verdict row; C9's
      within-split centring depends on it and C8's RESULTS.md sec.1 disagrees with its own CSV)
KA-1  305's 23 published rows recomputed from `sums305_<split>.npz`      -> max |diff| <= 1e-9
KA-2a 301's selection reproduced from its saved hits                     -> 4/13, 4/21, 4/31; 8; 83
KA-2b the inherited 209,070-spacer pool is the CORRECTED one             -> 0 escapes, 0 own-target
      (Horlbeck `="SYMBOL"` Excel-escaping bug)                             violations, rates exact
KA-3  the 31 NTC constructs' sums identical across the three splits      -> exact
KA-7  293's on-target `kd_lfc` for all 34 constructs                     -> max |diff| <= 1e-9

    PYTHONDONTWRITEBYTECODE=1 uv run python ka309.py
"""
from __future__ import annotations

import collections
import csv
import json
import pickle
import time

import numpy as np

import lib309 as X

t0 = time.time()
log: list[str] = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    log.append(s)


out: dict = {}
SP = {s: X.Split(s) for s in X.SPLITS}
cons_rows = X.load_constructs()
TARG = [r for r in cons_rows if r["unit"] == "targeting"]
P("loaded 305's sums: " + ", ".join(
    f"{s} {len(SP[s].constructs)} constructs ({len(SP[s].targeting)} targeting + "
    f"{len(SP[s].ntc)} NTC)" for s in X.SPLITS))

# ---------------------------------------------------------------- KA-0 split composition
comp = collections.Counter(r["split"] for r in TARG)
grp = collections.Counter(r["group"] for r in TARG)
P(f"KA-0 the {len(TARG)} targeting constructs by split: {dict(comp)}; by group {dict(grp)}; "
  f"wrong recipients available per construct: "
  f"{ {s: comp[s] - 1 for s in comp} }")
P("     C8 RESULTS.md sec.1 says 'train 21 / validation 7 / test 6'; C8's own constructs307.csv "
  "says the above.  The CSV is the object every C8 computation used (SPLIT_OF is read from it).")
out["KA0"] = dict(by_split=dict(comp), by_group=dict(grp),
                  c8_results_md_claim={"train": 21, "validation": 7, "test": 6})

# ---------------------------------------------------------------- KA-1
r305 = json.load(open(X.E305 / "results305.json"))["rows"]
d_pool, d_h0, d_h1, d_ctrl = [], [], [], []
for r in r305:
    sp = SP[r["split"]]
    c, g = r["construct"], r["candidate_gene"]
    loo = bool(r["loo_control"])
    gi = sp.gix[g]
    e = sp.effect(c, loo_control=loo)[gi]
    e0, e1 = (v[gi] for v in sp.half_effects(c, loo_control=loo))
    ref = sp.ctrl_pooled - (sp.S[c].sum(0) if loo else 0)
    d_pool.append(abs(e - r["e_pool"]))
    d_h0.append(abs(e0 - r["e_half0"]))
    d_h1.append(abs(e1 - r["e_half1"]))
    d_ctrl.append(abs(50000.0 * ref[gi] / ref.sum() - r["ctrl_cp50k"]))
ka1 = max(max(d_pool), max(d_h0), max(d_h1), max(d_ctrl))
P(f"KA-1 305's {len(r305)} published rows reproduced from its saved sums: max |diff| "
  f"pooled {max(d_pool):.3e}, half0 {max(d_h0):.3e}, half1 {max(d_h1):.3e}, "
  f"ctrl CP50k {max(d_ctrl):.3e}  -> {'PASS' if ka1 <= 1e-9 else 'FAIL'}")
out["KA1"] = dict(n_rows=len(r305), max_abs_diff=float(ka1), passed=bool(ka1 <= 1e-9))

# ---------------------------------------------------------------- KA-2a  301's selection
lookup = pickle.load(open(X.E307 / "seq307_lookup.pkl", "rb"))
CAND_L = lookup["cand_L"]
cnt = collections.Counter()
tot = collections.Counter()
l12 = []
for r in cons_rows:
    c = r["construct"]
    key = r["group"] if r["unit"] == "targeting" else "non-targeting"
    tot[key] += 1
    if len(CAND_L[c][12]):
        cnt[key] += 1
        if r["unit"] == "targeting":
            l12.append(c)
lstar = collections.Counter(int(r["L_star"]) for r in TARG)
npairs = sum(int(r["k"]) for r in TARG)
ka2a = (cnt["failing"] == 4 and tot["failing"] == 13 and cnt["working"] == 4 and tot["working"] == 21
        and cnt["non-targeting"] == 4 and tot["non-targeting"] == 31 and len(l12) == 8
        and npairs == 83 and len(TARG) == 34)
P(f"KA-2a 301's rule from its saved hits: failing {cnt['failing']}/{tot['failing']}, "
  f"working {cnt['working']}/{tot['working']}, NTC {cnt['non-targeting']}/{tot['non-targeting']}; "
  f"{len(l12)} L=12 targeting constructs; L* {dict(sorted(lstar.items(), reverse=True))}; "
  f"{npairs} L* pairs over {len(TARG)} constructs -> {'PASS' if ka2a else 'FAIL'}")
out["KA2a"] = dict(counts={k: [cnt[k], tot[k]] for k in ("failing", "working", "non-targeting")},
                   n_L12_targeting=len(l12), L12_constructs=sorted(l12),
                   L_star={str(k): v for k, v in lstar.items()}, n_pairs=npairs, passed=bool(ka2a))

# ---------------------------------------------------------------- KA-2b  the CORRECTED pool
D = pickle.load(open(X.E307 / "poolhits307.pkl", "rb"))
pool, hits = D["pool"], D["hits"]
lmer = lookup["lmer_genes"]
n_escaped = sum(1 for _, g in pool if str(g).startswith('="'))
viol = 0
for L in X.LS:
    for (sp, tgt), h in zip(pool, hits[L]):
        if tgt and tgt in h:
            viol += 1
rate_pool = {L: float(np.mean([1.0 if h else 0.0 for h in hits[L]])) for L in X.LS}
# the bug-equivalent rate: the own-target / distal filters never fire (an escaped label matches
# nothing), so every promoter-corpus seed match counts as an off-target
bug12 = float(np.mean([1.0 if lmer[12].get(sp[-12:]) else 0.0 for sp, _ in pool]))
meta = json.load(open(X.E307 / "pool307_meta.json"))
probe = json.load(open(X.E307 / "probe307.json"))
# the 68 real 2025 targeting spacers, located in the pool exactly as C8's probe did
sp2ix = collections.defaultdict(list)
for i, (s_, g_) in enumerate(pool):
    sp2ix[s_].append(i)
spac = {r["guide_id"]: r for r in csv.DictReader(open(X.E301 / "vcc2025_spacers.csv"))
        if r["route"] != "UNRESOLVED"}
real = []
for a in csv.DictReader(open(X.E298 / "seed298.csv")):
    key = a["construct"] + "_A|" + a["construct"] + "_B"
    r = spac.get(key) or next(v for k, v in spac.items() if k.startswith(a["construct"] + "_A|"))
    real += [r["spacer_A"].upper(), r["spacer_B"].upper()]
real_ix = [sp2ix[s][0] for s in real if s in sp2ix]
rate_real = {L: float(np.mean([1.0 if hits[L][i] else 0.0 for i in real_ix])) for L in X.LS}
d_pool12 = abs(rate_pool[12] - meta["frac_spacers_with_hit"]["12"])
d_real12 = abs(rate_real[12] - probe["real_rate"]["12"])
ka2b = (n_escaped == 0 and viol == 0 and len(real_ix) == 68
        and d_pool12 <= 1e-6 and d_real12 <= 1e-6 and bug12 > 0.5)
P(f"KA-2b the inherited pool: {len(pool):,} spacers, {n_escaped} Excel-escaped gene labels, "
  f"{viol} own-target-in-own-hit-set violations over all four L; "
  f"{len(real_ix)}/68 real 2025 spacers found verbatim")
P(f"      per-spacer P(>=1 off-target hit) at L=12: pool {rate_pool[12]:.6f} "
  f"(saved {meta['frac_spacers_with_hit']['12']:.6f}, |diff| {d_pool12:.1e}), "
  f"real {rate_real[12]:.6f} (saved {probe['real_rate']['12']:.6f}, |diff| {d_real12:.1e}); "
  f"BUG-EQUIVALENT (no own-target/distal filter) {bug12:.3f} -> {'PASS' if ka2b else 'FAIL'}")
out["KA2b"] = dict(n_pool=len(pool), n_escaped_labels=n_escaped, own_target_violations=viol,
                   n_real_found=len(real_ix),
                   rate_pool={str(L): rate_pool[L] for L in X.LS},
                   rate_real={str(L): rate_real[L] for L in X.LS},
                   bug_equivalent_rate_L12=bug12, diff_pool12=d_pool12, diff_real12=d_real12,
                   passed=bool(ka2b))

# ---------------------------------------------------------------- KA-3
ntc_sets = {s: set(SP[s].ntc) for s in X.SPLITS}
same_names = ntc_sets["train"] == ntc_sets["validation"] == ntc_sets["test"]
maxdiff = cellsdiff = 0
for c in sorted(ntc_sets["train"]):
    a = SP["train"].S[c]
    for s in ("validation", "test"):
        maxdiff = max(maxdiff, int(np.abs(a - SP[s].S[c]).max()))
        cellsdiff = max(cellsdiff, int(np.abs(SP["train"].cells[c] - SP[s].cells[c]).max()))
ntc_cells = int(sum(SP["train"].cells[c].sum() for c in ntc_sets["train"]))
ka3 = bool(same_names and maxdiff == 0 and cellsdiff == 0 and ntc_cells == 38176)
P(f"KA-3 the 31 NTC constructs across the three splits: identical names {same_names}, "
  f"max |sum diff| {maxdiff}, max |cell-count diff| {cellsdiff}, {ntc_cells:,} control cells "
  f"-> {'PASS' if ka3 else 'FAIL'}")
out["KA3"] = dict(same_names=bool(same_names), max_abs_sum_diff=maxdiff,
                  max_abs_cell_diff=cellsdiff, ntc_cells=ntc_cells, passed=ka3)

# ---------------------------------------------------------------- KA-7
arms = list(csv.DictReader(open(X.E298 / "seed298.csv")))
d_kd, worst = [], ""
for a in arms:
    c, tgt = a["construct"], a["target"]
    sp = next(SP[s] for s in X.SPLITS if c in SP[s].constructs)
    e = sp.effect(c)[sp.gix[tgt]]
    d = abs(e - float(a["kd_lfc"]))
    if not d_kd or d > max(d_kd):
        worst = f"{c} ({e:+.6f} vs {float(a['kd_lfc']):+.6f})"
    d_kd.append(d)
ka7 = max(d_kd)
P(f"KA-7 293's on-target kd_lfc for all {len(arms)} targeting constructs: max |diff| {ka7:.3e} "
  f"(worst {worst}) -> {'PASS' if ka7 <= 1e-9 else 'FAIL'}")
out["KA7"] = dict(n=len(arms), max_abs_diff=float(ka7), worst=worst, passed=bool(ka7 <= 1e-9))

out["peak_rss_gib"] = round(X.peak_gib(), 3)
out["seconds"] = round(time.time() - t0, 1)
(X.HERE / "ka309.json").write_text(json.dumps(out, indent=1, default=float))
(X.HERE / "ka309.log").write_text("\n".join(log) + "\n")
P(f"done {out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB")
