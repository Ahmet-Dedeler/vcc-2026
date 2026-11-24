"""307 step 3 -- the known-answer rows that touch expression, run BEFORE any candidate outcome.

KA-1  305's 15 published pooled effects and both half effects, recomputed from `sums305_<split>.npz`
      by code written here                                       -> max abs diff <= 1e-9
KA-3  305's finding that the 31 NTC constructs' sums are identical across the three splits
      -> exact integer equality of all 31 constructs' per-half sums
KA-7  293's on-target `kd_lfc` for all 34 targeting constructs recomputed from 305's sums
      -> 34/34 within 1e-9 of `seed298.csv`

    PYTHONDONTWRITEBYTECODE=1 uv run python ka307.py
"""
from __future__ import annotations

import csv
import json
import time
from pathlib import Path

import numpy as np

import lib307 as X

t0 = time.time()
log: list[str] = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    log.append(s)


SP = {s: X.Split(s) for s in X.SPLITS}
P(f"loaded 305's sums: " + ", ".join(
    f"{s} {len(SP[s].constructs)} constructs ({len(SP[s].targeting)} targeting + {len(SP[s].ntc)} NTC)"
    for s in X.SPLITS))

out: dict = {}

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
    # 305's `ctrl_cp50k` is RAW CP50k (50000*S/sum S), not its log1p, and it uses the same
    # leave-that-construct-out reference as the effect for a tested NTC construct.
    ref = sp.ctrl_pooled - (sp.S[c].sum(0) if loo else 0)
    d_pool.append(abs(e - r["e_pool"]))
    d_h0.append(abs(e0 - r["e_half0"]))
    d_h1.append(abs(e1 - r["e_half1"]))
    d_ctrl.append(abs(50000.0 * ref[gi] / ref.sum() - r["ctrl_cp50k"]))
ka1 = max(max(d_pool), max(d_h0), max(d_h1), max(d_ctrl))
P(f"KA-1 305's {len(r305)} published rows reproduced from its saved sums: "
  f"max |diff| pooled {max(d_pool):.3e}, half0 {max(d_h0):.3e}, half1 {max(d_h1):.3e}, "
  f"ctrl CP50k {max(d_ctrl):.3e}  -> {'PASS' if ka1 <= 1e-9 else 'FAIL'}")
out["KA1"] = dict(n_rows=len(r305), max_abs_diff=float(ka1), passed=bool(ka1 <= 1e-9),
                  max_pool=float(max(d_pool)), max_half0=float(max(d_h0)),
                  max_half1=float(max(d_h1)), max_ctrl=float(max(d_ctrl)))

# ---------------------------------------------------------------- KA-3
ntc_sets = {s: set(SP[s].ntc) for s in X.SPLITS}
same_names = ntc_sets["train"] == ntc_sets["validation"] == ntc_sets["test"]
maxdiff, cellsdiff = 0, 0
for c in sorted(ntc_sets["train"]):
    a = SP["train"].S[c]
    for s in ("validation", "test"):
        maxdiff = max(maxdiff, int(np.abs(a - SP[s].S[c]).max()))
        cellsdiff = max(cellsdiff, int(np.abs(SP["train"].cells[c] - SP[s].cells[c]).max()))
ntc_cells = int(sum(SP["train"].cells[c].sum() for c in ntc_sets["train"]))
ka3 = same_names and maxdiff == 0 and cellsdiff == 0
P(f"KA-3 the 31 NTC constructs across the three splits: identical names {same_names}, "
  f"max |sum diff| {maxdiff}, max |cell-count diff| {cellsdiff}, {ntc_cells:,} control cells "
  f"-> {'PASS' if ka3 else 'FAIL'}")
out["KA3"] = dict(same_names=bool(same_names), max_abs_sum_diff=maxdiff,
                  max_abs_cell_diff=cellsdiff, ntc_cells=ntc_cells, passed=bool(ka3))

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
(X.HERE / "ka307.json").write_text(json.dumps(out, indent=1))
(X.HERE / "ka307.log").write_text("\n".join(log) + "\n")
P(f"done {out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB")
