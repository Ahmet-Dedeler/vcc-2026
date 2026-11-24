"""263 -- extra diagnostics from the arrays on disk (no fold load): (1) the one positive single-fold reading, the
three-donor MAJORITY rule on kolf2 (per-target structure, retention, controls); (2) the abstaining oracle (how often
the mean IS the best of the three options); (3) every cell with a positive 3-fold mean on the cancer folds.
    PYTHONDONTWRITEBYTECODE=1 uv run python extra263.py  -> work/extra263.{json,log}"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
FOLDS = CANCER + ("kolf2",)
BASE = "cd4svd5"
SIGMA = 0.001494
rows, pt, picks, cov, meta = {}, {}, {}, {}, {}
for h in FOLDS:
    d = json.loads((WORK / f"sweep263-{h}.json").read_text())
    meta[h] = d["meta"]
    for r in d["rows"]:
        rows[(h, r["arm"])] = r
    pt[h] = dict(np.load(WORK / f"pertarget263-{h}.npz"))
    picks[h] = dict(np.load(WORK / f"picks263-{h}.npz"))
    cov[h] = dict(np.load(WORK / f"cov263-{h}.npz"))
lines, out = [], {}


def P(s=""):
    print(s, flush=True)
    lines.append(s)


def dp(h, a):
    return pt[h][f"{a}|pds"] - pt[h][f"{BASE}|pds"]


# ------------------------------------------------------------------ 1. the majority rule on kolf2
P("=== 1. kolf2 (apart): the three-donor MAJORITY rule (keep the two donors whose un-whitened rows agree, drop the odd one)")
h = "kolf2"
t = np.asarray([str(x) for x in meta[h]["targets"]])
b262 = cov[h]["base262_cos_truth"]
groups = {"(e) base·truth <= 0": b262 <= 0, "marginal (0, 0.05]": (b262 > 0) & (b262 <= 0.05), "saturated > 0.05": b262 > 0.05}
maj = {}
for a in ("majority:raw:all:nm", "majority:raw:all", "majority:raw:tau0:nm", "majority:wht:all:nm", "majority:raw:all:rev:nm"):
    d = dp(h, a)
    p = picks[h][a]
    sel = p >= 0
    order = np.argsort(-d)
    keep = np.ones(d.size, bool)
    keep[order[:3]] = False
    ret = d[keep].mean() / d.mean() if d.mean() != 0 else float("nan")
    by_group = {g: dict(n=int(m.sum()), sum=float(d[m].sum()), mean=float(d[m].mean()), selected=int((sel & m).sum())) for g, m in groups.items()}
    dropped = {int(j): int((p == j).sum()) for j in range(3)}
    maj[a] = dict(d_pds=float(d.mean()), n_selected=int(sel.sum()), gainers=int((d > 1e-12).sum()), losers=int((d < -1e-12).sum()),
                  retention_lt3o=float(ret), top3=[f"{t[i]} {d[i]:+.3f}" for i in order[:3]], bottom3=[f"{t[i]} {d[i]:+.3f}" for i in order[::-1][:3]],
                  by_group=by_group, dropped_donor_counts=dropped, sigma_units=float(d.mean()) / SIGMA,
                  rel_cos_agg=100 * (rows[(h, a)]["cos_agg"] / rows[(h, BASE)]["cos_agg"] - 1), rel_macro=100 * (rows[(h, a)]["a23_cos"] / rows[(h, BASE)]["a23_cos"] - 1),
                  d_nmae89=rows[(h, a)]["nmae_wald89"] - rows[(h, BASE)]["nmae_wald89"], delivered_ratio=rows[(h, a)]["delivered_norm"] / rows[(h, BASE)]["delivered_norm"])
    m = maj[a]
    P(f"  {a:28s} d pds {m['d_pds']:+.5f} ({m['sigma_units']:+.1f} sigma-units, ONE fold) on {m['n_selected']} targets; gainers {m['gainers']} / losers {m['losers']}; "
      f"LT3O retention {100*m['retention_lt3o']:+.0f} %; top3 {m['top3']}; bottom3 {m['bottom3']}")
    P(f"      by group: " + "; ".join(f"{g}: n {v['n']} sel {v['selected']} sum {v['sum']:+.3f} mean {v['mean']:+.4f}" for g, v in by_group.items()) +
      f"; dropped donor counts {dropped} (0 = {meta[h]['donors'][0]}, 1 = {meta[h]['donors'][1]}, 2 = {meta[h]['donors'][2]})")
    P(f"      rel cos_agg {m['rel_cos_agg']:+.2f} %, rel macro {m['rel_macro']:+.2f} %, d nmae89 {m['d_nmae89']:+.5f}, delivered ratio {m['delivered_ratio']:.4f}")
# which donor does the majority drop, and does dropping it agree with the oracle (the oracle's pick is one of the KEPT two)?
o = picks[h]["oracle_wht:all:nm"]
p = picks[h]["majority:raw:all:nm"]
m = (p >= 0) & (o >= 0)
kept_contains_oracle = float((p[m] != o[m]).mean())
P(f"  majority keeps the oracle's donor among its two on {100*kept_contains_oracle:.0f} % of {int(m.sum())} targets (chance 67 %)")
# per-donor: pds of "always drop donor j" is not measured; the majority's dropped-donor histogram tells if it is a disguised global weight
out["kolf2_majority"] = dict(cells=maj, keeps_oracle_donor=kept_contains_oracle, n=int(m.sum()))
# the same cell's split by the min pairwise donor-donor cosine
dd = cov[h]["dd_raw"]
d = dp(h, "majority:raw:all:nm")
P("  majority:raw:all:nm by min pairwise raw donor-donor cosine bin: " + "; ".join(
    f"{lo:+.2f}..{hi:+.2f}: n {int(((dd > lo) & (dd <= hi)).sum())} mean {d[(dd > lo) & (dd <= hi)].mean():+.4f}" for lo, hi in ((-1, -0.05), (-0.05, 0), (0, 0.05), (0.05, 0.1), (0.1, 1))))

# ------------------------------------------------------------------ 2. the abstaining oracle
P("\n=== 2. The abstaining oracle (best of {donor rows, the mean row} by whitened cosine to the truth): how often the mean wins")
ab = {}
for h in FOLDS:
    o3 = picks[h]["oracle3_wht:all:nm"]
    o = picks[h]["oracle_wht:all:nm"]
    two = cov[h]["avail"].sum(0) >= 2
    n2 = int(two.sum())
    mean_best = int(((o3 < 0) & two).sum())
    ab[h] = dict(n_two_donors=n2, mean_is_best=mean_best, frac=mean_best / n2, d_oracle=float(dp(h, "oracle_wht:all:nm").mean()),
                 d_oracle3=float(dp(h, "oracle3_wht:all:nm").mean()))
    P(f"  {h:14s} two-donor targets {n2}: the mean row is the best of the three on {mean_best} ({100*mean_best/n2:.0f} %); "
      f"oracle {ab[h]['d_oracle']:+.5f} -> abstaining oracle {ab[h]['d_oracle3']:+.5f}")
out["abstaining_oracle"] = ab

# ------------------------------------------------------------------ 3. positive cells on the cancer folds
P("\n=== 3. Every non-oracle cell with a POSITIVE 3-fold mean d pds on the cancer folds")
arms = sorted({a for (hh, a) in rows if hh == CANCER[0] and all((h2, a) in rows for h2 in CANCER)})
pos = []
for a in arms:
    if a.startswith(("oracle", BASE, "mean_explicit")):
        continue
    v = [float(dp(h2, a).mean()) for h2 in CANCER]
    if np.mean(v) > 0:
        pos.append((float(np.mean(v)), a, v))
pos.sort(reverse=True)
for mval, a, v in pos:
    P(f"  {a:30s} {mval:+.5f} ({mval/SIGMA:+.1f} sigma) per fold {['%+.4f' % x for x in v]}  {sum(x > 0 for x in v)}/3")
if not pos:
    P("  none")
out["positive_cells_cancer"] = [dict(arm=a, d_pds=mval, per_fold=v) for mval, a, v in pos]
(WORK / "extra263.json").write_text(json.dumps(out, indent=1, default=float))
(WORK / "extra263.log").write_text("\n".join(lines) + "\n")
