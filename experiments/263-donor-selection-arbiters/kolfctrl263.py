"""263 -- kolf2 control for the majority rule: is its gain a per-target choice or a disguised GLOBAL down-weight of one donor?
Rows on the kolf2 fold (three production donors, production operator), vs cd4svd5:
  drop:<donor>[:nm]        the two OTHER donors averaged on every three-donor target (a global weight; 010/019's closed axis)
  majority_k562only[:nm]   the majority rule applied only where the odd one out is k562_gwps_pm50
  majority_notk562[:nm]    the majority rule applied only where the odd one out is an Orion line
Plus, from the arrays on disk: the majority arm's per-target d pds split by the donor it dropped.
    PYTHONDONTWRITEBYTECODE=1 uv run python kolfctrl263.py -> work/kolfctrl263.{json,log}"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import lib263 as Y  # noqa: E402

held = "kolf2"
t0 = time.time()
sw = json.loads((Y.WORK / f"sweep263-{held}.json").read_text())
ref = {r["arm"]: r for r in sw["rows"]}
pt = dict(np.load(Y.WORK / f"pertarget263-{held}.npz"))
picks = dict(np.load(Y.WORK / f"picks263-{held}.npz"))
lines, out = [], {}


def P(s=""):
    print(s, flush=True)
    lines.append(s)


# ---- from the arrays: the majority arm's per-target d pds by the donor it dropped
d = pt["majority:raw:all:nm|pds"] - pt["cd4svd5|pds"]
p = picks["majority:raw:all:nm"]
donors = sw["meta"]["donors"]
P("majority:raw:all:nm per-target d pds by the DROPPED donor:")
bydrop = {}
for j in range(3):
    m = p == j
    bydrop[donors[j]] = dict(n=int(m.sum()), sum=float(d[m].sum()), mean=float(d[m].mean()))
    P(f"  dropped {donors[j]:16s}: n {int(m.sum()):3d}  sum {d[m].sum():+.3f}  mean {d[m].mean():+.4f}")
out["majority_by_dropped"] = bydrop

st = Y.stage(held)
ident = Y.X.Identity(st)
n, J = st["n"], 3
avail = st["avail"]
three = avail.all(0)
rows = []


def run(name, W, extra=None, norm_match=False):
    t1 = time.time()
    pre_sel = Y.combine(st, W, norm_match=norm_match)
    wh = Y.base_from(st, pre_sel)
    e, m = Y.evaluate(st, wh, ident)
    b = ref["cd4svd5"]
    r = dict(arm=name, pds=e["pds"], cos_agg=e["cos_agg"], a23_cos=e["a23_cos"], nmae_wald89=m["nmae_wald89"], delivered_norm=m["delivered_norm"],
             d_pds=e["pds"] - b["pds"], rel_agg=100 * (e["cos_agg"] / b["cos_agg"] - 1), rel_macro=100 * (e["a23_cos"] / b["a23_cos"] - 1),
             d_nmae89=m["nmae_wald89"] - b["nmae_wald89"], delivered_ratio=m["delivered_norm"] / b["delivered_norm"],
             n_changed=int((np.abs(pre_sel - st["pre"]).max(1) > 0).sum()), audit=dict(**m["audit"], **(extra or {})), seconds=round(time.time() - t1, 1))
    rows.append(r)
    P(f"  [{held}/{name}] d pds {r['d_pds']:+.5f}  rel cos_agg {r['rel_agg']:+.2f} %  rel macro {r['rel_macro']:+.2f} %  d nmae89 {r['d_nmae89']:+.5f}  "
      f"deliv {r['delivered_ratio']:.4f}  changed {r['n_changed']} ({r['seconds']}s)")
    return r


# KA: the staged base reproduces the sweep's cd4svd5 row
e0, m0 = Y.evaluate(st, st["whitened"], ident)
ka = max(abs(e0[k] - ref["cd4svd5"][k]) for k in ("pds", "cos_agg", "a23_cos"))
P(f"KA staged base == sweep263 cd4svd5: max|d| {ka:.1e} -> {'PASS' if ka == 0.0 else 'FAIL'}")
out["ka"] = ka
P("\nglobal drop-one-donor weights on every three-donor target (n = %d):" % int(three.sum()))
for j in range(J):
    W = np.full((n, J), np.nan)
    W[three] = 1.0
    W[three, j] = 0.0
    for nm in (False, True):
        run(f"drop:{donors[j]}{':nm' if nm else ''}", W, dict(control="global_drop", dropped=donors[j]), norm_match=nm)
P("\nthe majority rule restricted by which donor it drops:")
for label, mask in (("majority_k562only", p == 0), ("majority_notk562", (p == 1) | (p == 2))):
    W = np.full((n, J), np.nan)
    for t in np.flatnonzero(mask):
        W[t] = 1.0
        W[t, p[t]] = 0.0
    for nm in (False, True):
        run(f"{label}{':nm' if nm else ''}", W, dict(control="majority_restricted", n_targets=int(mask.sum())), norm_match=nm)
# reproduce the sweep's majority row through this script's path (a second KA)
W = np.full((n, J), np.nan)
for t in np.flatnonzero(p >= 0):
    W[t] = 1.0
    W[t, p[t]] = 0.0
r = run("majority:raw:all:nm (re-run)", W, norm_match=True)
ka2 = abs(r["pds"] - ref["majority:raw:all:nm"]["pds"])
P(f"KA majority:raw:all:nm re-run == sweep row: |d pds| {ka2:.1e} -> {'PASS' if ka2 == 0.0 else 'FAIL'}")
out["ka_majority_rerun"] = ka2
out["rows"] = rows
out["seconds"] = round(time.time() - t0, 1)
out["peak_rss_gib"] = Y.peak_gib()
(Y.WORK / "kolfctrl263.json").write_text(json.dumps(out, indent=1, default=float))
(Y.WORK / "kolfctrl263.log").write_text("\n".join(lines) + "\n")
P(f"done {out['seconds']}s peak {out['peak_rss_gib']:.2f} GiB")
