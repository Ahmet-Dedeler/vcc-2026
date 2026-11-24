"""450h (2026-09-26, after 446bja's null): is a line's shared response programme a property of the LINE or of the LAB /
protocol? Per-gene P_down (def_sameline's 396 rule, the A/B/C panel's rows excluded, official gene axis) for eight banks
from three labs, then pairwise agreement grouped by (same line?, same lab?).
  Weissman: K562-GWPS, K562-essential, Nadig Jurkat, Nadig HepG2, RPE1.  Xaira: Orion HCT116, Orion HEK293T.  Satija:
  VIPerturb K562.
Agreement: Pearson r of (P_down - 0.5) on genes both measure; sign agreement on genes where both |P_down - 0.5| >= 0.05.
usage: uv run python experiments/450-same-line-programme/labline450.py"""
import csv
import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import def_sameline as S  # noqa: E402

S.BANK_DIRS.setdefault("orion_hek293t", ROOT / "data/donors/orion_hek293t")
S.BANK_DIRS.setdefault("rpe1", ROOT / "data/donors/rpe1")
META = {"k562_gwps": ("K562", "Weissman"), "k562_essential": ("K562", "Weissman"), "viperturb_k562": ("K562", "Satija"),
        "nadig_jurkat": ("Jurkat", "Weissman"), "nadig_hepg2": ("HepG2", "Weissman"), "rpe1": ("RPE1", "Weissman"),
        "orion_hct116": ("HCT116", "Xaira"), "orion_hek293t": ("HEK293T", "Xaira")}
official = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
panel = {r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))}
P = {}
for b in META:
    delta, perts, genes, _ = S.load_bank(b)
    psign, pdown, nrows, _ = S.programme(delta, perts, genes, official, panel)
    P[b] = pdown
    print(f"{b:15s} rows {nrows:6d}  genes with P_down {int(np.isfinite(pdown).sum())}  mean P_down {np.nanmean(pdown):.3f}", flush=True)
rows = []
for a, b in itertools.combinations(META, 2):
    x, y = P[a] - 0.5, P[b] - 0.5
    m = np.isfinite(x) & np.isfinite(y)
    r = float(np.corrcoef(x[m], y[m])[0, 1])
    s = m & (np.abs(x) >= 0.05) & (np.abs(y) >= 0.05)
    agree = float(np.mean(np.sign(x[s]) == np.sign(y[s])))
    same_line, same_lab = META[a][0] == META[b][0], META[a][1] == META[b][1]
    rows.append(dict(a=a, b=b, same_line=same_line, same_lab=same_lab, n=int(m.sum()), r=r, n_sign=int(s.sum()), sign=agree))
print(f"\n{'pair':34s} {'line':5s} {'lab':5s} {'n':>6s} {'r':>7s} {'n_sgn':>6s} {'sign':>6s}")
for x in sorted(rows, key=lambda x: (not x["same_line"], not x["same_lab"], -x["r"])):
    print(f"{x['a'] + ' / ' + x['b']:34s} {'same' if x['same_line'] else 'diff':5s} {'same' if x['same_lab'] else 'diff':5s} "
          f"{x['n']:6d} {x['r']:+7.3f} {x['n_sign']:6d} {x['sign']:6.3f}")
print("\ngroup means (r / sign):")
for sl, sb in ((True, True), (True, False), (False, True), (False, False)):
    g = [x for x in rows if x["same_line"] == sl and x["same_lab"] == sb]
    if g:
        print(f"  line {'same' if sl else 'diff'}, lab {'same' if sb else 'diff'}: n pairs {len(g):2d}  r {np.mean([x['r'] for x in g]):+.3f}"
              f"  sign {np.mean([x['sign'] for x in g]):.3f}")
(Path(__file__).resolve().parent / "work/labline450.json").write_text(json.dumps(rows, indent=1))
