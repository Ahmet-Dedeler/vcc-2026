"""450h v3: the lab-vs-line 2x2 with the PERTURBATION SET controlled. For each bank pair, both programmes (P_down, 396's
rule) are computed over only the perturbations BOTH banks screened (A/B/C panel excluded), so a shared essential-gene set
cannot masquerade as a shared lab. Pairs sharing < 100 perturbations are skipped.
usage: uv run python experiments/450-same-line-programme/labline450_v3.py"""
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
T225 = ROOT / "experiments/225-source-sweep/work/tables"
for n in ("jiang_K562", "jiang_A549", "jiang_MCF7", "jiang_HT29", "jiang_HAP1", "jiang_BXPC3", "xucao_hek293"):
    S.BANK_DIRS.setdefault(n, T225 / n)
META = {"k562_gwps": ("K562", "W"), "k562_essential": ("K562", "W"), "viperturb_k562": ("K562", "S"),
        "nadig_jurkat": ("Jurkat", "W"), "nadig_hepg2": ("HepG2", "W"), "rpe1": ("RPE1", "W"),
        "orion_hct116": ("HCT116", "X"), "orion_hek293t": ("HEK293", "X"), "jiang_K562": ("K562", "J"),
        "jiang_A549": ("A549", "J"), "jiang_MCF7": ("MCF7", "J"), "jiang_HT29": ("HT29", "J"), "jiang_HAP1": ("HAP1", "J"),
        "jiang_BXPC3": ("BXPC3", "J"), "xucao_hek293": ("HEK293", "C")}
official = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
panel = {r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))}
B = {b: S.load_bank(b) for b in META}
rows = []
for a, b in itertools.combinations(META, 2):
    common = (set(B[a][1]) & set(B[b][1])) - panel
    if len(common) < 100:
        continue
    excl_a = set(B[a][1]) - common; excl_b = set(B[b][1]) - common
    _, pa, na, _ = S.programme(B[a][0], B[a][1], B[a][2], official, panel | excl_a)
    _, pb, nb, _ = S.programme(B[b][0], B[b][1], B[b][2], official, panel | excl_b)
    x, y = pa - 0.5, pb - 0.5
    m = np.isfinite(x) & np.isfinite(y)
    r = float(np.corrcoef(x[m], y[m])[0, 1])
    rows.append(dict(a=a, b=b, same_line=META[a][0] == META[b][0], same_lab=META[a][1] == META[b][1], n_perts=len(common),
                     rows_a=na, rows_b=nb, n_genes=int(m.sum()), r=r))
for x in sorted(rows, key=lambda x: (not x["same_line"], not x["same_lab"], -x["r"])):
    print(f"{x['a'] + ' / ' + x['b']:32s} {'same' if x['same_line'] else 'diff':5s} {'same' if x['same_lab'] else 'diff':5s} "
          f"shared perts {x['n_perts']:5d}  genes {x['n_genes']:6d}  r {x['r']:+.3f}")
print("\ngroup means over pairs sharing >= 100 perturbations:")
for sl, sb in ((True, True), (True, False), (False, True), (False, False)):
    g = [x for x in rows if x["same_line"] == sl and x["same_lab"] == sb]
    if g:
        print(f"  line {'same' if sl else 'diff'}, lab {'same' if sb else 'diff'}: pairs {len(g):2d}  mean r {np.mean([x['r'] for x in g]):+.3f}"
              f"  (median shared perts {int(np.median([x['n_perts'] for x in g]))})")
(Path(__file__).resolve().parent / "work/labline450_v3.json").write_text(json.dumps(rows, indent=1))
