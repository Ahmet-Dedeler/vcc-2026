"""450k: the same 2x2 as labline450, for PER-TARGET rows instead of the programme. For each bank pair, on the targets both
screened (the A/B/C panel excluded) and the genes both measure: centre each bank's rows by its own mean row, then the mean
per-target cosine, over all shared targets and over the strongest quartile (by the product of the two row norms).
usage: uv run python experiments/450-same-line-programme/labline_rows450.py"""
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import csv

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import def_sameline as S  # noqa: E402

S.BANK_DIRS.setdefault("orion_hek293t", ROOT / "data/donors/orion_hek293t")
S.BANK_DIRS.setdefault("rpe1", ROOT / "data/donors/rpe1")
META = {"k562_gwps": ("K562", "Weissman"), "k562_essential": ("K562", "Weissman"), "viperturb_k562": ("K562", "Satija"),
        "nadig_jurkat": ("Jurkat", "Weissman"), "nadig_hepg2": ("HepG2", "Weissman"), "rpe1": ("RPE1", "Weissman"),
        "orion_hct116": ("HCT116", "Xaira"), "orion_hek293t": ("HEK293T", "Xaira")}
panel = {r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))}
B = {}
for b in META:
    delta, perts, genes, _ = S.load_bank(b)
    keep = [i for i, p in enumerate(perts) if p not in panel]
    mean = np.zeros(len(genes))
    for s in range(0, len(keep), 1024):
        mean += np.asarray(delta[keep[s:s + 1024]], dtype=np.float64).sum(0)
    mean /= max(len(keep), 1)
    B[b] = dict(delta=delta, perts={p: i for i, p in enumerate(perts) if p not in panel}, genes={g: j for j, g in enumerate(genes)},
                mean=mean)
    print(f"{b:15s} rows {len(keep)}", flush=True)
COMMON = "--common" in sys.argv
T_ALL = set.intersection(*[set(v["perts"]) for v in B.values()]); G_ALL = set.intersection(*[set(v["genes"]) for v in B.values()])
print(f"common to all 8 banks: {len(T_ALL)} targets, {len(G_ALL)} genes; mode {'COMMON' if COMMON else 'pairwise'}", flush=True)
rows = []
for a, b in itertools.combinations(META, 2):
    A_, B_ = B[a], B[b]
    tg = sorted(T_ALL) if COMMON else sorted(set(A_["perts"]) & set(B_["perts"]))
    gg = sorted(G_ALL) if COMMON else sorted(set(A_["genes"]) & set(B_["genes"]))
    if len(tg) < 30:
        continue
    ga = np.array([A_["genes"][g] for g in gg]); gb = np.array([B_["genes"][g] for g in gg])
    Xa = np.asarray(A_["delta"][[A_["perts"][t] for t in tg]], dtype=np.float64)[:, ga] - A_["mean"][ga]
    Xb = np.asarray(B_["delta"][[B_["perts"][t] for t in tg]], dtype=np.float64)[:, gb] - B_["mean"][gb]
    Xa = np.nan_to_num(Xa); Xb = np.nan_to_num(Xb)
    na, nb = np.linalg.norm(Xa, axis=1), np.linalg.norm(Xb, axis=1)
    ok = (na > 0) & (nb > 0)
    cos = (Xa[ok] * Xb[ok]).sum(1) / (na[ok] * nb[ok])
    strength = (na * nb)[ok]
    top = cos[strength >= np.quantile(strength, 0.75)]
    rows.append(dict(a=a, b=b, same_line=META[a][0] == META[b][0], same_lab=META[a][1] == META[b][1], n_targets=int(ok.sum()),
                     n_genes=len(gg), cos_mean=float(cos.mean()), cos_top=float(top.mean())))
print(f"\n{'pair':34s} {'line':5s} {'lab':5s} {'targets':>7s} {'genes':>6s} {'cos all':>8s} {'cos top-q':>9s}")
for x in sorted(rows, key=lambda x: (not x["same_line"], not x["same_lab"], -x["cos_top"])):
    print(f"{x['a'] + ' / ' + x['b']:34s} {'same' if x['same_line'] else 'diff':5s} {'same' if x['same_lab'] else 'diff':5s} "
          f"{x['n_targets']:7d} {x['n_genes']:6d} {x['cos_mean']:+8.3f} {x['cos_top']:+9.3f}")
print("\ngroup means (cos all / cos strongest quartile):")
for sl, sb in ((True, True), (True, False), (False, True), (False, False)):
    g = [x for x in rows if x["same_line"] == sl and x["same_lab"] == sb]
    if g:
        print(f"  line {'same' if sl else 'diff'}, lab {'same' if sb else 'diff'}: pairs {len(g):2d}  "
              f"{np.mean([x['cos_mean'] for x in g]):+.3f} / {np.mean([x['cos_top'] for x in g]):+.3f}")
(Path(__file__).resolve().parent / ("work/labline_rows450" + ("_common" if COMMON else "") + ".json")).write_text(json.dumps(rows, indent=1))
