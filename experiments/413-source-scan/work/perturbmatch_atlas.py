"""413: LLM-built GEO bulk perturbation atlas (Soul & Young 2026; perturbMatchData, Zenodo 21792179) — panel coverage overall,
in challenge-context lines (HeLa, Jurkat, HNSCC) and top lines; knockdown/knockout only (overexpression excluded)."""
import json, pandas as pd
from overlap import to_current, PC, report
d = pd.read_csv("dl/perturbmatch_atlas_metadata.csv")
d = d[d.db == d.db.iloc[0]]  # one row per signature per method: keep a single method table
print("signatures", len(d), "types", d.pert_type.value_counts().to_dict())
lof = d[d.pert_type.isin(["KD", "KO"]) | d.pert_type.str.contains("KD|KO|knock", case=False, na=False)].copy()
lof["cur"] = lof.pert_symbol.astype(str).map(to_current); lof["panel"] = lof.cur.isin(set(PC))
out = {"n_signatures": int(len(d)), "n_lof": int(len(lof)), "genes_lof": int(lof.cur.nunique()),
       "panel_any_line": int(lof.loc[lof.panel, "cur"].nunique())}
lab = (lof.cell_line_efo_label.fillna(lof.cell_line).astype(str)).str.lower()
ctx = {"HeLa": r"^hela", "Jurkat": r"jurkat", "HNSCC": r"cal[- ]?(27|33)|fadu|detroit|scc[- ]?\d|hn\d|upci|^sas$|hsc-?\d|tu ?177|ca9-?22|oral|head and neck"}
for k, pat in ctx.items():
    sub = lof[lab.str.contains(pat, regex=True)]
    out[k] = {"signatures": int(len(sub)), "genes": int(sub.cur.nunique()), "panel_targets": sorted(sub.loc[sub.panel, "cur"].unique().tolist())}
    print(k, out[k])
top = lof[lof.panel].assign(line=lab[lof.panel]).groupby("line").cur.nunique().sort_values(ascending=False).head(12)
print("panel targets per line:\n", top.to_string()); out["top_lines"] = top.to_dict()
print("panel targets with any LoF signature:", out["panel_any_line"])
json.dump(out, open("perturbmatch_atlas_overlap.json", "w"), indent=1)
