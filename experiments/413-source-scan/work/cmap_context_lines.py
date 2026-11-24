"""413 (criterion b): CMap/LINCS 2020 genetic signatures (on disk) — panel targets per cell line, focus on HeLa/Jurkat/HNSCC and epithelial carcinomas."""
import json, pandas as pd
from overlap import to_current, PC
g = pd.read_parquet("../../../data/external/cmap/genetic_sigs.parquet")
print(g.shape, g.columns.tolist()[:20])
gene_col = [c for c in g.columns if c in ("cmap_name", "pert_iname", "gene")][0]
g["cur"] = g[gene_col].astype(str).map(to_current); g["panel"] = g.cur.isin(set(PC))
ci = pd.read_csv("../../../data/external/cmap/cellinfo_beta.txt", sep="\t")
tab = g[g.panel].groupby(["cell_iname", "pert_type"]).cur.nunique().unstack(fill_value=0)
tab["any"] = g[g.panel].groupby("cell_iname").cur.nunique()
tab = tab.join(ci.set_index("cell_iname")[["primary_disease", "cell_lineage"]], how="left").sort_values("any", ascending=False)
print(tab.head(40).to_string())
lines = sorted(g.cell_iname.unique()); print(len(lines), "lines with genetic sigs")
ctx = [l for l in lines if any(k in l.upper() for k in ["HELA", "JURKAT", "CAL", "FADU", "SCC", "DETROIT", "HN"])]
print("context-adjacent lines present:", ctx)
json.dump({"per_line": tab.reset_index().to_dict(orient="records"), "context_lines": ctx}, open("cmap_context_lines.json", "w"), indent=1, default=str)
