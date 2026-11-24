"""413 (criterion b): panel targets with ANY public bulk knockdown/knockout RNA-seq group in the challenge-context lines
(HeLa, Jurkat, CAL33/other HNSCC) and other epithelial carcinoma lines, from GPSAdb 2.0 metadata saved by exp 270."""
import json, re, pandas as pd
from overlap import to_current, PC
m = pd.read_csv("../../270-data-reconnaissance/work/gpsadb_metadata_merged.csv")
m["cur"] = m.gene.astype(str).map(to_current)
m["panel"] = m.cur.isin(set(PC))
ctx = {"HeLa": r"^hela", "Jurkat": r"jurkat", "HNSCC(CAL27/CAL33/FaDu/SCC/Detroit)": r"cal[- ]?(27|33)|fadu|detroit|scc[- ]?\d|hn\d|upci|sas$|hsc[- ]?\d|oecm|tu[- ]?177|ca9[- ]?22"}
out = {}
for name, pat in ctx.items():
    sub = m[m.cell_line.astype(str).str.lower().str.contains(pat, regex=True)]
    out[name] = {"groups": int(len(sub)), "lines": sorted(sub.cell_line.unique().tolist())[:12], "genes": int(sub.cur.nunique()),
                 "panel_targets": sorted(sub.loc[sub.panel, "cur"].unique().tolist()), "methods": sub.loc[sub.panel, "method"].value_counts().to_dict()}
    print(name, out[name]["groups"], "groups;", out[name]["genes"], "genes;", len(out[name]["panel_targets"]), "panel:", out[name]["panel_targets"], out[name]["methods"], out[name]["lines"])
# epithelial carcinoma lines ranked by panel coverage
pl = m[m.panel].groupby("cell_line").cur.nunique().sort_values(ascending=False)
print("\nTop lines by #panel targets with any GPSAdb group:\n", pl.head(15).to_string())
out["top_lines_by_panel"] = pl.head(25).to_dict()
json.dump(out, open("gpsadb_context_lines.json", "w"), indent=1)
