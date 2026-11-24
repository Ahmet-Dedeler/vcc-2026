"""413: GSE345058 (CZ Biohub A549 1,000-gene CRISPR-KO multimodal atlas) panel overlap from library and from cells."""
import gzip, json, pandas as pd
from overlap import report
lib = pd.read_csv("dl/GSE345058_sgRNA.csv.gz")
print(lib.design.value_counts().to_dict(), lib.source.value_counts().head().to_dict())
tg = lib.loc[lib.design == "targeting", "gene_symbol"].dropna().tolist()
cells = pd.read_csv("dl/GSE345058_cell_annotations.csv.gz", usecols=["Batch", "Perturbation", "Condition", "total_counts", "n_genes_by_counts"])
print(cells.shape, cells.Condition.value_counts().to_dict(), cells.Batch.value_counts().to_dict())
print("median UMI", cells.total_counts.median(), "median genes", cells.n_genes_by_counts.median())
vc = cells.Perturbation.value_counts()
print("top labels", vc.head(8).to_dict())
out = {"library": report("GSE345058 library", tg),
       "cells_any": report("GSE345058 cells>=1", vc.index.tolist()),
       "cells_ge20": report("GSE345058 cells>=20", vc[vc >= 20].index.tolist()),
       "cells_ge50": report("GSE345058 cells>=50", vc[vc >= 50].index.tolist())}
hits = out["cells_ge20"]["panel_hits"]
out["cells_per_panel_target"] = {g: int(vc.get(g, 0)) for g in hits}
out["median_cells_per_target_all"] = float(vc[vc.index.isin(set(tg))].median())
for k, v in out.items():
    if isinstance(v, dict) and "panel_exact" in v: print(k, v["n_unique"], v["panel_exact"], v["panel_resolved"])
print("panel hits", hits)
print("median cells/target", out["median_cells_per_target_all"], "median cells/panel target", pd.Series(out["cells_per_panel_target"]).median())
json.dump(out, open("gse345058_overlap.json", "w"), indent=1)
