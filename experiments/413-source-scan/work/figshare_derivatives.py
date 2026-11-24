"""413: panel overlap of (a) the third-party VIPerturb-seq mean-profile export (figshare 33719593/33732421) and
(b) Perturb-ME A375 (figshare 33264225): library (MAGeCK gene summary) and regression-model targets (elastic-net beta)."""
import json, gzip, pandas as pd
from overlap import report
out = {}
m1 = json.load(open("dl/vip_export_meta_v1.json"))
cells = pd.Series(m1["total_cells_available"], index=m1["targets"])
out["vip_export_v1"] = report("VIPerturb export v1 (>=20 cells)", m1["targets"]) | {"n_genes": m1["n_genes"],
    "panel_cells_median": float(cells[cells.index.isin(report("x", m1["targets"])["panel_hits"])].median())}
t3 = [l.strip() for l in open("dl/vip_export_targets_v3.txt") if l.strip()]
m3 = json.load(open("dl/vip_export_meta_v3.json"))
out["vip_export_v3_targetlist"] = report("VIPerturb export v3 supplied target list", t3)
out["vip_export_v3_fulldepth"] = report("VIPerturb export v3 full-depth", m3.get("targets_full", m3.get("full_targets", []))) if any(k for k in m3 if "full" in k) else None
print("v3 meta keys:", [k for k in m3 if not isinstance(m3[k], list)], [k for k in m3 if isinstance(m3[k], list)])
g = pd.read_csv("dl/perturbme_mageck_gene_summary.txt", sep="\t")
out["perturbme_library"] = report("Perturb-ME library (MAGeCK)", g["id"].astype(str).tolist())
b = pd.read_csv("dl/perturbme_beta.csv.gz", nrows=5)
print("beta cols", b.columns[:6].tolist(), b.shape)
bi = pd.read_csv("dl/perturbme_beta.csv.gz", usecols=[0])
out["perturbme_model_targets_rows"] = report("Perturb-ME beta rows", bi.iloc[:, 0].astype(str).tolist())
out["perturbme_model_targets_cols"] = report("Perturb-ME beta cols", [c.split("_")[0] for c in pd.read_csv("dl/perturbme_beta.csv.gz", nrows=0).columns[1:]])
for k, v in out.items():
    if v: print(k, v["n_unique"], v["panel_exact"], v["panel_resolved"], {kk: vv for kk, vv in v.items() if kk not in ("panel_hits", "source", "n_unique", "panel_exact", "panel_resolved")})
json.dump(out, open("figshare_derivatives_overlap.json", "w"), indent=1)
