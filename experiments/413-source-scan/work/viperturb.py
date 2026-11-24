"""413: VIPerturb-seq (Bradu et al. 2026, Zenodo 18460279) genome-wide K562 CRISPRi — panel overlap from the manifest."""
import json, pandas as pd
from overlap import report, P, PC, to_current
m = pd.read_csv("dl/viperturb_genome_wide_manifest.txt", sep=r"\s+")
print(m.shape, m.columns.tolist()); print(m.file.value_counts().to_dict(), m.is_in_filtered_obj.value_counts().to_dict())
print("cells/gene median", m.cell_count.median(), "total cells", m.cell_count.sum())
m["cur"] = m.gene.map(to_current)
pan = m[m.cur.isin(set(PC))]
out = {"manifest_all": report("VIPerturb all", m.gene.tolist()),
       "ge20": report("VIPerturb >=20 cells", m.loc[m.cell_count >= 20, "gene"].tolist()),
       "ge50": report("VIPerturb >=50 cells", m.loc[m.cell_count >= 50, "gene"].tolist()),
       "filtered_obj": report("VIPerturb filtered obj", m.loc[m.is_in_filtered_obj == True, "gene"].tolist()),
       "panel_cells_median": float(pan.cell_count.median()), "panel_in_filtered": int(pan.is_in_filtered_obj.sum())}
for k, v in out.items():
    print(k, (v["n_unique"], v["panel_exact"], v["panel_resolved"]) if isinstance(v, dict) else v)
k562 = set(l.strip() for l in open("../../../data/donors/k562_gwps/perts.txt"))
print("panel targets in VIPerturb>=20 but NOT in k562_gwps donor:", sorted(set(out["ge20"]["panel_hits"]) - k562))
out["panel_ge20_not_in_k562gwps"] = sorted(set(out["ge20"]["panel_hits"]) - k562)
json.dump(out, open("viperturb_overlap.json", "w"), indent=1)
