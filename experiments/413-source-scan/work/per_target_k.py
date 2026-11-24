"""413: per-panel-target count of held donor lines (K562-GWPS, Orion HCT116, Orion HEK293T, CD4 any condition), and the
increment each new candidate would add (a new LINE adds +1 on its panel hits; a second K562 screen adds a line only where
K562-GWPS is absent)."""
import json, pandas as pd
from overlap import PANEL, to_current, PC, ROOT
held = {d: {to_current(l.strip()) for l in open(ROOT / f"data/donors/{d}/perts.txt")} for d in ["k562_gwps", "orion_hct116", "orion_hek293t"]}
cd4 = pd.read_csv(ROOT / "data/external/cd4t/DE_stats.suppl_table.csv", usecols=["target_contrast_gene_name", "n_cells_target"])
held["cd4"] = {to_current(g) for g in cd4.loc[cd4.n_cells_target >= 20, "target_contrast_gene_name"].astype(str)}
pc = {to_current(p): p for p in PANEL}
k = pd.Series({p: sum(c in s for s in held.values()) for c, p in pc.items()})
print("held-line count per panel target:", k.value_counts().sort_index().to_dict())
cands = {
 "VIPerturb K562 (>=20 cells)": json.load(open("viperturb_overlap.json"))["ge20"]["panel_hits"],
 "GSE345058 A549 KO": json.load(open("gse345058_overlap.json"))["cells_ge20"]["panel_hits"],
 "DLD-1 (held, reference)": json.load(open("dld1_arms_overlap.json"))["Low1"]["panel_hits"],
}
out = {"held_k_hist": {int(a): int(b) for a, b in k.value_counts().sort_index().items()}}
for name, hits in cands.items():
    hits = [h for h in hits if h in k.index]
    if name.startswith("VIPerturb"):
        new_line_for = [h for h in hits if to_current(h) not in held["k562_gwps"]]
    else:
        new_line_for = hits
    out[name] = {"panel_hits": len(hits), "adds_a_line_for": len(new_line_for),
                 "of_which_currently_k<=2": int((k[new_line_for] <= 2).sum()), "of_which_currently_k<=1": int((k[new_line_for] <= 1).sum()),
                 "targets_k_le_2": sorted(k[new_line_for][k[new_line_for] <= 2].index.tolist())}
    print(name, out[name])
low = sorted(k[k <= 2].index.tolist()); print("panel targets with <=2 held lines:", len(low), low)
out["panel_targets_k_le_2"] = low
json.dump(out, open("per_target_k.json", "w"), indent=1)
