"""357 step 1: inventory JURKAT/HELA genetic signatures, intersect with panels."""
import json
import numpy as np
import pandas as pd

OUT = "experiments/357-jurkat-shrna-truth-proxy/work/"

g = pd.read_parquet("data/external/cmap/genetic_sigs.parquet")
print("genetic_sigs shape:", g.shape)
print("pert_types:", g["pert_type"].value_counts().to_dict())
print("cells (trt_sh):", sorted(g.loc[g.pert_type == "trt_sh", "cell_iname"].unique().tolist()))
print("cells (trt_xpr):", sorted(g.loc[g.pert_type == "trt_xpr", "cell_iname"].unique().tolist()))

# full siginfo for quality columns on the JURKAT/HELA genetic sigs of interest
need = g[g.cell_iname.isin(["JURKAT", "HELA"]) & g.pert_type.isin(["trt_sh", "trt_xpr"])]
need_ids = set(need["sig_id"])
print("JURKAT/HELA genetic sig rows:", len(need))

# panel targets (single panel file; check manifest for per-context panels)
panel = pd.read_csv("data/challenge/pert_counts.csv")["target_gene"].tolist()
print("panel n:", len(panel))
h1 = np.load("experiments/323-current-family-validation/work/reference.npz",
             allow_pickle=True)["targets"].tolist()
print("H1 eval n:", len(h1))

for (cell, pt), sub in need.groupby(["cell_iname", "pert_type"]):
    genes = sorted(sub["cmap_name"].dropna().unique().tolist())
    print(f"\n=== {cell} {pt}: {len(sub)} sigs, {len(genes)} distinct genes ===")
    print("pert_time:", sub["pert_time"].value_counts().to_dict())
    print("is_hiq rate:", float(sub["is_hiq"].mean()), " n_hiq:", int(sub["is_hiq"].sum()))
    print("tas median/IQR:", float(sub["tas"].median()),
          float(sub["tas"].quantile(0.25)), float(sub["tas"].quantile(0.75)))
    for name, universe in [("panel", panel), ("H1-85", h1)]:
        hit = sorted(set(genes) & set(universe))
        print(f"  intersect {name}: {len(hit)}: {hit}")
    # duplicates: multiple sigs per gene?
    dup = sub.groupby("cmap_name")["sig_id"].nunique()
    print("  genes with >1 sig:", int((dup > 1).sum()),
          ("max sigs/gene: %d" % int(dup.max())) if len(dup) else "")

need.to_parquet(OUT + "jurkat_hela_genetic_sigs.parquet", index=False)
print("\nsaved", OUT + "jurkat_hela_genetic_sigs.parquet")
