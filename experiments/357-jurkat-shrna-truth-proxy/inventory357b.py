"""357 step 1b: full inventory from siginfo_beta.txt (chunked, 1.2M rows)."""
import numpy as np
import pandas as pd

OUT = "experiments/357-jurkat-shrna-truth-proxy/work/"
usecols = ["sig_id", "pert_id", "pert_type", "cell_iname", "cmap_name",
           "pert_time", "pert_dose", "pert_itime", "pert_idose",
           "tas", "cc_q75", "ss_ngene", "is_hiq", "qc_pass", "is_exemplar_sig"]

panel = pd.read_csv("data/challenge/pert_counts.csv")["target_gene"].tolist()
panelset = set(panel)
h1 = set(np.load("experiments/323-current-family-validation/work/reference.npz",
                 allow_pickle=True)["targets"].tolist())

frames = []
for ch in pd.read_csv("data/external/cmap/siginfo_beta.txt", sep="\t",
                       usecols=usecols, chunksize=500000):
    j = ch[ch.cell_iname.isin(["JURKAT", "HELA"]) & ch.pert_type.isin(["trt_sh", "trt_xpr"])]
    if len(j):
        frames.append(j)
jk = pd.concat(frames, ignore_index=True)
print("JURKAT/HELA trt_sh/trt_xpr rows:", len(jk))
print(jk.groupby(["cell_iname", "pert_type"]).size())

for (cell, pt), sub in jk.groupby(["cell_iname", "pert_type"]):
    print(f"\n=== {cell} {pt}: {len(sub)} sigs ===")
    print("cmap_name sample:", sub["cmap_name"].dropna().unique()[:8].tolist())
    genes = sorted(sub["cmap_name"].dropna().unique().tolist())
    print("distinct cmap_name:", len(genes))
    print("pert_id distinct:", sub["pert_id"].nunique())
    print("pert_time:", sub["pert_time"].value_counts().to_dict())
    print("pert_dose:", sub["pert_dose"].value_counts().to_dict())
    print("is_hiq:", int(sub["is_hiq"].sum()), "/", len(sub),
          " qc_pass:", int(sub["qc_pass"].sum()), "/", len(sub))
    for c in ["tas", "cc_q75", "ss_ngene"]:
        print(f"  {c}: median {sub[c].median():.4f} "
              f"q25 {sub[c].quantile(0.25):.4f} q75 {sub[c].quantile(0.75):.4f}")
    hp = sorted(set(genes) & panelset)
    hh = sorted(set(genes) & h1)
    print(f"  intersect PANEL-A/B/C ({len(hp)}): {hp}")
    print(f"  intersect H1-85 ({len(hh)}): {hh}")

jk.to_parquet(OUT + "jurkat_hela_siginfo.parquet", index=False)
print("\nsaved", OUT + "jurkat_hela_siginfo.parquet")
