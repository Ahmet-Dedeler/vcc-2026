"""425: are genes near a CRISPRi target's TSS systematically repressed in the target's cells (H1 2025 reference DE)?

Reads the h1-test reference DE table (323's reference_de.parquet: target, feature, p_value, p_adj, log2_fold_change,
the scorer's Wilcoxon on arithmetic-mean CPM) and 425's TSS table. For every (target, feature) pair on the same
chromosome, bins |TSS distance| and reports: n pairs, share significant (p_adj < 0.05), mean LFC, share of significant
calls that are DOWN, and the share of the nmae denominator (sum |lfc| over the gate, gate = significant, target removed,
targets with >= 10 gate genes) carried by each bin.
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
de = pd.read_parquet(ROOT / "experiments/323-current-family-validation/work/reference_de.parquet")
tss = pd.read_csv(HERE / "work/tss_v44.tsv", sep="\t").set_index("gene_name")

de = de[de["target"] != de["feature"]].copy()
t = tss.reindex(de["target"].to_numpy())
g = tss.reindex(de["feature"].to_numpy())
same = (t["chrom"].to_numpy() == g["chrom"].to_numpy())
dist = np.where(same, np.abs(t["tss"].to_numpy(dtype=float) - g["tss"].to_numpy(dtype=float)), np.inf)
de["dist"] = dist
de["sig"] = de["p_adj"] < 0.05
bins = [0, 1e3, 5e3, 2e4, 1e5, 1e6, np.inf]
labels = ["<1kb", "1-5kb", "5-20kb", "20-100kb", "0.1-1Mb", ">1Mb/other chr"]
de["bin"] = pd.cut(de["dist"], bins=bins, labels=labels, right=False, include_lowest=True)
print(f"targets {de['target'].nunique()}  pairs {len(de)}  TSS-resolved targets "
      f"{t['chrom'].notna().groupby(de['target'].to_numpy()).any().sum()}")

gate = de[de["sig"] & np.isfinite(de["log2_fold_change"])]
gsize = gate.groupby("target").size()
keep = gsize[gsize >= 10].index
gate = gate[gate["target"].isin(keep)]
den = gate.groupby("target")["log2_fold_change"].apply(lambda s: s.abs().sum())

rows = []
for b in labels:
    x = de[de["bin"] == b]
    s = x[x["sig"]]
    gb = gate[gate["bin"] == b]
    share = (gb.groupby("target")["log2_fold_change"].apply(lambda s: s.abs().sum()).reindex(keep).fillna(0) / den)
    rows.append(dict(bin=b, pairs=len(x), sig_rate=x["sig"].mean(), mean_lfc=x["log2_fold_change"].mean(),
                     sig_down=(s["log2_fold_change"] < 0).mean() if len(s) else np.nan,
                     mean_lfc_sig=s["log2_fold_change"].mean() if len(s) else np.nan,
                     nmae_den_share_mean=share.mean(), nmae_den_share_max=share.max()))
out = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print(out.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

near = de[(de["dist"] < 2e4) & de["sig"]].sort_values("dist")
print("\nsignificant pairs within 20 kb:")
print(near[["target", "feature", "dist", "log2_fold_change", "p_adj"]].to_string(index=False))
out.to_csv(HERE / "work/h1cis425_bins.csv", index=False)
