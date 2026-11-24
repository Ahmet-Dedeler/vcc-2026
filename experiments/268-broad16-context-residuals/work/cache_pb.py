"""Parse the Broad-16 pseudobulk CSV once into an npz cache (in work/), and verify 111's counts."""
import json, time, numpy as np, pandas as pd
from pathlib import Path
t0 = time.time()
SRC = Path("~/Code/virtual-cell-challenge-2026/data/external/depmap26q1/extracted")
pb = pd.read_csv(SRC / "processed/pseudobulk_sum_by_ko.csv", header=[0, 1], index_col=0, low_memory=False)
print("pseudobulk shape", pb.shape, f"{time.time()-t0:.0f}s")
lines = [c[0] for c in pb.columns]; targets = [c[1] for c in pb.columns]
genes = np.array(pb.index.astype(str))
X = pb.to_numpy(dtype=np.float32)
print("genes unique", len(set(genes)) == len(genes), "min", float(np.nanmin(X)), "nan", int(np.isnan(X).sum()),
      "integral", bool(np.all(np.mod(X[~np.isnan(X)], 1) == 0)))
meta = pd.read_csv(SRC / "metadata/cell_line_metadata.csv")
cq = pd.read_csv(SRC / "processed/cell_quality_covariates.csv")
ct = pd.read_csv(SRC / "processed/crispr_table.csv")
np.savez_compressed("pb_cache.npz", X=X, genes=genes, lines=np.array(lines), targets=np.array(targets))
uniq_lines = sorted(set(lines)); uniq_t = sorted(set(targets))
ors = [t for t in uniq_t if t.startswith("OR")]
print("lines", len(uniq_lines), "targets", len(uniq_t), "OR controls", len(ors), "columns", X.shape[1])
present = {(l, t) for l, t in zip(lines, targets)}
missing = [(l, t) for l in uniq_lines for t in uniq_t if (l, t) not in present]
print("missing pairs", missing)
complete = [t for t in uniq_t if t not in ors and all((l, t) in present for l in uniq_lines)]
print("complete non-control targets", len(complete))
depth = X.sum(0)
print("column depth (total UMIs): min %.3g median %.3g max %.3g" % (depth.min(), np.median(depth), depth.max()))
panel = pd.read_csv("~/Code/virtual-cell-challenge-2026/data/challenge/pert_counts.csv")["target_gene"].tolist()
print("VCC panel", len(panel), "overlap with Broad targets", sorted(set(panel) & set(uniq_t)), "panel genes on response axis", len(set(panel) & set(genes)))
print(meta.to_string()); print(cq.to_string()); print(ct.head(3).to_string(), len(ct))
print("guides per line in crispr_table:", ct.groupby("cell_line").size().to_dict())
json.dump({"lines": uniq_lines, "targets": uniq_t, "controls": ors, "complete": complete, "missing": missing,
           "n_genes": int(len(genes)), "depth_by_column": {f"{l}|{t}": float(d) for l, t, d in zip(lines, targets, depth)}},
          open("inventory.json", "w"), indent=1)
print(f"done {time.time()-t0:.0f}s")
