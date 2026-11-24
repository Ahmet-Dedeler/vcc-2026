"""Build the promoter-accessibility matrix and its matched basal-expression matrix.

Downloads: ENCODE DNase-seq GRCh38 `bed narrowPeak` peak calls, up to TWO independent
experiments per cell line (the second is used only for the split-experiment reliability
control). Everything lands in this experiment's `work/peaks/`. Nothing goes into `data/`.

Frozen construction (declared before any correlation was computed):
  * gene set: Ensembl GRCh38 protein-coding genes with a unique HGNC symbol on chr1-22/X/Y,
    from `E301/work/biomart_genes.tsv` (read-only);
  * TSS = gene start on the + strand, gene end on the - strand; promoter window = TSS +/- 1000 bp;
  * a[line, gene] = max narrowPeak signalValue over peaks intersecting the window, 0 if none;
    stored as log1p(a);
  * e[line, gene] = log2(1 + CCLE TPM) from `ccle_cache.npz`.
"""
from __future__ import annotations

import gzip
import json
import time
import urllib.request
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
WORK = Path(__file__).resolve().parent
PEAKS = WORK / "peaks"
PEAKS.mkdir(exist_ok=True)
API = "https://www.encodeproject.org"
WIN = 1000

# ------------------------------------------------------------------ 1. download the peak calls
sel = json.loads((WORK / "line_selection.json").read_text())
jobs = []
for row in sel:
    for k, exp in enumerate(row["experiments"][:2]):
        f = exp["peak"]
        if not f.get("file_size"):  # HuH-7's only peak file is 0 bytes
            continue
        jobs.append({"line": row["ccle_sample"], "encode_term": row["encode_term_name"],
                     "rep": k, "experiment": exp["experiment"], "file": f["accession"],
                     "href": f["href"], "bytes": f["file_size"],
                     "output_type": f["output_type"]})

total = 0
for j in jobs:
    dest = PEAKS / f"{j['file']}.bed.gz"
    if not dest.exists() or dest.stat().st_size == 0:
        urllib.request.urlretrieve(f"{API}{j['href']}", dest)
        time.sleep(0.05)  # ENCODE asks for <= 10 GET/s
    total += dest.stat().st_size
print(f"[download] {len(jobs)} peak files, {total/1e6:.1f} MB in {PEAKS}")

# ------------------------------------------------------------------------------ 2. the TSS set
bm = pd.read_csv(ROOT / "experiments/301-metadata-acquisitions/work/biomart_genes.tsv", sep="\t")
bm = bm[bm["Gene type"] == "protein_coding"]
bm = bm[bm["Chromosome/scaffold name"].astype(str).isin(
    [str(i) for i in range(1, 23)] + ["X", "Y"])]
bm = bm.dropna(subset=["Gene name"])
bm = bm[~bm["Gene name"].duplicated(keep=False)]  # unique symbol only; no arbitration
tss = np.where(bm["Strand"].to_numpy() > 0,
               bm["Gene start (bp)"].to_numpy(), bm["Gene end (bp)"].to_numpy())
genes_all = bm["Gene name"].astype(str).to_numpy()
chrom_all = ("chr" + bm["Chromosome/scaffold name"].astype(str)).to_numpy()
print(f"[tss] {len(genes_all)} protein-coding genes with a unique symbol")

# --------------------------------------------------------------- 3. CCLE basal expression axis
z = np.load(ROOT / "data/external/ccle_broad_2025/ccle_cache.npz", allow_pickle=False)
X, ccle_genes, ccle_samples = z["X"], [str(g) for g in z["genes"]], [str(s) for s in z["samples"]]
gi = {g: i for i, g in enumerate(ccle_genes)}   # last occurrence wins, as in def_identity
si = {s: i for i, s in enumerate(ccle_samples)}
print(f"[ccle] X {X.shape}; value range {X.min():.3f} .. {X.max():.1f} (raw TPM if max >> 20)")

lines = sorted({j["line"] for j in jobs if j["line"] in si})
keep = np.array([g in gi for g in genes_all])
genes, chrom, tss = genes_all[keep], chrom_all[keep], tss[keep]
print(f"[axis] {len(lines)} lines x {len(genes)} genes on the shared axis")

E = np.log2(1.0 + np.asarray([X[si[l]][[gi[g] for g in genes]] for l in lines], dtype=np.float64))

# ----------------------------------------------------- 4. promoter accessibility from the peaks
by_chrom: dict[str, list[int]] = defaultdict(list)
for idx, c in enumerate(chrom):
    by_chrom[c].append(idx)


def promoter_signal(path: Path) -> np.ndarray:
    """max narrowPeak signalValue over TSS +/- WIN; 0 where no peak intersects."""
    out = np.zeros(len(genes), dtype=np.float64)
    iv: dict[str, list[tuple[int, int, float]]] = defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 7:
                continue
            iv[p[0]].append((int(p[1]), int(p[2]), float(p[6])))
    for c, idxs in by_chrom.items():
        seg = iv.get(c)
        if not seg:
            continue
        seg.sort()
        starts = [s for s, _, _ in seg]
        for i in idxs:
            lo, hi = tss[i] - WIN, tss[i] + WIN
            k = bisect_left(starts, lo)
            # peaks are short; step back a bounded window to catch one that starts before `lo`
            best = 0.0
            for j in range(max(0, k - 40), len(seg)):
                s, e, v = seg[j]
                if s >= hi:
                    break
                if e > lo and v > best:
                    best = v
            out[i] = best
    return out


A = {}
for j in jobs:
    if j["line"] not in si:
        continue
    A[(j["line"], j["rep"])] = promoter_signal(PEAKS / f"{j['file']}.bed.gz")
    n_open = int((A[(j['line'], j['rep'])] > 0).sum())
    print(f"  {j['line'][:34]:34s} rep{j['rep']} {j['file']} open promoters {n_open:6d} "
          f"({100*n_open/len(genes):4.1f}%)")

A0 = np.log1p(np.asarray([A[(l, 0)] for l in lines]))
rep1_lines = [l for l in lines if (l, 1) in A]
A1 = np.log1p(np.asarray([A[(l, 1)] for l in rep1_lines])) if rep1_lines else None

np.savez_compressed(
    WORK / "matrix315.npz", A0=A0.astype(np.float32), E=E.astype(np.float32),
    A1=(A1.astype(np.float32) if A1 is not None else np.zeros((0, len(genes)), np.float32)),
    lines=np.asarray(lines, dtype="U"), rep1_lines=np.asarray(rep1_lines, dtype="U"),
    genes=genes.astype("U"),
)
meta = {"n_lines": len(lines), "n_genes": int(len(genes)), "window_bp": WIN,
        "n_lines_with_second_experiment": len(rep1_lines),
        "peak_files": [{k: j[k] for k in ("line", "encode_term", "rep", "experiment", "file",
                                          "bytes", "output_type")} for j in jobs],
        "download_bytes": total}
(WORK / "matrix315.meta.json").write_text(json.dumps(meta, indent=2) + "\n")
print(f"\nA0 {A0.shape}  E {E.shape}  second experiments for {len(rep1_lines)} lines")
print(f"wrote {WORK / 'matrix315.npz'}")
