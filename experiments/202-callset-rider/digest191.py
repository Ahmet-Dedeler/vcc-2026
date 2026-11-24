"""191 falsifiers 1 and 2 — one streaming pass over each of the two 19 GB artifacts.

FALSIFIER 1 (pseudobulk identity).  `pds` (BRIEF §1) and the profile half of `mse` are read
off each (context, perturbation) block's PER-GENE PANEL TOTAL.  `rc_hybrid` claims those are
`a15_local`'s count for count, so `Delta pds = 0` and `Delta mse_profile = 0` are IDENTITIES,
not estimates.  One sha256 of the delivered `int64` totals per block, on the artifacts.

Non-vacuity, which is the whole reason 156's first a23 check was worthless: the NONZERO
STRUCTURE digest must DIFFER on essentially every block.  If it matched, the two files would
simply be the same file and the totals identity would prove nothing.

FALSIFIER 2 (`nmae`).  BRIEF §3 reads `lfc` off the ARITHMETIC mean of per-cell CPM, so the
two arms' predicted fold changes differ exactly where their per-cell library sizes do -- the
quasi law leaves them scattered (sd ~519 counts), `rc_hybrid` welds them to `T = S // 400`.
The reference cancels out of the difference:

    lfc_rc - lfc_a15 = log2( (m_rc + 1e-9) / (m_a15 + 1e-9) ),   m = mean_i 1e6 c_ig / L_i

so this is computable with no ground truth, exactly, on all 900 blocks.  `nmae_p` is
`sum_{S_p} |lfc_hat - lfc| / sum_{S_p} |lfc|` over the REFERENCE's DE set, so
`|Delta nmae_p| <= max_g |Delta lfc_hat_g| / mean_{S_p} |lfc|`, and the denominator is
measured on h1-test by `score191.py`.

Also checked here, because the pass is already paying for the I/O: the per-cell library
totals (the weld), the maximum cell total against the 1,000,000 format cap, and the panel
total sums.

    uv run python digest191.py <rc_hybrid.h5ad> <a15_local.h5ad>
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

import h5py
import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = Path(__file__).resolve().parent
BUNDLE = ROOT / "data/challenge"
N_GENES = 18533
SCORER_MIN_CPM = 5.0


def read_obs_blocks(path: Path):
    """(context, target_gene) per row block, straight off the categorical encoding."""
    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        out = {}
        for name in ("context", "target_gene"):
            group = obs[name]
            cats = np.asarray(group["categories"][:]).astype(str)
            codes = np.asarray(group["codes"][:])
            out[name] = cats[codes]
        indptr = np.asarray(handle["X"]["indptr"][:])
        shape = tuple(handle["X"].attrs["shape"])
    labels = list(zip(out["context"], out["target_gene"]))
    blocks, start = [], 0
    for i in range(1, len(labels) + 1):
        if i == len(labels) or labels[i] != labels[start]:
            blocks.append((labels[start], start, i))
            start = i
    return blocks, indptr, shape


def scan(path: Path, blocks, indptr, tag: str):
    """Per block: totals digest, structure digest, total sum, mean per-cell CPM, cell stats."""
    totals_sha, struct_sha, totals_sum, cpm_mean = {}, {}, {}, {}
    cell_stats = {}
    t0 = time.time()
    with h5py.File(path, "r") as handle:
        data = handle["X"]["data"]
        indices = handle["X"]["indices"]
        for k, (key, lo_row, hi_row) in enumerate(blocks):
            lo, hi = int(indptr[lo_row]), int(indptr[hi_row])
            cols = indices[lo:hi].astype(np.int64)
            vals = data[lo:hi].astype(np.float64)
            row_ptr = indptr[lo_row:hi_row + 1] - lo
            libraries = np.add.reduceat(vals, row_ptr[:-1])
            libraries[row_ptr[1:] == row_ptr[:-1]] = 0.0
            totals = np.bincount(cols, weights=vals, minlength=N_GENES).astype(np.int64)
            totals_sha[key] = hashlib.sha256(totals.tobytes()).hexdigest()[:16]
            totals_sum[key] = int(totals.sum())
            struct_sha[key] = hashlib.sha256(
                np.bincount(cols, minlength=N_GENES).astype(np.int64).tobytes()
            ).hexdigest()[:16]
            # mean over cells of 1e6 * count / library -- the object BRIEF §3's lfc reads
            safe = np.where(libraries > 0, libraries, 1.0)
            weights = vals * (1e6 / np.repeat(safe, np.diff(row_ptr)))
            n_cells = hi_row - lo_row
            cpm_mean[key] = (
                np.bincount(cols, weights=weights, minlength=N_GENES) / n_cells
            )
            cell_stats[key] = (
                float(libraries.min()), float(libraries.max()), float(libraries.std())
            )
            if (k + 1) % 150 == 0:
                print(f"  [{tag}] {k + 1}/{len(blocks)} [{time.time() - t0:.0f}s]", flush=True)
    return totals_sha, struct_sha, totals_sum, cpm_mean, cell_stats


def de_universe(context: str) -> np.ndarray:
    """The scorer's own DE gate for one context: mean per-cell CPM > 5 on its controls."""
    import anndata as ad
    import scipy.sparse as sp

    controls = ad.read_h5ad(BUNDLE / f"context_{context}.h5ad")
    matrix = controls.X.tocsr() if not sp.isspmatrix_csr(controls.X) else controls.X
    n_cells = matrix.shape[0]
    cpm = np.zeros(matrix.shape[1], dtype=np.float64)
    for lo in range(0, n_cells, 2000):
        block = np.asarray(matrix[lo:min(lo + 2000, n_cells)].todense(), dtype=np.float64)
        sums = block.sum(1)
        sums[sums == 0] = 1.0
        cpm += (block * (1e6 / sums[:, None])).sum(0)
    cpm /= n_cells
    del controls, matrix
    return cpm > SCORER_MIN_CPM


def main() -> None:
    built, incumbent = Path(sys.argv[1]), Path(sys.argv[2])
    b_blocks, b_indptr, b_shape = read_obs_blocks(built)
    a_blocks, a_indptr, a_shape = read_obs_blocks(incumbent)
    print(f"rc_hybrid  {built}  shape {b_shape}  {len(b_blocks)} blocks")
    print(f"a15_local  {incumbent}  shape {a_shape}  {len(a_blocks)} blocks")
    b_keys = [k for k, _, _ in b_blocks]
    a_keys = [k for k, _, _ in a_blocks]
    if b_keys != a_keys:
        raise SystemExit("block (context, target) ORDER differs -- comparison meaningless")
    sizes = {hi - lo for _, lo, hi in b_blocks} | {hi - lo for _, lo, hi in a_blocks}
    print(f"block sizes {sizes} (must be {{400}})")

    b_tot, b_str, b_sum, b_cpm, b_cell = scan(built, b_blocks, b_indptr, "rc_hybrid")
    a_tot, a_str, a_sum, a_cpm, a_cell = scan(incumbent, a_blocks, a_indptr, "a15_local")

    same_tot = sum(int(b_tot[k] == a_tot[k]) for k in b_keys)
    same_str = sum(int(b_str[k] == a_str[k]) for k in b_keys)
    same_sum = sum(int(b_sum[k] == a_sum[k]) for k in b_keys)
    print()
    print(f"per-gene PANEL TOTAL digest identical on   {same_tot}/{len(b_keys)}")
    print(f"panel total SUM identical on               {same_sum}/{len(b_keys)}")
    print(f"nonzero-STRUCTURE digest identical on      {same_str}/{len(b_keys)}"
          "   (must be far below 900, or the two files are the same file)")
    bad = [k for k in b_keys if b_tot[k] != a_tot[k]]
    if bad:
        print(f"FIRST MISMATCHES: {bad[:5]}")
    print("\nsample rows:")
    for k in b_keys[:3] + b_keys[-2:]:
        print(f"  {k}  totals {b_tot[k]} vs {a_tot[k]}   struct {b_str[k]} vs {a_str[k]}")
    print(f"\nFALSIFIER 1 {'PASS' if same_tot == len(b_keys) else 'FAIL'}: "
          f"{same_tot}/{len(b_keys)} delivered pseudobulk digests identical")

    print("\nPER-CELL LIBRARY TOTALS (the weld) and the 1,000,000 format cap")
    print(f"  {'arm':<12}{'min':>10}{'max':>10}{'sd median':>12}{'sd max':>10}")
    for tag, stats in (("a15_local", a_cell), ("rc_hybrid", b_cell)):
        lows = np.array([v[0] for v in stats.values()])
        highs = np.array([v[1] for v in stats.values()])
        sds = np.array([v[2] for v in stats.values()])
        print(f"  {tag:<12}{lows.min():>10.0f}{highs.max():>10.0f}"
              f"{np.median(sds):>12.3f}{sds.max():>10.3f}")
    cap_ok = max(v[1] for v in b_cell.values()) <= 1_000_000

    print("\nFALSIFIER 2 — the predicted log2 fold change, all 900 blocks, reference-free")
    gates = {ctx: de_universe(ctx) for ctx in ("A", "B", "C")}
    for ctx in ("A", "B", "C"):
        print(f"  context {ctx}: {int(gates[ctx].sum())}/{N_GENES} genes past the "
              f"{SCORER_MIN_CPM}-CPM scorer gate")
    per_block = []
    for key in b_keys:
        gate = gates[key[0]]
        delta = np.log2((b_cpm[key][gate] + 1e-9) / (a_cpm[key][gate] + 1e-9))
        per_block.append((float(np.abs(delta).max()), float(np.abs(delta).mean()),
                          float(np.abs(delta).sum()), float(delta.mean())))
    per_block = np.asarray(per_block)
    print(f"  |d lfc| over the DE universe: max {per_block[:, 0].max():.6f}   "
          f"mean {per_block[:, 1].mean():.6f}   median-of-block-means "
          f"{np.median(per_block[:, 1]):.6f}")
    print(f"  signed mean d lfc per block: median {np.median(per_block[:, 3]):+.6f}   "
          f"range [{per_block[:, 3].min():+.6f}, {per_block[:, 3].max():+.6f}]")

    report = dict(
        blocks=len(b_keys), totals_identical=same_tot, sums_identical=same_sum,
        structure_identical=same_str, cap_ok=bool(cap_ok),
        dlfc_max=float(per_block[:, 0].max()),
        dlfc_mean=float(per_block[:, 1].mean()),
        dlfc_block_mean_median=float(np.median(per_block[:, 1])),
        dlfc_signed_block_mean_median=float(np.median(per_block[:, 3])),
        cell_totals={
            "a15_local": [float(min(v[0] for v in a_cell.values())),
                          float(max(v[1] for v in a_cell.values())),
                          float(np.median([v[2] for v in a_cell.values()]))],
            "rc_hybrid": [float(min(v[0] for v in b_cell.values())),
                          float(max(v[1] for v in b_cell.values())),
                          float(np.median([v[2] for v in b_cell.values()]))],
        },
    )
    with open(HERE / "work" / "digest191.json", "w") as handle:
        json.dump(report, handle, indent=2)
    print(f"\nwrote {HERE / 'work' / 'digest191.json'}")


if __name__ == "__main__":
    main()
