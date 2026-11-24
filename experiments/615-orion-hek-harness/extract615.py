"""615 extractor: one streaming pass over Orion RAW HEK293T (223 batch parquets) -> work/hek-raw.h5.

Keeps, on the official 18,533-gene axis (Orion tokens mapped by symbol, as build_orion_harness.py does):
  * every pass-guide cell of the 300 A/B/C panel targets (kind 1), with its batch and official-axis UMI;
  * a seeded uniform sample of up to 300 pass-guide controls per batch (kind 0; seed 615 + batch), from which
    the harness's single control pool is later drawn with the pool-matched batch mixture (PREREGISTRATION-blob.md P');
  * per-batch control count sums over ALL controls (ctrl_sums, 223 x 18,533 float64) and control counts, which
    give the per-batch control means for the P / P' / BM / G profiles of the blob test.
Controls are gene_target == "Non-Targeting" exactly. Counts are checked integral and stored int32.
Memory: streams 4,096-cell record batches; peak ~2-3 GB. Nothing outside work/ is written.

    uv run python extract615.py [--line HEK293T] [--out work/hek-raw.h5]
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
CTRL = "Non-Targeting"
PER_BATCH_CTRL = 300


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--line", default="HEK293T")
    ap.add_argument("--out", default=str(HERE / "work" / "hek-raw.h5"))
    ap.add_argument("--max-batches", type=int, default=0, help="smoke test only")
    a = ap.parse_args()

    genes = pd.read_csv(REPO / "data/challenge/gene_names.csv")["gene_name"].astype(str).to_numpy()
    gidx = {g: i for i, g in enumerate(genes)}
    G = len(genes)
    panel = pd.read_csv(REPO / "data/challenge/pert_counts.csv")["target_gene"].astype(str).to_numpy()
    gm = pq.read_table(REPO / "data/external/orion/metadata/gene_metadata.parquet").to_pandas()
    n_tok = int(gm["gene_token_id"].max()) + 1
    tok2col = np.full(n_tok, -1, dtype=np.int64)
    for t, nm in zip(gm["gene_token_id"].to_numpy(), gm["gene_name"].astype(str).to_numpy()):
        c = gidx.get(nm)
        if c is not None:
            tok2col[int(t)] = c
    print(f"gene axis: {int((tok2col >= 0).sum()):,} Orion tokens -> official ({G:,})", flush=True)

    files = sorted((REPO / "data/external/orion/data").glob(f"{a.line}_Batch*.parquet"),
                   key=lambda p: int(re.search(r"Batch(\d+)", p.name)[1]))
    if a.max_batches:
        files = files[:a.max_batches]
    bnum = np.array([int(re.search(r"Batch(\d+)", p.name)[1]) for p in files])
    print(f"{len(files)} batches", flush=True)

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".h5.part")
    h = h5py.File(tmp, "w")
    dsd = h.create_dataset("X/data", (0,), maxshape=(None,), dtype="i4", chunks=(1 << 20,), compression="lzf")
    dsi = h.create_dataset("X/indices", (0,), maxshape=(None,), dtype="i4", chunks=(1 << 20,), compression="lzf")
    indptr = [0]
    lab_all, bat_all, umi_all, kind_all = [], [], [], []
    ctrl_sums = np.zeros((len(files), G), dtype=np.float64)
    n_ctrl = np.zeros(len(files), dtype=np.int64)
    n_cells = np.zeros(len(files), dtype=np.int64)
    ctrl_umi_med = np.zeros(len(files), dtype=np.float64)
    panel_set = set(panel)
    t0 = time.time()

    def append(rows):
        """rows: list of (cols int32, vals int32) -> append to the CSR datasets."""
        if not rows:
            return
        cols = np.concatenate([r[0] for r in rows])
        vals = np.concatenate([r[1] for r in rows])
        n0 = dsd.shape[0]
        dsd.resize((n0 + len(vals),))
        dsi.resize((n0 + len(cols),))
        dsd[n0:] = vals
        dsi[n0:] = cols
        for r in rows:
            indptr.append(indptr[-1] + len(r[0]))

    for bi, f in enumerate(files):
        pf = pq.ParquetFile(f)
        ctrl_rows, ctrl_umi, pan_rows, pan_lab, pan_umi = [], [], [], [], []
        for rb in pf.iter_batches(batch_size=4096, columns=["gene_token_id", "gene_expression",
                                                            "gene_target", "pass_guide_filter"]):
            tok, val = rb.column(0), rb.column(1)
            off = tok.offsets.to_numpy()
            flat = tok.values.to_numpy()
            v = val.values.to_numpy()
            assert np.array_equal(off, val.offsets.to_numpy())
            tgt = rb.column(2).to_numpy(zero_copy_only=False).astype(str)
            pas = rb.column(3).to_numpy()
            n_cells[bi] += int((pas == 1).sum())
            for j in range(rb.num_rows):
                if pas[j] != 1:
                    continue
                name = tgt[j]
                isc = name == CTRL
                if not isc and name not in panel_set:
                    continue
                s, e = off[j], off[j + 1]
                col = tok2col[flat[s:e]]
                m = col >= 0
                cv = v[s:e][m]
                ci = col[m]
                iv = np.rint(cv).astype(np.int64)
                assert np.array_equal(iv, cv), "non-integral counts"
                # Several Orion tokens can carry one official symbol; their counts are summed (what a CSR
                # pseudobulk with duplicate columns would do). np.unique sorts, bincount aligns with it.
                uc, inv = np.unique(ci, return_inverse=True)
                agg = np.bincount(inv, weights=iv, minlength=len(uc))
                ci, iv = uc.astype(np.int32), np.rint(agg).astype(np.int32)
                assert np.all(np.diff(ci) > 0), "duplicate gene after mapping"
                if isc:
                    ctrl_sums[bi, ci] += iv
                    n_ctrl[bi] += 1
                    ctrl_rows.append((ci, iv))
                    ctrl_umi.append(int(iv.sum()))
                else:
                    pan_rows.append((ci, iv))
                    pan_lab.append(name)
                    pan_umi.append(int(iv.sum()))
        ctrl_umi_med[bi] = float(np.median(ctrl_umi)) if ctrl_umi else np.nan
        rng = np.random.default_rng(615 + int(bnum[bi]))
        keep = np.sort(rng.choice(len(ctrl_rows), size=min(PER_BATCH_CTRL, len(ctrl_rows)), replace=False))
        rows = pan_rows + [ctrl_rows[k] for k in keep]
        append(rows)
        lab_all += pan_lab + [CTRL] * len(keep)
        bat_all += [int(bnum[bi])] * len(rows)
        umi_all += pan_umi + [ctrl_umi[k] for k in keep]
        kind_all += [1] * len(pan_rows) + [0] * len(keep)
        if bi % 10 == 0 or bi == len(files) - 1:
            print(f"  {bi + 1}/{len(files)} batch {bnum[bi]}: {n_cells[bi]:,} cells, {n_ctrl[bi]:,} ctrl, "
                  f"{len(pan_rows):,} panel | rows {len(lab_all):,} nnz {indptr[-1]:,} | {time.time() - t0:.0f}s",
                  flush=True)

    h.create_dataset("X/indptr", data=np.asarray(indptr, dtype=np.int64))
    h["X"].attrs["shape"] = (len(lab_all), G)
    h.create_dataset("obs/target_gene", data=np.asarray(lab_all, dtype="S"))
    h.create_dataset("obs/batch", data=np.asarray(bat_all, dtype=np.int16))
    h.create_dataset("obs/umi", data=np.asarray(umi_all, dtype=np.int32))
    h.create_dataset("obs/kind", data=np.asarray(kind_all, dtype=np.int8))
    h.create_dataset("var/gene_name", data=genes.astype("S"))
    h.create_dataset("batches/number", data=bnum)
    h.create_dataset("batches/ctrl_sums", data=ctrl_sums, compression="lzf")
    h.create_dataset("batches/n_ctrl", data=n_ctrl)
    h.create_dataset("batches/n_cells", data=n_cells)
    h.create_dataset("batches/ctrl_umi_median", data=ctrl_umi_med)
    h.close()
    os.replace(tmp, out)
    lab = np.asarray(lab_all)
    kind = np.asarray(kind_all)
    vc = pd.Series(lab[kind == 1]).value_counts()
    summ = dict(line=a.line, batches=len(files), rows=len(lab_all), nnz=int(indptr[-1]),
                controls_total=int(n_ctrl.sum()), controls_kept=int((kind == 0).sum()),
                panel_targets_present=int(len(vc)), panel_targets_ge100=int((vc >= 100).sum()),
                panel_cells=int((kind == 1).sum()), panel_cells_per_target_median=float(vc.median()),
                panel_umi_median=float(np.median(np.asarray(umi_all)[kind == 1])),
                ctrl_umi_median_all=float(np.nanmedian(ctrl_umi_med)), seconds=round(time.time() - t0))
    with open(out.with_suffix(".json"), "w") as fh:
        json.dump(summ, fh, indent=1)
    print(json.dumps(summ), flush=True)


if __name__ == "__main__":
    sys.exit(main())
