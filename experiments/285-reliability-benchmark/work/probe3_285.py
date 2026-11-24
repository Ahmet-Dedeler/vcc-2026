"""285 -- CD4 probe, part 3: build `indptr` from `obs/n_genes_by_counts` and VALIDATE it, then the byte plan.

Reading the whole `indptr` through h5py's chunk index took > 20 min (2,861 chunks, non-linear layout), and a
single byte-range over the chunk region returns garbage because the chunks are NOT laid out linearly.  But X
stores no explicit zeros, so `indptr = [0, cumsum(n_genes_by_counts)]`.  That construction is validated two
ways before it is used: (i) its total equals the listed `X/data` length exactly; (ii) it equals the real
`indptr` at 40 scattered indices read through h5py (a handful of chunks -- seconds).

    PYTHONDONTWRITEBYTECODE=1 uv run python probe3_285.py <OBJECT>
"""
from __future__ import annotations

import csv
import json
import sys
import time
from collections import Counter
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/experiments/280-gwcd4i-donor-replicates/work")
from remote_h5 import RemoteFile  # noqa: E402

ROOT = Path("~/Code/virtual-cell-challenge-2026")
WORK = ROOT / "experiments/285-reliability-benchmark/work"
B = "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data"
NAME = sys.argv[1] if len(sys.argv) > 1 else "D2_Stim48hr.assigned_guide.h5ad"
STEM = NAME.split(".")[0]
LOG: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def main() -> None:
    t0 = time.time()
    with open(ROOT / "data/challenge/pert_counts.csv") as fh:
        PAN = {r["target_gene"] for r in csv.DictReader(fh)}
    z = np.load(WORK / f"probe2_{STEM}.npz", allow_pickle=True)
    tgt, guide, ggrp, gtyp = (z[k].astype(str) for k in ("tgt", "guide", "ggrp", "gtyp"))
    lane = z["lane"].astype(str)
    lowq = z["lowq"].astype(bool)
    ngc = np.asarray(z["ngc"], dtype=np.int64)
    tot = np.asarray(z["total"], dtype=np.float64)
    ip = np.zeros(ngc.size + 1, dtype=np.int64)
    np.cumsum(ngc, out=ip[1:])
    P(f"285 probe3 {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} -- {NAME}; constructed indptr from "
      f"n_genes_by_counts: total {ip[-1]:,}")

    url = f"{B}/{NAME}"
    rf = RemoteFile(url, block=1 << 22, max_blocks=64)
    f = h5py.File(rf, "r")
    X = f["X"]
    listed = int(X["data"].shape[0])
    itemsize = X["data"].dtype.itemsize + X["indices"].dtype.itemsize
    rng = np.random.default_rng(20260906)
    probe_idx = np.sort(rng.choice(ngc.size, size=40, replace=False))
    real = np.asarray([int(X["indptr"][int(i)]) for i in probe_idx])
    d = int(np.abs(real - ip[probe_idx]).max())
    ka = dict(row="KA-probe indptr == [0, cumsum(n_genes_by_counts)]",
              total_constructed=int(ip[-1]), total_listed=listed, total_ok=bool(int(ip[-1]) == listed),
              max_abs_diff_at_40_scattered_indices=d, ok=bool(int(ip[-1]) == listed and d == 0))
    P(f"KA-probe: constructed total nnz {ip[-1]:,} vs listed {listed:,}; max|d| at 40 scattered rows {d} "
      f"-> {'PASS' if ka['ok'] else 'FAIL'}  [{time.time()-t0:.0f}s, {rf.n_requests} req, {rf.n_bytes/1e6:.0f} MB]")

    ct = Counter(tgt.tolist())
    pan_hit = sorted(set(ct) & PAN)
    single = (ggrp == "targeting single sgRNA") & (~lowq)
    onpan = np.isin(tgt, pan_hit)
    out: dict = dict(object=NAME, bytes=int(rf.size), shape=[int(x) for x in X.attrs["shape"]], ka=ka,
                     itemsize=int(itemsize), n_panel_present=len(pan_hit),
                     median_nnz_per_cell=float(np.median(ngc)), median_umi=float(np.median(tot)),
                     n_lanes=int(len(set(lane.tolist()))),
                     guide_group=dict(Counter(ggrp.tolist())), guide_type=dict(Counter(gtyp.tolist())))
    for tag, mask in (("panel_all", onpan), ("panel_single_qc", onpan & single),
                      ("ntc_single_qc", (gtyp == "non-targeting") & (~lowq) & (ggrp != "multi sgRNA"))):
        idx = np.flatnonzero(mask)
        nnz = int(ngc[idx].sum())
        runs = np.split(idx, np.flatnonzero(np.diff(idx) != 1) + 1) if idx.size else []
        span = int(sum(int(ip[r[-1] + 1] - ip[r[0]]) for r in runs))
        out[f"plan_{tag}"] = dict(n_rows=int(idx.size), nnz=nnz, payload_gb=nnz * itemsize / 1e9,
                                  span_gb=span * itemsize / 1e9, n_runs=len(runs),
                                  bank_i32_gb=nnz * 8 / 1e9)
        P(f"  PLAN[{tag}]: {idx.size:,} rows, {len(runs):,} runs, payload {nnz*itemsize/1e9:.2f} GB, "
          f"span {span*itemsize/1e9:.2f} GB, kept sparse-i32 bank {nnz*8/1e9:.2f} GB")
    ntc_idx = np.flatnonzero((gtyp == "non-targeting") & (~lowq) & (ggrp != "multi sgRNA"))
    keep_ntc = np.sort(rng.choice(ntc_idx, size=min(15000, ntc_idx.size), replace=False)) if ntc_idx.size else np.asarray([], dtype=np.int64)
    plan_idx = np.sort(np.concatenate([np.flatnonzero(onpan & single), keep_ntc]).astype(np.int64))
    nnz = int(ngc[plan_idx].sum())
    out["plan_FINAL"] = dict(n_rows=int(plan_idx.size), n_panel=int((onpan & single).sum()), n_ntc=int(keep_ntc.size),
                             nnz=nnz, payload_gb=nnz * itemsize / 1e9, bank_i32_gb=nnz * 8 / 1e9,
                             admissible=bool(nnz * itemsize <= 25e9))
    P(f"  PLAN[FINAL = panel single-guide QC-pass + {keep_ntc.size:,} NTC]: {plan_idx.size:,} rows, "
      f"payload {nnz*itemsize/1e9:.2f} GB, kept bank {nnz*8/1e9:.2f} GB, ADMISSIBLE "
      f"{out['plan_FINAL']['admissible']}")
    np.save(WORK / f"plan_{STEM}.npy", plan_idx)
    np.save(WORK / f"indptr_{STEM}.npy", ip)
    (WORK / f"probe3_{STEM}.json").write_text(json.dumps(out, indent=1, default=str))
    (WORK / f"probe3_{STEM}.log").write_text("\n".join(LOG) + "\n")
    P(f"done {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
