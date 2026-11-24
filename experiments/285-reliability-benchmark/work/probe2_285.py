"""285 -- CD4 probe, part 2: the byte plan WITHOUT letting h5py walk the chunk index.

`obs/n_genes_by_counts` is the per-cell nonzero count, so `indptr = [0, cumsum(n_genes_by_counts)]` if X stores
no explicit zeros.  That is checked against the real `indptr`, read in ONE range request from the dataset's own
chunk byte offsets (the dataset is uncompressed and its chunks are laid out linearly -- verified here on 7
chunks, 280's `layout280.py` check).  Reading `indptr` through h5py's own chunk index took > 20 min and was
abandoned; this route takes seconds.

    PYTHONDONTWRITEBYTECODE=1 uv run python probe2_285.py <OBJECT>
"""
from __future__ import annotations

import csv
import json
import sys
import time
import urllib.request
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


def cat_col(h, key):
    o = h[f"obs/{key}"]
    if isinstance(o, h5py.Group) and "categories" in o:
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in o["categories"][:]])
        return cats[o["codes"][:]]
    v = o[:]
    if v.dtype.kind in "SO":
        return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])
    return v


def get_range(url: str, lo: int, hi: int) -> bytes:
    req = urllib.request.Request(url, headers={"Range": f"bytes={lo}-{hi}"})
    with urllib.request.urlopen(req, timeout=300) as r:
        return r.read()


def main() -> None:
    t0 = time.time()
    with open(ROOT / "data/challenge/pert_counts.csv") as fh:
        PAN = {r["target_gene"] for r in csv.DictReader(fh)}
    url = f"{B}/{NAME}"
    rf = RemoteFile(url, block=1 << 22, max_blocks=64)
    f = h5py.File(rf, "r")
    X = f["X"]
    shape = [int(x) for x in X.attrs["shape"]]
    P(f"285 probe2 {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} -- {NAME}: {rf.size:,} B, X {shape}, "
      f"data {X['data'].dtype}/{X['data'].compression}, indices {X['indices'].dtype}/{X['indices'].compression}")

    tgt = cat_col(f, "perturbed_gene_name")
    guide = cat_col(f, "guide_id")
    ggrp = cat_col(f, "guide_group")
    gtyp = cat_col(f, "guide_type")
    lane = cat_col(f, "lane_id")
    lowq = np.asarray(f["obs/low_quality"][:])
    ngc = np.asarray(f["obs/n_genes_by_counts"][:], dtype=np.int64)
    tot = np.asarray(f["obs/total_counts"][:], dtype=np.float64)
    P(f"[obs read: {time.time()-t0:.0f}s, {rf.n_requests} req, {rf.n_bytes/1e6:.0f} MB]")
    P(f"  guide_group: {dict(Counter(ggrp.tolist()))}")
    P(f"  guide_type : {dict(Counter(gtyp.tolist()))}")
    P(f"  low_quality: {dict(Counter(lowq.tolist()))}")

    # the real indptr, one range request from the chunk byte offsets
    dsid = X["indptr"].id
    nch = dsid.get_num_chunks()
    infos = [dsid.get_chunk_info(i) for i in (0, 1, 2, nch // 2, nch - 2, nch - 1)]
    base, cb = infos[0].byte_offset, X["indptr"].chunks[0] * X["indptr"].dtype.itemsize
    linear = all(inf.byte_offset == base + i * cb for i, inf in zip((0, 1, 2, nch // 2, nch - 2, nch - 1), infos))
    P(f"  indptr: {X['indptr'].shape[0]:,} int64 in {nch:,} chunks of {X['indptr'].chunks[0]}, base {base:,}, "
      f"linear layout {linear}")
    lo = base
    hi = base + nch * cb - 1
    raw = get_range(url, lo, hi)
    ip = np.frombuffer(raw, dtype="<i8")[:X["indptr"].shape[0]].copy()
    P(f"[indptr via one {len(raw)/1e6:.0f} MB range request: {time.time()-t0:.0f}s]")
    nnz_row = np.diff(ip)
    ka = dict(row="KA-probe indptr == cumsum(n_genes_by_counts)",
              max_abs_diff=float(np.abs(nnz_row - ngc).max()), monotone=bool((nnz_row >= 0).all()),
              total_nnz=int(ip[-1]), listed_nnz=int(X["data"].shape[0]))
    ka["ok"] = bool(ka["max_abs_diff"] == 0 and ka["total_nnz"] == ka["listed_nnz"])
    P(f"KA-probe: max|diff(indptr) - n_genes_by_counts| = {ka['max_abs_diff']:.0f}; total nnz {ka['total_nnz']:,} "
      f"vs listed {ka['listed_nnz']:,}  {'PASS' if ka['ok'] else 'FAIL'}")

    ct = Counter(tgt.tolist())
    pan_hit = sorted(set(ct) & PAN)
    itemsize = X["data"].dtype.itemsize + X["indices"].dtype.itemsize
    rng = np.random.default_rng(20260906)
    out: dict = dict(object=NAME, bytes=int(rf.size), shape=shape, itemsize=int(itemsize), ka=ka,
                     n_panel_present=len(pan_hit), n_panel_cells=int(sum(ct[t] for t in pan_hit)),
                     guide_group=dict(Counter(ggrp.tolist())), guide_type=dict(Counter(gtyp.tolist())),
                     low_quality_true=int(lowq.sum()), n_lanes=int(len(set(lane.tolist()))),
                     median_nnz_per_cell=float(np.median(ngc)), median_umi=float(np.median(tot)))

    # the QC-clean single-guide panel cells (the extraction cohort)
    single = (ggrp == "targeting single sgRNA") & (~lowq)
    onpan = np.isin(tgt, pan_hit)
    P(f"  cells: {len(tgt):,}; panel {int(onpan.sum()):,}; panel & single-guide & QC-pass {int((onpan & single).sum()):,}")
    gp: dict[str, Counter] = {}
    for t, g, s in zip(tgt[onpan], guide[onpan], single[onpan]):
        if s:
            gp.setdefault(t, Counter())[g] += 1
    nper = Counter(len(v) for v in gp.values())
    cells_per_target = {t: sum(v.values()) for t, v in gp.items()}
    out["single_guide"] = dict(n_targets=len(gp), guides_per_target_hist=dict(nper),
                               cells_per_target=dict(median=float(np.median(list(cells_per_target.values()))),
                                                     min=int(min(cells_per_target.values())),
                                                     max=int(max(cells_per_target.values())),
                                                     ge50=int(sum(v >= 50 for v in cells_per_target.values())),
                                                     ge100=int(sum(v >= 100 for v in cells_per_target.values())),
                                                     ge200=int(sum(v >= 200 for v in cells_per_target.values()))),
                               min_guide_cells=dict(median=float(np.median([min(v.values()) for v in gp.values()]))))
    P(f"  single-guide QC-pass panel targets {len(gp)}; guides/target {dict(nper)}; cells/target median "
      f"{out['single_guide']['cells_per_target']['median']:.0f} (>=50 {out['single_guide']['cells_per_target']['ge50']}, "
      f">=100 {out['single_guide']['cells_per_target']['ge100']}); median smaller-guide cells "
      f"{out['single_guide']['min_guide_cells']['median']:.0f}")

    out["plan"] = {}
    for tag, mask, cap in (("all_panel_cells", onpan, None),
                           ("panel_single_qc", onpan & single, None),
                           ("panel_single_qc_cap900", onpan & single, 900)):
        idx = np.flatnonzero(mask)
        if cap:
            keep = []
            for t in pan_hit:
                ix = idx[tgt[idx] == t]
                keep.append(ix if ix.size <= cap else np.sort(rng.choice(ix, cap, replace=False)))
            idx = np.sort(np.concatenate(keep)) if keep else idx
        payload = int(nnz_row[idx].sum()) * itemsize
        runs = np.split(idx, np.flatnonzero(np.diff(idx) != 1) + 1)
        span = int(sum(ip[r[-1] + 1] - ip[r[0]] for r in runs)) * itemsize
        out["plan"][tag] = dict(n_rows=int(idx.size), total_nnz=int(nnz_row[idx].sum()), payload_bytes=payload,
                                contiguous_span_bytes=span, n_runs=len(runs),
                                sparse_i32_bank_bytes=int(nnz_row[idx].sum()) * 8,
                                admissible=bool(span <= 25e9))
        P(f"  PLAN[{tag}]: {idx.size:,} rows, {len(runs):,} runs, payload {payload/1e9:.2f} GB, span "
          f"{span/1e9:.2f} GB, sparse i32 bank {payload*8/12/1e9:.2f} GB, admissible {out['plan'][tag]['admissible']}")

    out["requests"], out["bytes_fetched"], out["seconds"] = int(rf.n_requests), int(rf.n_bytes) + len(raw), round(time.time() - t0, 1)
    np.savez_compressed(WORK / f"probe2_{STEM}.npz", indptr=ip, tgt=tgt, guide=guide, ggrp=ggrp, gtyp=gtyp,
                        lane=lane, lowq=lowq, ngc=ngc, total=tot)
    (WORK / f"probe2_{STEM}.json").write_text(json.dumps(out, indent=1, default=str))
    (WORK / f"probe2_{STEM}.log").write_text("\n".join(LOG) + "\n")
    P(f"\ndone {out['seconds']}s, {rf.n_requests} h5py range requests + 1 bulk, {out['bytes_fetched']/1e6:.0f} MB fetched")


if __name__ == "__main__":
    main()
