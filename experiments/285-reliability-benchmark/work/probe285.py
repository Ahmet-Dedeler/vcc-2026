"""285 -- the BOUNDED CD4 cell-level extraction PROBE (metadata / layout only, no expression value).

Probes ONE object over HTTP range requests with 280's `remote_h5.RemoteFile`:
    D2_Stim48hr.assigned_guide.h5ad   (137.6 GB; Stim48hr = the run-clean condition, R280 section 1.2)

Reads: the group/dataset tree (shapes, dtypes, chunks, compression), `obs` (perturbed_gene_name, guide_id,
lane_id, total_counts, low_quality), and `X`'s encoding attrs + `indptr` if CSR.  Then plans the rows of every
cell whose target is one of the 300 panel targets plus a non-targeting sample of <= 2x those cells, and
estimates the bytes from `indptr`.

    PYTHONDONTWRITEBYTECODE=1 uv run python probe285.py [OBJECT]

Writes work/probe285-<object>.json + the log on stdout.  Nothing is downloaded beyond the metadata.
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
LOG: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def walk(g, prefix="", depth=0, maxdepth=3, out=None):
    for k, v in g.items():
        p = f"{prefix}/{k}"
        if isinstance(v, h5py.Dataset):
            out.append(dict(path=p, kind="D", shape=list(v.shape), dtype=str(v.dtype), chunks=str(v.chunks),
                            compression=str(v.compression), nbytes=int(np.prod(v.shape) * v.dtype.itemsize)))
        else:
            out.append(dict(path=p, kind="G", enc=str(dict(v.attrs).get("encoding-type", "")),
                            attrs={kk: str(vv) for kk, vv in v.attrs.items()}))
            if depth < maxdepth:
                walk(v, p, depth + 1, maxdepth, out)


def cat_col(h, key):
    o = h[f"obs/{key}"]
    if isinstance(o, h5py.Group) and "categories" in o:
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in o["categories"][:]])
        return cats[o["codes"][:]], len(cats)
    v = o[:]
    if v.dtype.kind in "SO":
        return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v]), None
    return v, None


def main() -> None:
    t0 = time.time()
    with open(ROOT / "data/challenge/pert_counts.csv") as fh:
        PAN = {r["target_gene"] for r in csv.DictReader(fh)}
    rf = RemoteFile(f"{B}/{NAME}", block=1 << 22, max_blocks=128)
    P(f"285 probe285.py {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} -- {NAME}: {rf.size:,} bytes")
    f = h5py.File(rf, "r")
    tree: list[dict] = []
    walk(f, maxdepth=3, out=tree)
    for t in tree:
        if t["kind"] == "D":
            P(f"  D {t['path']:<58} {str(t['shape']):<16} {t['dtype']:<10} chunks={t['chunks']:<18} {t['compression']}")
        else:
            P(f"  G {t['path']:<58} {t['enc']}")
    P(f"[tree: {time.time() - t0:.0f}s, {rf.n_requests} req, {rf.n_bytes/1e6:.1f} MB]")

    X = f["X"]
    xattrs = {k: (v.tolist() if hasattr(v, "tolist") else str(v)) for k, v in X.attrs.items()} if isinstance(X, h5py.Group) else {}
    out: dict = dict(object=NAME, bytes=int(rf.size), tree=tree, x_attrs=xattrs)
    is_csr = isinstance(X, h5py.Group) and str(xattrs.get("encoding-type", "")).startswith("csr")
    uncompressed = bool(is_csr and X["data"].compression is None and X["indices"].compression is None)
    out["is_csr"] = is_csr
    out["uncompressed"] = uncompressed
    P(f"X: csr={is_csr} uncompressed={uncompressed} "
      f"(data {X['data'].dtype if is_csr else '-'} comp {X['data'].compression if is_csr else '-'}; "
      f"indices {X['indices'].dtype if is_csr else '-'} comp {X['indices'].compression if is_csr else '-'})")

    # ---------------------------------------------------------------- obs
    obs_keys = list(f["obs"].keys())
    out["obs_keys"] = obs_keys
    cols: dict[str, object] = {}
    tgt = guide = lane = lowq = None
    for k in ("perturbed_gene_name", "guide_id", "lane_id", "low_quality", "guide_type", "total_counts",
              "n_genes_by_counts", "guide_group"):
        if k not in obs_keys:
            continue
        v, ncat = cat_col(f, k)
        cols[k] = dict(n=int(len(v)), n_unique=int(len(set(v.tolist()))) if v.dtype.kind in "US" else int(np.unique(v).size),
                       n_categories=ncat, sample=[str(x) for x in v[:4]],
                       dtype=str(v.dtype))
        if v.dtype.kind in "iuf":
            cols[k].update(median=float(np.median(v)), q25=float(np.percentile(v, 25)), q75=float(np.percentile(v, 75)))
        elif v.dtype.kind == "b":
            cols[k].update(n_true=int(v.sum()), frac_true=float(v.mean()))
        if k == "perturbed_gene_name":
            tgt = v
        elif k == "guide_id":
            guide = v
        elif k == "lane_id":
            lane = v
        elif k == "low_quality":
            lowq = v
        P(f"  obs/{k}: n {len(v):,} unique {cols[k]['n_unique']} sample {cols[k]['sample']}")
    out["obs_cols"] = cols
    P(f"[obs: {time.time() - t0:.0f}s, {rf.n_requests} req, {rf.n_bytes/1e6:.1f} MB]")

    if tgt is None:
        out["verdict"] = "no perturbed_gene_name column"
        (WORK / f"probe285-{NAME.split('.')[0]}.json").write_text(json.dumps(out, indent=1, default=str))
        return

    ct = Counter(tgt.tolist())
    ntc_labels = [k for k in ct if "non" in k.lower() and "target" in k.lower()]
    pan_present = sorted(set(ct) - set(ntc_labels))
    pan_hit = sorted(set(pan_present) & PAN)
    n_pan_cells = int(sum(ct[t] for t in pan_hit))
    n_ntc = int(sum(ct[t] for t in ntc_labels))
    out["targets"] = dict(n_labels=len(ct), ntc_labels=ntc_labels, n_ntc_cells=n_ntc,
                          n_panel_targets_present=len(pan_hit), n_panel_cells=n_pan_cells,
                          cells_per_panel_target=dict(
                              median=float(np.median([ct[t] for t in pan_hit])) if pan_hit else None,
                              min=int(min(ct[t] for t in pan_hit)) if pan_hit else None,
                              max=int(max(ct[t] for t in pan_hit)) if pan_hit else None,
                              ge100=int(sum(ct[t] >= 100 for t in pan_hit)), ge200=int(sum(ct[t] >= 200 for t in pan_hit)),
                              ge400=int(sum(ct[t] >= 400 for t in pan_hit))))
    if lowq is not None:
        lq = np.asarray([str(x) for x in lowq])
        out["targets"]["low_quality_values"] = dict(Counter(lq.tolist()))
    P(f"targets: {len(ct)} labels, {len(pan_hit)} of the 300 panel present, {n_pan_cells:,} panel cells "
      f"(median/target {out['targets']['cells_per_panel_target']['median']}), {n_ntc:,} non-targeting cells")
    if guide is not None:
        gp = Counter()
        for t, g in zip(tgt, guide):
            if t in set(pan_hit):
                gp[(t, g)] += 1
        per = Counter()
        for (t, g) in gp:
            per[t] += 1
        out["targets"]["guides_per_panel_target_hist"] = dict(Counter(per.values()))
        P(f"  guides per panel target: {dict(Counter(per.values()))}")
    if lane is not None:
        out["targets"]["n_lanes"] = int(len(set(lane.tolist())))
        P(f"  lanes: {out['targets']['n_lanes']}")

    # ---------------------------------------------------------------- byte plan from indptr
    if is_csr:
        ip = np.asarray(X["indptr"][:], dtype=np.int64)
        P(f"[indptr read: {time.time() - t0:.0f}s, {rf.n_requests} req, {rf.n_bytes/1e6:.1f} MB]")
        nnz_per = np.diff(ip)
        itemsize = X["data"].dtype.itemsize + X["indices"].dtype.itemsize
        rng = np.random.default_rng(20260906)
        keep = np.isin(tgt, pan_hit)
        ntc_idx = np.flatnonzero(np.isin(tgt, ntc_labels))
        out["plan"] = {}
        for tag, n_ntc_take, cap in (("mandate_2x_ntc", 2 * int(keep.sum()), None),
                                     ("panel_only", 0, None),
                                     ("panel_cap900_plus_15k_ntc", 15000, 900)):
            if cap is None:
                pk = np.flatnonzero(keep)
            else:
                pk = []
                for t in pan_hit:
                    ix = np.flatnonzero(tgt == t)
                    pk.append(ix if ix.size <= cap else rng.choice(ix, size=cap, replace=False))
                pk = np.sort(np.concatenate(pk)) if pk else np.asarray([], dtype=np.int64)
            n_take = min(len(ntc_idx), n_ntc_take)
            ntc_take = rng.choice(ntc_idx, size=n_take, replace=False) if n_take else np.asarray([], dtype=np.int64)
            rows = np.sort(np.concatenate([pk, ntc_take]).astype(np.int64))
            payload = int(nnz_per[rows].sum()) * itemsize
            runs = np.split(rows, np.flatnonzero(np.diff(rows) != 1) + 1)
            span = int(sum(ip[r[-1] + 1] - ip[r[0]] for r in runs)) * itemsize
            out["plan"][tag] = dict(n_rows=int(rows.size), n_panel_rows=int(len(pk)), n_ntc_rows=int(n_take),
                                    total_nnz=int(nnz_per[rows].sum()), itemsize=int(itemsize),
                                    payload_bytes=payload, contiguous_span_bytes=span, n_runs=len(runs),
                                    sparse_i32_bank_bytes=int(nnz_per[rows].sum()) * 8,
                                    dense_f32_bytes=int(rows.size) * int(X.attrs["shape"][1]) * 4,
                                    admissible=bool(uncompressed and span <= 25e9))
            P(f"PLAN[{tag}]: {rows.size:,} rows ({len(pk):,} panel + {n_take:,} NTC) in {len(runs):,} runs; "
              f"payload {payload/1e9:.2f} GB, span {span/1e9:.2f} GB, sparse-i32 bank "
              f"{out['plan'][tag]['sparse_i32_bank_bytes']/1e9:.2f} GB, dense f32 "
              f"{out['plan'][tag]['dense_f32_bytes']/1e9:.2f} GB; admissible {out['plan'][tag]['admissible']}")
        P(f"median nnz/cell {np.median(nnz_per):.0f}; uncompressed CSR {uncompressed}")
    else:
        out["plan"] = dict(admissible=False, reason="X is not CSR / not range-addressable")
        P("NOT ADMISSIBLE: X is not an uncompressed CSR")

    out["requests"] = int(rf.n_requests)
    out["bytes_fetched"] = int(rf.n_bytes)
    out["seconds"] = round(time.time() - t0, 1)
    (WORK / f"probe285-{NAME.split('.')[0]}.json").write_text(json.dumps(out, indent=1, default=str))
    (WORK / f"probe285-{NAME.split('.')[0]}.log").write_text("\n".join(LOG) + "\n")
    P(f"\ndone {out['seconds']}s, {rf.n_requests} range requests, {rf.n_bytes/1e6:.1f} MB fetched")


if __name__ == "__main__":
    main()
