#!/usr/bin/env python3
"""614 stage 3 (generic): a harness reference in evaluate323's exact format, built with evaluate323's own functions
(labels, block_stats, de_frame; rc_emit.control_block_from_counts, exact_de323.control_block, cell-eval2's
jackknife_correction). Only H1's frozen-axis / 400-cell asserts are dropped: the counts per target are recorded instead.
Any CSR-h5 harness truth (obs/target_gene, var/_index, X/{data,indices,indptr} + shape) plus a target list.
    uv run python experiments/614-k562-real-emitter/ref614.py <truth.h5> <targets.txt> <outdir>
Writes <outdir>/reference.npz, reference_de.parquet, reference.json exactly as evaluate323.reference() does, so that
evaluate323.load_reference() and score() work with run323.WORK = <outdir>."""
from __future__ import annotations

import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/323-current-family-validation"))
import h5py  # noqa: E402
import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import evaluate323 as E  # noqa: E402  (read-only)
R = E.R


def main() -> None:
    real, tfile, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    with h5py.File(real, "r") as h:
        lab = E.labels(h, "obs/target_gene")
        genes = E.labels(h, "var/_index")
        x = sp.csr_matrix((h["X/data"][:], h["X/indices"][:], h["X/indptr"][:]), shape=tuple(h["X"].attrs["shape"]))
    want = [ln.strip() for ln in open(tfile) if ln.strip()]
    targets = np.array(sorted(want))
    missing = set(targets) - set(lab)
    if missing:
        raise ValueError(f"targets absent from the truth: {sorted(missing)[:5]}")
    xc = x[lab == E.CTRL]
    emission_control = E.rc_emit.control_block_from_counts(xc)
    c = E.exact_de323.control_block(xc)
    if not np.array_equal(c.keep, emission_control.keep):
        raise AssertionError("emitter and scorer CPM gate differ")
    cp = np.asarray(xc.astype(np.float64).sum(0)).ravel().astype(np.int64)
    cj = E.jackknife_correction(xc, np.zeros(xc.shape[0], dtype=np.int64), 1, E.TS, chunk=32)[0]
    sums, jks, ps, qs, lfcs, counts = [], [], [], [], [], []
    for i, target in enumerate(targets):
        rows = x[lab == target]
        s, j, p, q, l = E.block_stats(rows, c)
        sums.append(s); jks.append(j); ps.append(p); qs.append(q); lfcs.append(l)
        counts.append(rows.shape[0])
        if (i + 1) % 10 == 0:
            print(f"reference {i + 1}/{len(targets)} ({time.monotonic() - t0:.1f}s)", flush=True)
    np.savez(out / "reference.npz", genes=genes, targets=targets,
             perts=np.append(targets, E.CTRL), sums=np.vstack([*sums, cp]),
             counts=np.array([*counts, xc.shape[0]]), jk=np.array([*jks, cj]),
             p=ps, q=qs, lfc=lfcs, keep=c.keep, cs=c.cs, tie=c.tie,
             ref_lin=c.ref_lin, emit_cs=emission_control.cs,
             emit_tie=emission_control.tie, emit_ref_lin=emission_control.ref_lin,
             libraries=np.asarray(xc.sum(1)).ravel().astype(float),
             b_ctrl=np.log1p(E.TS * cp / cp.sum()))
    frame = E.de_frame(targets, genes[c.keep], ps, qs, lfcs)
    frame.write_parquet(out / "reference_de.parquet")
    own = np.equal.outer(targets, genes[c.keep])
    n_real = ((np.array(qs) < .05) & ~own).sum(1)
    R.write_json(out / "reference.json", {
        "reference_sha256": R.sha(real),
        "statistics_sha256": R.sha(out / "reference.npz"),
        "de_sha256": R.sha(out / "reference_de.parquet"),
        "n_targets": len(targets), "n_controls": xc.shape[0],
        "cells_per_target": {"min": int(min(counts)), "median": float(np.median(counts)), "max": int(max(counts))},
        "n_keep_genes": int(np.asarray(c.keep).sum()), "median_control_library": float(np.median(np.asarray(xc.sum(1)).ravel())),
        "n_real_quantiles": np.percentile(n_real, [0, 10, 25, 50, 75, 90, 100]).tolist(),
        "n_real_by_target": dict(zip(map(str, targets), map(int, n_real))),
        "seconds": time.monotonic() - t0,
    })
    print(f"wrote {out} ({len(targets)} targets, keep {int(np.asarray(c.keep).sum())}, {time.monotonic() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
