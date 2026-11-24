#!/usr/bin/env python3
"""453: map a D/E/F artifact built on 180's `def-bundle` (the A/B/C controls and panel, relabelled D->A... in reverse:
the bundle maps A->D, B->E, C->F) back onto A/B/C labels, so the board can score the October path.

It APFS-clones the source (never touches the build artifact, which `def_day.sh` hard-links), then rewrites only
obs['context']'s categories by NAME (D->A, E->B, F->C). X, obs index, target_gene and var are untouched, and it verifies
that. Also writes the clone's meta.json with a `relabel453` record.

usage: uv run python experiments/453-abc-dress-rehearsal/relabel453.py <src.h5ad> <dst.h5ad>
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

MAP = {"D": "A", "E": "B", "F": "C"}
src, dst = Path(sys.argv[1]), Path(sys.argv[2])
if dst.exists():
    raise SystemExit(f"{dst} exists; refusing to overwrite")
subprocess.run(["cp", "-c", str(src), str(dst)], check=True)
with h5py.File(dst, "r+") as h:
    node = h["obs"]["context"]
    cats = [x.decode() if isinstance(x, bytes) else str(x) for x in node["categories"][:]]
    if sorted(cats) != ["D", "E", "F"]:
        raise SystemExit(f"expected categories D/E/F, got {cats}")
    codes = np.asarray(node["codes"][:])
    counts_before = {c: int((codes == i).sum()) for i, c in enumerate(cats)}
    attrs = dict(node["categories"].attrs)
    new = [MAP[c] for c in cats]
    del node["categories"]
    node.create_dataset("categories", data=np.asarray(new, dtype=object), dtype=h5py.string_dtype())
    for k, v in attrs.items():
        node["categories"].attrs[k] = v
with h5py.File(dst, "r") as h, h5py.File(src, "r") as g:
    newcats = [x.decode() if isinstance(x, bytes) else str(x) for x in h["obs/context/categories"][:]]
    assert np.array_equal(h["obs/context/codes"][:], g["obs/context/codes"][:])
    for k in ("X/indptr",):
        assert np.array_equal(h[k][:], g[k][:]), k
    n = h["X/data"].shape[0]
    for s in (0, n // 2, max(n - 100000, 0)):
        assert np.array_equal(h["X/data"][s:s + 100000], g["X/data"][s:s + 100000])
        assert np.array_equal(h["X/indices"][s:s + 100000], g["X/indices"][s:s + 100000])
counts_after = {MAP[c]: v for c, v in counts_before.items()}
print(f"relabelled {cats} -> {newcats}; cells per context {counts_after}")
meta_src = Path(str(src) + ".meta.json")
if meta_src.exists():
    meta = json.loads(meta_src.read_text())
    meta["relabel453"] = {"source": str(src), "map": MAP, "categories_before": cats, "categories_after": newcats,
                          "cells": counts_after}
    # the novelty check reads output.sha256 as the prediction's identity: make it this file's, keep the source's
    if isinstance(meta.get("output"), dict) and "sha256" in meta["output"]:
        h = hashlib.sha256()
        with open(dst, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 24), b""):
                h.update(chunk)
        meta["relabel453"]["source_output_sha256"] = meta["output"]["sha256"]
        meta["output"]["sha256"] = h.hexdigest()
        meta["output"]["bytes"] = dst.stat().st_size if "bytes" in meta["output"] else meta["output"].get("bytes")
    Path(str(dst) + ".meta.json").write_text(json.dumps(meta, indent=1))
print("wrote", dst)
