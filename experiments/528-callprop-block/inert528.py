#!/usr/bin/env python3
"""528 — the INERT block-key artifact for KA-2: key_<c> = the context's control CPM computed with exactly the arithmetic of
`rc_emit.control_block_from_counts`' first pass (dense 2,000-row chunks, float64, count * (1e6 / row sum), summed over rows,
divided by n). On the 5-CPM universe that IS the emitter's `rc_control.ref_lin`, so `VCC_RC_CALLSET_BLOCK=prop:<this>:<sha>`
must rebuild the unset build byte for byte: it runs the whole keyed path (parse, sha, per-context load, keep-mapping, the
emitter's block_from_key) with a key whose order is the CPM order.

    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/528-callprop-block/inert528.py [--bundle data/challenge] [--contexts A,B,C]
Self-check: on h1-test's controls the same arithmetic reproduces 323's reference `emit_ref_lin` bit for bit.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import h5py  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def labels(f):
    o = f["obs/target_gene"]
    if isinstance(o, h5py.Group):
        return o["categories"][:].astype(str)[o["codes"][:]] if "categories" in o else o["values"][:].astype(str)
    return o[:].astype(str)


def rc_ref_lin(path, ctrl_only=True, chunk=2000):
    """control_block_from_counts' first pass on the file's control rows (file order), without materialising the matrix."""
    f = h5py.File(path, "r")
    ipd = f["X/indptr"][:].astype(np.int64)
    G = int(f["X"].attrs["shape"][1])
    rows = np.flatnonzero(labels(f) == "non-targeting") if ctrl_only else np.arange(ipd.size - 1)
    n = rows.size
    cpm = np.zeros(G, dtype=np.float64)
    for lo in range(0, n, chunk):
        r = rows[lo:lo + chunk]
        block = np.zeros((r.size, G), dtype=np.float64)
        for j, i in enumerate(r):
            a, b = int(ipd[i]), int(ipd[i + 1])
            block[j, f["X/indices"][a:b]] = f["X/data"][a:b]
        sums = block.sum(1)
        sums[sums == 0] = 1.0
        cpm += (block * (1e6 / sums[:, None])).sum(0)
    f.close()
    return cpm / n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", default=str(ROOT / "data/challenge"))
    ap.add_argument("--contexts", default="A,B,C")
    ap.add_argument("--out", default=str(HERE / "work/callblock528_inert_abc.npz"))
    ap.add_argument("--selfcheck", action="store_true")
    a = ap.parse_args()
    if a.selfcheck:
        z = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
        rl = rc_ref_lin(ROOT / "data/harness/h1-test.h5ad")
        keep = rl > 5.0
        ok = bool(np.array_equal(keep, np.asarray(z["keep"], bool)) and np.array_equal(rl[keep], np.asarray(z["emit_ref_lin"], float)))
        print(f"self-check on h1-test: keep and emit_ref_lin bit-equal: {ok}")
        sys.exit(0 if ok else 1)
    bundle = Path(a.bundle)
    official = [r["gene_name"] for r in csv.DictReader(open(bundle / "gene_names.csv"))]
    out = {"genes": np.array(official)}
    for c in [s.strip() for s in a.contexts.split(",") if s.strip()]:
        out[f"key_{c}"] = rc_ref_lin(bundle / f"context_{c}.h5ad")
        print(c, "universe", int((out[f"key_{c}"] > 5).sum()), flush=True)
    out["meta"] = np.array(json.dumps({"inert": True, "key": "control CPM (rc_emit.control_block_from_counts arithmetic)"}))
    np.savez(a.out, **out)
    sha = hashlib.sha256(Path(a.out).read_bytes()).hexdigest()
    Path(a.out + ".json").write_text(json.dumps({"sha256": sha, "path": str(a.out)}, indent=1))
    print(f"wrote {a.out} sha256 {sha}")


if __name__ == "__main__":
    main()
