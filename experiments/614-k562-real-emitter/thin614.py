#!/usr/bin/env python3
"""614 (descriptive control, the coordinator's request): h1-test binomially thinned to K562's median control library, so
every 614 reading can separate "cancer line vs H1" from "11k vs 20k depth". Each count c becomes Binomial(c, p) with
p = 11,327.5 / 20,190.5 (the median control libraries of the K562 truth and of h1-test), seed 614, applied to every cell,
controls and perturbed alike. Explicit zeros are dropped. Written in the minimal CSR-h5 layout that ref614 reads.
    uv run python experiments/614-k562-real-emitter/thin614.py   -> work/h1-thin11k.h5"""
import json
import sys
from pathlib import Path
import h5py
import numpy as np
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/323-current-family-validation"))
import evaluate323 as E  # noqa: E402  (labels)
K562_MED = json.loads((HERE / "work/ref_k562/reference.json").read_text())["median_control_library"]
src = ROOT / "data/harness/h1-test.h5ad"
with h5py.File(src, "r") as h:
    lab = E.labels(h, "obs/target_gene")
    genes = E.labels(h, "var/_index")
    shape = tuple(int(v) for v in h["X"].attrs["shape"])
    data = h["X/data"][:]
    indices = h["X/indices"][:]
    indptr = h["X/indptr"][:].astype(np.int64)
assert np.all(data == np.floor(data))
lib = np.add.reduceat(data.astype(np.float64), indptr[:-1]) if data.size else np.zeros(shape[0])
H1_MED = float(np.median(lib[lab == "non-targeting"]))
p = K562_MED / H1_MED
rng = np.random.default_rng(614)
new = rng.binomial(data.astype(np.int64), p)
del data
keep = new > 0
rowid = np.repeat(np.arange(shape[0]), np.diff(indptr))
cnt = np.bincount(rowid[keep], minlength=shape[0])
ip2 = np.concatenate([[0], np.cumsum(cnt)]).astype(np.int64)
d2 = new[keep].astype(np.float32)
i2 = indices[keep].astype(np.int32)
del new, indices, rowid
out = HERE / "work/h1-thin11k.h5"
with h5py.File(out, "w") as o:
    o.create_dataset("obs/target_gene", data=lab.astype("S"))
    o.create_dataset("var/_index", data=genes.astype("S"))
    x = o.create_group("X")
    x.attrs["shape"] = shape
    x.create_dataset("data", data=d2)
    x.create_dataset("indices", data=i2)
    x.create_dataset("indptr", data=ip2)
lib2 = np.add.reduceat(d2.astype(np.float64), ip2[:-1])
info = {"p": p, "h1_median_control_library": H1_MED, "target_median": K562_MED,
        "thinned_median_control_library": float(np.median(lib2[lab == "non-targeting"])), "nnz": int(d2.size)}
(HERE / "work/h1-thin11k.json").write_text(json.dumps(info, indent=1))
print(json.dumps(info))
