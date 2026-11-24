#!/usr/bin/env python3
"""616: a CLEAN version of an Orion thinned harness truth, from 615's raw extraction.

Same design as `data/harness/orion-<line>-thinned.h5ad` except the cell selection:
  * targets: exactly the thinned harness's 132 (its .meta.json `perts`), so Δ(clean − thinned) on any fold instrument
    isolates the truth's construction (the key and the call set depend on the targets only, not on the truth);
  * perturbed cells: a seeded random 150 of each target's raw cells over ALL batches (rng [616, panel index]), where
    the thinned build took the first 150 in batch order;
  * controls: 18,400 drawn with 615's preregistered P' mixture (the selected perturbed cells' batch mixture) from the
    per-batch kept controls, seed 616, where the thinned build took the first batches' controls and thinned them.
    615 measured this pool's split-half blob share at 0.060 (h1-test 0.049; thinned 0.85).
Output: work/truth-<line>.h5 (obs/target_gene, var/_index, X/{data float32, indices, indptr} + shape; gates480 reads
it), work/truth-<line>.json.

    uv run python truth616.py HEK293T|HCT116
"""
from __future__ import annotations

import json
import os
import sys

import h5py
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
P615 = os.path.join(REPO, "experiments/615-orion-hek-harness/work")
N_CELL, N_CTRL, SEED = 150, 18400, 616


def main() -> None:
    line = sys.argv[1]
    raw = os.path.join(P615, {"HEK293T": "hek-raw.h5", "HCT116": "hct-raw.h5"}[line])
    meta = json.load(open(os.path.join(REPO, f"data/harness/orion-{line.lower()}-thinned.h5ad.meta.json")))
    targets = sorted(meta["perts"])
    assert len(targets) == meta["n_perturbations"] == 132, len(targets)
    panel = pd.read_csv(os.path.join(REPO, "data/challenge/pert_counts.csv"))["target_gene"].astype(str).tolist()
    pidx = {g: i for i, g in enumerate(panel)}
    with h5py.File(raw, "r") as h:
        lab = h["obs/target_gene"][:].astype(str)
        bat = h["obs/batch"][:].astype(np.int64)
        kind = h["obs/kind"][:]
        indptr = h["X/indptr"][:]
        genes = h["var/gene_name"][:]
        bnum = h["batches/number"][:]
        G = int(h["X"].attrs["shape"][1])
    pan = np.flatnonzero(kind == 1)
    rows = {}
    for t in targets:
        r = pan[lab[pan] == t]
        assert len(r) >= N_CELL, (t, len(r))
        rng = np.random.default_rng([SEED, pidx[t]])
        rows[t] = np.sort(rng.choice(r, size=N_CELL, replace=False))
    allp = np.concatenate(list(rows.values()))
    bpos = {int(x): i for i, x in enumerate(bnum)}
    w = np.bincount([bpos[int(x)] for x in bat[allp]], minlength=len(bnum)).astype(np.float64)
    q = N_CTRL * w / w.sum()
    nb = np.floor(q).astype(np.int64)
    nb[np.argsort(-(q - nb), kind="stable")[:N_CTRL - int(nb.sum())]] += 1
    rng = np.random.default_rng(SEED)
    ctrl = []
    for i, x in enumerate(bnum):
        avail = np.flatnonzero((kind == 0) & (bat == int(x)))
        assert nb[i] <= len(avail), (int(x), int(nb[i]), len(avail))
        if nb[i]:
            ctrl.append(np.sort(rng.choice(avail, size=int(nb[i]), replace=False)))
    ctrl = np.concatenate(ctrl)
    order = np.concatenate([rows[t] for t in targets] + [ctrl])
    labels = np.array(sum(([t] * N_CELL for t in targets), []) + ["non-targeting"] * N_CTRL)
    out = os.path.join(HERE, "work", f"truth-{line}.h5")
    lens = indptr[order + 1] - indptr[order]
    ip = np.concatenate([[0], np.cumsum(lens)]).astype(np.int64)
    lib = np.zeros(len(order))
    with h5py.File(raw, "r") as h, h5py.File(out + ".part", "w") as o:
        dd, di = h["X/data"], h["X/indices"]
        od = o.create_dataset("X/data", (int(ip[-1]),), dtype="f4")
        oi = o.create_dataset("X/indices", (int(ip[-1]),), dtype="i4")
        srt = np.argsort(order, kind="stable")
        for s0 in range(0, len(order), 2000):
            ks = srt[s0:s0 + 2000]
            r0, r1 = int(order[ks].min()), int(order[ks].max()) + 1
            a0 = int(indptr[r0])
            bd, bi = dd[a0:int(indptr[r1])], di[a0:int(indptr[r1])]
            for k in ks:
                r = int(order[k])
                s, e = int(indptr[r]) - a0, int(indptr[r + 1]) - a0
                od[ip[k]:ip[k + 1]] = bd[s:e].astype(np.float32)
                oi[ip[k]:ip[k + 1]] = bi[s:e]
                lib[k] = float(bd[s:e].sum())
        o.create_dataset("X/indptr", data=ip)
        o["X"].attrs["shape"] = np.array([len(order), G])
        o.create_dataset("obs/target_gene", data=labels.astype("S"))
        o.create_dataset("var/_index", data=genes)
    os.replace(out + ".part", out)
    isc = labels == "non-targeting"
    info = dict(line=line, targets=len(targets), cells_per_target=N_CELL, controls=int(isc.sum()), nnz=int(ip[-1]),
                library_median_pert=float(np.median(lib[~isc])), library_median_ctrl=float(np.median(lib[isc])),
                ctrl_batches_used=int((nb > 0).sum()), source=os.path.relpath(raw, REPO),
                design="thinned harness's 132 targets; random 150 cells over all batches; P' control pool")
    json.dump(info, open(os.path.join(HERE, "work", f"truth-{line}.json"), "w"), indent=1)
    print(json.dumps(info))


if __name__ == "__main__":
    sys.exit(main())
