"""615 stage 1: the harness truth (CSR h5 in 614's format) from work/hek-raw.h5, ONLY after the blob test PASSES.

  * perturbed: the 100 harness targets' exact cells (work/harness615_rows.npz, written by blob615.py: each target's
    seeded <= 250-cell subset, both halves), so the truth is the object the blob test read;
  * controls: 18,400 drawn from the per-batch control samples (kind 0) with the PREREGISTERED P' mixture, n_b
    proportional to the harness's perturbed cells in batch b (largest remainder), seed 615. Asserts n_b <= available.
Output: work/hek-truth.h5 (X/{data float32, indices int32, indptr int64} + shape attr, obs/target_gene, var/_index),
work/truth615.json. Refuses to run unless work/blob615.json says PASS.

    uv run python truth615.py
"""
from __future__ import annotations

import json
import os
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
N_CTRL, SEED = 18400, 615


def main() -> None:
    v = json.load(open(os.path.join(WORK, "blob615.json")))
    if v["verdict"] != "PASS":
        raise SystemExit(f"blob test verdict is {v['verdict']!r}: no truth is built (PREREGISTRATION-blob.md)")
    z = np.load(os.path.join(WORK, "harness615_rows.npz"), allow_pickle=False)
    targets = [str(t) for t in z["targets"]]
    with open(os.path.join(WORK, "targets615.txt")) as f:
        assert [ln.strip() for ln in f if ln.strip()] == targets, "targets615.txt differs from the saved harness rows"
    src = os.path.join(WORK, "hek-raw.h5")
    with h5py.File(src, "r") as h:
        lab = h["obs/target_gene"][:].astype(str)
        bat = h["obs/batch"][:].astype(np.int64)
        kind = h["obs/kind"][:]
        indptr = h["X/indptr"][:]
        genes = h["var/gene_name"][:]
        bnum = h["batches/number"][:]
        G = int(h["X"].attrs["shape"][1])
    pert_rows = {t: np.sort(np.concatenate([z[f"{t}__0"], z[f"{t}__1"]])) for t in targets}
    for t, r in pert_rows.items():
        assert np.all(lab[r] == t) and np.all(kind[r] == 1), t
    allp = np.concatenate(list(pert_rows.values()))
    # P' allocation of the control pool over batches
    bpos = {int(x): i for i, x in enumerate(bnum)}
    w = np.bincount([bpos[int(x)] for x in bat[allp]], minlength=len(bnum)).astype(np.float64)
    q = N_CTRL * w / w.sum()
    nb = np.floor(q).astype(np.int64)
    rem = N_CTRL - int(nb.sum())
    nb[np.argsort(-(q - nb), kind="stable")[:rem]] += 1
    rng = np.random.default_rng(SEED)
    ctrl_rows = []
    short = []
    for i, x in enumerate(bnum):
        avail = np.flatnonzero((kind == 0) & (bat == int(x)))
        if nb[i] > len(avail):
            short.append((int(x), int(nb[i]), len(avail)))
            continue
        if nb[i]:
            ctrl_rows.append(np.sort(rng.choice(avail, size=int(nb[i]), replace=False)))
    if short:
        raise SystemExit(f"P' needs more controls than kept in batches {short[:5]}")
    ctrl_rows = np.concatenate(ctrl_rows)
    assert len(ctrl_rows) == N_CTRL
    order = np.concatenate([pert_rows[t] for t in targets] + [ctrl_rows])
    labels = np.concatenate([[t] * len(pert_rows[t]) for t in targets] + [["non-targeting"] * N_CTRL])

    out = os.path.join(WORK, "hek-truth.h5")
    tmp = out + ".part"
    with h5py.File(src, "r") as h, h5py.File(tmp, "w") as o:
        dd, di = h["X/data"], h["X/indices"]
        lens = indptr[order + 1] - indptr[order]
        ip = np.concatenate([[0], np.cumsum(lens)]).astype(np.int64)
        od = o.create_dataset("X/data", (int(ip[-1]),), dtype="f4")
        oi = o.create_dataset("X/indices", (int(ip[-1]),), dtype="i4")
        # read in row order of the source (sorted) for locality, write to the output positions
        src_sorted = np.argsort(order, kind="stable")
        step = 2000
        lib = np.zeros(len(order))
        for s0 in range(0, len(order), step):
            ks = src_sorted[s0:s0 + step]
            r0, r1 = int(order[ks].min()), int(order[ks].max()) + 1
            a0, a1 = int(indptr[r0]), int(indptr[r1])
            bd, bi = dd[a0:a1], di[a0:a1]
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
    os.replace(tmp, out)
    isc = labels == "non-targeting"
    info = dict(targets=len(targets), cells=int(len(order)), controls=int(isc.sum()),
                cells_per_target_median=float(np.median([len(pert_rows[t]) for t in targets])),
                cells_per_target_min=int(min(len(pert_rows[t]) for t in targets)),
                library_median_pert=float(np.median(lib[~isc])), library_median_ctrl=float(np.median(lib[isc])),
                library_q25_q50_q75_ctrl=[float(np.quantile(lib[isc], q_)) for q_ in (0.25, 0.5, 0.75)],
                ctrl_batches_used=int((nb > 0).sum()), ctrl_per_batch_max=int(nb.max()), nnz=int(ip[-1]),
                source="work/hek-raw.h5", rows="work/harness615_rows.npz", mixture="P' (PREREGISTRATION-blob.md)")
    with open(os.path.join(WORK, "truth615.json"), "w") as f:
        json.dump(info, f, indent=1)
    print(json.dumps(info))


if __name__ == "__main__":
    sys.exit(main())
