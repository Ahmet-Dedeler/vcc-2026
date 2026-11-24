"""184 Phase 1 -- build the (treatment x cell line x gene) pseudobulk tensor + per-line basal.

Sources (data/external/tahoe/, downloaded 2026-09-03 from tahoebio/tahoe-de-rhaister):
  cell_eval/plate_plate{1..14}.parquet  one row per (cell_line, treatment, plate), 2,000 gene
                                        columns holding the pseudobulk delta
                                        mean(treated) - mean(plate-matched DMSO), library-size
                                        normalised (verified: ||delta|| is uncorrelated with the
                                        line's control library size, r = -0.36, while ||basal||
                                        is; see RESULTS.md).
  zeroshot/control_expression.parquet   per (cell_line, plate, gene) DMSO-control mean expression
                                        on the RAW count scale (sum over all 62,710 genes is
                                        0.69 x mean library size), plus library size and control
                                        cell count.  Converted here to a library fraction.

156 treatments appear on more than one plate, giving 7,641 (line, treatment) pairs with 2-3
independent measurements and 50 pairs with 14.  Those are kept separately: they are the only
same-line, same-perturbation replicate structure in the dataset and they carry the reliability
anchor.

Writes work/tensor.npz and work/reps.npz only.
"""
from __future__ import annotations

import hashlib
import json
import os
from collections import defaultdict

import numpy as np
import polars as pl

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
TAHOE = os.path.join(REPO, "data/external/tahoe")
WORK = os.path.join(HERE, "work")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    genes = None
    parts = []
    for p in range(1, 15):
        df = pl.read_parquet(os.path.join(TAHOE, f"cell_eval/plate_plate{p}.parquet"))
        g = df.columns[2:]
        if genes is None:
            genes = g
        assert g == genes, f"gene axis differs on plate {p}"
        parts.append((p, df))
        print(f"plate {p:2d}: {df.height:5d} rows, {df['cell_line'].n_unique():3d} lines, "
              f"{df['treatment'].n_unique():4d} treatments")
    genes = np.asarray(genes)

    lines = np.asarray(sorted(set().union(*[set(d["cell_line"].to_list()) for _, d in parts])))
    treats = np.asarray(sorted(set().union(*[set(d["treatment"].to_list()) for _, d in parts])))
    li = {l: i for i, l in enumerate(lines)}
    ti = {t: i for i, t in enumerate(treats)}

    X = np.zeros((len(treats), len(lines), len(genes)), dtype=np.float64)
    C = np.zeros((len(treats), len(lines)), dtype=np.int16)
    # group key for control-noise centring: the set of plates a (treatment, line) measurement
    # averages over.  Singletons get their own plate id; replicated pairs get a synthetic id
    # shared by every pair averaging the same plate set, because those share the same controls.
    GK = np.zeros((len(treats), len(lines)), dtype=np.int32)
    plate_sets = defaultdict(set)
    plate_seen = defaultdict(set)
    rep_store = defaultdict(list)          # (t_idx, l_idx) -> list of (plate, vector)
    for p, df in parts:
        ri = np.asarray([ti[t] for t in df["treatment"].to_list()])
        ci = np.asarray([li[l] for l in df["cell_line"].to_list()])
        V = df.select(genes.tolist()).to_numpy()
        np.add.at(X, (ri, ci), V)
        np.add.at(C, (ri, ci), 1)
        for a, b in zip(ri, ci):
            plate_sets[(int(a), int(b))].add(p)
        for t in df["treatment"].unique().to_list():
            plate_seen[t].add(p)
        del V, df
    # second pass, only for the replicated treatments (cheap, keeps peak RSS down)
    rep_treats = sorted([t for t, ps in plate_seen.items() if len(ps) > 1])
    rep_set = set(rep_treats)
    for p in range(1, 15):
        df = pl.read_parquet(os.path.join(TAHOE, f"cell_eval/plate_plate{p}.parquet"))
        df = df.filter(pl.col("treatment").is_in(pl.Series(list(rep_set)).implode()))
        if df.height == 0:
            continue
        V = df.select(genes.tolist()).to_numpy().astype(np.float32)
        for k, (t, l) in enumerate(zip(df["treatment"].to_list(), df["cell_line"].to_list())):
            rep_store[(ti[t], li[l])].append((p, V[k]))
        del V, df

    # bit p of GK[t, l] is set iff that measurement averages plate p (plates are 1..14)
    for (a, b), ps in plate_sets.items():
        GK[a, b] = int(sum(1 << p for p in ps))
    print(f"\nplate bitmask: {len(set(GK.ravel().tolist()))} distinct plate sets")

    n_pop = int((C > 0).sum())
    print(f"\n{len(treats)} treatments x {len(lines)} lines x {len(genes)} genes")
    print(f"populated cells {n_pop} / {len(treats)*len(lines)} = {n_pop/(len(treats)*len(lines)):.4f}")
    print(f"replicate counts: {np.bincount(C.ravel())}")
    X = np.where(C[:, :, None] > 0, X / np.maximum(C, 1)[:, :, None], np.nan).astype(np.float32)

    # ---- basal, per (cell_line, plate) and per line, as a fraction of the control library
    ce = pl.read_parquet(os.path.join(TAHOE, "zeroshot/control_expression.parquet"),
                         columns=["feature", "ref_mean", "plate", "cell_line", "ref_n_cells"])
    tot = (ce.group_by(["cell_line", "plate"]).agg(pl.col("ref_mean").sum().alias("tot")))
    ce = ce.filter(pl.col("feature").is_in(pl.Series(genes.tolist()).implode())).join(
        tot, on=["cell_line", "plate"], how="left")
    gi = {g: i for i, g in enumerate(genes)}
    B = np.full((len(lines), 15, len(genes)), np.nan, dtype=np.float32)
    N = np.zeros((len(lines), 15), dtype=np.int64)
    bl = np.asarray([li.get(l, -1) for l in ce["cell_line"].to_list()])
    bp = ce["plate"].to_numpy()
    bg = np.asarray([gi[g] for g in ce["feature"].to_list()])
    keep = bl >= 0
    frac = (ce["ref_mean"].to_numpy() / ce["tot"].to_numpy())[keep]
    B[bl[keep], bp[keep], bg[keep]] = frac
    N[bl[keep], bp[keep]] = ce["ref_n_cells"].to_numpy()[keep]
    w = N.astype(np.float64)
    Bl = np.nansum(np.nan_to_num(B) * w[:, :, None], axis=1) / np.maximum(w.sum(1), 1)[:, None]
    print(f"line basal (fraction of control library on the 2,000-gene axis): "
          f"row sums {Bl.sum(1).min():.4f}..{Bl.sum(1).max():.4f}, "
          f"n plates per line {int((N > 0).sum(1).min())}..{int((N > 0).sum(1).max())}")

    out = os.path.join(WORK, "tensor.npz")
    np.savez_compressed(out, delta=X, nrep=C, genes=genes, lines=lines, treatments=treats,
                        basal_plate=B, basal_line=Bl.astype(np.float32), ctrl_cells=N,
                        group=GK)
    np.save(os.path.join(WORK, "plate_of.npy"), GK)

    # the per-(line, plate) mean delta over that plate's ~94 drugs: the estimate of the
    # plate-matched-DMSO control-noise blob that every delta on that plate carries.
    PM = np.zeros((15, len(lines), len(genes)), dtype=np.float32)
    for p in range(1, 15):
        df = pl.read_parquet(os.path.join(TAHOE, f"cell_eval/plate_plate{p}.parquet"))
        V = df.select(genes.tolist()).to_numpy().astype(np.float32)
        ci = np.asarray([li[l] for l in df["cell_line"].to_list()])
        for l in np.unique(ci):
            PM[p, l] = V[ci == l].mean(0)
        del V, df
    np.save(os.path.join(WORK, "platemean.npy"), PM)

    keys = sorted(rep_store)
    RT = np.asarray([k[0] for k in keys], dtype=np.int32)
    RL = np.asarray([k[1] for k in keys], dtype=np.int32)
    nmax = max(len(v) for v in rep_store.values())
    RV = np.full((len(keys), nmax, len(genes)), np.nan, dtype=np.float32)
    RP = np.zeros((len(keys), nmax), dtype=np.int16)
    for i, k in enumerate(keys):
        for j, (p, v) in enumerate(sorted(rep_store[k])):
            RV[i, j] = v
            RP[i, j] = p
    rout = os.path.join(WORK, "reps.npz")
    np.savez_compressed(rout, t_idx=RT, l_idx=RL, delta=RV, plate=RP)
    print(f"replicate store: {len(keys)} (treatment, line) pairs, up to {nmax} plates each")

    meta = {
        "n_treatments": int(len(treats)), "n_lines": int(len(lines)), "n_genes": int(len(genes)),
        "n_populated": n_pop, "grid_fill": float(n_pop / (len(treats) * len(lines))),
        "n_replicated_pairs": int(len(keys)), "max_replicates": int(nmax),
        "tensor_sha256": sha256(out), "tensor_bytes": os.path.getsize(out),
        "reps_sha256": sha256(rout), "reps_bytes": os.path.getsize(rout),
        "source": "tahoebio/tahoe-de-rhaister cell_eval/* + zeroshot/control_expression.parquet",
    }
    with open(os.path.join(WORK, "tensor.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
