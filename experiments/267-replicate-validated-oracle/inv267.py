"""267 -- item 1, the INVENTORY (metadata only; no expression is read).

Per held line: per-target cell counts, number of distinct guide constructs per target, batch / gemgroup structure,
control cell counts and their batches.  Writes work/inv/<line>.json and work/inv/<line>-targets.csv, plus the per-cell
tables the heavy pass (`half267.py`) reuses so that the split assignment is fixed HERE, once:
  work/inv/<line>-cells.npz : per cell (in file order): pert code, batch code, control flag, half assignment for the
                              batch split (parity of the batch number) and for the random-cell split (seed 20260905).

Orion: scalar parquet columns only (cell_barcode, sample, guide_target, gene_target, total_counts, pass_guide_filter).
K562-GWPS: the h5ad obs (perturbation, guide_id, batch = gemgroup).  Nothing written outside this experiment.
"""
from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/267-replicate-validated-oracle"
INV = HERE / "work/inv"
INV.mkdir(parents=True, exist_ok=True)
SEED = 20260905
CONTROL_PATTERNS = (r"^non.?targeting", r"^ntc", r"^control", r"^safe.?harbou?r", r"^neg")   # build_orion_donor.py
PANEL = [l.split(",")[0] for l in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:] if l.strip()]
assert len(PANEL) == 300


def is_control(name: str) -> bool:
    n = str(name).lower()
    return any(re.match(p, n) for p in CONTROL_PATTERNS)


def summarise(name: str, cells: pd.DataFrame, bank_perts: set[str], guide_col: str) -> dict:
    """cells: one row per cell with columns pert, guide, batch, ctrl, half_batch, half_rand."""
    ctrl = cells[cells.ctrl]
    pert = cells[~cells.ctrl]
    per_batch_ctrl = ctrl.groupby("batch").size()
    rows = []
    for t in PANEL:
        sub = pert[pert.pert == t]
        if not len(sub):
            continue
        g = sub.groupby("guide").size().sort_values(ascending=False)
        rows.append(dict(
            target=t, cells=int(len(sub)), in_bank=bool(t in bank_perts), guides=int(len(g)),
            top_guide_share=float(g.iloc[0] / len(sub)), guide_ids="|".join(map(str, g.index.tolist())),
            batches=int(sub.batch.nunique()), max_cells_per_batch=int(sub.groupby("batch").size().max()),
            cells_half_batch_A=int((sub.half_batch == 0).sum()), cells_half_batch_B=int((sub.half_batch == 1).sum()),
            cells_half_rand_A=int((sub.half_rand == 0).sum()), cells_half_rand_B=int((sub.half_rand == 1).sum()),
        ))
    df = pd.DataFrame(rows)
    df.to_csv(INV / f"{name}-targets.csv", index=False)
    inb = df[df.in_bank]
    out = dict(
        line=name, n_cells=int(len(cells)), n_control_cells=int(len(ctrl)), n_batches=int(cells.batch.nunique()),
        n_perturbations_seen=int(pert.pert.nunique()), n_bank_perturbations=len(bank_perts),
        panel_targets_with_cells=int(len(df)), panel_targets_in_bank=int(inb.shape[0]),
        guides_per_target=inb.guides.value_counts().sort_index().to_dict(),
        targets_with_ge2_guides=inb[inb.guides >= 2].target.tolist(),
        cells_per_target=dict(min=int(inb.cells.min()), q25=float(inb.cells.quantile(.25)), median=float(inb.cells.median()),
                              q75=float(inb.cells.quantile(.75)), max=int(inb.cells.max())),
        batches_per_target=dict(min=int(inb.batches.min()), median=float(inb.batches.median()), max=int(inb.batches.max())),
        max_cells_per_batch_per_target=dict(median=float(inb.max_cells_per_batch.median()), max=int(inb.max_cells_per_batch.max())),
        control_cells_per_batch=dict(min=int(per_batch_ctrl.min()), median=float(per_batch_ctrl.median()), max=int(per_batch_ctrl.max()),
                                     n_batches_with_controls=int(len(per_batch_ctrl))),
        control_guides=int(ctrl.guide.nunique()),
        control_share_of_cells=float(len(ctrl) / len(cells)),
        half_batch=dict(cells_A=int((cells.half_batch == 0).sum()), cells_B=int((cells.half_batch == 1).sum()),
                        controls_A=int((ctrl.half_batch == 0).sum()), controls_B=int((ctrl.half_batch == 1).sum()),
                        batches_A=int(cells[cells.half_batch == 0].batch.nunique()), batches_B=int(cells[cells.half_batch == 1].batch.nunique()),
                        panel_min_cells_in_a_half=int(min(inb.cells_half_batch_A.min(), inb.cells_half_batch_B.min())),
                        panel_median_cells_A=float(inb.cells_half_batch_A.median()), panel_median_cells_B=float(inb.cells_half_batch_B.median())),
        half_rand=dict(cells_A=int((cells.half_rand == 0).sum()), controls_A=int((ctrl.half_rand == 0).sum()),
                       panel_min_cells_in_a_half=int(min(inb.cells_half_rand_A.min(), inb.cells_half_rand_B.min())),
                       panel_median_cells_A=float(inb.cells_half_rand_A.median())),
        split_unit={"(a) by guide": int((inb.guides >= 2).sum()), "(b) by batch": int((inb.guides < 2).sum()),
                    "(c) random cells": "comparison split on every target"},
        seed=SEED,
    )
    (INV / f"{name}.json").write_text(json.dumps(out, indent=1))
    return out


def orion(line: str) -> None:
    t0 = time.time()
    files = sorted((ROOT / "data/external/orion/data").glob(f"{line}_Batch*.parquet"),
                   key=lambda p: int(re.search(r"Batch(\d+)", p.name)[1]))
    bank = set((ROOT / f"data/donors/orion_{line.lower()}/perts.txt").read_text().split("\n"))
    parts = []
    rng = np.random.default_rng(SEED)
    per_file = []
    for i, f in enumerate(files, 1):
        bno = int(re.search(r"Batch(\d+)", f.name)[1])
        t = pq.read_table(f, columns=["cell_barcode", "sample", "guide_target", "gene_target", "total_counts",
                                      "pass_guide_filter"]).to_pandas()
        n_all = len(t)
        t = t[t["pass_guide_filter"] == 1].reset_index(drop=True)
        assert (t["sample"] == f"{line}_Batch{bno}").all(), f"{f.name}: sample column != file name"
        ctrl = np.fromiter((is_control(x) for x in t["gene_target"].to_numpy()), bool, len(t))
        half_rand = rng.integers(0, 2, size=len(t)).astype(np.int8)          # one draw per cell, in file order
        parts.append(pd.DataFrame(dict(pert=t["gene_target"].to_numpy(), guide=t["guide_target"].to_numpy(),
                                       batch=np.full(len(t), bno, np.int32), ctrl=ctrl,
                                       half_batch=np.full(len(t), bno % 2, np.int8),   # odd batch -> 1 (B), even -> 0 (A)
                                       half_rand=half_rand, total=t["total_counts"].to_numpy())))
        per_file.append(dict(file=f.name, batch=bno, rows=n_all, kept=len(t), controls=int(ctrl.sum()),
                             median_total_counts=float(np.median(t["total_counts"]))))
        if i % 25 == 0 or i == len(files):
            print(f"  [{line}] {i}/{len(files)} files, {sum(len(p) for p in parts):,} cells, {time.time() - t0:.0f}s", flush=True)
    cells = pd.concat(parts, ignore_index=True)
    # NOTE the parity convention: batch number even -> half 0 = 'A', odd -> half 1 = 'B'.
    codes_p, pert_names = pd.factorize(cells.pert, sort=True)
    np.savez_compressed(INV / f"orion_{line.lower()}-cells.npz", pert_code=codes_p.astype(np.int32),
                        pert_names=np.asarray(pert_names, dtype=object), batch=cells.batch.to_numpy(),
                        ctrl=cells.ctrl.to_numpy(), half_batch=cells.half_batch.to_numpy(), half_rand=cells.half_rand.to_numpy(),
                        total=cells.total.to_numpy(), file_order=np.asarray([p["file"] for p in per_file], dtype=object),
                        file_kept=np.asarray([p["kept"] for p in per_file]))
    out = summarise(f"orion_{line.lower()}", cells, bank, "guide")
    out["files"] = per_file
    out["seconds"] = time.time() - t0
    (INV / f"orion_{line.lower()}.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k != "files"}, indent=1), flush=True)


def k562() -> None:
    import h5py
    t0 = time.time()
    P = ROOT / "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad"
    bank = set((ROOT / "data/donors/k562_gwps/perts.txt").read_text().split("\n"))
    with h5py.File(P, "r") as f:
        pc = f["obs/perturbation/codes"][:]
        pcat = f["obs/perturbation/categories"][:].astype(str)
        gc = f["obs/guide_id/codes"][:]
        gcat = f["obs/guide_id/categories"][:].astype(str)
        batch = f["obs/batch"][:].astype(np.int32)
        ncounts = f["obs/ncounts"][:]
    rng = np.random.default_rng(SEED)
    half_rand = rng.integers(0, 2, size=len(pc)).astype(np.int8)
    ctrl = pcat[pc] == "control"
    cells = pd.DataFrame(dict(pert=pcat[pc], guide=gcat[gc], batch=batch, ctrl=ctrl,
                              half_batch=(batch % 2).astype(np.int8), half_rand=half_rand, total=ncounts))
    # guide split for the two-construct targets: the construct whose name sorts first -> A (deterministic)
    np.savez_compressed(INV / "k562_gwps-cells.npz", pert_code=pc.astype(np.int32), pert_names=np.asarray(pcat, dtype=object),
                        guide_code=gc.astype(np.int32), guide_names=np.asarray(gcat, dtype=object), batch=batch, ctrl=ctrl,
                        half_batch=cells.half_batch.to_numpy(), half_rand=half_rand, total=ncounts)
    out = summarise("k562_gwps", cells, bank, "guide")
    out["seconds"] = time.time() - t0
    (INV / "k562_gwps.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1), flush=True)


if __name__ == "__main__":
    which = sys.argv[1:] or ["k562_gwps", "HCT116", "HEK293T"]
    for w in which:
        if w == "k562_gwps":
            k562()
        else:
            orion(w)
