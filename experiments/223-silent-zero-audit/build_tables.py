"""223 Phase 2a — materialise donor-shaped tables for the partial-coverage sources that
do not already have one, INSIDE this experiment's work dir only.

Two of the sources named in the 21 hold their effect matrix in a non-standard layout, so
218's loader cannot read them.  Nothing outside `experiments/223-silent-zero-audit/` is
written, and no source array is modified -- each table is a re-layout of bytes already on
disk, with NaN (unmeasured) mapped to exact 0.0, which is how every other partial-coverage
source in the registry already stores an unmeasured coordinate.

  dld1     116's `transfer_lfc_{low1,low2}.npy` are (18533 genes x 4715 targets) and split
           by library.  116 measured Low1 vs Low2 split retrieval PDS 0.947 on the 67 VCC
           targets, so the two splits are averaged coordinatewise over their measured masks
           -- the same estimator 116 itself validated.  All 300 panel response coordinates
           are already exact zero in both splits (116's own check) and that is re-asserted.
  altair   119's `log2_fold_change.npy` is (1000 x 18533) with NaN at every unmeasured
           coordinate and three all-NaN rows (COQ10B, HSD17B12, NEK1), which are dropped.

    uv run python build_tables.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "work/tables"
OUT.mkdir(parents=True, exist_ok=True)

PANEL = [line.split(",")[0]
         for line in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:]
         if line.strip()]
PANEL_SET = set(PANEL)


def lines(p: Path) -> list[str]:
    return [x.strip() for x in p.read_text().splitlines() if x.strip()]


def write_table(name: str, delta: np.ndarray, perts: list[str], genes: list[str]) -> dict:
    d = OUT / name
    d.mkdir(parents=True, exist_ok=True)
    delta = np.ascontiguousarray(np.nan_to_num(delta, nan=0.0, posinf=0.0, neginf=0.0),
                                 dtype=np.float32)
    assert delta.shape == (len(perts), len(genes))
    np.save(d / "delta.npy", delta)
    (d / "perts.txt").write_text("\n".join(perts) + "\n")
    (d / "genes.txt").write_text("\n".join(genes) + "\n")
    # a flat control profile: the residual construction centres on the bank's own row mean
    # and never reads b_ctrl, but 218's Bank.load wants the path to exist.
    np.save(d / "b_ctrl.npy", np.zeros(len(genes), dtype=np.float32))
    nz = int(np.count_nonzero(delta))
    row_nz = (delta != 0).sum(1)
    meta = dict(name=name, shape=list(delta.shape), nonzero=nz,
                panel_targets=len(set(perts) & PANEL_SET),
                median_nonzero_genes_per_row=int(np.median(row_nz)),
                rows_all_zero=int((row_nz == 0).sum()))
    (d / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"[{name}] {delta.shape} panel {meta['panel_targets']}/300  "
          f"median nz genes/row {meta['median_nonzero_genes_per_row']}  "
          f"all-zero rows {meta['rows_all_zero']}")
    return meta


def build_dld1() -> dict:
    src = ROOT / "experiments/116-dld1-lowmoi/work/compact-source-v1"
    genes = lines(src / "canonical_genes.txt")
    perts = lines(src / "target_symbols.txt")
    lo1 = np.load(src / "transfer_lfc_low1.npy").T.astype(np.float64)
    lo2 = np.load(src / "transfer_lfc_low2.npy").T.astype(np.float64)
    m1 = np.load(src / "transfer_measured_low1.npy").T
    m2 = np.load(src / "transfer_measured_low2.npy").T
    tot = np.where(m1, np.nan_to_num(lo1), 0.0) + np.where(m2, np.nan_to_num(lo2), 0.0)
    cnt = m1.astype(np.int8) + m2.astype(np.int8)
    delta = np.zeros_like(tot)
    live = cnt > 0
    delta[live] = tot[live] / cnt[live]
    # 116 zeroes every panel response coordinate; assert it survived the re-layout.
    pos = {g: i for i, g in enumerate(genes)}
    pcols = np.asarray([pos[g] for g in PANEL if g in pos])
    assert np.all(delta[:, pcols] == 0.0), "116's panel-coordinate zeroing did not survive"
    del lo1, lo2, m1, m2, tot, cnt
    return write_table("dld1_lowmoi", delta, perts, genes)


def build_altair() -> dict:
    src = ROOT / "experiments/119-altair-a549/work/compact-source"
    genes = lines(src / "genes.txt")
    perts = lines(src / "perts.txt")
    lfc = np.load(src / "log2_fold_change.npy").astype(np.float64)
    keep = ~np.all(np.isnan(lfc), axis=1)
    dropped = [p for p, k in zip(perts, keep) if not k]
    lfc = lfc[keep]
    perts = [p for p, k in zip(perts, keep) if k]
    # panel response coordinates are zeroed by the residual builder anyway, but do it here
    # too so the table is self-consistent with dld1's and telohaec's.
    pos = {g: i for i, g in enumerate(genes)}
    pcols = np.asarray([pos[g] for g in PANEL if g in pos])
    lfc = np.nan_to_num(lfc, nan=0.0)
    lfc[:, pcols] = 0.0
    m = write_table("altair_a549", lfc, perts, genes)
    m["dropped_all_nan_rows"] = dropped
    (OUT / "altair_a549/meta.json").write_text(json.dumps(m, indent=2))
    return m


def main() -> None:
    out = dict(dld1_lowmoi=build_dld1(), altair_a549=build_altair())
    (HERE / "work/build_tables.json").write_text(json.dumps(out, indent=2))
    print(f"-> {HERE / 'work/build_tables.json'}")


if __name__ == "__main__":
    main()
