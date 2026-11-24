"""225 Phase 1 — the census.  Panel overlap for EVERY perturbation asset reachable from
this machine, computed from each asset's own file, never from an abstract.

166 built the same census against a >=100-target bar, because that was the bar for a fifth
*panel-covering donor line*.  223 then priced `dld1_lowmoi` at 67/300 and found it live, so
that bar is wrong for the residual channel and the sub-100 class has to be enumerated too.
This census therefore reports EVERY asset with >=1 panel target, and the count is computed
at the >=20-cells-per-group threshold the donor builder actually uses — not at the deposit's
advertised library size, which is always larger.

    uv run python census225.py
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"

PANEL = {line.split(",")[0]
         for line in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:]
         if line.strip()}
assert len(PANEL) == 300
MIN_CELLS = 20


def _dec(a):
    if a.dtype.kind in "SO":
        return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in a])
    return np.asarray(a)


def read_col(h, name):
    g = h["obs"][name]
    if isinstance(g, h5py.Group):
        cats = _dec(g["categories"][:])
        return cats[np.clip(g["codes"][:], 0, len(cats) - 1)]
    return _dec(g[:])


def h5ad_row(p: Path) -> dict:
    with h5py.File(p, "r") as h:
        obs = list(h["obs"].keys())
        pc = next((c for c in ("perturbation", "target_gene", "gene") if c in obs), None)
        if pc is None:
            return dict(file=str(p.relative_to(ROOT)), status="no perturbation column",
                        obs_columns=obs[:12])
        lab = read_col(h, pc)
        u, c = np.unique(lab, return_counts=True)
        cnt = dict(zip(u.tolist(), c.tolist()))
        usable = [x for x in u.tolist()
                  if cnt[x] >= MIN_CELLS and x not in ("control", "non-targeting")]
        hit = sorted(PANEL & set(usable))
        line = None
        for cc in ("cell_line", "celltype", "tissue_type"):
            if cc in obs:
                vals = sorted(set(read_col(h, cc).tolist()))
                line = ",".join(vals[:3]) + ("…" if len(vals) > 3 else "")
                break
        X = h["X"]
        enc = (X.attrs.get("encoding-type", "?") if isinstance(X, h5py.Group) else "dense")
        return dict(file=str(p.relative_to(ROOT)), status="ok", pert_column=pc,
                    gb=round(p.stat().st_size / 2**30, 3), n_cells=int(len(lab)),
                    n_labels=int(len(u)), n_usable=len(usable),
                    panel_overlap=len(hit), panel_targets=hit, cell_line=line,
                    control_cells=int(cnt.get("control", cnt.get("non-targeting", 0))),
                    storage=str(enc))


def table_row(d: Path) -> dict:
    perts = [x for x in (d / "perts.txt").read_text().split("\n") if x.strip()]
    genes = [x for x in (d / "genes.txt").read_text().split("\n") if x.strip()]
    delta = np.load(d / "delta.npy", mmap_mode="r")
    hit = sorted(PANEL & set(perts))
    meta = {}
    if (d / "meta.json").exists():
        meta = json.loads((d / "meta.json").read_text())
    return dict(file=str(d.relative_to(ROOT)), status="donor-table",
                shape=list(delta.shape), n_targets=len(perts), n_genes=len(genes),
                panel_overlap=len(hit), panel_targets=hit,
                median_nz_genes_per_panel_row=meta.get(
                    "median_nonzero_genes_per_PANEL_row"),
                effect_norm_median_panel=meta.get("effect_norm_median_panel"))


def main() -> None:
    rows = []
    for p in sorted((ROOT / "data/external/scperturb").glob("*.h5ad")):
        rows.append(h5ad_row(p))
        print(f"{p.name:60s} panel {rows[-1].get('panel_overlap')}", flush=True)
    for rel in ("data/external/vcc2025/train/adata_Training.h5ad",
                "data/external/vcc2025/validation/adata_Validation.h5ad",
                "data/external/vcc2025/test/adata_Test.h5ad"):
        p = ROOT / rel
        if p.exists():
            rows.append(h5ad_row(p))
            print(f"{rel:60s} panel {rows[-1].get('panel_overlap')}", flush=True)
    for d in sorted(list((HERE / "work/tables").glob("*"))
                    + list((ROOT / "experiments/223-silent-zero-audit/work/tables")
                           .glob("*"))
                    + list((ROOT / "data/donors").glob("*"))):
        if d.is_dir() and (d / "delta.npy").exists():
            rows.append(table_row(d))
            print(f"{d.name:60s} TABLE panel {rows[-1]['panel_overlap']}", flush=True)

    # flat DE deposits that are not h5ad and not yet a table
    extra = []
    p = ROOT / "data/external/igvf-cardiomyocyte/IGVFFI0830FXFI.csv.gz"
    if p.exists():
        import csv
        t = set()
        with gzip.open(p, "rt") as f:
            for r in csv.DictReader(f):
                if r["intended_target_name"].strip():
                    t.add(r["intended_target_name"].strip())
        extra.append(dict(file=str(p.relative_to(ROOT)), status="flat-DE",
                          n_targets=len(t), panel_overlap=len(PANEL & t)))
    rows += extra

    (WORK / "census225.json").write_text(json.dumps(rows, indent=2))
    hits = [r for r in rows if (r.get("panel_overlap") or 0) >= 1]
    hits.sort(key=lambda r: -(r.get("panel_overlap") or 0))
    print(f"\n{len(rows)} assets, {len(hits)} with >=1 panel target")
    for r in hits:
        print(f"  {r.get('panel_overlap'):>3d}  {r['file']}")
    print(f"-> {WORK / 'census225.json'}")


if __name__ == "__main__":
    main()
