"""278 item 1 -- OBS-ONLY depth scan of every genetic-perturbation h5ad on disk.

Reads ONLY the `obs` group of each file (h5py, categorical codes + categories); never
touches `X`.  For each file it writes, per perturbation label, the number of cells (and,
where the file has a cell-line / donor / batch / guide column, the same broken down).

Output: work/h5ad_counts/<key>.json  -- {"file", "n_cells", "label_col", "counts": {label: n},
                                        "by_line": {line: {label: n}}, "n_batches", "n_guides"}

    PYTHONDONTWRITEBYTECODE=1 uv run python scan278_h5ad.py
"""
from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

import h5py
import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/278-depth-inventory"
WORK = HERE / "work"
OUT = WORK / "h5ad_counts"
OUT.mkdir(parents=True, exist_ok=True)

PANEL = [l.strip() for l in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:] if l.strip()]
PANEL_SET = set(PANEL)

# candidate obs columns, in priority order
LABEL_CANDIDATES = ["perturbation", "gene", "gene_target", "target_gene", "target",
                    "guide_identity", "sgRNA_gene", "grna_target", "knockout", "gene_symbol",
                    "sgrna_symbol"]
LINE_CANDIDATES = ["cell_line", "cell_type", "celltype", "line", "Cell_Line", "donor", "sample_line"]
BATCH_CANDIDATES = ["batch", "gemgroup", "gem_group", "sample", "channel", "lane", "replicate", "Batch"]
GUIDE_CANDIDATES = ["guide_id", "sgRNA", "guide", "sgID_AB", "guide_identity", "Guide_Call", "grna", "protospacer"]


def read_col(obs: h5py.Group, name: str) -> np.ndarray | None:
    """Return an obs column as a string array, or None."""
    if name not in obs:
        return None
    node = obs[name]
    if isinstance(node, h5py.Group):  # categorical
        if "categories" not in node or "codes" not in node:
            return None
        cats = node["categories"][:]
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in cats], dtype=object)
        codes = node["codes"][:]
        out = np.empty(len(codes), dtype=object)
        ok = codes >= 0
        out[ok] = cats[codes[ok]]
        out[~ok] = ""
        return out
    arr = node[:]
    if arr.dtype.kind in "SO":
        return np.asarray([a.decode() if isinstance(a, bytes) else str(a) for a in arr], dtype=object)
    return np.asarray([str(a) for a in arr], dtype=object)


def pick(obs: h5py.Group, cands: list[str]) -> str | None:
    keys = list(obs.keys())
    lower = {k.lower(): k for k in keys}
    for c in cands:
        if c in keys:
            return c
        if c.lower() in lower:
            return lower[c.lower()]
    return None


def scan(key: str, path: Path) -> dict:
    t0 = time.time()
    rec: dict = {"key": key, "file": str(path.relative_to(ROOT)), "exists": path.exists()}
    if not path.exists():
        return rec
    with h5py.File(path, "r") as f:
        if "obs" not in f:
            rec["error"] = "no obs group"
            return rec
        obs = f["obs"]
        rec["obs_cols"] = sorted(obs.keys())
        lab_c = pick(obs, LABEL_CANDIDATES)
        rec["label_col"] = lab_c
        if lab_c is None:
            return rec
        labels = read_col(obs, lab_c)
        if labels is None:
            rec["error"] = f"cannot read {lab_c}"
            return rec
        rec["n_cells"] = int(len(labels))
        uniq, counts = np.unique(labels.astype(str), return_counts=True)
        rec["n_labels"] = int(len(uniq))
        rec["counts"] = {str(u): int(c) for u, c in zip(uniq, counts)}
        # panel restriction (exact symbol match)
        rec["panel_counts"] = {str(u): int(c) for u, c in zip(uniq, counts) if str(u) in PANEL_SET}
        rec["n_panel"] = len(rec["panel_counts"])
        # aux columns
        for tag, cands in (("line", LINE_CANDIDATES), ("batch", BATCH_CANDIDATES), ("guide", GUIDE_CANDIDATES)):
            col = pick(obs, cands)
            rec[f"{tag}_col"] = col
            if col is None or col == lab_c:
                continue
            vals = read_col(obs, col)
            if vals is None:
                continue
            vu = np.unique(vals.astype(str))
            rec[f"n_{tag}s"] = int(len(vu))
            rec[f"{tag}_vals"] = sorted(vu.tolist())[:40]
            if tag == "line" and len(vu) > 1 and rec["n_panel"] > 0:
                by = {}
                sv = vals.astype(str)
                sl = labels.astype(str)
                for v in vu:
                    m = sv == v
                    lu, lc = np.unique(sl[m], return_counts=True)
                    d = {str(a): int(b) for a, b in zip(lu, lc) if str(a) in PANEL_SET}
                    if d:
                        by[str(v)] = d
                rec["panel_by_line"] = by
            if tag == "guide" and rec["n_panel"] > 0:
                sv = vals.astype(str)
                sl = labels.astype(str)
                g = {}
                for t in rec["panel_counts"]:
                    g[t] = int(len(np.unique(sv[sl == t])))
                rec["guides_per_panel_target"] = g
    rec["seconds"] = round(time.time() - t0, 1)
    return rec


FILES: list[tuple[str, str]] = []
sc = ROOT / "data/external/scperturb"
for p in sorted(sc.glob("*.h5ad")):
    FILES.append((f"scperturb::{p.stem}", str(p)))
FILES += [
    ("nadig_jurkat", "data/external/nadig2025/jurkat_raw_singlecell.h5ad"),
    ("nadig_hepg2", "data/external/nadig2025/hepg2_raw_singlecell.h5ad"),
    ("altair_a549", "data/external/altair-a549/altair-scrna.h5ad"),
    ("myllia_a549_ch1", "data/external/myllia-a549/GSE337804_A549_ch_1.RDS.h5ad"),
    ("myllia_a549_ch2", "data/external/myllia-a549/GSE337804_A549_ch_2.RDS.h5ad"),
    ("igvf_telohaec", "data/external/igvf-telohaec/IGVFFI8645YMGX.h5ad"),
]
for extra in sorted((ROOT / "data/external/hcasmc-igvf").rglob("*.h5ad")):
    FILES.append((f"hcasmc::{extra.stem}", str(extra)))


def main() -> None:
    only = sys.argv[1] if len(sys.argv) > 1 else None
    index = []
    for key, rel in FILES:
        if only and only not in key:
            continue
        path = Path(rel) if Path(rel).is_absolute() else ROOT / rel
        safe = key.replace("::", "__").replace("/", "_")
        dest = OUT / f"{safe}.json"
        if dest.exists():
            rec = json.loads(dest.read_text())
        else:
            try:
                rec = scan(key, path)
            except Exception as exc:  # noqa: BLE001
                rec = {"key": key, "file": rel, "error": f"{type(exc).__name__}: {exc}"}
            dest.write_text(json.dumps(rec))
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30
        print(f"{key:52s} cells={rec.get('n_cells', '?'):>10} panel={rec.get('n_panel', '?'):>4} "
              f"lab={rec.get('label_col')} rss={rss:.2f}GiB", flush=True)
        index.append({k: rec.get(k) for k in
                      ("key", "file", "n_cells", "n_labels", "n_panel", "label_col", "line_col",
                       "batch_col", "guide_col", "n_lines", "n_batchs", "n_guides", "error")})
    (WORK / "h5ad_index.json").write_text(json.dumps(index, indent=1))
    print("peak RSS GiB", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30)


if __name__ == "__main__":
    main()
