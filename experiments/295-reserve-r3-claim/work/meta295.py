"""295 item 1 — METADATA ONLY, before the preregistration is hashed.

Reads: the R3 manifest, the eight files' bytes (sha256 over the whole file — a digest, explicitly
free under the R3 access rule), 278's obs-only label counts, 292's frozen label grammars, 279's
held map, and the HDF5 *structure* of each reserve file (group names, `encoding-type` /
`encoding-version` / `shape` attributes, dataset shapes and dtypes, and the `var` index).

NO expression value is read.  `X/data`, `X/indices`, `X/indptr` and every `layers/*` dataset are
never sliced; only their attributes and shapes are inspected, which is the same class of read as
`obs` / `var` / `uns` and is declared as such in PREREGISTRATION.md §0.

Writes work/meta295.json, work/reserve_genes295.csv.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/295-reserve-r3-claim"
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/292-reserve-declaration/work"))
import parse292 as P  # noqa: E402  (frozen grammars + control pattern)

MANIFEST = ROOT / "experiments/292-reserve-declaration/work/reserve292.json"
C278 = ROOT / "experiments/278-depth-inventory/work/h5ad_counts"
HELD279 = ROOT / "experiments/279-znf-state-dependence/work/held279.npz"

MIN_CELLS = 50

# the reserve file -> its cell line, and the 28-bank keys whose LINE equals it (excluded from the
# consensus).  Frozen here from RESERVE.md §2 and 279's `lines` array; no outcome involved.
FILE_LINE = {
    "AdamsonWeissman2016_GSM2406681_10X010": "K562",
    "DatlingerBock2017": "Jurkat",
    "DatlingerBock2021": "Jurkat",
    "DixitRegev2016": "K562",
    "PapalexiSatija2021_eccite_RNA": "THP-1",
    "PapalexiSatija2021_eccite_arrayed_RNA": "THP-1",
    "ShifrutMarson2018": "primary CD8 T",
    "TianKampmann2019_day7neuron": "iPSC-neuron",
}
# extra exclusions beyond the line match, stated per file
EXTRA_EXCLUDE = {
    # same STUDY as the map bank tian2019_ipsc (Kampmann lab), and the differentiated state of the
    # same iPSC line as tian_ipsc_neuron -> both banks are excluded for this file.
    "TianKampmann2019_day7neuron": ["tian2019_ipsc", "tian_ipsc_neuron"],
}
MODALITY = {
    "AdamsonWeissman2016_GSM2406681_10X010": "CRISPRi",
    "DatlingerBock2017": "Cas9-KO",
    "DatlingerBock2021": "Cas9-KO",
    "DixitRegev2016": "CRISPRi",
    "PapalexiSatija2021_eccite_RNA": "Cas9-KO",
    "PapalexiSatija2021_eccite_arrayed_RNA": "Cas9-KO",
    "ShifrutMarson2018": "Cas9-KO",
    "TianKampmann2019_day7neuron": "CRISPRi",
}


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


def h5_structure(p: Path) -> dict:
    """Group/dataset NAMES, attributes, shapes and dtypes only.  No dataset is sliced."""
    out = {}
    with h5py.File(p, "r") as f:
        out["root_keys"] = sorted(f.keys())
        x = f["X"]
        if isinstance(x, h5py.Group):
            out["X"] = dict(kind="group", keys=sorted(x.keys()),
                            attrs={k: (v.tolist() if hasattr(v, "tolist") else str(v))
                                   for k, v in x.attrs.items()},
                            data_dtype=str(x["data"].dtype), data_len=int(x["data"].shape[0]),
                            indices_dtype=str(x["indices"].dtype),
                            indptr_len=int(x["indptr"].shape[0]))
        else:
            out["X"] = dict(kind="dataset", shape=list(x.shape), dtype=str(x.dtype),
                            attrs={k: str(v) for k, v in x.attrs.items()})
        out["layers"] = sorted(f["layers"].keys()) if "layers" in f else []
        out["obsm"] = sorted(f["obsm"].keys()) if "obsm" in f else []
        v = f["var"]
        out["var_keys"] = sorted(v.keys())
        out["var_attrs"] = {k: (v.attrs[k].tolist() if hasattr(v.attrs[k], "tolist")
                                else str(v.attrs[k])) for k in v.attrs}
        o = f["obs"]
        out["obs_keys"] = sorted(o.keys())
        out["obs_attrs"] = {k: (o.attrs[k].tolist() if hasattr(o.attrs[k], "tolist")
                                else str(o.attrs[k])) for k in o.attrs}
    return out


def var_symbols(p: Path) -> list[str]:
    """The var axis symbols (metadata)."""
    with h5py.File(p, "r") as f:
        v = f["var"]
        idx = v.attrs.get("_index", "_index")
        idx = idx if isinstance(idx, str) else idx.decode()
        d = v[idx]
        if isinstance(d, h5py.Group):          # categorical
            cats = np.asarray(d["categories"]).astype(str)
            return cats[np.asarray(d["codes"])].tolist()
        return np.asarray(d).astype(str).tolist()


def main() -> None:
    man = json.loads(MANIFEST.read_text())
    files = man["files"]
    z = np.load(HELD279, allow_pickle=True)
    held = set(np.asarray(z["genes"]).astype(str).tolist())
    axis = [str(g) for g in z["axis"]]
    keys = [str(k) for k in z["keys"]]
    lines = [str(l) for l in z["lines"]]

    out = {"files": {}, "n_axis": len(axis), "n_held": len(held)}
    rows = []
    for stem, rec in files.items():
        p = ROOT / rec["path"]
        got = sha256_file(p)
        ok = got == rec["sha256"]
        print(f"KA-1 sha256 {stem:<40} {'MATCH' if ok else 'MISMATCH'}")
        if not ok:
            raise SystemExit(f"digest mismatch on {stem}: {got} != {rec['sha256']}")

        counts = json.loads((C278 / f"scperturb__{stem}.json").read_text())["counts"]
        res = P.resolve(stem, counts)
        # rebuild the per-gene cell counts with the frozen grammar
        desc, fn = P.GRAMMAR.get(stem, P.DEFAULT)
        per_gene: dict[str, dict[str, int]] = {}
        ctrl_labels, ctrl_cells = [], 0
        for lab, n in counts.items():
            s = str(lab).strip()
            if P._ctrl(s):
                ctrl_labels.append(s); ctrl_cells += n
                continue
            g = fn(s)
            if g is None or g not in P.UNIVERSE:
                continue
            per_gene.setdefault(g, {})[s] = n
        cells = {g: sum(v.values()) for g, v in per_gene.items()}
        sel = sorted(g for g, n in cells.items() if n >= MIN_CELLS and g in held)

        vs = var_symbols(p)
        cov = sorted(set(vs) & set(axis), key=axis.index)
        struct = h5_structure(p)
        line = FILE_LINE[stem]
        excl = sorted({k for k, l in zip(keys, lines) if l == line}
                      | set(EXTRA_EXCLUDE.get(stem, [])))
        out["files"][stem] = dict(
            path=rec["path"], sha256=got, size_bytes=rec["size_bytes"], sha_ok=ok,
            line=line, modality=MODALITY[stem], grammar=desc,
            n_labels=len(counts), n_control_cells=ctrl_cells,
            n_control_labels=len(ctrl_labels), control_labels=sorted(ctrl_labels),
            n_genes_resolved=len(cells), n_genes_ge50=int(sum(1 for n in cells.values() if n >= MIN_CELLS)),
            n_genes_ge50_in_held=len(sel), genes=sel,
            cells_per_gene={g: int(cells[g]) for g in sel},
            guides_per_gene={g: len(per_gene[g]) for g in sel},
            labels_per_gene={g: sorted(per_gene[g]) for g in sel},
            n_var=len(vs), n_axis_covered=len(cov),
            excluded_banks=excl, n_banks_used=len(keys) - len(excl),
            structure=struct)
        for g in sel:
            rows.append(dict(file=stem, gene=g, cells=int(cells[g]),
                             guides=len(per_gene[g]), line=line, modality=MODALITY[stem]))
        print(f"    {stem:<40} genes>=50&held {len(sel):3d}  axis cov {len(cov):5d}/{len(axis)}  "
              f"ctrl {ctrl_cells:6d}  excl banks {excl}")

    pairs = len(rows)
    genes = sorted({r['gene'] for r in rows})
    dup = sorted({g for g in genes if sum(1 for r in rows if r['gene'] == g) > 1})
    out["n_pairs"] = pairs
    out["n_distinct_genes"] = len(genes)
    out["duplicated_genes"] = dup
    print(f"\n(file, gene) rows {pairs}; distinct genes {len(genes)}; genes in >1 file {len(dup)}")

    import csv as _csv
    with open(WORK / "reserve_genes295.csv", "w", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=["file", "gene", "cells", "guides", "line", "modality"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    (WORK / "meta295.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
