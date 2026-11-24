"""590 — build control pseudobulk profiles of every labelled cell line on disk, on the official gene axis.

One row per (dataset, line): summed raw counts over the line's control cells, on `data/challenge/gene_names.csv`'s
18,533-gene axis (NaN where the dataset does not carry the gene). The matcher (`match590.py`) turns a row into log2 CPM
over its present genes, the way `def_identity.py` treats a D/E/F context (whose file carries only the official genes).

Matrices are read with h5py in row or column chunks (never materialised): every process stays well under 4 GB.
Writes `work/profiles/<dataset>.npz` (counts float64 [n_rows, 18533], names, meta JSON) and prints a one-line summary
per row. Usage: `uv run python profiles590.py [dataset ...]` (default: all).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUT = HERE / "work" / "profiles"
SCP = ROOT / "data/external/scperturb"
OFFICIAL = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).tolist()
OFF_POS = {g: i for i, g in enumerate(OFFICIAL)}
MIN_CELLS = 50
CHUNK_NNZ = 12_000_000


def _str(a) -> list[str]:
    return [x.decode() if isinstance(x, bytes) else str(x) for x in a]


def obs_column(h: h5py.File, key: str) -> np.ndarray:
    """An obs column as an array of python strings (categorical or plain)."""
    o = h["obs"][key]
    if isinstance(o, h5py.Group) and "categories" in o:
        cats = np.asarray(_str(o["categories"][:]) + ["__NA__"], dtype=object)
        codes = o["codes"][:].astype(np.int64)
        codes[codes < 0] = len(cats) - 1
        return cats[codes]
    return np.asarray(_str(o[:]), dtype=object)


def var_names(h: h5py.File, prefer: tuple[str, ...] = ()) -> list[str]:
    v = h["var"]
    for k in prefer + ("_index", "gene_symbol", "genes", "features", "gene_name"):
        if k in v:
            o = v[k]
            if isinstance(o, h5py.Group) and "categories" in o:
                cats = _str(o["categories"][:])
                return [cats[c] for c in o["codes"][:]]
            return _str(o[:])
    raise KeyError("no var names")


def group_sums(h: h5py.File, groups: np.ndarray, n_groups: int) -> np.ndarray:
    """Sum X's rows by group id (−1 = drop) without materialising X. Returns [n_groups, n_var] float64."""
    X = h["X"]
    if isinstance(X, h5py.Dataset):  # dense
        n_obs, n_var = X.shape
        out = np.zeros((n_groups, n_var))
        keep = np.flatnonzero(groups >= 0)
        step = max(1, int(2e8 // max(n_var, 1)))
        for s in range(0, keep.size, step):
            idx = keep[s:s + step]
            block = X[idx, :] if idx.size > 1 else X[idx[0]:idx[0] + 1, :]
            np.add.at(out, groups[idx], np.asarray(block, dtype=np.float64))
        return out
    enc = X.attrs["encoding-type"]
    shape = tuple(int(x) for x in X.attrs["shape"])
    indptr = X["indptr"][:].astype(np.int64)
    data, indices = X["data"], X["indices"]
    if enc == "csr_matrix":
        n_obs, n_var = shape
        out = np.zeros(n_groups * n_var)
        r0 = 0
        while r0 < n_obs:
            r1 = int(np.searchsorted(indptr, indptr[r0] + CHUNK_NNZ, side="right")) - 1
            r1 = min(max(r1, r0 + 1), n_obs)
            g = groups[r0:r1]
            if (g >= 0).any():
                a, b = indptr[r0], indptr[r1]
                d = data[a:b].astype(np.float64)
                j = indices[a:b].astype(np.int64)
                rg = np.repeat(g, np.diff(indptr[r0:r1 + 1]))
                m = rg >= 0
                out += np.bincount(rg[m] * n_var + j[m], weights=d[m], minlength=n_groups * n_var)
            r0 = r1
        return out.reshape(n_groups, n_var)
    if enc == "csc_matrix":
        n_obs, n_var = shape
        out = np.zeros((n_groups, n_var))
        c0 = 0
        while c0 < n_var:
            c1 = int(np.searchsorted(indptr, indptr[c0] + CHUNK_NNZ, side="right")) - 1
            c1 = min(max(c1, c0 + 1), n_var)
            a, b = indptr[c0], indptr[c1]
            d = data[a:b].astype(np.float64)
            rows = indices[a:b].astype(np.int64)
            cols = np.repeat(np.arange(c0, c1), np.diff(indptr[c0:c1 + 1]))
            g = groups[rows]
            m = g >= 0
            w = c1 - c0
            blk = np.bincount(g[m] * w + (cols[m] - c0), weights=d[m], minlength=n_groups * w)
            out[:, c0:c1] += blk.reshape(n_groups, w)
            c0 = c1
        return out
    raise ValueError(enc)


def to_official(sums: np.ndarray, genes: list[str]) -> np.ndarray:
    """[n, n_var] dataset-axis sums -> [n, 18533] official-axis sums, NaN where the dataset lacks the gene."""
    out = np.full((sums.shape[0], len(OFFICIAL)), np.nan)
    src, dst = [], []
    seen = set()
    for j, g in enumerate(genes):
        i = OFF_POS.get(g)
        if i is not None and i not in seen:  # first occurrence of a duplicated symbol wins
            seen.add(i)
            src.append(j)
            dst.append(i)
    out[:, dst] = sums[:, src]
    return out


def save(name: str, rows: list[dict], counts: np.ndarray) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / f"{name}.npz", counts=counts, meta=json.dumps(rows))
    for r, c in zip(rows, counts):
        present = np.isfinite(c)
        print(f"  [{name}] {r['query']:<44} cells {r['n_cells']:>7}  umi/cell {np.nansum(c) / max(r['n_cells'], 1):>9.0f}"
              f"  official genes {int(present.sum()):>5}  true {r['true']}", flush=True)


def scperturb_lines(name: str, fname: str, line_key: str | None, ctrl: set[str] | None, chem: str,
                    truth: dict[str, str | None], fixed_line: str | None = None, pert_key: str = "perturbation",
                    extra_filter=None, min_cells: int = MIN_CELLS, prefer_var: tuple[str, ...] = ()) -> None:
    t0 = time.time()
    with h5py.File(SCP / fname if not fname.startswith("/") else fname, "r") as h:
        X = h["X"]
        n_obs = int(X.shape[0]) if isinstance(X, h5py.Dataset) else int(X.attrs["shape"][0])
        lines = obs_column(h, line_key) if line_key else np.asarray([fixed_line] * n_obs, dtype=object)
        sel = np.ones(n_obs, bool)
        if ctrl is not None:
            pert = obs_column(h, pert_key)
            sel &= np.isin(pert, list(ctrl))
        if extra_filter is not None:
            sel &= extra_filter(h)
        labs, inv = np.unique(lines, return_inverse=True)
        counts_per = np.bincount(inv[sel], minlength=len(labs))
        keep_lab = [i for i, l in enumerate(labs) if counts_per[i] >= min_cells and l in truth]
        gid = np.full(len(labs), -1)
        for k, i in enumerate(keep_lab):
            gid[i] = k
        groups = np.where(sel, gid[inv], -1)
        if not keep_lab:
            print(f"[{name}] no line with >= {min_cells} control cells", flush=True)
            return
        sums = group_sums(h, groups, len(keep_lab))
        genes = var_names(h, prefer_var)
    rows = [{"query": f"{name}:{labs[i]}", "dataset": name, "line": str(labs[i]), "chemistry": chem,
             "true": truth[labs[i]], "n_cells": int(counts_per[i])} for i in keep_lab]
    save(name, rows, to_official(sums, genes))
    print(f"[{name}] {len(rows)} rows in {time.time() - t0:.0f}s", flush=True)


CTRL = {"control", "Vehicle", "vehicle", "DMSO", "non-targeting", "NTC", "sgLACZ", "sgOR2J2", "safe-targeting"}


def mcfarland() -> None:
    clin = pd.read_csv(ROOT / "data/external/ccle_broad_2025/data_clinical_sample.txt", sep="\t", comment="#")
    m = clin.set_index("DEPMAP_ID")["SAMPLE_ID"].to_dict()
    with h5py.File(SCP / "McFarlandTsherniak2020.h5ad", "r") as h:
        ids = sorted(set(obs_column(h, "DepMap_ID")))
        q = obs_column(h, "cell_quality") if "cell_quality" in h["obs"] else None
        print("[mcfarland] cell_quality:", None if q is None else dict(zip(*np.unique(q, return_counts=True))), flush=True)
    truth = {i: m.get(i) for i in ids}

    def normal(h):
        if "cell_quality" not in h["obs"]:
            return np.ones(h["obs"]["DepMap_ID"]["codes"].shape[0], bool)
        return obs_column(h, "cell_quality") == "normal"

    scperturb_lines("mcfarland_ctrl", "McFarlandTsherniak2020.h5ad", "DepMap_ID", {"control", "sgLACZ", "sgOR2J2"},
                    "10x 3' v2 (MIX-seq)", truth, extra_filter=normal)
    scperturb_lines("mcfarland_all", "McFarlandTsherniak2020.h5ad", "DepMap_ID", None,
                    "10x 3' v2 (MIX-seq, all conditions)", truth, extra_filter=normal)


def others() -> None:
    A549, MCF7, K562 = "A549_LUNG", "MCF7_BREAST", "K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"
    THP1, JURKAT = "THP1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"
    jobs = [
        ("sciplex3", "SrivatsanTrapnell2020_sciplex3.h5ad", "cell_line", {"control"}, "sci-RNA-seq3",
         {"A549": A549, "MCF7": MCF7, "K562": K562}),
        ("sciplex4", "SrivatsanTrapnell2020_sciplex4.h5ad", "cell_line", {"control"}, "sci-RNA-seq3",
         {"A549": A549, "MCF7": MCF7}),
        ("sciplex2", "SrivatsanTrapnell2020_sciplex2.h5ad", "cell_line", None, "sci-RNA-seq3 (all conditions)",
         {"A549": A549}),
        ("lotfollahi", "LotfollahiTheis2023.h5ad", "cell_line", {"control"}, "sci-RNA-seq3", {"A549": A549}),
        ("sunshine", "SunshineHein2023.h5ad", "cell_line", {"control"}, "10x 3'", {"Calu-3": "CALU3_LUNG"}),
        ("papalexi", "PapalexiSatija2021_eccite_RNA.h5ad", "cell_line", {"control"}, "ECCITE (10x 5')", {"THP-1": THP1}),
        ("papalexi_arr", "PapalexiSatija2021_eccite_arrayed_RNA.h5ad", "cell_line", {"control"}, "ECCITE (10x 5')",
         {"THP-1": THP1}),
        ("wessels", "WesselsSatija2023.h5ad", "cell_line", None, "CaRPool-seq (10x 3', all conditions)", {"THP-1": THP1}),
        ("datlinger17", "DatlingerBock2017.h5ad", "cell_line", {"control"}, "CROP-seq", {"Jurkat cells": JURKAT}),
        ("datlinger21", "DatlingerBock2021.h5ad", "cell_line", {"control"}, "scifi-RNA-seq", {"Jurkat cells": JURKAT}),
        ("norman", "NormanWeissman2019_filtered.h5ad", "cell_line", {"control"}, "10x 3'", {"K562": K562}),
        ("gasperini_low", "GasperiniShendure2019_lowMOI.h5ad", "cell_line", {"control"}, "10x 3'", {"K562": K562}),
        ("schraivogel8", "SchraivogelSteinmetz2020_TAP_SCREEN__chromosome_8_screen.h5ad", "cell_line", {"control"},
         "TAP-seq (targeted)", {"K562": K562}),
        ("xie", "XieHon2017.h5ad", "cell_line", {"control"}, "10x 3'", {"K562": K562}),
        ("adamson5", "AdamsonWeissman2016_GSM2406677_10X005.h5ad", "cell_line", None, "10x 3' (all conditions)",
         {"K562": K562}),
        ("dixit", "DixitRegev2016.h5ad", "cell_line", None, "Perturb-seq (all conditions)", {"K562": K562}),
        ("xu", "XuCao2023.h5ad", "cell_line", {"control"}, "10x 3'", {"HEK293": None}),
        ("aissa", "AissaBenevolenskaya2021.h5ad", "cell_line", {"control"}, "10x 3' (xenograft)", {"PC9_xenograft": "PC9_LUNG"}),
        ("frangieh", "FrangiehIzar2021_RNA.h5ad", None, {"control"}, "Perturb-CITE (10x 3')", {"melanoma_patient": None}),
        ("zhao", "ZhaoSims2021.h5ad", None, {"control"}, "10x 3' (patient GBM slices)", {"gbm_patient": None}),
        ("tian_ipsc", "TianKampmann2019_iPSC.h5ad", None, {"control"}, "10x 3'", {"iPSC": None}),
        ("tian_neuron", "TianKampmann2019_day7neuron.h5ad", None, {"control"}, "10x 3'", {"iPSC_neuron": None}),
        ("shifrut", "ShifrutMarson2018.h5ad", None, {"control"}, "10x 3'", {"primary_T": None}),
        ("liang", "LiangWang2023.h5ad", None, {"control"}, "10x 3'", {"cholangiocyte_organoid": None}),
        ("gehring", "GehringPachter2019.h5ad", None, None, "10x 3' (BMP4)", {"NSC": None}),
    ]
    for name, fname, key, ctrl, chem, truth in jobs:
        fixed = next(iter(truth)) if key is None else None
        try:
            scperturb_lines(name, fname, key, ctrl, chem, truth, fixed_line=fixed)
        except Exception as e:  # one unreadable set must not stop the rest
            print(f"[{name}] FAILED: {type(e).__name__}: {e}", flush=True)
    # Cui: primary immune cell types, controls per type
    try:
        with h5py.File(SCP / "CuiHacohen2023.h5ad", "r") as h:
            types = sorted(set(obs_column(h, "celltype")))
        scperturb_lines("cui", "CuiHacohen2023.h5ad", "celltype", {"control"}, "10x 3' (primary immune)",
                        {t: None for t in types}, min_cells=30)
    except Exception as e:
        print(f"[cui] FAILED: {type(e).__name__}: {e}", flush=True)


def myllia() -> None:
    """Myllia A549 (10x 5' HT, protein-coding only): non-targeting + gene-desert guides, two channels."""
    A549 = "A549_LUNG"
    try:
        sums, n = None, 0
        for ch in (1, 2):
            f = ROOT / f"data/external/myllia-a549/GSE337804_A549_ch_{ch}.RDS.h5ad"
            with h5py.File(f, "r") as h:
                sg = obs_column(h, "sgrna_symbol")
                g = np.where(np.isin(sg, ["non-targeting", "gene_desert"]), 0, -1)
                s = to_official(group_sums(h, g, 1), var_names(h))
            # the two channels carry different gene axes: sum on the official axis (NaN only where both lack a gene)
            sums = s if sums is None else np.where(np.isnan(sums), 0, sums) + np.where(np.isnan(s), 0, s) \
                + np.where(np.isnan(sums) & np.isnan(s), np.nan, 0)
            n += int((g >= 0).sum())
        save("myllia", [{"query": "myllia:A549", "dataset": "myllia", "line": "A549", "chemistry": "10x 5' HT",
                         "true": A549, "n_cells": n}], sums)
    except Exception as e:
        print(f"[myllia] FAILED: {type(e).__name__}: {e}", flush=True)


if __name__ == "__main__":
    which = sys.argv[1:] or ["mcfarland", "others"]
    for w in which:
        {"mcfarland": mcfarland, "others": others, "myllia": myllia}[w]()
