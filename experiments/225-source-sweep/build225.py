"""225 Phase 1b — materialise a donor-shaped table for every fetchable source in the census
that does not already have one.  Writes ONLY into `experiments/225-source-sweep/work/tables/`.

Three input shapes, three builders:

  scPerturb h5ad   `.X` raw counts, `obs['perturbation']`, control label `control`.
                   166's `build_extra.py` construction verbatim: TS = 5e4,
                   `b = log1p(TS * P / sum(P))`, `delta = b_pert - b_ctrl`, >=20 cells/group.
                   Handles CSC (gene-block walk) and CSR/dense (row-block walk).

  IGVF global-DE   a per-(element, gene) SIGNIFICANT-HIT table: `intended_target_name`,
                   `gene_names`, `log2fc`.  Sparse by construction (median 6 genes per panel
                   target) -- every unlisted coordinate is an exact 0.0, which is how every
                   other partial-coverage source in the registry stores an unmeasured
                   coordinate.  Values are scaled by ln2 so the units match CD4's
                   `log2FC x ln2` semantics (218's registry).

  Jiang 2025 zip   per-(target, stimulus) DE text files with one `log2FC_<LINE>` column per
                   cell line.  Averaged over the six lines and the five stimuli on the
                   coordinates each carries -- the 6-line consensus, which is the only thing
                   in Jiang that is not a single-line duplicate of something we hold.

CRISPRa sources (Tian 2021 CRISPRa, Norman 2019) are built AS MEASURED and a sign-flipped
twin is written next to them (`*_flip`), because activation is the opposite intervention to
the CRISPRi knockdown every recipient fold measures.  Which one (if either) helps is then a
measurement, not an assumption.

    uv run python build225.py [--only NAME ...]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import zipfile
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "work/tables"
OUT.mkdir(parents=True, exist_ok=True)
TS = 5e4
LN2 = float(np.log(2.0))

PANEL = [line.split(",")[0]
         for line in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:]
         if line.strip()]
PANEL_SET = set(PANEL)

SCP = "data/external/scperturb"
H5AD = {
    # name                      file                              stratum
    "tian_ipsc_neuron": (f"{SCP}/TianKampmann2021_CRISPRi.h5ad", None),
    "tian_ipsc_neuron_crispra": (f"{SCP}/TianKampmann2021_CRISPRa.h5ad", None),
    "xucao_hek293": (f"{SCP}/XuCao2023.h5ad", None),
    "sunshine_calu3": (f"{SCP}/SunshineHein2023.h5ad", None),
    "norman_k562": (f"{SCP}/NormanWeissman2019_filtered.h5ad", None),
    "wessels_thp1": (f"{SCP}/WesselsSatija2023.h5ad", None),
    "gasperini_k562": (f"{SCP}/GasperiniShendure2019_lowMOI.h5ad", None),
    "tian2019_ipsc": (f"{SCP}/TianKampmann2019_iPSC.h5ad", None),
}
FLIP_TWINS = ("tian_ipsc_neuron_crispra", "norman_k562")


# ------------------------------------------------------------------ shared table writer
def write_table(name: str, delta: np.ndarray, perts: list[str], genes: list[str],
                extra: dict | None = None, b_ctrl: np.ndarray | None = None) -> dict:
    d = OUT / name
    d.mkdir(parents=True, exist_ok=True)
    delta = np.ascontiguousarray(np.nan_to_num(delta, nan=0.0, posinf=0.0, neginf=0.0),
                                 dtype=np.float32)
    assert delta.shape == (len(perts), len(genes)), (delta.shape, len(perts), len(genes))
    # 218's `Bank.load` reads these back with `splitlines()` and asserts the shape, so a
    # symbol carrying embedded whitespace silently changes the row/column count.  Sanitise
    # and round-trip-assert here rather than discovering it three scripts downstream.
    def _clean(xs):
        return ["_".join(str(x).split()) or "UNNAMED" for x in xs]

    perts, genes = _clean(perts), _clean(genes)
    np.save(d / "delta.npy", delta)
    (d / "perts.txt").write_text("\n".join(perts) + "\n")
    (d / "genes.txt").write_text("\n".join(genes) + "\n")
    for f, n in (("perts.txt", delta.shape[0]), ("genes.txt", delta.shape[1])):
        back = [x for x in (d / f).read_text().split("\n") if x.strip()]
        assert len(back) == n, f"{name}/{f}: {len(back)} lines != {n} in delta"
    np.save(d / "b_ctrl.npy",
            np.zeros(len(genes), np.float32) if b_ctrl is None else b_ctrl.astype(np.float32))
    row_nz = (delta != 0).sum(1)
    norms = np.linalg.norm(delta.astype(np.float64), axis=1)
    panel_rows = [p for p in perts if p in PANEL_SET]
    pidx = [i for i, p in enumerate(perts) if p in PANEL_SET]
    meta = dict(name=name, shape=list(delta.shape), panel_targets=len(panel_rows),
                nonzero=int(np.count_nonzero(delta)),
                median_nonzero_genes_per_row=int(np.median(row_nz)),
                median_nonzero_genes_per_PANEL_row=(int(np.median(row_nz[pidx]))
                                                    if pidx else 0),
                rows_all_zero=int((row_nz == 0).sum()),
                effect_norm_median=float(np.median(norms)),
                effect_norm_median_panel=(float(np.median(norms[pidx])) if pidx else 0.0),
                **(extra or {}))
    (d / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"[{name}] {delta.shape} panel {meta['panel_targets']}/300  "
          f"nz genes/panel-row median {meta['median_nonzero_genes_per_PANEL_row']}  "
          f"norm(panel) {meta['effect_norm_median_panel']:.3f}", flush=True)
    return meta


# ----------------------------------------------------------------------- scPerturb h5ad
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


def var_names(h):
    v = h["var"]
    key = v.attrs.get("_index", "_index")
    return _dec(v[key.decode() if isinstance(key, bytes) else key][:])


def group_sums(h, codes, n_groups, n_genes, block=2000):
    """Group-sum `.X` regardless of storage: CSC walks gene blocks, CSR/dense walk rows."""
    X = h["X"]
    keep = codes >= 0
    ind = sp.csr_matrix((np.ones(int(keep.sum()), np.float32),
                         (codes[keep], np.flatnonzero(keep))),
                        shape=(n_groups, len(codes)))
    sums = np.zeros((n_groups, n_genes), np.float64)
    if isinstance(X, h5py.Group) and X.attrs.get("encoding-type") == "csc_matrix":
        indptr = X["indptr"][:]
        for g0 in range(0, n_genes, block):
            g1 = min(g0 + block, n_genes)
            a, b = int(indptr[g0]), int(indptr[g1])
            blk = sp.csc_matrix((X["data"][a:b], X["indices"][a:b], indptr[g0:g1 + 1] - a),
                                shape=(len(codes), g1 - g0))
            sums[:, g0:g1] = (ind @ blk).toarray()
        return sums
    if isinstance(X, h5py.Group):                        # csr_matrix
        indptr = X["indptr"][:]
        for r0 in range(0, len(codes), 20000):
            r1 = min(r0 + 20000, len(codes))
            a, b = int(indptr[r0]), int(indptr[r1])
            blk = sp.csr_matrix((X["data"][a:b], X["indices"][a:b], indptr[r0:r1 + 1] - a),
                                shape=(r1 - r0, n_genes))
            sums += (ind[:, r0:r1] @ blk).toarray()
        return sums
    for r0 in range(0, len(codes), 20000):               # dense
        r1 = min(r0 + 20000, len(codes))
        sums += ind[:, r0:r1] @ np.asarray(X[r0:r1], np.float64)
    return sums


def build_h5ad(name: str, rel: str, stratum, min_cells: int = 20) -> dict:
    path = ROOT / rel
    with h5py.File(path, "r") as h:
        labels = read_col(h, "perturbation")
        genes = var_names(h)
        n_genes, n_cells = len(genes), len(labels)
        in_s = np.ones(n_cells, bool)
        if stratum:
            in_s = read_col(h, stratum[0]) == stratum[1]
        lab = np.where(in_s, labels, "")
        uniq, counts = np.unique(lab[in_s], return_counts=True)
        cnt = dict(zip(uniq.tolist(), counts.tolist()))
        assert "control" in cnt, f"{name}: no control cells, saw {list(cnt)[:8]}"
        keep = [n for n in uniq.tolist() if n == "control" or cnt[n] >= min_cells]
        index = {n: i for i, n in enumerate(keep)}
        codes = np.fromiter((index.get(x, -1) for x in lab), np.int64, count=n_cells)
        sums = group_sums(h, codes, len(keep), n_genes)

    uniq_g, inv = np.unique(genes.astype(str), return_inverse=True)
    if len(uniq_g) != len(genes):
        mp = sp.csr_matrix((np.ones(len(genes)), (np.arange(len(genes)), inv)),
                           shape=(len(genes), len(uniq_g)))
        sums = np.asarray(sums @ mp)
        genes, n_genes = uniq_g.astype(object), len(uniq_g)

    ci = index["control"]
    Pc = sums[ci]
    b_ctrl = np.log1p(TS * Pc / Pc.sum()).astype(np.float32)
    tot = sums.sum(1, keepdims=True)
    good = tot.ravel() > 0
    good[ci] = False
    prof = np.zeros_like(sums, np.float32)
    prof[good] = np.log1p(TS * sums[good] / tot[good]).astype(np.float32)
    delta = prof[good] - b_ctrl
    names = [n for i, n in enumerate(keep) if good[i]]
    return write_table(name, delta, names, list(map(str, genes)), b_ctrl=b_ctrl,
                       extra=dict(source=rel, builder="scperturb-h5ad", target_sum=TS,
                                  min_cells=min_cells, n_control_cells=int(cnt["control"]),
                                  median_cells_per_pert=float(np.median(
                                      [cnt[n] for n in names]))))


# ------------------------------------------------------------------- IGVF global-DE csv
def build_igvf_cardiomyocyte() -> dict:
    """IGVFFI0830FXFI — WTC11-derived cardiac muscle cell, CRISPRi Perturb-seq (Gary Hon).

    A SIGNIFICANT-HIT table: only (element, gene) pairs that passed the deposit's own
    significance filter appear.  Every other coordinate is an exact 0.0.  081 measured this
    same asset through a retrieval-gallery adapter in Aug 2026 and closed it; it has never
    been through the post-whitening residual arm, which is the 208 rider rule.
    """
    src = ROOT / "data/external/igvf-cardiomyocyte/IGVFFI0830FXFI.csv.gz"
    rows = defaultdict(dict)
    genes: set[str] = set()
    with gzip.open(src, "rt") as f:
        for r in csv.DictReader(f):
            t, g = r["intended_target_name"].strip(), r["gene_names"]
            # the deposit carries a blank `intended_target_name` for its non-targeting and
            # region-only elements; those are not gene perturbations and must not be rows
            if not t or t.upper() in {"NA", "NAN", "NONE", "NON-TARGETING"} or not g:
                continue
            try:
                v = float(r["log2fc"])
            except (TypeError, ValueError):
                continue
            if not np.isfinite(v):
                continue
            rows[t][g] = v
            genes.add(g)
    genes_l = sorted(genes)
    gpos = {g: i for i, g in enumerate(genes_l)}
    perts = sorted(rows)
    delta = np.zeros((len(perts), len(genes_l)), np.float64)
    for i, t in enumerate(perts):
        for g, v in rows[t].items():
            delta[i, gpos[g]] = v * LN2          # log2FC -> natural-log units, as CD4
    # the panel's own response coordinates are zeroed by every source in the registry
    pcols = np.asarray([gpos[g] for g in PANEL if g in gpos], dtype=np.int64)
    if pcols.size:
        delta[:, pcols] = 0.0
    return write_table("igvf_cardiomyocyte", delta, perts, genes_l,
                       extra=dict(source=str(src.relative_to(ROOT)),
                                  builder="igvf-global-de",
                                  accession="IGVFFI0830FXFI",
                                  analysis_set="IGVFDS6332VCTO",
                                  semantics="published log2FC x ln2, significant hits only",
                                  panel_response_coords_zeroed=int(pcols.size)))


# -------------------------------------------------------------------------- Jiang 2025
def build_jiang() -> dict:
    """Jiang 2025 — 6 cancer lines x 5 stimuli, per-target DE text files.

    Averaged over every (line, stimulus) column that carries a coordinate: the 6-line
    CONSENSUS.  Per-line arms would each duplicate a line we already hold (K562) or add a
    single line at 9/300, so the consensus is the only construction worth a table.
    """
    z = zipfile.ZipFile(ROOT / "data/external/jiang2025/DE_results_all_pathway.zip")
    names = [n for n in z.namelist()
             if n.endswith("_DE_results.txt") and not n.startswith("__")]
    # match on the BASENAME only — the archive nests <root>/<stimulus>/<TARGET>_<STIM>_...
    pat = re.compile(r"^(.+?)_([A-Z0-9]+)_pathway_DE_results\.txt$")
    tot: dict[str, dict[str, list]] = defaultdict(lambda: defaultdict(list))
    genes: set[str] = set()
    for n in names:
        m = pat.match(n.rsplit("/", 1)[-1])
        if not m:
            continue
        target = m.group(1)
        with z.open(n) as f:
            # R `write.table` default: single-space separated, not tab
            rd = csv.DictReader((ln.decode() for ln in f), delimiter=" ")
            lcols = [c for c in (rd.fieldnames or []) if c.startswith("log2FC_")]
            assert lcols, f"{n}: no log2FC_<LINE> columns in {rd.fieldnames}"
            for r in rd:
                g = r["gene_ID"]
                vs = []
                for c in lcols:
                    try:
                        v = float(r[c])
                    except (TypeError, ValueError):
                        continue
                    if np.isfinite(v):
                        vs.append(v)
                if vs:
                    tot[target][g].extend(vs)
                    genes.add(g)
    genes_l = sorted(genes)
    gpos = {g: i for i, g in enumerate(genes_l)}
    perts = sorted(tot)
    delta = np.zeros((len(perts), len(genes_l)), np.float64)
    for i, t in enumerate(perts):
        for g, vs in tot[t].items():
            delta[i, gpos[g]] = float(np.mean(vs)) * LN2
    pcols = np.asarray([gpos[g] for g in PANEL if g in gpos], dtype=np.int64)
    if pcols.size:
        delta[:, pcols] = 0.0
    return write_table("jiang_6line", delta, perts, genes_l,
                       extra=dict(source="data/external/jiang2025/DE_results_all_pathway.zip",
                                  builder="jiang-de-zip", n_files=len(names),
                                  semantics="mean log2FC over 6 lines x 5 stimuli, x ln2",
                                  panel_response_coords_zeroed=int(pcols.size)))


def flip(name: str) -> dict:
    """Sign-flipped twin of a CRISPRa table — activation is the opposite intervention."""
    s = OUT / name
    d = OUT / f"{name}_flip"
    d.mkdir(parents=True, exist_ok=True)
    np.save(d / "delta.npy", -np.load(s / "delta.npy"))
    for f in ("perts.txt", "genes.txt"):
        (d / f).write_text((s / f).read_text())
    np.save(d / "b_ctrl.npy", np.load(s / "b_ctrl.npy"))
    m = json.loads((s / "meta.json").read_text())
    m["name"] = f"{name}_flip"
    m["sign_flipped_twin_of"] = name
    (d / "meta.json").write_text(json.dumps(m, indent=2))
    print(f"[{name}_flip] sign-flipped twin written")
    return m


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args()
    todo = args.only
    out = {}
    if todo is None or "igvf_cardiomyocyte" in todo:
        out["igvf_cardiomyocyte"] = build_igvf_cardiomyocyte()
    if todo is None or "jiang_6line" in todo:
        out["jiang_6line"] = build_jiang()
    for name, (rel, stratum) in H5AD.items():
        if todo is not None and name not in todo:
            continue
        if not (ROOT / rel).exists():
            print(f"[{name}] SKIP — {rel} absent")
            continue
        out[name] = build_h5ad(name, rel, stratum)
        if name in FLIP_TWINS:
            out[f"{name}_flip"] = flip(name)
    p = HERE / "work/build225.json"
    prev = json.loads(p.read_text()) if p.exists() else {}
    prev.update(out)
    p.write_text(json.dumps(prev, indent=2))
    print(f"-> {p}")


if __name__ == "__main__":
    main()
