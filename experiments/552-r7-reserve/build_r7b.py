"""552 -- build R7b's cell-level surface (RESERVE.md section 3 + amendment 1): CONSTRUCTION ONLY, no statistic of any kind.

Source: data/external/scperturb/SunshineHein2023.h5ad (sha256 asserted). Streams the CSC matrix by gene chunks (h5py
slices; peak ~3 GB) and writes, for the selected cells only, raw integer counts on the official 18,533-gene axis:

  cells    single-target labels (the obs `perturbation` holds one symbol: 1-2 guides of one target) of the ELIGIBLE targets,
           and every `control` cell (the non-targeting guides). Unassigned (nperts 0) and joined labels ("A_B") are dropped.
  targets  >= 100 such cells; carried by >= 2 of k562_gwps / orion_hct116 / orion_hek293t (data/donors/<name>/perts.txt);
           not in data/donors/k562_essential/perts.txt; on the official gene axis (data/challenge/gene_names.csv).
  genes    the official axis, in gene_names.csv order, mapped by symbol (var `gene_symbol`); an official gene with no source
           column is all-zero; a symbol present twice in the source keeps its first column (both recorded).
  no downsampling of cells; per-cell depth as measured.

Writes experiments/552-r7-reserve/work/r7b/ (gitignored matrices, tracked record):
  truth_r7b.h5ad        every selected cell (targets + controls); obs: target_gene ("non-targeting" for controls), guide_id,
                        gem_group, ncounts_src; var index = the official axis. CSR, integer-valued float32.
  bundle/context_R.h5ad the control cells only (obs: target_gene = non-targeting, context = R, ntc_id = guide_id) -- the
                        prediction side's only view of the line; pert_counts.csv (target_gene), gene_names.csv, manifest.json
  build_r7b.json        counts, target list, the gene-map audit, sha256 of every output. NO per-gene or per-target
                        expression statistic is computed or written (RESERVE.md section 3: sealed until the access).

    uv run python experiments/552-r7-reserve/build_r7b.py

Amendment 3 (RESERVE.md section 7): `--floor F --exclude-targets T1,T2 --min-cells-after M --out <dir name>` keeps a cell
(control or target) iff its raw counts over the 18,533 official-gene columns are >= F, drops the excluded targets' cells,
then drops targets left with < M cells; cell names are those of the unfloored selection. `--controls-must-equal <h5ad>`
asserts the floored controls equal that file cell for cell and count for count (553's prediction-side control file).
Target cells' library sizes are used for the threshold only; only per-target cell COUNTS are recorded. Defaults
(F 0, no exclusions, M 0, out r7b) are the original build.
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard552  # noqa: E402,F401
import collections  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import shutil  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import argparse  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SRC = ROOT / "data/external/scperturb/SunshineHein2023.h5ad"
SRC_SHA = "4dd14fb68dc16057e1cf4a05d0bcd2efd26c1a9e4ef8c7d21afccaf2cfc84270"
OUT = HERE / "work" / "r7b"
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
MIN_CELLS = 100
CHUNK = 512                                   # genes per streamed slice


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def strings(node) -> np.ndarray:
    if isinstance(node, h5py.Group):
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in node["categories"][:]], dtype=object)
        codes = node["codes"][:]
        out = np.asarray(["NA"] * len(codes), dtype=object)
        ok = codes >= 0
        out[ok] = cats[codes[ok]]
        return out
    v = node[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v], dtype=object) if v.dtype.kind in "SO" else v


def perts(name: str) -> set[str]:
    return {ln.strip() for ln in (ROOT / f"data/donors/{name}/perts.txt").read_text().splitlines() if ln.strip()}


def write_h5ad(path: Path, X: sp.csr_matrix, obs: pd.DataFrame, genes: list[str]) -> None:
    import anndata as ad
    a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name=None)))
    a.write_h5ad(path)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--floor", type=float, default=0.0)
    ap.add_argument("--exclude-targets", default="")
    ap.add_argument("--min-cells-after", type=int, default=0)
    ap.add_argument("--out", default="r7b")
    ap.add_argument("--controls-must-equal", default=None)
    args = ap.parse_args()
    out_dir = HERE / "work" / args.out
    exclude = [t for t in args.exclude_targets.split(",") if t]
    t0 = time.time()
    assert sha256(SRC) == SRC_SHA, "source sha256 differs from the declaration"
    official = [str(g) for g in pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"]]
    assert len(official) == 18533 and len(set(official)) == 18533
    with h5py.File(SRC, "r") as h:
        pert = strings(h["obs"]["perturbation"])
        guide = strings(h["obs"]["guide_id"])
        gem = strings(h["obs"]["gem_group"])
        ncounts = strings(h["obs"]["ncounts"])
        sym = strings(h["var"]["gene_symbol"])
        shape = tuple(int(x) for x in h["X"].attrs["shape"])
        assert h["X"].attrs["encoding-type"] == "csc_matrix"
    n_cells, n_genes = shape
    # ---------------------------------------------------------------- the eligible targets (metadata only)
    ess = perts("k562_essential")
    cover = [perts(d) for d in DONORS]
    offset = set(official)
    cnt = collections.Counter(pert)
    single = {k: v for k, v in cnt.items() if k not in ("control", "NA") and "_" not in k and "+" not in k}
    audit = {}
    eligible = []
    for t, v in sorted(single.items()):
        why = []
        if v < MIN_CELLS:
            why.append(f"cells {v} < {MIN_CELLS}")
        if sum(t in c for c in cover) < 2:
            why.append("carried by < 2 production donors")
        if t in ess:
            why.append("K562-essential")
        if t not in offset:
            why.append("off the official axis")
        audit[t] = dict(cells=int(v), excluded=why)
        if not why:
            eligible.append(t)
    eligible = sorted(eligible)
    keep_cell = np.isin(pert, eligible) | (pert == "control")
    rows = np.flatnonzero(keep_cell)
    print(f"{len(single)} single-target labels; eligible {len(eligible)}; cells kept {len(rows)} "
          f"({int((pert[rows] == 'control').sum())} controls)", flush=True)
    # ---------------------------------------------------------------- the gene map (symbol -> first source column)
    first = {}
    dup = collections.Counter(sym)
    for j, s in enumerate(sym):
        first.setdefault(s, j)
    src_col = np.asarray([first.get(g, -1) for g in official], dtype=np.int64)
    mapped = int((src_col >= 0).sum())
    print(f"official genes mapped {mapped} / 18,533; source duplicate symbols {sum(1 for v in dup.values() if v > 1)}", flush=True)
    # ---------------------------------------------------------------- stream the kept cells' counts (CSC by gene chunks)
    new_row = -np.ones(n_cells, dtype=np.int64)
    new_row[rows] = np.arange(len(rows))
    want_cols = np.flatnonzero(src_col >= 0)
    col_to_official = {int(src_col[i]): int(i) for i in want_cols}
    src_cols_sorted = np.asarray(sorted(col_to_official))
    R_, C_, V_ = [], [], []
    with h5py.File(SRC, "r") as h:
        indptr = h["X"]["indptr"][:]
        D, I = h["X"]["data"], h["X"]["indices"]
        for s0 in range(0, len(src_cols_sorted), CHUNK):
            cols = src_cols_sorted[s0:s0 + CHUNK]
            lo, hi = int(indptr[cols[0]]), int(indptr[cols[-1] + 1])
            dd, ii = D[lo:hi], I[lo:hi]
            for c in cols:
                a, b = int(indptr[c]) - lo, int(indptr[c + 1]) - lo
                r = new_row[ii[a:b]]
                m = r >= 0
                if m.any():
                    R_.append(r[m]); C_.append(np.full(int(m.sum()), col_to_official[int(c)], np.int32)); V_.append(dd[a:b][m])
    R_ = np.concatenate(R_); C_ = np.concatenate(C_); V_ = np.concatenate(V_)
    assert np.all(V_ >= 0) and np.all(V_ == np.floor(V_)), "non-integer counts in the source"
    X = sp.csr_matrix((V_.astype(np.float32), (R_, C_)), shape=(len(rows), len(official)))
    X.sum_duplicates()
    X.eliminate_zeros()
    del R_, C_, V_
    assert X.max() <= 1e6
    labels = np.where(pert[rows] == "control", "non-targeting", pert[rows]).astype(object)
    names = np.asarray([f"r7b_{i}" for i in range(len(rows))], dtype=object)
    floor_rec = None
    if args.floor > 0 or exclude or args.min_cells_after > 0:          # amendment 3 (RESERVE.md section 7)
        keep = np.ones(len(rows), bool)
        n_before = {t: int((labels == t).sum()) for t in eligible}
        n_ctrl_before = int((labels == "non-targeting").sum())
        if args.floor > 0:
            lib = np.asarray(X.sum(1)).ravel()                         # official-axis library: used for the threshold only
            keep &= lib >= args.floor
            ctl_lib = lib[(labels == "non-targeting") & keep]          # controls' quantiles are allowed (RESERVE section 6)
            del lib
        assert set(exclude) <= set(eligible), f"excluded targets not in the panel: {set(exclude) - set(eligible)}"
        keep &= ~np.isin(labels, exclude)
        after = collections.Counter(labels[keep])
        dropped = sorted(t for t in eligible if t not in exclude and after.get(t, 0) < args.min_cells_after)
        keep &= ~np.isin(labels, dropped)
        eligible = sorted(t for t in eligible if t not in exclude and t not in dropped)
        X, labels, names, rows = X[np.flatnonzero(keep)], labels[keep], names[keep], rows[keep]
        floor_rec = dict(floor=args.floor, excluded=exclude, min_cells_after=args.min_cells_after, dropped_after_floor=dropped,
                         controls_before=n_ctrl_before, controls_after=int((labels == "non-targeting").sum()),
                         cells_per_target_before=n_before,
                         cells_per_target_after={t: int((labels == t).sum()) for t in eligible},
                         control_library_quantiles_after=[float(v) for v in np.quantile(ctl_lib, [0, .1, .25, .5, .75, .9, 1])]
                         if args.floor > 0 else None,
                         rule="keep a cell iff its raw counts over the 18,533 official-gene columns >= floor (controls and "
                              "targets alike); excluded targets' cells dropped; then targets with < min_cells_after cells "
                              "dropped (RESERVE.md section 7)")
        print(f"floor {args.floor:g}: controls {n_ctrl_before} -> {floor_rec['controls_after']}; targets {len(n_before)} -> "
              f"{len(eligible)} (excluded {exclude}, dropped {dropped}); cells {int(keep.size)} -> {int(keep.sum())}", flush=True)
    obs = pd.DataFrame({"target_gene": pd.Categorical(labels), "guide_id": guide[rows].astype(str),
                        "gem_group": gem[rows].astype(str), "ncounts_src": np.asarray(ncounts[rows], float)},
                       index=pd.Index(list(names)))
    OUT = out_dir
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "bundle").mkdir(exist_ok=True)
    write_h5ad(OUT / "truth_r7b.h5ad", X, obs, official)
    ctl = np.flatnonzero(labels == "non-targeting")
    cobs = pd.DataFrame({"target_gene": pd.Categorical(["non-targeting"] * len(ctl)), "context": pd.Categorical(["R"] * len(ctl)),
                         "ntc_id": guide[rows][ctl].astype(str)}, index=obs.index[ctl])
    write_h5ad(OUT / "bundle" / "context_R.h5ad", X[ctl], cobs, official)
    if args.controls_must_equal:                                      # the prediction side's floored control file
        import anndata as ad
        ref = ad.read_h5ad(args.controls_must_equal)
        assert list(ref.obs_names) == list(obs.index[ctl]), "floored controls: cell names differ from the prediction side's"
        assert list(ref.var_names) == official, "floored controls: gene order differs"
        R_ = sp.csr_matrix(ref.X)
        d = (R_ != X[ctl]).nnz
        assert d == 0, f"floored controls: {d} count entries differ from the prediction side's"
        print(f"floored controls == {args.controls_must_equal}: {len(ctl)} cells, every count equal", flush=True)
    pd.DataFrame({"target_gene": eligible}).to_csv(OUT / "bundle" / "pert_counts.csv", index=False)
    shutil.copyfile(ROOT / "data/challenge/gene_names.csv", OUT / "bundle" / "gene_names.csv")
    (OUT / "bundle" / "manifest.json").write_text(json.dumps(dict(partition="reserve-r7b", panel_id=f"r7b-sunshine-calu3{'' if args.out == 'r7b' else '-' + args.out}",
                                                                 contexts=["R"], source=str(SRC.relative_to(ROOT))), indent=2) + "\n")
    rec = dict(written=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), source=str(SRC.relative_to(ROOT)), source_sha256=SRC_SHA,
               rules=dict(min_cells=MIN_CELLS, donors=DONORS, min_donors=2, essential="data/donors/k562_essential/perts.txt",
                          gene_map="official symbol -> first source column with that gene_symbol"),
               n_single_labels=len(single), n_eligible=len(eligible), targets=eligible,
               cells_per_target={t: int(single[t]) for t in eligible}, n_controls=int(len(ctl)), n_cells=int(len(labels)),
               floor=floor_rec, controls_checked_against=args.controls_must_equal,
               genes=dict(official=18533, mapped=mapped, unmapped=18533 - mapped), target_audit=audit,
               outputs={p: sha256(OUT / p) for p in ("truth_r7b.h5ad", "bundle/context_R.h5ad", "bundle/pert_counts.csv",
                                                     "bundle/gene_names.csv", "bundle/manifest.json")},
               note="construction only: no per-gene or per-target statistic was computed (RESERVE.md section 3)",
               seconds=round(time.time() - t0, 1))
    (OUT / "build_r7b.json").write_text(json.dumps(rec, indent=1))
    print(f"-> {OUT.relative_to(ROOT)}: {len(eligible)} targets, {len(rows)} cells, {mapped} mapped genes ({rec['seconds']} s)")


if __name__ == "__main__":
    main()
