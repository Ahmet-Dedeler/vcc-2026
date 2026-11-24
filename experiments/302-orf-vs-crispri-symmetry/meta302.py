"""302 — metadata pass (FREE under the R4 access rule: obs / var / uns only, `X` NEVER opened).

Establishes the cohort, the gene axis, the isoform-pooling arithmetic and every known-answer row
that does not need an expression value, so the preregistration can name real numbers before the
single access.

    PYTHONDONTWRITEBYTECODE=1 uv run python meta302.py

Writes work/meta302.json + work/cohort302.csv.  Read-only outside experiments/302-*/.
"""
from __future__ import annotations

import csv
import hashlib
import json
import resource
import sys
from collections import Counter
from pathlib import Path

import h5py
import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/302-orf-vs-crispri-symmetry"
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)

ATLAS = ROOT / "data/external/scperturb/JoungZhang2023_atlas.h5ad"
ORFMAP = ROOT / "experiments/301-metadata-acquisitions/work/joung_orf_map.csv"
H1BANK = ROOT / "experiments/115-gse295214-h1/work/compact-donor"
HELD279 = ROOT / "experiments/279-znf-state-dependence/work/held279.npz"
ZNF279 = ROOT / "experiments/279-znf-state-dependence/work/znf279.csv"
PERSRC290 = ROOT / "experiments/290-remainder-response-features/work/persrc290.npz"

MIN_CELLS_SECONDARY = 50
MIN_CELLS_PRIMARY = 100
LOG = []


def P(msg: str = "") -> None:
    print(msg, flush=True)
    LOG.append(msg)


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


def main() -> None:
    out: dict = {}

    # ---------------------------------------------------------------- the atlas, obs/var only
    with h5py.File(ATLAS, "r") as f:
        assert "X" in f
        shape = [int(x) for x in f["X"].attrs["shape"]]
        cats = np.asarray(f["obs/perturbation/categories"][:]).astype(str)
        codes = f["obs/perturbation/codes"][:]
        ncounts = f["obs/ncounts"][:]
        batch_cats = np.asarray(f["obs/batch/categories"][:]).astype(str)
        batch = f["obs/batch/codes"][:]
        var_sym = np.asarray(f["var/gene_symbol"][:]).astype(str)
    P(f"atlas {ATLAS.name}: X shape {shape}, {len(cats)} perturbation categories, "
      f"{len(var_sym)} var symbols, batches {list(batch_cats)}")
    out["atlas"] = dict(shape=shape, n_categories=len(cats), n_var=len(var_sym),
                        batches=list(batch_cats),
                        cells_per_batch=[int((batch == i).sum()) for i in range(len(batch_cats))])

    # categories are the raw code strings '0'..'3368'
    cat_int = np.asarray([int(c) for c in cats])
    assert np.array_equal(cat_int, np.arange(len(cats))), "categories are not 0..n-1"
    cells_per_code = np.bincount(codes.astype(np.int64), minlength=len(cats))

    dup = Counter(var_sym)
    n_dup_sym = sum(1 for s, c in dup.items() if c > 1)
    P(f"var symbols: {len(var_sym)} rows, {len(dup)} distinct, {n_dup_sym} symbols duplicated")
    out["var_duplicate_symbols"] = n_dup_sym

    # ---------------------------------------------------------------- 301's ORF id map
    code2sym, code2orf, code2iso = {}, {}, {}
    with open(ORFMAP) as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        sym = (r["gene_symbol"] or "").strip()
        if not sym:
            continue
        c = int(r["perturbation_code"])
        code2sym[c] = sym
        code2orf[c] = r["tforf_id"]
        code2iso[c] = r["isoform_id"]
    P(f"301 map: {len(rows)} rows, {len(code2sym)} codes with a symbol "
      f"(sha256 {sha256(ORFMAP)[:16]}…)")
    out["orf_map"] = dict(path=str(ORFMAP.relative_to(ROOT)), sha256=sha256(ORFMAP),
                          n_rows=len(rows), n_codes_with_symbol=len(code2sym))

    # KA-A: the two controls
    ctrl_codes = {}
    for c, s in code2sym.items():
        if s in ("GFP", "mCherry"):
            ctrl_codes[s] = c
    P(f"KA-A controls: GFP code {ctrl_codes.get('GFP')} n={cells_per_code[ctrl_codes['GFP']]} "
      f"(want 49,539); mCherry code {ctrl_codes.get('mCherry')} "
      f"n={cells_per_code[ctrl_codes['mCherry']]} (want 78,274)")
    ka_a = dict(gfp_code=int(ctrl_codes["GFP"]), gfp_cells=int(cells_per_code[ctrl_codes["GFP"]]),
                mcherry_code=int(ctrl_codes["mCherry"]),
                mcherry_cells=int(cells_per_code[ctrl_codes["mCherry"]]))
    ka_a["pass_"] = bool(ka_a["gfp_cells"] == 49539 and ka_a["mcherry_cells"] == 78274
                         and ka_a["gfp_code"] == 3367 and ka_a["mcherry_code"] == 3368)
    out["KA_A_controls"] = ka_a

    # ---------------------------------------------------------------- pooled per-symbol counts
    sym_cells: dict[str, int] = {}
    sym_codes: dict[str, list[int]] = {}
    for c, s in code2sym.items():
        if s in ("GFP", "mCherry"):
            continue
        sym_cells[s] = sym_cells.get(s, 0) + int(cells_per_code[c])
        sym_codes.setdefault(s, []).append(c)
    n_cov = sum(int(cells_per_code[c]) for c in code2sym)
    P(f"resolved symbols {len(sym_cells)} (want 1,647 incl. controls -> {len(sym_cells) + 2}); "
      f"cells covered {n_cov} / {len(codes)} = {n_cov / len(codes):.4%} (want 99.48 %)")
    out["symbols"] = dict(n_tf_symbols=len(sym_cells), cells_covered=n_cov,
                          frac_covered=float(n_cov / len(codes)),
                          n_ge50=int(sum(v >= 50 for v in sym_cells.values())))

    # ---------------------------------------------------------------- the H1 CRISPRi bank
    h1_perts = [x.strip() for x in (H1BANK / "perts.txt").read_text().splitlines() if x.strip()]
    h1_genes = [x.strip() for x in (H1BANK / "genes.txt").read_text().splitlines() if x.strip()]
    P(f"H1 bank gse295214_h1: {len(h1_perts)} perturbations, {len(h1_genes)} genes")
    out["h1_bank"] = dict(path=str(H1BANK.relative_to(ROOT)), n_perts=len(h1_perts),
                          n_genes=len(h1_genes))

    h1_set = set(h1_perts)
    overlap = sorted(s for s in sym_cells if s in h1_set)
    P(f"TF symbols in BOTH the ORF atlas (resolved) and the H1 bank: {len(overlap)} (want 249)")

    held = np.load(HELD279, allow_pickle=True)
    held_genes = np.asarray(held["genes"]).astype(str)
    held_set = set(held_genes.tolist())
    n_in_held = sum(1 for s in overlap if s in held_set)
    P(f"  of those, in 279's 17,106-gene held map: {n_in_held} (want 233)")

    coh50 = [s for s in overlap if sym_cells[s] >= MIN_CELLS_SECONDARY]
    coh100 = [s for s in overlap if sym_cells[s] >= MIN_CELLS_PRIMARY]
    P(f"  >= {MIN_CELLS_SECONDARY} ORF cells: {len(coh50)};  "
      f">= {MIN_CELLS_PRIMARY} ORF cells: {len(coh100)} (want 148)")
    ka_b = dict(n_overlap=len(overlap), n_in_held=n_in_held, n_ge50=len(coh50),
                n_ge100=len(coh100))
    ka_b["pass_"] = bool(len(overlap) == 249 and n_in_held == 233 and len(coh100) == 148)
    out["KA_B_cohort"] = ka_b

    # ---------------------------------------------------------------- the gene axis
    h1_gset = set(h1_genes)
    keep_sym = {s for s, c in dup.items() if c == 1}      # duplicated symbols dropped, stated rule
    axis = [g for g in h1_genes if g in keep_sym and g in set(var_sym.tolist())]
    P(f"axis = H1-bank genes ∩ atlas var (unique symbols only): {len(axis)} genes "
      f"(H1 {len(h1_genes)}, atlas {len(var_sym)}, duplicated atlas symbols dropped: {n_dup_sym})")
    out["axis"] = dict(n=len(axis), rule="H1 bank genes.txt order; atlas var symbols that are "
                                          "unique in var; intersection")

    # how many cohort TFs have their own column on the axis (it is REMOVED per the prereg)
    on_axis = sum(1 for s in coh100 if s in set(axis))
    P(f"  cohort TFs whose own gene is on the axis (column removed): {on_axis} / {len(coh100)}")
    out["axis"]["cohort_on_axis"] = on_axis

    # ---------------------------------------------------------------- 279 / 290 feature columns
    z = {}
    with open(ZNF279) as fh:
        for r in csv.DictReader(fh):
            z[r["gene"]] = r
    have_C = [s for s in coh100 if s in z and z[s]["C_held"] not in ("", "nan")]
    P(f"cohort TFs with a 279 C_held value: {len(have_C)} / {len(coh100)}")

    d290 = np.load(PERSRC290, allow_pickle=True)
    keys290 = np.asarray(d290["keys"]).astype(str).tolist()
    g290 = np.asarray(d290["genes"]).astype(str).tolist()
    qi = keys290.index("gse295214_h1")
    gpos = {g: i for i, g in enumerate(g290)}
    pc20 = d290["pc20"][qi]
    have_pc = [s for s in coh100 if s in gpos and np.isfinite(pc20[gpos[s]])]
    P(f"cohort TFs with a 290 gse295214_h1 gen_pc20 value: {len(have_pc)} / {len(coh100)}")
    out["features"] = dict(n_C_held=len(have_C), n_gen_pc20=len(have_pc),
                           src290_index=int(qi), n_sources290=len(keys290))

    # ---------------------------------------------------------------- isoforms + depth
    iso_counts = {s: len(sym_codes[s]) for s in coh100}
    iso_ge50 = {s: sum(1 for c in sym_codes[s] if cells_per_code[c] >= MIN_CELLS_SECONDARY)
                for s in coh100}
    n_multi = sum(1 for s in coh100 if iso_ge50[s] >= 2)
    cells100 = np.asarray([sym_cells[s] for s in coh100])
    P(f"isoforms per cohort TF: median {np.median(list(iso_counts.values())):.0f}, "
      f"max {max(iso_counts.values())}; TFs with >= 2 isoforms at >= 50 cells: {n_multi}")
    P(f"cells per cohort TF: median {np.median(cells100):.0f}, q10 {np.percentile(cells100, 10):.0f}, "
      f"q90 {np.percentile(cells100, 90):.0f}, min {cells100.min()}, max {cells100.max()}")

    # UMI depth per cell for the cohort's cells and the controls (obs['ncounts'], free)
    code_of_sym = {s: np.asarray(sym_codes[s]) for s in coh100}
    sel = np.zeros(len(codes), bool)
    for s in coh100:
        sel |= np.isin(codes, code_of_sym[s])
    ctrl_sel = np.isin(codes, [ka_a["gfp_code"], ka_a["mcherry_code"]])
    P(f"UMI/cell: cohort median {np.median(ncounts[sel]):.0f}, controls median "
      f"{np.median(ncounts[ctrl_sel]):.0f}; cohort cells {int(sel.sum())}, "
      f"control cells {int(ctrl_sel.sum())}")
    out["depth"] = dict(cohort_cells=int(sel.sum()), control_cells=int(ctrl_sel.sum()),
                        umi_median_cohort=float(np.median(ncounts[sel])),
                        umi_median_ctrl=float(np.median(ncounts[ctrl_sel])),
                        cells_per_tf_median=float(np.median(cells100)),
                        cells_per_tf_q10=float(np.percentile(cells100, 10)),
                        n_multi_isoform=int(n_multi))

    # batch composition (a declared covariate, never conditioned on)
    for i, b in enumerate(batch_cats):
        P(f"  batch {b}: cohort cells {int((batch[sel] == i).sum())}, "
          f"control cells {int((batch[ctrl_sel] == i).sum())}")

    # ---------------------------------------------------------------- cohort table
    with open(WORK / "cohort302.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["symbol", "orf_cells", "n_isoforms", "n_isoforms_ge50", "codes",
                    "in_held_map", "cohort"])
        for s in sorted(overlap, key=lambda x: -sym_cells[x]):
            w.writerow([s, sym_cells[s], len(sym_codes[s]),
                        sum(1 for c in sym_codes[s] if cells_per_code[c] >= 50),
                        " ".join(str(c) for c in sorted(sym_codes[s])),
                        int(s in held_set),
                        "primary" if sym_cells[s] >= MIN_CELLS_PRIMARY
                        else ("secondary" if sym_cells[s] >= MIN_CELLS_SECONDARY else "excluded")])

    out["cohort_primary"] = sorted(coh100)
    out["cohort_secondary"] = sorted(set(coh50) - set(coh100))
    out["rss_gib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30
    out["log"] = LOG
    (WORK / "meta302.json").write_text(json.dumps(out, indent=1))
    P(f"\npeak RSS {out['rss_gib']:.2f} GiB -> work/meta302.json, work/cohort302.csv")


if __name__ == "__main__":
    sys.exit(main())
