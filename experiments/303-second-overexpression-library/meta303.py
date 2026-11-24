"""303 — the metadata pass.  FREE: obs / var / uns / measure_mask only.

Reads no `X` and no `layers/control` for any dataset.  Writes
`work/meta303.json`, `work/R5-ROWS.tsv`, `work/AXIS.txt`, `work/COHORT.tsv`,
`work/INPUT-MANIFEST.json`.

    PYTHONDONTWRITEBYTECODE=1 uv run python meta303.py
"""
from __future__ import annotations

import csv
import json
import resource
import sys
import time
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib303 as L  # noqa: E402


def read_symmetry302() -> dict[str, dict]:
    out = {}
    if not L.SYM302.exists():
        return out
    with open(L.SYM302) as fh:
        for r in csv.DictReader(fh):
            out[r.get("tf") or r.get("symbol") or r.get("gene")] = r
    return out


def main() -> None:
    t0 = time.time()
    L.WORK.mkdir(exist_ok=True)
    man: dict = {}

    # ---------------------------------------------------------------- source
    size = L.DB.stat().st_size
    sha_path = L.WORK / "source.sha256"
    sha = sha_path.read_text().split()[0] if sha_path.exists() else L.sha256_file(L.DB)
    man["source"] = dict(path=str(L.DB.relative_to(L.ROOT)), size_bytes=size, sha256=sha,
                         matches_115_meta=bool(sha == L.DB_SHA and size == L.DB_BYTES))
    print(f"source sha256 {sha} ({size} B) — matches 115's pin: {man['source']['matches_115_meta']}")

    m = L.read_db_meta()
    canon = L.canonical_genes()
    bank = L.read_bank()
    proj = L.build_projection(m.var_symbol, canon)

    # ------------------------------------------------------------ uns is empty
    with h5py.File(L.DB, "r") as f:
        uns_keys = sorted(f["uns"].keys())
        obsm_keys = sorted(f["obsm"].keys())
        x_dtype = str(f["X"].dtype)
        x_chunks = f["X"].chunks
        layers = {k: dict(shape=list(f[f"layers/{k}"].shape), dtype=str(f[f"layers/{k}"].dtype),
                          chunks=(list(f[f"layers/{k}"].chunks) if f[f"layers/{k}"].chunks else None))
                  for k in sorted(f["layers"].keys())}
    man["file"] = dict(n_obs=m.n_obs, n_var=m.n_var, X_dtype=x_dtype,
                       X_chunks=(list(x_chunks) if x_chunks else None), layers=layers,
                       uns_keys=uns_keys, obsm_keys=obsm_keys,
                       obs_columns=["condition", "perturbed_gene", "perturbation_sign",
                                    "perturbed_gene_name", "dataset"],
                       var_index="human_id", var_columns=["gene_name"])
    print(f"file {m.n_obs} x {m.n_var}; uns keys {uns_keys}; layers {sorted(layers)}")

    # --------------------------------------------------------- dataset census
    ds, counts = np.unique(m.dataset, return_counts=True)
    man["datasets"] = {str(a): int(b) for a, b in zip(ds, counts)}

    r5 = np.flatnonzero(m.dataset == L.RESERVE_DATASET).astype(np.int64)
    bnk = np.flatnonzero(m.dataset == L.BANK_DATASET).astype(np.int64)
    jng = np.flatnonzero(m.dataset == L.JOUNG_DATASET).astype(np.int64)
    r5.sort(); bnk.sort(); jng.sort()

    # -------------------------------------------------------- R5 row identity
    r5_sym = m.pert_symbol[r5]
    r5_ens = m.pert_gene[r5]
    r5_id = m.obs_id[r5]
    r5_sign = m.sign[r5]
    r5_cond = np.unique(m.condition[r5]).tolist()
    dup_sym = sorted({s for s in r5_sym if int((r5_sym == s).sum()) > 1})
    id_ok = all(i == f"hESC_TF_screen_{s}" for i, s in zip(r5_id, r5_sym))
    man["R5"] = dict(
        dataset=L.RESERVE_DATASET, n_rows=int(r5.size),
        contiguous=bool(np.array_equal(r5, np.arange(r5[0], r5[0] + r5.size))),
        first_row=int(r5[0]), last_row=int(r5[-1]),
        row_index_sha256=L.sha256_array(r5),
        symbols_sha256=L.sha256_strings(r5_sym.tolist()),
        obs_id_sha256=L.sha256_strings(r5_id.tolist()),
        duplicate_symbols=dup_sym, n_unique_symbols=int(np.unique(r5_sym).size),
        signs=sorted({int(x) for x in r5_sign}),
        conditions=r5_cond,
        obs_id_pattern_ok=bool(id_ok),
        ensembl_sha256=L.sha256_strings(r5_ens.tolist()),
    )
    print(f"R5 rows {r5.size} [{r5[0]}..{r5[-1]}] contiguous={man['R5']['contiguous']} "
          f"unique symbols {man['R5']['n_unique_symbols']} signs {man['R5']['signs']}")

    # ------------------------------------------------------------- cohorts
    r5_set = set(r5_sym.tolist())
    bank_set = set(bank["perts"])
    cohort100 = sorted(r5_set & bank_set)
    s302 = read_symmetry302()
    s302_primary = {k for k, v in s302.items() if v.get("cohort") == "primary"}
    c302 = sorted(r5_set & s302_primary)
    c302_all = sorted(r5_set & set(s302))
    joung_sym = set(m.pert_symbol[jng].tolist())
    c13 = sorted(r5_set & joung_sym)
    man["cohorts"] = dict(
        nakatake_x_bank=dict(n=len(cohort100), claimed=100, sha256=L.sha256_strings(cohort100)),
        nakatake_x_302_cohort=dict(n=len(c302), claimed=58, n_302_rows=len(s302),
                                   n_302_primary=len(s302_primary),
                                   n_with_302_scored_any=len(c302_all),
                                   sha256=L.sha256_strings(c302)),
        nakatake_x_db_joung=dict(n=len(c13), claimed=13, n_joung_rows=int(jng.size),
                                 sha256=L.sha256_strings(c13), members=c13),
    )
    print(f"cohorts: Nakatake x bank {len(cohort100)} (claim 100) | "
          f"x 302 cohort {len(c302)} of {len(s302)} (claim 58) | "
          f"x db-Joung {len(c13)} of {jng.size} (claim 13)")

    # a TF is a ROW on both sides; it is also a COLUMN when it is on the axis
    on_axis = [t for t in cohort100 if t in bank["col_of"]]
    man["cohorts"]["own_column_present"] = dict(n=len(on_axis))

    # ------------------------------------------------- measure_mask (FREE)
    def mask_rows(rows: np.ndarray) -> np.ndarray:
        out = []
        with h5py.File(L.DB, "r") as f:
            d = f["layers/measure_mask"]
            for lo in range(0, rows.size, L.BLOCK):
                out.append(np.asarray(d[rows[lo:lo + L.BLOCK], :], np.float32))
        return np.concatenate(out, axis=0)

    mm_r5 = mask_rows(r5)
    if np.any((mm_r5 != 0) & (mm_r5 != 1)):
        raise L.Refuse("R5 measure_mask is not binary")
    mm_r5b = mm_r5.astype(bool)
    cov_src = mm_r5b.sum(axis=1)
    r5_row_const = bool(np.all(mm_r5b == mm_r5b[0]))
    can_r5 = L.collapse_mask_only(mm_r5b, proj)
    cov_can = can_r5.sum(axis=1)
    man["R5"]["mask"] = dict(
        source_coverage_min=int(cov_src.min()), source_coverage_max=int(cov_src.max()),
        source_coverage_median=float(np.median(cov_src)),
        source_coverage_fraction_median=float(np.median(cov_src) / m.n_var),
        row_constant=r5_row_const,
        canonical_coverage_min=int(cov_can.min()), canonical_coverage_max=int(cov_can.max()),
        canonical_coverage_median=float(np.median(cov_can)),
        n_rows_all_measured=int((cov_src == m.n_var).sum()),
    )
    print(f"R5 mask: source coverage {cov_src.min()}..{cov_src.max()} "
          f"(median frac {np.median(cov_src)/m.n_var:.4f}), row-constant {r5_row_const}; "
          f"canonical {cov_can.min()}..{cov_can.max()}")
    del mm_r5

    mm_bank = mask_rows(bnk).astype(bool)
    can_bank = L.collapse_mask_only(mm_bank, proj)
    bank_const = bool(np.all(can_bank == can_bank[0]))
    man["bank_source_rows"] = dict(
        n=int(bnk.size), first_row=int(bnk[0]), last_row=int(bnk[-1]),
        row_index_sha256=L.sha256_array(bnk),
        mask_row_constant_canonical=bank_const,
        canonical_coverage_median=float(np.median(can_bank.sum(axis=1))),
    )
    del mm_bank

    # ------------------------------------------------------------ gene axis
    stored = np.load(L.H1BANK / "measured_mask.npy", mmap_mode="r")
    bank_mask0 = np.asarray(stored[0]).astype(bool)
    bank_mask_const = bool(np.array_equal(np.asarray(stored).all(axis=0),
                                          np.asarray(stored).any(axis=0)))
    del stored

    # cohort rows on each side, in the FROZEN cohort order
    r5_row_of = {s: int(i) for s, i in zip(r5_sym, r5)}
    cohort = cohort100                       # frozen, sorted
    r5_cohort_rows = np.asarray([r5_row_of[s] for s in cohort], np.int64)
    kd_cohort_rows = np.asarray([bank["row_of"][s] for s in cohort], np.int64)

    r5_can_cohort = can_r5[np.searchsorted(r5, r5_cohort_rows)]
    kd_can_cohort = np.asarray(
        np.load(L.H1BANK / "measured_mask.npy", mmap_mode="r")[kd_cohort_rows]).astype(bool)

    axis_mask = r5_can_cohort.all(axis=0) & kd_can_cohort.all(axis=0)
    axis = [g for g, k in zip(canon, axis_mask) if k]
    man["axis"] = dict(
        n=len(axis), sha256=L.sha256_strings(axis),
        canonical_genes=len(canon),
        bank_mask_row_constant=bank_mask_const,
        bank_measured=int(bank_mask0.sum()),
        source_covered_canonical=int(proj.covered.sum()),
        duplicate_collapse_groups=len(proj.dup_groups),
        n_source_features=m.n_var, n_unique_source_symbols=int(np.unique(m.var_symbol).size),
        definition=("canonical genes measured in EVERY cohort Nakatake row (after 115's "
                    "duplicate-collapse mask AND) and in EVERY cohort gse295214_h1 bank row"),
    )
    print(f"axis: {len(axis)} genes of {len(canon)} canonical "
          f"(bank measured {int(bank_mask0.sum())}, bank mask row-constant {bank_mask_const})")

    # own-column presence on the axis
    apos = {g: i for i, g in enumerate(axis)}
    own_on_axis = [t for t in cohort if t in apos]
    man["cohorts"]["own_column_on_axis"] = dict(n=len(own_on_axis))
    print(f"own target column present on the axis for {len(own_on_axis)}/{len(cohort)} cohort TFs")

    # ------------------------------------------------------------- write out
    (L.WORK / "AXIS.txt").write_text("\n".join(axis) + "\n")
    with open(L.WORK / "R5-ROWS.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["order", "source_row_index", "obs_id", "symbol", "ensembl",
                    "sign", "source_mask_coverage", "canonical_mask_coverage",
                    "in_bank_cohort"])
        cs = set(cohort)
        for k in range(r5.size):
            w.writerow([k, int(r5[k]), r5_id[k], r5_sym[k], r5_ens[k], int(r5_sign[k]),
                        int(cov_src[k]), int(cov_can[k]), int(r5_sym[k] in cs)])
    with open(L.WORK / "COHORT.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["order", "symbol", "r5_source_row", "bank_row", "own_column_on_axis",
                    "in_302_cohort", "in_db_joung"])
        c302s, c13s = set(c302), set(c13)
        for k, s in enumerate(cohort):
            w.writerow([k, s, r5_row_of[s], bank["row_of"][s], int(s in apos),
                        int(s in c302s), int(s in c13s)])
    man["artifacts"] = {
        "work/AXIS.txt": L.sha256_file(L.WORK / "AXIS.txt"),
        "work/R5-ROWS.tsv": L.sha256_file(L.WORK / "R5-ROWS.tsv"),
        "work/COHORT.tsv": L.sha256_file(L.WORK / "COHORT.tsv"),
    }
    man["inputs"] = {}
    for p in (L.H1BANK / "delta.npy", L.H1BANK / "b_ctrl.npy", L.H1BANK / "measured_mask.npy",
              L.H1BANK / "genes.txt", L.H1BANK / "perts.txt", L.H1BANK / "meta.json",
              L.HELD279, L.ZNF279, L.PERSRC290, L.SYM302, L.SUMS302, L.CANON):
        man["inputs"][str(p.relative_to(L.ROOT))] = (
            dict(exists=True, size_bytes=p.stat().st_size,
                 sha256=(L.sha256_file(p) if p.stat().st_size < 300_000_000 else "SKIPPED_LARGE"))
            if p.exists() else dict(exists=False))
    man["runtime"] = dict(seconds=round(time.time() - t0, 1),
                          peak_rss_gib=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30, 2))
    man["read_contract"] = dict(
        X_rows_read=0, control_rows_read=0,
        measure_mask_rows_read=int(r5.size + bnk.size),
        note="metadata pass: obs/var/uns + layers/measure_mask only")
    L.dump_json(L.WORK / "meta303.json", man)
    print(f"\npeak RSS {man['runtime']['peak_rss_gib']:.2f} GiB, {man['runtime']['seconds']}s")


if __name__ == "__main__":
    sys.exit(main())
