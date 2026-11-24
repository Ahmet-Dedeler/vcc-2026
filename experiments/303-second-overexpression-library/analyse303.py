"""Saved-array C7 driver. Never imports access303/lib303 or opens any atlas.

--preflight reads only metadata and USED KD/feature artifacts; prints input
hashes and coverage to stdout, writes nothing. --analyse requires the completed
hash-bound access303.json/marker/release and reads R5 ONLY from nak303.npz.
The coordinator must include required_inputs() and CODE_NAMES in manifest.files,
and analyse303.py/test_analysis303.py in manifest.analysis_code before release.
No release, marker, preregistration or source product is created by this driver.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import platform
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy
import stats303 as S

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
BANK = ROOT / "experiments/115-gse295214-h1/work/compact-donor"
HELD = ROOT / "experiments/279-znf-state-dependence/work/held279.npz"
ZNF = ROOT / "experiments/279-znf-state-dependence/work/znf279.csv"
GEOM = ROOT / "experiments/290-remainder-response-features/work/persrc290.npz"
SYM302 = ROOT / "experiments/302-orf-vs-crispri-symmetry/work/symmetry302.csv"
SUMS302 = ROOT / "experiments/302-orf-vs-crispri-symmetry/work/sums302.npz"
SOURCE_KEY = "data/external/gse295214/pseudo_bulk_perturbation_database.h5ad"
KD_SOURCE = "gse295214_h1"
CODE_NAMES = ("analyse303.py", "stats303.py", "test_analysis303.py")
LIMITATIONS = [
    "Cross-study/cross-line observed alignment, conditional on target-label exchangeability; not a causal null.",
    "No OE replicate structure: no disattenuation, latent-equivalence bound or reliability plug-in.",
    "HGNC family-block bootstrap NOT implemented; no family-adjusted generalization claim.",
    "Joung OE-vector rebuild from sums302 NOT implemented: 247 measured-axis genes are absent there; no padding or raw-atlas fallback.",
    "Derangements are available as a pure helper but not reported by this driver.",
]


class Refuse(ValueError):
    """Integrity/contract failure; no outcomes may be written."""


def require(ok, message):
    if not ok:
        raise Refuse(message)


def sha_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(4 << 20), b""):
            h.update(b)
    return h.hexdigest()


_LOADED_CODE_HASHES = {str(p.relative_to(ROOT)): sha_file(p)
                       for p in (Path(__file__).resolve(), Path(S.__file__).resolve())}


def sha_array(a):
    a = np.ascontiguousarray(a)
    return hashlib.sha256(a.tobytes() + f"|{a.dtype.str}|{a.shape}".encode()).hexdigest()


def sha_strings(x):
    return hashlib.sha256(("\n".join(x) + "\n").encode()).hexdigest()


def manifest_sha256(manifest):
    return hashlib.sha256(json.dumps(manifest, sort_keys=True, separators=(",", ":"),
                                    ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def read_json(path):
    def pairs(items):
        out = {}
        for k, v in items:
            require(k not in out, f"duplicate JSON key {k}")
            out[k] = v
        return out

    def invalid(s):
        raise Refuse(f"nonfinite JSON token {s}")
    value = json.loads(Path(path).read_text(), object_pairs_hook=pairs, parse_constant=invalid)
    require(type(value) is dict, "JSON must be an object")
    return value


def labels(values, name):
    a = np.asarray(values)
    require(a.ndim == 1 and a.dtype.kind in "UO", f"{name}: expected string labels")
    out = a.tolist()
    require(all(type(x) is str and x and x.strip() == x for x in out), f"{name}: invalid label")
    require(len(set(out)) == len(out), f"{name}: duplicate label")
    return tuple(out)


def text_labels(path):
    return labels(np.asarray(Path(path).read_text().splitlines()), str(path))


def read_table(path, delimiter=","):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def integer(text):
    require(isinstance(text, str) and text.isdecimal(), f"non-integer metadata: {text!r}")
    return int(text)


@dataclass(frozen=True)
class Design:
    axis: tuple[str, ...]
    cohort: tuple[str, ...]
    r5_symbols: tuple[str, ...]
    r5_ids: tuple[str, ...]
    r5_rows: np.ndarray
    oe_rows: np.ndarray
    kd_rows: np.ndarray
    own: np.ndarray
    db_joung: np.ndarray


def load_design(work=WORK) -> Design:
    """Exact, ordered manifest alignment; no data-dependent cohort changes."""
    axis = text_labels(work / "AXIS.txt")
    c = read_table(work / "COHORT.tsv", "\t")
    r = read_table(work / "R5-ROWS.tsv", "\t")
    require(len(c) == 100 and len(r) == 510 and len(axis) == 15309,
            "expected metadata-fixed 100 TFs, 510 R5 rows, 15309 measured genes")
    for table in (c, r):
        require([integer(x["order"]) for x in table] == list(range(len(table))),
                "metadata row order mismatch")
    cohort = labels(np.asarray([x["symbol"] for x in c]), "cohort")
    require(cohort == tuple(sorted(cohort)), "cohort is not in frozen sorted order")
    rs = labels(np.asarray([x["symbol"] for x in r]), "R5 symbols")
    ids = labels(np.asarray([x["obs_id"] for x in r]), "R5 obs IDs")
    rows = np.asarray([integer(x["source_row_index"]) for x in r], np.int64)
    require(np.all(np.diff(rows) > 0), "R5 row indices not strictly increasing")
    pos, apos = {s: i for i, s in enumerate(rs)}, {g: i for i, g in enumerate(axis)}
    require(set(cohort) <= set(rs), "cohort member absent from R5 metadata")
    oe_rows = np.asarray([pos[s] for s in cohort], np.int64)
    require(np.array_equal(rows[oe_rows], [integer(x["r5_source_row"]) for x in c]),
            "COHORT/R5 source-row mapping disagrees")
    own = np.asarray([apos.get(s, -1) for s in cohort], np.int64)
    require(int((own >= 0).sum()) == 98, "expected 98 cohort own-target coordinates")
    for x, present in zip(c, own >= 0):
        require(x["own_column_on_axis"] == str(int(present)), "own-column flag mismatch")
        require(x["in_db_joung"] in ("0", "1"), "invalid db-Joung flag")
    meta = read_json(work / "meta303.json")
    for name in ("AXIS.txt", "COHORT.tsv", "R5-ROWS.tsv"):
        require(meta["artifacts"][f"work/{name}"] == sha_file(work / name),
                f"metadata artifact pin mismatch: {name}")
    require(meta["cohorts"]["nakatake_x_bank"]["sha256"] == sha_strings(cohort),
            "cohort symbol content hash mismatch (distinct from TSV file hash)")
    require(meta["axis"]["sha256"] == sha_strings(axis), "axis content hash mismatch")
    require(meta["R5"]["row_index_sha256"] == sha_array(rows), "R5 row content hash mismatch")
    require(meta["R5"]["symbols_sha256"] == sha_strings(rs), "R5 symbol content hash mismatch")
    return Design(axis, cohort, rs, ids, rows, oe_rows,
                  np.asarray([integer(x["bank_row"]) for x in c], np.int64), own,
                  np.asarray([x["in_db_joung"] == "1" for x in c], bool))


def required_inputs():
    """All readout data dependencies; exact repo-relative paths, no raw source."""
    return (
        *(WORK / n for n in ("AXIS.txt", "COHORT.tsv", "R5-ROWS.tsv", "meta303.json")),
        *(BANK / n for n in ("delta.npy", "measured_mask.npy", "genes.txt", "perts.txt")),
        HELD, GEOM, ZNF, SYM302,
    )


def input_hashes():
    return {p.relative_to(ROOT).as_posix(): sha_file(p) for p in
            (*required_inputs(), *(HERE / n for n in CODE_NAMES))}


def verify_manifest(manifest):
    """Freshly bind USED inputs only: NEVER open/hash the raw R5 source here."""
    require(type(manifest) is dict and type(manifest.get("files")) is dict, "missing manifest.files")
    pins = input_hashes()
    for name, digest in pins.items():
        require(manifest["files"].get(name) == digest, f"unbound/changed readout input: {name}")
    for name, digest in _LOADED_CODE_HASHES.items():
        require(pins.get(name) == digest, f"readout code changed since import: {name}")
    additions = {str((HERE / n).relative_to(ROOT)) for n in ("analyse303.py", "test_analysis303.py")}
    require(additions <= set(manifest.get("analysis_code", [])), "analysis_code omits readout/tests")
    return pins


def load_used(design: Design):
    """Read the already-used KD bank and frozen per-target feature artifacts."""
    genes, perts = text_labels(BANK / "genes.txt"), text_labels(BANK / "perts.txt")
    require(len(perts) == 2354, "KD bank row count changed")
    rows = S.index_vector(design.kd_rows, len(design.cohort), upper=len(perts))
    require(tuple(perts[i] for i in rows) == design.cohort, "KD row labels do not match COHORT.tsv")
    axis_set = set(design.axis)
    require(tuple(g for g in genes if g in axis_set) == design.axis, "KD gene-axis order mismatch")
    col_of = {g: i for i, g in enumerate(genes)}
    cols = np.asarray([col_of[g] for g in design.axis])
    D = np.load(BANK / "delta.npy", mmap_mode="r", allow_pickle=False)
    M = np.load(BANK / "measured_mask.npy", mmap_mode="r", allow_pickle=False)
    require(D.shape == M.shape == (len(perts), len(genes)), "KD shape mismatch")
    require(D.dtype == np.float32 and M.dtype == np.bool_, "KD float32/boolean types required")
    total = np.zeros(len(cols), np.float64)
    for lo in range(0, len(perts), 128):
        part = S.numeric(D[lo:lo+128, :], "KD bank", 2)
        require(np.all(M[lo:lo+128, :][:, cols]), "unmeasured KD value on fixed axis")
        total += part[:, cols].sum(0)
    KD = np.asarray(D[rows][:, cols], np.float64)
    mean = total / len(perts)
    del D, M
    features, provenance = load_features(design.cohort)
    sym = read_table(SYM302)
    syms = labels(np.asarray([x["symbol"] for x in sym]), "302 scored symbols")
    smap = dict(zip(syms, sym))
    sg = np.asarray([float(smap[s]["S_g"]) if s in smap else np.nan for s in design.cohort])
    primary302 = np.asarray([s in smap and smap[s]["cohort"] == "primary" for s in design.cohort])
    all302 = np.asarray([s in smap for s in design.cohort])
    require(int(all302.sum()) == 58 and int(primary302.sum()) == 45,
            "302 metadata-fixed overlaps changed (58 scored / 45 primary)")
    require(np.isfinite(sg[all302]).all(), "nonfinite scored S_g in SYM302")
    return dict(KD=KD, bank_mean=mean, features=features, feature_provenance=provenance,
                Sg302=sg, all302=all302, primary302=primary302)


def _map_values(cohort, genes, values):
    pos = {g: i for i, g in enumerate(genes)}
    require(np.asarray(values).shape == (len(genes),), "feature axis shape mismatch")
    v = S.numeric(values, "feature", 1, missing=True)
    return np.asarray([v[pos[s]] if s in pos else np.nan for s in cohort])


def load_features(cohort):
    # Legacy USED npz archives store string labels as object arrays. These
    # specific local artifacts are release-hashed before post-access loading.
    with np.load(HELD, allow_pickle=True) as z:
        genes, keys = labels(z["genes"], "held genes"), labels(z["keys"], "held sources")
        pairs = np.asarray(z["pairs_cross"])
        require(pairs.ndim == 2 and pairs.shape[1] == 2 and pairs.dtype.kind in "iu",
                "held pair indices invalid")
        require(np.all((pairs >= 0) & (pairs < len(keys))), "held pair index out of bounds")
        require(len({tuple(sorted(x)) for x in pairs}) == len(pairs), "duplicate held pairs")
        q = keys.index(KD_SOURCE)
        use = np.all(pairs != q, axis=1)
        require(int(use.sum()) == 338, "expected 338 source-excluded held pairs")
        cross = S.numeric(z["Ccross"], "Ccross", 2, missing=True)
        require(cross.shape == (len(pairs), len(genes)), "held feature alignment mismatch")
        finite = np.isfinite(cross[use])
        count = finite.sum(0)
        excl = np.full(len(genes), np.nan)
        np.divide(np.where(finite, cross[use], 0).sum(0), count, out=excl, where=count > 0)
        f = dict(C_excl=_map_values(cohort, genes, excl),
                 C_held=_map_values(cohort, genes, z["c_held"]))
        counts = _map_values(cohort, genes, count)
        provenance = dict(excluded_source=KD_SOURCE, n_pairs=int(use.sum()),
                          pairs=[[keys[i], keys[j]] for i, j in pairs[use]],
                          n_pairs_per_tf=[int(v) if np.isfinite(v) else 0 for v in counts])
    require(int(np.isfinite(f["C_excl"]).sum()) == 97, "C_excl finite coverage changed from 97")
    with np.load(GEOM, allow_pickle=True) as z:
        genes, keys = labels(z["genes"], "290 genes"), labels(z["keys"], "290 sources")
        q = keys.index(KD_SOURCE)
        for src, dst in (("pc20", "gen_pc20"), ("breadth", "breadth"), ("mag", "mag290")):
            a = S.numeric(z[src], src, 2, missing=True)
            require(a.shape == (len(keys), len(genes)), "290 feature alignment mismatch")
            f[dst] = _map_values(cohort, genes, a[q])
    records = read_table(ZNF)
    names = labels(np.asarray([x["gene"] for x in records]), "279 CSV genes")
    rr = np.asarray([float(x["r_bar"]) if x["r_bar"] else np.nan for x in records])
    f["r_bar"] = _map_values(cohort, names, rr)
    return f, provenance


def validate_product(z, d: Design, access: dict):
    """Pure schema validation, with metadata-exact row/column alignment."""
    required = {"rows", "symbols", "obs_id", "axis", "delta", "mask", "b_ctrl",
                "control_cpm_sum", "perturbed_cpm_sum", "control_rows_identical"}
    require(required <= set(z), "incomplete saved product schema")
    require(labels(z["axis"], "product axis") == d.axis, "product AXIS order mismatch")
    require(labels(z["symbols"], "product symbols") == d.r5_symbols, "product R5 row order mismatch")
    require(labels(z["obs_id"], "product IDs") == d.r5_ids, "product R5 ID order mismatch")
    rows = S.index_vector(z["rows"], len(d.r5_rows))
    require(np.array_equal(rows, d.r5_rows), "product source-row mapping mismatch")
    shape = (len(d.r5_rows), len(d.axis))
    delta, mask = np.asarray(z["delta"]), np.asarray(z["mask"])
    require(delta.shape == mask.shape == shape, "product delta/mask shape mismatch")
    require(delta.dtype == np.float32 and mask.dtype == np.bool_, "product must use float32 delta / boolean mask")
    require(np.all(mask), "unmeasured R5 coordinates on frozen axis; no imputation")
    S.numeric(delta, "product delta", 2)
    b = S.numeric(z["b_ctrl"], "b_ctrl", 1)
    require(b.shape == (shape[1],) and np.all(b >= 0), "invalid b_ctrl")
    valid = True
    for name in ("control_cpm_sum", "perturbed_cpm_sum"):
        v = S.numeric(z[name], name, 1)
        require(v.shape == (shape[0],) and np.all(v >= 0), "invalid CPM diagnostic")
        valid &= bool(np.all(np.abs(v - 10000) <= 1))
    identical = np.asarray(z["control_rows_identical"])
    require(identical.dtype == np.bool_ and identical.shape == (1,), "shared-control flag must be boolean [1]")
    for name, actual in (("cpm_gate_pass", valid), ("control_rows_identical", bool(identical[0])),
                         ("mask_row_constant", True), ("all_finite", True)):
        require(type(access.get(name)) is bool and access[name] == actual, f"access diagnostic mismatch: {name}")
    require(access.get("delta_sha256") == sha_array(delta), "delta content digest mismatch")
    require(access.get("mask_sha256") == sha_array(mask), "mask content digest mismatch")
    require(type(access.get("n_rows")) is int and access["n_rows"] == shape[0], "access row count mismatch")
    require(type(access.get("n_axis")) is int and access["n_axis"] == shape[1], "access axis count mismatch")
    return dict(cpm_gate_pass=valid, validity="VALID" if valid else "VOID",
                shared_control_identical=bool(identical[0]),
                shared_control_interpretation=("retained" if identical[0] else "WITHDRAWN"),
                cpm_rule="all 510 source-axis CPM sums within 1 of 10000")


def load_completed_product(d: Design, work=WORK):
    """Check completed access and all readout pins before decoding any R5 array."""
    paths = {k: work / v for k, v in dict(release="RELEASE.json", marker="ACCESS-MARKER.json",
                                         access="access303.json", product="nak303.npz").items()}
    require(all(p.is_file() and not p.is_symlink() for p in paths.values()),
            "completed release/marker/access/product required; R5 source is never a fallback")
    rel, marker, acc = (read_json(paths[k]) for k in ("release", "marker", "access"))
    require(type(rel.get("schema_version")) is int and rel["schema_version"] == 1, "release schema mismatch")
    require(rel.get("claim") == "C7" and rel.get("reserve") == "R5"
            and rel.get("authorization") == "ONE_SHOT_R5_EXPRESSION", "wrong release authority")
    man = rel["manifest"]
    pins = verify_manifest(man)
    core = manifest_sha256(man)
    release_hash = sha_file(paths["release"])
    require(acc.get("manifest_sha256") == marker.get("live", {}).get("manifest_sha256") == core,
            "access/marker manifest mismatch")
    require(acc.get("release_sha256") == marker.get("release_sha256") == release_hash,
            "access/marker release mismatch")
    require(marker.get("release") == rel and marker.get("claim") == "C7", "marker release mismatch")
    require(marker.get("live", {}).get("source_sha256") == man["files"].get(SOURCE_KEY),
            "source provenance mismatch")
    final = HERE / "PREREGISTRATION.md"
    require(man["files"].get(str(final.relative_to(ROOT))) == sha_file(final), "final preregistration changed")
    # Hash the SAME byte buffer decoded below; no hash-then-reopen race.
    raw = paths["product"].read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    require(acc.get("product_sha256") == digest, "access303.json product hash mismatch")
    with np.load(io.BytesIO(raw), allow_pickle=False) as z:
        product = {k: z[k] for k in z.files}
    validity = validate_product(product, d, acc)
    provenance = dict(manifest_sha256=core, release_sha256=release_hash,
                      product_sha256=digest, access_sha256=sha_file(paths["access"]),
                      access_started_at=marker.get("started_at"),
                      source_sha256=man["files"][SOURCE_KEY], input_sha256=pins,
                      extraction_code_sha256={
                          k: v for k, v in man["files"].items()
                          if k in (str((HERE / "lib303.py").relative_to(ROOT)),
                                   str((HERE / "access303.py").relative_to(ROOT)))},
                      readout_version=S.VERSION, python=platform.python_version(),
                      numpy=np.__version__, scipy=scipy.__version__, final_preregistration_sha256=sha_file(final))
    return product, validity, provenance


def association(x, y):
    x, y = np.asarray(x), np.asarray(y)
    keep = np.isfinite(x) & np.isfinite(y)
    return dict(n=int(keep.sum()), rho=float(S.rank_correlation(x[keep], y[keep])))


def analyse_arrays(product, d, used, validity, *, draws=S.N_NULL, boot_draws=S.N_BOOT):
    """Pure readout from supplied arrays; usable with entirely synthetic inputs."""
    OE = np.asarray(product["delta"][d.oe_rows], np.float64)
    KD, features = used["KD"], used["features"]
    O, K = S.cohort_residuals(OE, KD)
    keep = S.fixed_axis_mask(d.own, len(d.cohort), len(d.axis))
    C = S.cosine_matrix(O, K, d.own)
    primary = S.evaluate_c7(C, features["C_excl"], features["gen_pc20"],
                            draws=draws, boot_draws=boot_draws)
    for cell in primary["cells"].values():
        cell["statistical_verdict"] = cell["verdict"]
        if validity["validity"] == "VOID":
            cell["verdict"] = "VOID"      # preserve raw readings, never call this biological FAIL
    s = primary["S"]
    matrices = {
        "uncentered_fixed_axis": S.cosine_matrix(OE, KD, d.own),
        "bankwide_centered_fixed_axis": S.cosine_matrix(
            OE - np.asarray(product["delta"], float).mean(0), KD - used["bank_mean"], d.own),
        "original_pairwise_own_exclusion": S.cosine_matrix(O, K, d.own, exclusion="pairwise"),
    }
    secondary = {}
    for name, mat in matrices.items():
        sn = mat[np.arange(len(s)), primary["permutations"]]
        secondary[name] = dict(observed=S.cell_stats(np.diag(mat)),
                               null_median=S.summary(np.median(sn, axis=1)),
                               null_tail=S.summary(S._tail(sn)), verdict="DESCRIPTIVE ONLY")
    own_present = d.own >= 0
    own_oe, own_kd = np.full(len(s), np.nan), np.full(len(s), np.nan)
    ii = np.flatnonzero(own_present)
    own_oe[ii], own_kd[ii] = OE[ii, d.own[ii]], KD[ii, d.own[ii]]
    secondary["own_transcript"] = dict(
        n=int(own_present.sum()), oe_median=float(np.nanmedian(own_oe)),
        kd_median=float(np.nanmedian(own_kd)), oe_share_up=float(np.mean(own_oe[ii] > 0)),
        kd_share_down=float(np.mean(own_kd[ii] < 0)),
        S_on_own_up=S.cell_stats(s[own_oe > 0]),
        db_joung_cohort_overlap_n=int(d.db_joung.sum()),
        S_on_db_joung_cohort_overlap=S.cell_stats(s[d.db_joung]))
    secondary["302"] = {
        name: association(used["Sg302"][mask], s[mask])
        for name, mask in (("ALL58", used["all302"]), ("primary45", used["primary302"]))}
    secondary["302"]["scope"] = "saved SYM302 S_g versus C7 S; descriptive; C6 processing VOID; no reliability plug-in"
    kmag, omag = np.linalg.norm(K[:, keep], axis=1), np.linalg.norm(O[:, keep], axis=1)
    secondary["baselines"] = {
        name: dict(minusS=association(x, -s), absS=association(x, np.abs(s)))
        for name, x in {**features, "KD_magnitude": kmag, "OE_magnitude": omag}.items()}
    controls = np.column_stack((kmag, omag))
    secondary["partial_rho_given_magnitudes"] = dict(
        C7b=S.partial_spearman(features["C_excl"], -s, controls),
        C7c=S.partial_spearman(features["gen_pc20"], -s, controls),
        C7d=S.partial_spearman(features["gen_pc20"], np.abs(s), controls))
    rows = []
    for i, symbol in enumerate(d.cohort):
        rows.append(dict(symbol=symbol, S=float(s[i]), abs_S=float(abs(s[i])),
                         own_oe=float(own_oe[i]), own_kd=float(own_kd[i]),
                         KD_magnitude=float(kmag[i]), OE_magnitude=float(omag[i]),
                         Sg302=float(used["Sg302"][i]),
                         **{name: float(x[i]) for name, x in features.items()},
                         **{name: float(np.diag(mat)[i]) for name, mat in matrices.items()}))
    report = dict(readout_version=S.VERSION, validity=validity, cells=primary["cells"],
                  design=dict(n_tf=len(s), measured_axis=len(d.axis), scored_axis=int(keep.sum()),
                              dropped_own_columns=int((~keep).sum()),
                              axis_sha256=sha_strings([g for g, k in zip(d.axis, keep) if k]),
                              cohort_sha256=sha_strings(d.cohort), seed=S.SEED,
                              n_null=draws, n_boot=boot_draws, fixed_points=primary["fixed_points"],
                              permutation_scope=primary["permutation_scope"], ci_scope=S.CI_SCOPE),
                  observed=primary["observed_summary"], secondaries=secondary,
                  no_signal_rows=dict(OE=[s for s, norm in zip(d.cohort, omag) if norm == 0],
                                      KD=[s for s, norm in zip(d.cohort, kmag) if norm == 0]),
                  feature_provenance=used["feature_provenance"], limitations=LIMITATIONS)
    arrays = dict(cosine_matrix=C, S=s, null_S=primary["null_S"], permutations=primary["permutations"],
                  cohort=np.asarray(d.cohort), axis=np.asarray(d.axis)[keep],
                  **{f"null_{k}": v for k, v in primary["null_stats"].items()},
                  **{f"bootstrap_{k}": v for k, v in primary["bootstrap"].items()})
    return report, rows, arrays


def json_clean(x):
    if isinstance(x, dict):
        return {k: json_clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)):
        return [json_clean(v) for v in x]
    if isinstance(x, (float, np.floating)):
        return float(x) if np.isfinite(x) else None
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def render_results(report):
    """Full machine-readable readings embedded in the human report as well."""
    lines = ["# C7 — saved-array observed-response readout", "",
             f"Validity: **{report['validity']['validity']}**. Four-cell Holm FWER 0.05.", "",
             "C7a: observed median < 0 AND negative-tail excess > 0 AND Holm(max of two p).",
             "C7b/c/d: positive observed rho in the declared direction AND Holm.", "",
             "| Cell | Verdict | Raw p | Holm p | Estimable |",
             "|---|---|---:|---:|---|"]
    for name, c in report["cells"].items():
        lines.append(f"| {name} | {c['verdict']} | {c['p_raw']:.6g} | {c['p_holm']:.6g} | {c['estimable']} |")
    lines += ["", S.CI_SCOPE + ".", "",
              "Undefined expression/rank cells remain in the family at p=1; no zero-signal imputation.",
              "CPM premise failure voids every conclusion, regardless of nominal significance.", "",
              "## Full readout, intervals, null summaries, provenance and limitations", "",
              "```json", json.dumps(json_clean(report), indent=2, allow_nan=False), "```", ""]
    return "\n".join(lines)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--preflight", action="store_true")
    mode.add_argument("--analyse", action="store_true")
    args = p.parse_args(argv)
    try:
        # A missing completed access is rejected before even loading USED arrays.
        if args.analyse:
            require((WORK / "access303.json").is_file(), "no completed access; R5 remains locked")
        d = load_design()
        if args.preflight:
            used = load_used(d)
            out = dict(status="PREACCESS_ONLY_NO_R5_READ", input_sha256=input_hashes(),
                       cohort=len(d.cohort), axis=len(d.axis), primary_axis=15211,
                       features={k: int(np.isfinite(v).sum()) for k, v in used["features"].items()},
                       source_excluded_pairs=used["feature_provenance"]["n_pairs"],
                       overlap302=[int(used["all302"].sum()), int(used["primary302"].sum())],
                       limitations=LIMITATIONS)
            print(json.dumps(out, indent=2))
            return 0
        product, validity, provenance = load_completed_product(d)
        used = load_used(d)
        report, rows, arrays = analyse_arrays(product, d, used, validity)
        report["provenance"] = provenance
        # Refuse output if any readout input changed during the analysis.
        require(input_hashes() == provenance["input_sha256"], "inputs changed during analysis")
        (WORK / "stats303.json").write_text(json.dumps(json_clean(report), indent=2, allow_nan=False) + "\n")
        with open(WORK / "symmetry303.csv", "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        np.savez_compressed(WORK / "analysis303.npz", **arrays)
        (HERE / "RESULTS.md").write_text(render_results(report))
        print(" ".join(f"{k}={c['verdict']}" for k, c in report["cells"].items()))
        return 0
    except (Refuse, ValueError, KeyError, OSError) as exc:
        print(f"REFUSED: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
