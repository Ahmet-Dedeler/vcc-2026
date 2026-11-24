"""317 — build the OPERATING-POINT sweep as h1-test candidates, for the REAL scorer.

Derived verbatim from `experiments/221-cloud-scorer/build_h1_candidates.py` (all four of its
known-answer gates are kept). The ONE structural change: `CONFIGS` now carries the emission /
amplitude knobs per config, so a rung differs from `base` by exactly one of

    overdispersion · tail_spread · keep_top · keep_spread · effect_scale (delivered amplitude)

`base` is 221's `rc_callset` rebuilt bit-for-bit: same effect matrix, same seed, same emitter,
same knobs. Its sha256 must equal 221's recorded
`048b9286030ff0124514134b55d935fdb79a4d002cc04d590d83676b32553c77`, and a byte-identical
candidate through a deterministic scorer must reproduce 221's recorded six-member row. That is
this experiment's known-answer row and it is checked twice: on the digest here, and on the
scored `agg_results.csv` in the cloud.

GATE 4 IS RE-PURPOSED AS THE STRUCTURAL CLAIM. `emit_counts(mode="quasi")` repairs every gene's
panel total exactly, so an emission-knob arm must be panel-total digest-IDENTICAL to `base` on
85/85 blocks (=> `pds` and `mse` cannot move), while an amplitude arm must differ on 85/85
(=> it is the only knob in the sweep that can move `mse`). The gate asserts both directions.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

# a15's flags. Set BEFORE the builder is imported — it reads them at module scope.
os.environ.setdefault("VCC_ALPHA", "0.982")
os.environ.setdefault("VCC_MASS_FIX", "minnorm")
os.environ.setdefault("VCC_WHITEN_GAMMA", "0.65")
os.environ.setdefault("VCC_WHITEN_LAM", "0.0")
os.environ.setdefault("VCC_CD4_TARGET_F", "0.10")

sys.path.insert(0, str(REPO / "src"))

TS = 50_000.0            # bulk_target_sum, the scorer's own (095 run_params.yaml)
NCELL = 400
CTRL = "non-targeting"
PERT_COL = "target_gene"

# a15's shipped emission policy, read off `submissions/a15-cd4f010.vcc.meta.json`.
MIN_CPM = 10.0
SPREAD, TAIL_SPREAD, KEEP_TOP, KEEP_SPREAD = 0.25, 1.0, 8000, 0.25
SPREAD_MAP = "block"
OVERDISPERSION, DISPERSION = 4.0, 0.1
SEED = 0                                    # library_pick_seed; emission_seed = SEED + 1_000_003
RC_JITTER_SEED_BASE = 1000
RC_SCORER_MIN_CPM = 5.0
RC_CALLSET_HEAD, RC_CALLSET_CPM, RC_CALLSET_Z_CAP = 400, 2200, 1.90

# gate 1 — the two shipped board points on the operator path (170/gate.py's ANCHORS).
OPERATOR_SHA = "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"
ANCHORS = {
    (0.65, 0.0): {   # a15
        "before_operator_median_norm": 4.122892101713489,
        "after_operator_median_norm": 25.73792992728807,
        "operator_preserve_scale": 0.16018740098217024,
        "after_all_median_norm": 4.2292956044543,
        "cd4_f_per_unit": 1.055664536904655,
        "cd4_solved_fraction": 0.09472706196346517,
        "cd4_candidate_pre_preserve_median_norm": 4.166553356202115,
        "cd4_post_add_global_preserve_scale": 0.9895210139518232,
    },
    (0.80, 0.007): {  # a19
        "before_operator_median_norm": 4.122892101713494,
        "after_operator_median_norm": 22.030919431001,
        "operator_preserve_scale": 0.1871411728696139,
        "after_all_median_norm": 4.235393246076997,
        "cd4_f_per_unit": 1.0541447101386758,
        "cd4_solved_fraction": 0.094863635929876,
    },
}
RTOL = 1e-11

# The four rungs. `cd4_target_f` and `emitter` are the ONLY things that move.
# `a15_seed1` is not a rung: it is a RE-SEEDED TWIN of `a15`, whose only difference is the
# emission seed. Its distance from `a15` on every member is this surface's own emission-noise
# floor, and a rung's transfer factor is only readable where its local delta clears it.
# ---------------------------------------------------------------------------- 317's sweep
#
# `base` reproduces 221's `rc_callset` exactly (sha256 gate below). Every other row moves
# EXACTLY ONE knob off it. `effect_scale` is the delivered log-effect amplitude, applied in
# `profile_from_effect`, which that function's own docstring says "composes multiplicatively
# with the VCC_ALPHA env var, which scales the same" — so `effect_scale = s` is the shipped
# amplitude `VCC_ALPHA = 0.982` moved to `0.982 * s` at the point where the metric reads it.
BASE_KNOBS = {"tail_spread": TAIL_SPREAD, "keep_top": KEEP_TOP,
              "keep_spread": KEEP_SPREAD, "overdispersion": OVERDISPERSION,
              "effect_scale": 1.0}
RC_CALLSET_SHA_221 = "048b9286030ff0124514134b55d935fdb79a4d002cc04d590d83676b32553c77"

CONFIGS = [
    {"name": "base",    "cd4_target_f": 0.10, "emitter": "rc_callset", "seed": 0},
    # --- emission knobs: pseudobulk-preserving by construction, must be digest-identical ---
    {"name": "od1",     "cd4_target_f": 0.10, "emitter": "rc_callset", "seed": 0,
     "overdispersion": 1.0},
    {"name": "od8",     "cd4_target_f": 0.10, "emitter": "rc_callset", "seed": 0,
     "overdispersion": 8.0},
    {"name": "ks010",   "cd4_target_f": 0.10, "emitter": "rc_callset", "seed": 0,
     "keep_spread": 0.10},
    {"name": "ts070",   "cd4_target_f": 0.10, "emitter": "rc_callset", "seed": 0,
     "tail_spread": 0.70},
    {"name": "kt16000", "cd4_target_f": 0.10, "emitter": "rc_callset", "seed": 0,
     "keep_top": 16000},
    # --- amplitude: the ONLY knob in the sweep that can move `mse` ---
    {"name": "amp070",  "cd4_target_f": 0.10, "emitter": "rc_callset", "seed": 0,
     "effect_scale": 0.70},
    {"name": "amp085",  "cd4_target_f": 0.10, "emitter": "rc_callset", "seed": 0,
     "effect_scale": 0.85},
    {"name": "amp120",  "cd4_target_f": 0.10, "emitter": "rc_callset", "seed": 0,
     "effect_scale": 1.20},
    {"name": "amp145",  "cd4_target_f": 0.10, "emitter": "rc_callset", "seed": 0,
     "effect_scale": 1.45},
    {"name": "amp200",  "cd4_target_f": 0.10, "emitter": "rc_callset", "seed": 0,
     "effect_scale": 2.00},
]


def peak_rss_gib() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024 ** 3) if sys.platform == "darwin" else r / (1024 ** 2)


def builder(tag: str):
    """Import the PRODUCTION builder module (a fresh instance, so env is re-read)."""
    path = REPO / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location(f"prod_builder_221_{tag}", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[f"prod_builder_221_{tag}"] = mod
    spec.loader.exec_module(mod)
    return mod


def h5_labels(h, path):
    o = h[path]
    if hasattr(o, "keys") and "categories" in o:
        c = np.asarray([x.decode() if isinstance(x, bytes) else str(x)
                        for x in o["categories"][:]])
        return c[o["codes"][:]]
    v = o["values"][:] if hasattr(o, "keys") else o[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def to_h1_axis(effect_prod, prod_genes, h1_genes):
    """170's `to_h1_axis`: 18,533 production columns -> h1-test's own 18,080-gene axis."""
    pos = {g: i for i, g in enumerate(map(str, prod_genes))}
    idx = np.asarray([pos.get(str(g), -1) for g in h1_genes], dtype=np.int64)
    out = np.zeros((effect_prod.shape[0], idx.size), dtype=np.float64)
    have = idx >= 0
    out[:, have] = effect_prod[:, idx[have]]
    return out, int((~have).sum())


def block_digest(values, columns, pointers, n_genes):
    """Per-gene panel total of one emitted block, as a sha256 over the exact int totals.

    Rule 13's profile-preserving test. `pds`, `mse`'s profile half and `nmae` all read the
    emitted pseudobulk and nothing else, so two arms with equal per-gene totals on every
    block are pinned on those members BY CONSTRUCTION — but AGENTS.md is explicit that the
    identity must be proven with a digest, never asserted.
    """
    tot = np.zeros(n_genes, dtype=np.int64)
    np.add.at(tot, columns.astype(np.int64), values.astype(np.int64))
    return hashlib.sha256(np.ascontiguousarray(tot).tobytes()).hexdigest(), tot


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", required=True, help="h1-test.h5ad")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--pre", default=str(REPO / "experiments/170-gamma-mse/work/pre.npz"))
    ap.add_argument("--cd4", default=str(REPO / "experiments/170-gamma-mse/work/cd4.npz"))
    ap.add_argument("--h1-ka",
                    default=str(REPO / "experiments/170-gamma-mse/work/eff/h1_g0.65.npz"))
    ap.add_argument("--whitener", default=str(REPO / "data/donors/whitener_union_pruned.npz"))
    ap.add_argument("--only", default="", help="comma-separated subset of config names")
    ap.add_argument("--max-targets", type=int, default=0,
                    help="SMOKE TEST ONLY: emit just the first N targets. The resulting "
                         "h5ad is not a scoreable candidate; it exists to exercise the "
                         "emitters cheaply before a cloud run.")
    args = ap.parse_args()

    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    Path(args.report).parent.mkdir(parents=True, exist_ok=True)
    report: dict[str, object] = {"gates": {}, "configs": {}}
    t0 = time.time()

    import h5py
    import pandas as pd
    import scipy.sparse as sp
    import anndata as ad
    from vcc2026 import rc_emit
    from vcc2026.donor_transfer import emit_counts

    B = builder("f010")
    B.assert_frozen_contract()
    print(f"[221] ALPHA={B.ALPHA} MASS_FIX={B.MASS_FIX} GAMMA={B.WHITEN_GAMMA} "
          f"LAM={B.WHITEN_LAM} MIN_CPM={B.MIN_CPM} CD4_TARGET_F={B.CD4_TARGET_F}", flush=True)
    if (B.ALPHA, B.MASS_FIX, B.WHITEN_GAMMA, B.WHITEN_LAM, B.CD4_TARGET_F) != (
            0.982, "minnorm", 0.65, 0.0, 0.10):
        raise SystemExit("builder is not on a15's shipped flags")

    pre = np.load(args.pre, allow_pickle=True)
    cd4z = np.load(args.cd4, allow_pickle=True)
    cols = pre["operator_columns"].astype(np.int64)
    panel_pos = pre["panel_positions"].astype(np.int64)
    prod_genes = pre["genes"].astype(str)
    h1_targets = pre["h1_targets"].astype(str)
    h1_genes_staged = pre["h1_genes"].astype(str)
    mode = str(cd4z["mode"][0])
    cd4_chal = B.EXP062.CD4Prepared(mode=mode, mean_effect=cd4z["chal_mean"],
                                    gene_votes=cd4z["chal_votes"],
                                    target_condition_count=cd4z["chal_counts"],
                                    condition_details=[])
    cd4_h1 = B.EXP062.CD4Prepared(mode=mode, mean_effect=cd4z["h1_mean"],
                                  gene_votes=cd4z["h1_votes"],
                                  target_condition_count=cd4z["h1_counts"],
                                  condition_details=[])
    n_h1_carry = int(np.count_nonzero(np.abs(cd4z["h1_mean"]).sum(axis=1)))
    print(f"[221] CD4 support: challenge "
          f"{int(np.count_nonzero(np.abs(cd4z['chal_mean']).sum(axis=1)))}/300, "
          f"h1-test {n_h1_carry}/85 targets", flush=True)
    report["builder_sha256"] = B.file_sha256(REPO / "scripts/make_submission_pm_cd4.py")
    report["exp062_sha256"] = B.file_sha256(REPO / "experiments/062-cd4-post-whiten/run.py")
    report["cd4_support"] = {"challenge_of_300":
                             int(np.count_nonzero(np.abs(cd4z["chal_mean"]).sum(axis=1))),
                             "h1_of_85": n_h1_carry}

    sha = B.file_sha256(Path(args.whitener))
    if sha != OPERATOR_SHA:
        raise SystemExit(f"operator sha {sha} != shipped {OPERATOR_SHA}")
    whit = B.Whitener.load(Path(args.whitener))
    if whit.evecs.shape != (cols.size, cols.size):
        raise SystemExit(f"operator shape {whit.evecs.shape} != {(cols.size, cols.size)}")
    print(f"[221] operator {Path(args.whitener).name} sha {sha[:8]}… "
          f"shape {whit.evecs.shape}", flush=True)

    # ------------------------------------------------------------------ GATE 1: a15 + a19
    print("\n===== GATE 1: the two shipped board points on the operator path =====", flush=True)
    gate1 = {}
    ok1 = True
    for (g, lam), want in ANCHORS.items():
        whitened, waudit = B.apply_v9_whitening(pre["cand_chal"], whit, cols, gamma=g, lam=lam)
        # Index rather than unpack. `scripts/make_submission_pm_cd4.py` is a LIVE shared file
        # and its HEAD moved under this experiment on 2026-09-04 (commit 481cd09): 218's
        # `VCC_H1_UNION_G` second residual took `apply_frozen_cd4_residual` from a 2-tuple to
        # a 3-tuple, and a tuple-unpacking call site broke 20 s into a cloud run. The first
        # two elements are the contract this experiment uses; anything appended after them
        # belongs to another arm and must not break this one.
        caudit = B.apply_frozen_cd4_residual(whitened, cd4_chal, cols, panel_pos)[1]
        got = {
            "before_operator_median_norm": float(waudit["before_operator_median_norm"]),
            "after_operator_median_norm": float(waudit["after_operator_median_norm"]),
            "operator_preserve_scale": float(waudit["operator_preserve_scale"]),
            "after_all_median_norm": float(waudit["after_all_median_norm"]),
            "cd4_f_per_unit": float(caudit["residual_frame"]["f_per_unit_fraction"]),
            "cd4_solved_fraction": float(caudit["residual_frame"]["solved_fraction"]),
        }
        if "cd4_candidate_pre_preserve_median_norm" in want:
            got["cd4_candidate_pre_preserve_median_norm"] = float(
                caudit["post_add_preserve"]["candidate_pre_preserve_median_norm"])
            got["cd4_post_add_global_preserve_scale"] = float(
                caudit["post_add_preserve"]["post_add_global_preserve_scale"])
        label = "a15" if g == 0.65 else "a19"
        rows = {}
        for k, v in want.items():
            rel = abs(got[k] / v - 1.0) if v else abs(got[k])
            rows[k] = {"want": v, "got": got[k], "rel": rel, "pass": bool(rel <= RTOL)}
            ok1 &= rel <= RTOL
            print(f"  {label:4s} {k:<44} {'PASS' if rel <= RTOL else 'FAIL'}  "
                  f"{got[k]!r}  rel {rel:.2e}")
        gate1[label] = rows
        if g == 0.65:
            chal_fraction_010 = float(caudit["residual_frame"]["solved_fraction"])
            f_per_unit = float(caudit["residual_frame"]["f_per_unit_fraction"])
    report["gates"]["gate1_operator_anchors"] = {"pass": bool(ok1), "detail": gate1}
    if not ok1:
        raise SystemExit("GATE 1 FAILED — this is not a15/a19's operator path. STOP.")
    print("  ==> GATE 1 PASSED at both points, every scalar within 1e-11 relative.\n", flush=True)

    # ------------------------------------------------------ GATE 3: the f = 0.34 fraction
    os.environ["VCC_CD4_TARGET_F"] = "0.34"
    B34 = builder("f034")
    whitened_chal, _ = B34.apply_v9_whitening(pre["cand_chal"], whit, cols, gamma=0.65, lam=0.0)
    caudit34 = B34.apply_frozen_cd4_residual(whitened_chal, cd4_chal, cols, panel_pos)[1]
    chal_fraction_034 = float(caudit34["residual_frame"]["solved_fraction"])
    predicted_034 = chal_fraction_010 * 0.34 / 0.10
    rel34 = abs(chal_fraction_034 / predicted_034 - 1.0)
    print(f"[GATE 3] cd4 fraction at f=0.34: re-solved {chal_fraction_034!r} vs "
          f"fraction(0.10)x3.4 {predicted_034!r}  rel {rel34:.2e}")
    report["gates"]["gate3_cd4_fraction"] = {
        "f_per_unit": f_per_unit, "fraction_f010": chal_fraction_010,
        "fraction_f034_resolved": chal_fraction_034,
        "fraction_f034_scaled": predicted_034, "rel": rel34,
        "pass": bool(rel34 <= 1e-12)}
    if rel34 > 1e-12:
        raise SystemExit("GATE 3 FAILED — the CD4 fraction is not linear in f. STOP.")
    del whitened_chal, whitened
    os.environ["VCC_CD4_TARGET_F"] = "0.10"

    # ---------------------------------------------------- h1-test effects at both doses
    whitened_h1, wa_h1 = B.apply_v9_whitening(pre["cand_h1"], whit, cols, gamma=0.65, lam=0.0)
    del whit
    eff_prod = {}
    for f, frac in ((0.10, chal_fraction_010), (0.34, chal_fraction_034)):
        _r = B.EXP062.add_post_whiten_residual(whitened_h1, cd4_h1, frac, cols)
        out, preserve = _r[0], _r[1]
        eff_prod[f] = out
        print(f"[221] h1 effect at f={f}: fraction {frac!r}, preserve scale "
              f"{float(preserve['post_add_global_preserve_scale'])!r}, median row norm "
              f"{float(np.median(np.linalg.norm(out, axis=1))):.6f}", flush=True)

    # ------------------------------------------- GATE 2: reproduce 170's own h1 artifact
    #
    # NOT an ENTRYWISE relative test. `apply_v9_whitening` is a pair of 17,942 x 17,942 GEMMs
    # and BLAS reduction order is architecture-dependent, so the last bits move between the
    # arm64 machine that wrote 170's artifact and this x86_64 worker — 170's own METHOD.md
    # says the same thing about a15's `final_sha256`. On a matrix whose entries span many
    # orders of magnitude a 1e-14 ABSOLUTE wobble on a near-zero entry is a 1e-6 RELATIVE
    # one, and a first version of this gate failed on exactly that (run `ladder-20260904`,
    # max rel 6.9e-06 at max abs 1.3e-14). The scale-appropriate statements are the
    # Frobenius ratio and the per-row cosine, and both are far tighter than the entrywise
    # test ever was.
    ka = np.load(args.h1_ka)["effect"]
    diff = eff_prod[0.10] - ka
    fro = float(np.linalg.norm(diff) / np.linalg.norm(ka))
    max_abs = float(np.abs(diff).max())
    med_norm = float(np.median(np.linalg.norm(ka, axis=1)))
    rn_a = np.linalg.norm(eff_prod[0.10], axis=1)
    rn_b = np.linalg.norm(ka, axis=1)
    cos = (eff_prod[0.10] * ka).sum(1) / np.maximum(rn_a * rn_b, 1e-300)
    min_cos = float(cos.min())
    ok2 = fro <= 1e-12 and (1.0 - min_cos) <= 1e-12 and max_abs <= 1e-9 * med_norm
    print(f"[GATE 2] h1 effect at f=0.10 vs 170/work/eff/h1_g0.65.npz: "
          f"||d||F/||E||F {fro:.3e}, 1-min_row_cos {1.0 - min_cos:.3e}, "
          f"max|d| {max_abs:.3e} against a median row norm {med_norm:.6f} "
          f"-> {'PASS' if ok2 else 'FAIL'}")
    report["gates"]["gate2_h1_effect_vs_170"] = {
        "frobenius_ratio": fro, "one_minus_min_row_cosine": 1.0 - min_cos,
        "max_abs": max_abs, "median_row_norm": med_norm, "pass": bool(ok2)}
    if not ok2:
        raise SystemExit("GATE 2 FAILED — the h1 effect is not 170's. STOP.")
    del ka, diff, cos, rn_a, rn_b

    # -------------------------------------------------------------- h1-test's own control side
    with h5py.File(args.real, "r") as h:
        n, gcount = (int(x) for x in h["X"].attrs["shape"])
        lab = h5_labels(h, "obs/target_gene")
        h1_genes = h5_labels(h, "var/_index")
        Xall = sp.csr_matrix((h["X/data"][:], h["X/indices"][:], h["X/indptr"][:]),
                             shape=(n, gcount))
    if not np.array_equal(h1_genes, h1_genes_staged):
        raise SystemExit("h1-test's gene axis differs from the one 170 staged against")
    ctrl_rows = np.flatnonzero(lab == CTRL)
    Xc = Xall[ctrl_rows].tocsr()
    del Xall
    Pc = np.asarray(Xc.sum(0)).ravel().astype(np.float64)
    b_ctrl = np.log1p(TS * Pc / Pc.sum())
    libraries = np.asarray(Xc.sum(1)).ravel().astype(np.float64)
    surface_targets = np.asarray(sorted(set(map(str, lab)) - {CTRL}))
    if not np.array_equal(surface_targets, h1_targets):
        raise SystemExit("h1-test's target set differs from the one 170 staged")
    print(f"[221] h1-test: {Xc.shape[0]} controls, {surface_targets.size} targets, "
          f"{h1_genes.size} genes; peak RSS {peak_rss_gib():.2f} GiB", flush=True)

    eff_h1 = {}
    for f, E in eff_prod.items():
        eff_h1[f], missing = to_h1_axis(E, prod_genes, h1_genes)
        report.setdefault("axis", {})[str(f)] = {"h1_genes_absent_from_production": missing}
    del eff_prod

    rc_control = rc_emit.control_block_from_counts(Xc, min_cpm=RC_SCORER_MIN_CPM)
    print(f"[221] rc control block: {rc_control.n_controls} controls, "
          f"{int(rc_control.keep.sum())}/{h1_genes.size} genes past the "
          f"{RC_SCORER_MIN_CPM}-CPM scorer gate; peak RSS {peak_rss_gib():.2f} GiB", flush=True)

    # ------------------------------------------------------------------------- emit + write
    wanted = [c for c in CONFIGS if not args.only or c["name"] in args.only.split(",")]
    if args.max_targets:
        surface_targets = surface_targets[:args.max_targets]
        print(f"[221] SMOKE TEST: {surface_targets.size} targets only — NOT a candidate")
    digests: dict[str, list[str]] = {}
    for cfg in wanted:
        name = cfg["name"]
        E = eff_h1[cfg["cd4_target_f"]]
        seed = int(cfg.get("seed", SEED))
        knobs = {**BASE_KNOBS, **{k: cfg[k] for k in BASE_KNOBS if k in cfg}}
        moved = sorted(k for k in BASE_KNOBS if knobs[k] != BASE_KNOBS[k])
        if name != "base" and len(moved) != 1:
            raise SystemExit(f"{name}: one variable per run — moved {moved}")
        k_ts, k_kt = float(knobs["tail_spread"]), int(knobs["keep_top"])
        k_ks, k_od = float(knobs["keep_spread"]), float(knobs["overdispersion"])
        k_es = float(knobs["effect_scale"])
        print(f"[317] {name}: knobs {knobs}  moved={moved or ['(base)']}", flush=True)
        pick_rng, emission_rng = B.make_rngs(seed)
        data, indices, indptr, obs = [], [], [0], []
        digs, npred, delivered = [], [], []
        tstart = time.time()
        for i, target in enumerate(surface_targets):
            profile = (B.profile_from_effect(E[i], b_ctrl, min_cpm=MIN_CPM)
                       if k_es == 1.0 else
                       B.profile_from_effect(E[i], b_ctrl, min_cpm=MIN_CPM,
                                             effect_scale=k_es))
            picks = pick_rng.choice(len(libraries), NCELL,
                                    replace=len(libraries) < NCELL)
            # The spread vector is held at the BASE effect for every arm, including the
            # amplitude arms: the emission POLICY is not the variable there, the delivered
            # amplitude is. `emission_spread_vector` reads a rank axis, which a positive
            # global scale cannot reorder, so this is also the no-op choice.
            spread_vector = B.emission_spread_vector(
                E[i], b_ctrl, min_cpm=MIN_CPM, tail_spread=k_ts,
                keep_top=k_kt, keep_top_frac=None, keep_spread=k_ks,
                rank_score=None, spread_map=SPREAD_MAP, spread_lo=0.0, spread_hi=0.60)
            if cfg["emitter"] is None:
                v, col, ip = emit_counts(
                    profile, libraries[picks], DISPERSION, emission_rng, mode="quasi",
                    overdispersion=k_od, spread=spread_vector)
            else:
                gate = B.emission_rank_axis(E[i], b_ctrl, min_cpm=MIN_CPM)[0]
                key = np.where(gate, np.asarray(E[i], dtype=np.float64), 0.0)
                block, diag = rc_emit.rc_hybrid_block(
                    profile, libraries[picks], DISPERSION, emission_rng,
                    spread=spread_vector, overdispersion=k_od,
                    control=rc_control, key=key,
                    jitter_seed=RC_JITTER_SEED_BASE + i,
                    call_set_first=True,
                    call_set=(rc_emit.CallSet(RC_CALLSET_HEAD, RC_CALLSET_CPM,
                                              RC_CALLSET_Z_CAP)
                              if cfg["emitter"] == "rc_callset" else None),
                    emit_counts=emit_counts)
                v, col, ip = rc_emit.to_csr(block)
                npred.append(int(diag["n_pred_pass1"]))
                if "call_set_delivered" in diag:
                    delivered.append(int(diag["call_set_delivered"]))
                del block
            dig, _ = block_digest(v, col, ip, h1_genes.size)
            digs.append(dig)
            data.append(v); indices.append(col)
            indptr.extend(list(ip[1:] + indptr[-1]))
            obs += [str(target)] * NCELL
            if (i + 1) % 20 == 0:
                print(f"    {name}: {i + 1}/{surface_targets.size} "
                      f"({time.time() - tstart:.0f}s, peak RSS {peak_rss_gib():.2f} GiB)",
                      flush=True)
        # the real control block, verbatim, exactly as `experiments/031-somatic/run.py` does
        data.append(Xc.data); indices.append(Xc.indices)
        indptr.extend(list(Xc.indptr[1:] + indptr[-1]))
        obs += [CTRL] * Xc.shape[0]
        M = sp.csr_matrix((np.concatenate(data).astype(np.float32),
                           np.concatenate(indices).astype(np.int32),
                           np.asarray(indptr, dtype=np.int64)),
                          shape=(len(obs), h1_genes.size))
        del data, indices, indptr
        path = out_dir / f"{name}.h5ad"
        ad.AnnData(X=M, obs=pd.DataFrame({PERT_COL: obs},
                                         index=[f"c{i}" for i in range(len(obs))]),
                   var=pd.DataFrame(index=h1_genes)).write_h5ad(path, compression="lzf")
        digests[name] = digs
        report["configs"][name] = {
            "cd4_target_f": cfg["cd4_target_f"], "emitter": cfg["emitter"],
            "library_pick_seed": seed, "emission_seed": seed + 1_000_003,
            "path": str(path), "sha256": B.file_sha256(path),
            "cells": len(obs), "genes": int(h1_genes.size), "nnz": int(M.nnz),
            "seconds": round(time.time() - tstart, 1),
            "median_n_pred_pass1": (int(np.median(npred)) if npred else None),
            "median_call_set_delivered": (int(np.median(delivered)) if delivered else None),
            "panel_total_digest_sha256": hashlib.sha256(
                "".join(digs).encode()).hexdigest(),
        }
        print(f"[221] wrote {path} ({M.nnz} nnz, {time.time() - tstart:.0f}s, "
              f"peak RSS {peak_rss_gib():.2f} GiB)", flush=True)
        del M

    # ------------------------------------------- GATE 4: the profile-preserving digest test
    #
    # `emit_counts(mode="quasi")` repairs each gene's panel total EXACTLY and the rank-
    # controlled emitters re-key pass 1's own multiset without moving a count, so two arms
    # that share an EFFECT matrix AND a seed share their per-gene totals block for block
    # whatever their emitter. Two things break that and both are expected:
    #   * a different SEED. `pick_rng` draws each target's 400 control library sizes, and the
    #     panel total is `profile x sum(libraries)`. A first version of this gate expected
    #     `a15` and `a15_seed1` to agree and measured 0/85 — the gate was wrong, the build was
    #     not. It is recorded because it is the reason the twin is a valid noise floor on ALL
    #     SIX members rather than on the DE members alone.
    #   * `cd4f034`, whose exception is measurable: the frozen CD4 bank's `complete mean3`
    #     policy covers 4 of h1-test's 85 targets, so exactly 4 blocks differ and 81 must not.
    n_t = len(next(iter(digests.values()), []))
    by_name = {c["name"]: c for c in CONFIGS}
    expectations = []
    for cfg in wanted:
        name = cfg["name"]
        if name == "base":
            continue
        is_amp = float(cfg.get("effect_scale", 1.0)) != 1.0
        expectations.append(
            (name, "base", 0 if is_amp else n_t,
             "AMPLITUDE arm: the delivered pseudobulk MUST move on every block"
             if is_amp else
             "EMISSION arm: the quasi law repairs every gene total, so the delivered "
             "pseudobulk MUST be bit-identical => `pds` and `mse` are pinned"))
    g4 = {}
    for a, b, want_same, why in expectations:
        if a not in digests or b not in digests:
            continue
        same = sum(x == y for x, y in zip(digests[a], digests[b]))
        ok = same == want_same
        g4[f"{a}_vs_{b}"] = {"blocks_identical": same, "expected_identical": want_same,
                             "of": n_t, "why": why, "pass": bool(ok)}
        print(f"[GATE 4] {a:11s} vs {b:11s}: {same:3d}/{n_t} blocks digest-identical "
              f"(expected {want_same}) {'PASS' if ok else 'FAIL'}   {why}")
    report["gates"]["gate4_panel_total_digest"] = g4
    if g4 and not all(v["pass"] for v in g4.values()):
        raise SystemExit("GATE 4 FAILED — the delivered pseudobulk is not what the arm class "
                         "claims. STOP.")

    # -------------------------------------------- GATE 5: `base` IS 221's `rc_callset`
    #
    # A byte-identical candidate through a deterministic scorer reproduces a recorded row
    # exactly. This is the known-answer row this experiment reads its whole sweep against.
    if "base" in report["configs"]:
        got = report["configs"]["base"]["sha256"]
        ok5 = got == RC_CALLSET_SHA_221
        print(f"[GATE 5] base sha256 {got}\n         221 rc_callset {RC_CALLSET_SHA_221}"
              f"  -> {'PASS' if ok5 else 'FAIL'}")
        report["gates"]["gate5_base_is_221_rc_callset"] = {
            "got": got, "want": RC_CALLSET_SHA_221, "pass": bool(ok5)}
        if not ok5:
            raise SystemExit("GATE 5 FAILED — `base` is not 221's `rc_callset`. STOP.")

    report["peak_rss_gib"] = peak_rss_gib()
    report["seconds_total"] = round(time.time() - t0, 1)
    Path(args.report).write_text(json.dumps(report, indent=1))
    print(f"[221] report -> {args.report}; total {report['seconds_total']}s, "
          f"peak RSS {report['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main()
