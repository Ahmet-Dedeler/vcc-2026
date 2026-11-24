#!/usr/bin/env python3
"""Score actual emitted count blocks without repeatedly sorting 18,400 controls.

Production rc_emit supplies the pre-sorted-control Wilcoxon kernel; cell-eval2
supplies PDS, jackknife MSE, DE preparation and all four DE metrics. No
profile-only surrogate or fitted board conversion is used. A recorded emitted
panel digest AND six-metric row must reproduce before new arms are scored.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import subprocess
import sys
import time

# Parallel threads within ONE memory-guarded process, never parallel scorers.
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("POLARS_MAX_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import h5py
import numpy as np
import polars as pl
import scipy.sparse as sp

import run323 as R
import exact_de323
from cell_eval2.config import EvalConfig
from cell_eval2.de import prepare_de
from cell_eval2.metrics.discrimination import discrimination_score
from cell_eval2.moments import jackknife_correction
from cell_eval2.prep import bulk_lognorm_means
from cell_eval2.run import dispatch_de_metrics
from vcc2026 import rc_emit
from vcc2026.donor_transfer import emit_counts

MEMBERS = {
    "pds": "pds_cosine",
    "mse": "expr_mse_unbiased_capped_norm",
    "nmae": "de_wilcoxon_lfc_nmae",
    "fid": "de_wilcoxon_direction_fidelity_yield_raw",
    "reach": "de_wilcoxon_direction_reach_raw",
    "jac": "de_wilcoxon_sig_jaccard",
}
CTRL = "non-targeting"
TS = 50000.
HIST_EFFECT = R.ROOT / "experiments/170-gamma-mse/work/eff/h1_g0.65.npz"
HIST_BUILD = R.ROOT / "experiments/221-cloud-scorer/cloud-output/ladder5-20260904/prologue-out/build221.json"
HIST_RUN = R.ROOT / "experiments/221-cloud-scorer/cloud-output/ladder8-20260904/rc_callset/run"
TOLERANCE = {"pds": 1e-12, "mse": 1e-10, "nmae": 1e-6,
             "fid": 1e-6, "reach": 1e-6, "jac": 1e-6}


def freeze():
    path = R.HERE / "evaluation-inputs.json"
    if path.exists():
        raise FileExistsError(path)
    files = [Path(__file__), R.HERE / "EVALUATION-NOTE.md", HIST_EFFECT,
             HIST_BUILD, HIST_RUN / "agg_results.csv", HIST_RUN / "results.csv",
             R.ROOT / "experiments/096-mse-floor/msekit.py",
             R.WORK / "effects.npz", R.WORK / "compose.json"]
    import cell_eval2
    package = Path(cell_eval2.__file__).parent
    files += [package / f for f in ("moments.py", "prep.py", "run.py", "de.py",
              "config.py", "catalog.py", "configs/vcc2026.yaml",
              "metrics/delta.py", "metrics/discrimination.py",
              "metrics/de.py", "metrics/direction.py")]
    R.write_json(path, {
        "cell_eval2_version": version("cell-eval2"),
        "inputs_sha256": {str(p.relative_to(R.ROOT)): R.sha(p) for p in files},
        "known_answer_absolute_tolerances": TOLERANCE,
    })


def verify():
    R.verify_inputs()
    p = R.HERE / "evaluation-inputs.json"
    frozen = json.loads(p.read_text())
    for file, digest in frozen["inputs_sha256"].items():
        if R.sha(R.ROOT / file) != digest:
            raise ValueError(f"evaluation input drift: {file}")
    if version("cell-eval2") != frozen["cell_eval2_version"]:
        raise ValueError("scorer version drift")
    committed = subprocess.check_output(
        ["git", "log", "-1", "--format=%H", "--", str(p)],
        cwd=R.ROOT, text=True).strip()
    if not committed:
        raise ValueError("commit evaluation code and input manifest first")


def labels(h, path):
    obj = h[path]
    if isinstance(obj, h5py.Group):
        if "categories" in obj:
            cats = obj["categories"].asstr()[:]
            codes = obj["codes"][:]
            if np.any(codes < 0):
                raise ValueError(f"missing categorical label at {path}")
            return cats[codes].astype(str)
        if "values" in obj:
            if "mask" in obj and np.any(obj["mask"][:]):
                raise ValueError(f"missing nullable-string label at {path}")
            return obj["values"].asstr()[:].astype(str)
        raise ValueError(f"unsupported string encoding at {path}")
    return obj.asstr()[:].astype(str)


def read_real():
    with h5py.File(R.REAL, "r") as h:
        lab = labels(h, "obs/target_gene")
        genes = labels(h, "var/_index")
        x = sp.csr_matrix((h["X/data"][:], h["X/indices"][:], h["X/indptr"][:]),
                          shape=tuple(h["X"].attrs["shape"]))
    return lab, genes, x


def block_stats(block, control):
    dense = np.asarray(block.toarray() if sp.issparse(block) else block,
                       dtype=np.float64)
    if not np.isfinite(dense).all() or np.any(dense < 0) or np.any(dense != np.floor(dense)):
        raise ValueError("non-integer or invalid count block")
    if np.any(dense.sum(1) > 1e6):
        raise ValueError("counts per cell exceed contract")
    sums = dense.sum(0).astype(np.int64)
    jk = jackknife_correction(sp.csr_matrix(dense),
                              np.zeros(len(dense), dtype=np.int64), 1, TS, chunk=32)[0]
    _, p, q, lfc = exact_de323.de_from_counts(dense, control)
    return sums, float(jk), p, q, lfc


def de_frame(targets, genes, p, q, lfc):
    return pl.DataFrame({
        "target": np.repeat(targets, len(genes)),
        "feature": np.tile(genes, len(targets)),
        "p_value": np.asarray(p).ravel(),
        "p_adj": np.asarray(q).ravel(),
        "log2_fold_change": np.asarray(lfc).ravel(),
    })


def reference():
    verify()
    t0 = time.monotonic()
    lab, genes, x = read_real()
    targets = np.array(sorted(set(lab) - {CTRL}))
    expected = R.axes()
    if not np.array_equal(genes, expected[2]) or not np.array_equal(targets, expected[1]):
        raise ValueError("reference axes differ from frozen construction")
    xc = x[lab == CTRL]
    if xc.shape[0] != 18400:
        raise ValueError("unexpected control count")
    emission_control = rc_emit.control_block_from_counts(xc)
    c = exact_de323.control_block(xc)
    if not np.array_equal(c.keep, emission_control.keep):
        raise AssertionError("emitter and scorer CPM gate differ")
    cp = np.asarray(xc.astype(np.float64).sum(0)).ravel().astype(np.int64)
    cj = jackknife_correction(xc, np.zeros(xc.shape[0], dtype=np.int64), 1, TS, chunk=32)[0]
    sums, jks, ps, qs, lfcs, counts = [], [], [], [], [], []
    for i, target in enumerate(targets):
        rows = x[lab == target]
        if rows.shape[0] != 400:
            raise ValueError(f"{target}: not 400 reference cells")
        s, j, p, q, l = block_stats(rows, c)
        sums.append(s); jks.append(j); ps.append(p); qs.append(q); lfcs.append(l)
        counts.append(rows.shape[0])
        if (i + 1) % 10 == 0:
            print(f"reference {i + 1}/85 ({time.monotonic()-t0:.1f}s)", flush=True)
    np.savez(R.WORK / "reference.npz", genes=genes, targets=targets,
             perts=np.append(targets, CTRL), sums=np.vstack([*sums, cp]),
             counts=np.array([*counts, xc.shape[0]]), jk=np.array([*jks, cj]),
             p=ps, q=qs, lfc=lfcs, keep=c.keep, cs=c.cs, tie=c.tie,
             ref_lin=c.ref_lin, emit_cs=emission_control.cs,
             emit_tie=emission_control.tie, emit_ref_lin=emission_control.ref_lin,
             libraries=np.asarray(xc.sum(1)).ravel().astype(float),
             b_ctrl=np.log1p(TS * cp / cp.sum()))
    frame = de_frame(targets, genes[c.keep], ps, qs, lfcs)
    frame.write_parquet(R.WORK / "reference_de.parquet")
    own = np.equal.outer(targets, genes[c.keep])
    n_real = ((np.array(qs) < .05) & ~own).sum(1)
    R.write_json(R.WORK / "reference.json", {
        "reference_sha256": R.sha(R.REAL),
        "statistics_sha256": R.sha(R.WORK / "reference.npz"),
        "de_sha256": R.sha(R.WORK / "reference_de.parquet"),
        "n_targets": len(targets), "n_controls": xc.shape[0],
        "n_real_quantiles": np.percentile(n_real, [0, 10, 25, 50, 75, 90, 100]).tolist(),
        "n_real_by_target": dict(zip(map(str, targets), map(int, n_real))),
        "seconds": time.monotonic()-t0,
    })


def load_reference():
    info = json.loads((R.WORK / "reference.json").read_text())
    if R.sha(R.WORK / "reference.npz") != info["statistics_sha256"]:
        raise ValueError("reference statistics changed")
    if R.sha(R.WORK / "reference_de.parquet") != info["de_sha256"]:
        raise ValueError("reference DE changed")
    z = np.load(R.WORK / "reference.npz", allow_pickle=False)
    r = {k: z[k] for k in z.files}
    z.close()
    c = rc_emit.RcControlBlock(keep=r.pop("keep"), cs=r.pop("cs"),
                               tie=r.pop("tie"), ref_lin=r.pop("ref_lin"),
                               n_controls=int(r["counts"][-1]), min_cpm=5)
    emit_c = rc_emit.RcControlBlock(
        keep=c.keep, cs=r.pop("emit_cs"), tie=r.pop("emit_tie"),
        ref_lin=r.pop("emit_ref_lin"), n_controls=c.n_controls, min_cpm=5)
    return r, c, emit_c


def score(r, c, sums, jks, p, q, lfc):
    pred_sums = np.vstack([sums, r["sums"][-1]])
    pred_jks = np.append(jks, r["jk"][-1])
    pm = bulk_lognorm_means(pred_sums.astype(float), TS)
    rm = bulk_lognorm_means(r["sums"].astype(float), TS)
    pd = discrimination_score(
        pred_bulk=(r["perts"], pm), real_bulk=(r["perts"], rm),
        genes=r["genes"], pert_col="target_gene", control=CTRL,
        distance="cosine", rank_denominator="n-1", tie_policy="midrank",
        exclude_target_gene=True, exclusion_scope="panel", control_source="real")
    msekit = R.load_module("msekit323", R.ROOT / "experiments/096-mse-floor/msekit.py")
    mse, num, den = msekit.score(
        r["perts"], pred_sums, r["counts"], pred_jks,
        r["perts"], r["sums"], r["counts"], r["jk"], r["genes"])
    pred = de_frame(r["targets"], r["genes"][c.keep], p, q, lfc)
    real = pl.read_parquet(R.WORK / "reference_de.parquet")
    prepared = prepare_de(pred, real, control=CTRL, sort_by="abs_log2_fold_change",
                          p_adj_threshold=.05, nan_lfc_policy="mask")
    cfg = EvalConfig.from_preset("vcc2026")
    rows = dispatch_de_metrics([MEMBERS[m] for m in ("nmae", "fid", "reach", "jac")],
                               prepared, cfg)
    rows += [{"perturbation": t, "metric": MEMBERS["pds"], "value": v}
             for t, v in pd.items()]
    raw = {"mse": mse}
    for short, full in MEMBERS.items():
        if short != "mse":
            vals = [row["value"] for row in rows if row["metric"] == full]
            raw[short] = float(np.nanmean(vals))
    return raw, rows, {"numerator": num, "denominator": den}


def emit_and_score(name, r, c, emit_c, effects, *, historical=False):
    b = R.builder()
    if historical:
        # 221's historical adapter explicitly carried minnorm. This mode is only
        # a known-answer row; current arms retain production MASS_FIX=None.
        b.MASS_FIX = "minnorm"
    pick_rng, emit_rng = b.make_rngs(0)
    sums, jks, ps, qs, lfcs, digests = [], [], [], [], [], []
    counts_sha = hashlib.sha256()
    t0 = time.monotonic()
    for i, target in enumerate(r["targets"]):
        profile = b.profile_from_effect(effects[i], r["b_ctrl"], min_cpm=10.)
        picks = pick_rng.choice(len(r["libraries"]), 400, replace=False)
        spread = b.emission_spread_vector(
            effects[i], r["b_ctrl"], min_cpm=10., tail_spread=1., keep_top=8000,
            keep_top_frac=None, keep_spread=.25, rank_score=None, spread_map="block",
            spread_lo=0., spread_hi=.60)
        gate = b.emission_rank_axis(effects[i], r["b_ctrl"], min_cpm=10.)[0]
        block, _ = rc_emit.rc_hybrid_block(
            profile, r["libraries"][picks], .1, emit_rng, spread=spread,
            overdispersion=4., control=emit_c, key=np.where(gate, effects[i], 0.),
            jitter_seed=1000+i, call_set_first=True,
            call_set=rc_emit.CallSet(400, 2200, 1.9), emit_counts=emit_counts)
        s, j, p, q, l = block_stats(block, c)
        sums.append(s); jks.append(j); ps.append(p); qs.append(q); lfcs.append(l)
        digests.append(hashlib.sha256(np.ascontiguousarray(s).tobytes()).hexdigest())
        counts_sha.update(np.asarray(block, dtype="<u4", order="C").tobytes())
        if (i+1) % 10 == 0:
            print(f"{name} {i+1}/85 ({time.monotonic()-t0:.1f}s)", flush=True)
    panel_sha = hashlib.sha256("".join(digests).encode()).hexdigest()
    if historical:
        expected = json.loads(HIST_BUILD.read_text())["configs"]["rc_callset"]["panel_total_digest_sha256"]
        if panel_sha != expected:
            raise AssertionError(f"historical emitted panel totals mismatch: {panel_sha} != {expected}")
    raw, rows, legs = score(r, c, sums, jks, ps, qs, lfcs)
    np.savez(R.WORK / f"stats-{name}.npz", targets=r["targets"], sums=sums,
             jk=jks, p=ps, q=qs, lfc=lfcs)
    pl.DataFrame(rows).write_csv(R.WORK / f"per-target-{name}.csv")
    out = {"raw": raw, "panel_total_digest_sha256": panel_sha,
           "emitted_count_digest_sha256": counts_sha.hexdigest(),
           "mse_legs": legs, "seconds": time.monotonic()-t0,
           "statistics_sha256": R.sha(R.WORK / f"stats-{name}.npz")}
    if historical:
        expected = pl.read_csv(HIST_RUN / "agg_results.csv").filter(
            pl.col("statistic") == "mean").to_dicts()[0]
        differences = {m: abs(raw[m]-expected[k]) for m, k in MEMBERS.items()}
        out["known_answer"] = {
            "expected_raw": {m: expected[k] for m, k in MEMBERS.items()},
            "absolute_differences": differences, "tolerance": TOLERANCE,
            "pass": all(differences[m] <= TOLERANCE[m] for m in MEMBERS)}
    R.write_json(R.WORK / f"metrics-{name}.json", out)
    print(name, json.dumps(out["raw"]), flush=True)
    if historical and not out["known_answer"]["pass"]:
        raise AssertionError(f"known-answer metric gate failed: {out['known_answer']}")
    return out


def evaluate(only):
    verify()
    r, c, emit_c = load_reference()
    if only == "known":
        with np.load(HIST_EFFECT, allow_pickle=False) as z:
            effects = R.map_axis(z["effect"], R.axes()[0], r["genes"])
        emit_and_score("known", r, c, emit_c, effects, historical=True)
        return
    ka = json.loads((R.WORK / "metrics-known.json").read_text())
    if not ka["known_answer"]["pass"]:
        raise ValueError("known-answer gate has not passed")
    wanted = [*R.ARMS, "shuffled_cap15"] if only == "all" else [only]
    with np.load(R.WORK / "effects.npz", allow_pickle=False) as z:
        for name in wanted:
            if (R.WORK / f"metrics-{name}.json").exists():
                raise FileExistsError(f"refusing to rescore {name}")
            effect = R.map_axis(z[name], z["genes"], r["genes"])
            emit_and_score(name, r, c, emit_c, effect)
            gc.collect()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", choices=["freeze", "reference", "known", "all",
                                    *R.ARMS, "shuffled_cap15"])
    args = ap.parse_args()
    if args.mode in ("freeze", "reference"):
        globals()[args.mode]()
    else:
        evaluate(args.mode)


if __name__ == "__main__":
    main()
