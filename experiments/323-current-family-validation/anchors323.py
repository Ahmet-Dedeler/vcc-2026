#!/usr/bin/env python3
"""Reconstruct H1-local scaling anchors with official constructors and split rules.

These are local H1 anchors, never A/B/C anchors and never a server attestation.
Synthetic full-pipeline parity is required before reading H1 for this stage.
"""
from dataclasses import replace
import gc
import json
from pathlib import Path
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import scipy.sparse as sp

import evaluate323 as E
import exact_de323 as D
import run323 as R
from cell_eval2.anchor import _derive_seeds, compute_replicate_anchor
from cell_eval2.baseline import generic_response_profile, build_baseline_prediction, baseline_config
from cell_eval2.ceiling import _disjoint_halves
from cell_eval2.config import EvalConfig
from cell_eval2.de import prepare_de, resolve_target_genes
from cell_eval2.lfc_nmae_ref import _assert_disjoint_controls, _nmae_ref_from_tables
from cell_eval2.metrics.discrimination import discrimination_score
from cell_eval2.moments import jackknife_correction
from cell_eval2.prep import bulk_lognorm_means
from cell_eval2.run import compute_metrics, aggregate_metrics, dispatch_de_metrics


def stats(a, *, de_control=None):
    lab = a.obs.target_gene.to_numpy().astype(str)
    genes = np.asarray(a.var_names, dtype=str)
    perts = np.array(sorted(set(lab)), dtype=str)
    ctrl = D.control_block(a.X[lab == E.CTRL]) if de_control is None else de_control
    sums, jk, counts, ps, qs, lfcs, targets = [], [], [], [], [], [], []
    for target in perts:
        rows = a.X[lab == target].astype(np.float64)
        sums.append(np.asarray(rows.sum(0)).ravel())
        counts.append(rows.shape[0])
        jk.append(jackknife_correction(rows, np.zeros(rows.shape[0], dtype=np.int64),
                                      1, E.TS, chunk=32)[0])
        if target != E.CTRL:
            _, p, q, lfc = D.de_from_counts(rows, ctrl)
            ps.append(p); qs.append(q); lfcs.append(lfc); targets.append(target)
    de = E.de_frame(np.array(targets), genes[ctrl.keep], ps, qs, lfcs)
    return {"perts": perts, "genes": genes, "sums": np.array(sums),
            "counts": np.array(counts), "jk": np.array(jk), "de": de}


def fast_score(pred, real, *, own_controls=False):
    if not np.array_equal(pred["genes"], real["genes"]):
        raise ValueError("anchor gene-axis mismatch")
    pm = bulk_lognorm_means(pred["sums"], E.TS)
    rm = bulk_lognorm_means(real["sums"], E.TS)
    pds = discrimination_score(
        pred_bulk=(pred["perts"], pm), real_bulk=(real["perts"], rm),
        genes=real["genes"], pert_col="target_gene", control=E.CTRL,
        distance="cosine", rank_denominator="n-1", tie_policy="midrank",
        exclude_target_gene=True, exclusion_scope="panel",
        control_source="pred" if own_controls else "real")
    mk = R.load_module("msekit_anchor323", R.ROOT / "experiments/096-mse-floor/msekit.py")
    mse, _, _ = mk.score(
        pred["perts"], pred["sums"], pred["counts"], pred["jk"],
        real["perts"], real["sums"], real["counts"], real["jk"], real["genes"])
    prepared = prepare_de(pred["de"], real["de"], control=E.CTRL,
                          sort_by="abs_log2_fold_change", p_adj_threshold=.05)
    cfg = EvalConfig.from_preset("vcc2026")
    rows = dispatch_de_metrics([E.MEMBERS[m] for m in ("nmae", "fid", "reach", "jac")],
                               prepared, cfg)
    raw = {"mse": float(mse), "pds": float(np.mean(list(pds.values())))}
    for key in ("nmae", "fid", "reach", "jac"):
        raw[key] = float(np.nanmean([row["value"] for row in rows
                                    if row["metric"] == E.MEMBERS[key]]))
    return raw


def replicate(real_ad, full, *, n_splits=5):
    rows = []
    resolution = resolve_target_genes(
        full["de"], sorted(full["de"]["target"].unique().to_list()))
    for i, seed in enumerate(_derive_seeds(0, n_splits)):
        a, b = _disjoint_halves(real_ad, "target_gene", E.CTRL, seed)
        _assert_disjoint_controls(a, b, pert_col="target_gene", control=E.CTRL)
        sa, sb = stats(a), stats(b)
        raw = fast_score(sb, sa, own_controls=True)
        full_gate = _nmae_ref_from_tables(
            full["de"], sa["de"], sb["de"], p_adj_threshold=.05,
            min_gate_size=10, target_resolution=resolution)
        raw["nmae"] = float(full_gate["nmae_ref_raw"].mean())
        rows.append({"split": i, "seed": seed, "raw": raw,
                     "nmae_full_gate_cohort": full_gate.height})
        print(f"anchor split {i+1}/{n_splits}: {raw}", flush=True)
        del a, b, sa, sb
        gc.collect()
    return rows, {m: float(np.mean([r["raw"][m] for r in rows])) for m in E.MEMBERS}


def fixture():
    rng = np.random.default_rng(323)
    arrays = [rng.poisson(8, (80, 40)).astype(np.float32)]
    labels = [E.CTRL]*80
    for i in range(3):
        x = rng.poisson(8, (60, 40)).astype(np.float32)
        x[:, 3+8*i:3+8*(i+1)] *= 5
        arrays.append(x)
        labels += [f"G{i}"]*60
    return ad.AnnData(
        X=sp.csr_matrix(np.vstack(arrays)),
        obs=pd.DataFrame({"target_gene": labels}, index=[f"c{i}" for i in range(len(labels))]),
        var=pd.DataFrame(index=[f"G{i}" for i in range(40)]))


def selftest():
    a = fixture()
    cfg = replace(EvalConfig.from_preset("vcc2026"), pert_col="target_gene",
                  de=replace(EvalConfig.from_preset("vcc2026").de, backend="scanpy"))
    full = stats(a)
    p = build_baseline_prediction(generic_response_profile(
        a, pert_col="target_gene", control=E.CTRL),
        a, pert_col="target_gene", control=E.CTRL, emit="dispersed", seed=0)
    fast_base = fast_score(stats(p), full)
    tidy = compute_metrics(p, a, config=baseline_config(cfg))
    official = aggregate_metrics(tidy)
    expected = dict(zip(official["metric"], official["mean"]))
    bd = {m: abs(fast_base[m]-expected[k]) for m, k in E.MEMBERS.items()}
    _, ra = replicate(a, full)
    _, official_anchor = compute_replicate_anchor(a, config=cfg, base_seed=0, n_splits=5)
    expected_a = dict(zip(official_anchor["metric"], official_anchor["replicate"]))
    rd = {m: abs(ra[m]-expected_a[k]) for m, k in E.MEMBERS.items()}
    passed = all(bd[m] <= E.TOLERANCE[m] and rd[m] <= E.TOLERANCE[m] for m in E.MEMBERS)
    report = {"pass": passed, "baseline_abs_diff": bd, "replicate_abs_diff": rd,
              "tolerances": E.TOLERANCE, "scope": "synthetic full-pipeline parity"}
    R.write_json(R.WORK / "anchor-selftest.json", report)
    if not passed:
        raise AssertionError(report)


def freeze():
    files = [Path(__file__), R.HERE / "ANCHOR-NOTE.md",
             R.WORK / "anchor-selftest.json", R.WORK / "metrics-known.json"]
    import cell_eval2
    package = Path(cell_eval2.__file__).parent
    files += [package / n for n in ("anchor.py", "ceiling.py", "baseline.py",
                                    "lfc_nmae_ref.py")]
    R.write_json(R.HERE / "anchor-inputs.json",
                 {"inputs_sha256": {str(p.relative_to(R.ROOT)): R.sha(p) for p in files}})


def run():
    E.verify()
    for path, digest in json.loads((R.HERE / "anchor-inputs.json").read_text())["inputs_sha256"].items():
        if R.sha(R.ROOT / path) != digest:
            raise ValueError(f"anchor input drift: {path}")
    if not json.loads((R.WORK / "anchor-selftest.json").read_text())["pass"]:
        raise ValueError("synthetic anchor parity failed")
    t0 = time.monotonic()
    a = ad.read_h5ad(R.REAL)
    full = stats(a)
    p = build_baseline_prediction(generic_response_profile(
        a, pert_col="target_gene", control=E.CTRL),
        a, pert_col="target_gene", control=E.CTRL, emit="dispersed", seed=0)
    # The official baseline constructor copies real controls, so each side's
    # control is identical. Fractional counts are permitted ONLY for baseline.
    base = fast_score(stats(p), full)
    del p
    gc.collect()
    print(f"local H1 baseline: {base}", flush=True)
    splits, anchor = replicate(a, full)
    for m in E.MEMBERS:
        if not np.isfinite(base[m]) or not np.isfinite(anchor[m]) or abs(anchor[m]-base[m]) < 1e-12:
            raise ValueError(f"degenerate H1 scale for {m}")
    R.write_json(R.WORK / "h1-anchors.json", {
        "population": "H1 hESC, 85 targets, this reference only",
        "reference_sha256": R.sha(R.REAL), "baseline": base, "replicate": anchor,
        "replicate_splits": splits, "base_seed": 0, "n_splits": 5,
        "nmae_estimator": "full_gate_raw", "other_estimator": "split_half_raw",
        "baseline_constructor": "cell_eval2.baseline dispersed, seed 0",
        "classification": "locally reconstructed anchors; not A/B/C and not a server bundle",
        "seconds": time.monotonic()-t0,
    })


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode not in ("selftest", "freeze", "run"):
        raise ValueError(mode)
    globals()[mode]()
