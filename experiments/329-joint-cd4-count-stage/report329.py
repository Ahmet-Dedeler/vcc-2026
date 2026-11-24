"""Paired count-stage readout; panel MSE and its claim cap recomputed per draw."""

import json

import numpy as np
import pandas as pd
import polars as pl

import run329 as A
from estimator329 import geometry
import report323 as B
from cell_eval2.metrics.delta import _across_pert_budget
from cell_eval2.prep import bulk_lognorm_means


def main():
    A.verify()
    if (A.WORK / "results.json").exists():
        raise FileExistsError("terminal report already exists")
    if not json.loads((A.WORK / "known-answer.json").read_text())["pass"]:
        raise AssertionError("known-answer failure")
    with np.load(A.BASE / "work/reference.npz", allow_pickle=False) as z:
        r = {k: z[k] for k in z.files}
    targets, genes = r["targets"], r["genes"]
    truth_all = bulk_lognorm_means(r["sums"].astype(float), A.E.TS)
    truth, control = truth_all[:-1], truth_all[-1]
    lookup = {g: i for i, g in enumerate(genes)}
    own = np.array([lookup[t] for t in targets])
    mask, n = ~np.isin(genes, targets), len(targets)
    rng = np.random.default_rng(20260909)
    draws = rng.integers(0, n, (10000, n))
    counts = np.zeros((10000, n), dtype=np.int64)
    np.add.at(counts, (np.repeat(np.arange(10000), n), draws.ravel()), 1)
    counts = np.vstack((np.ones(n, dtype=np.int64), counts))
    members = list(A.E.MEMBERS)
    if members != ["pds", "mse", "nmae", "fid", "reach", "jac"]:
        raise ValueError("metric order changed")
    anchors = json.loads((A.BASE / "work/h1-anchors.json").read_text())
    base = np.array([anchors["baseline"][m] for m in members])
    replicate = np.array([anchors["replicate"][m] for m in members])
    raw_by_arm, total_by_arm, records, identity, emitted = {}, {}, [], {}, {}
    baseline_pb = None
    for name in A.NAMES:
        metrics = json.loads((A.WORK / f"metrics-{name}.json").read_text())
        stats_path = A.WORK / f"stats-{name}.npz"
        A.A.require_same("statistics hash", A.R.sha(stats_path), metrics["statistics_sha256"])
        with np.load(stats_path, allow_pickle=False) as z:
            if not np.array_equal(z["targets"], targets):
                raise ValueError("statistics target order")
            pm = bulk_lognorm_means(z["sums"].astype(float), A.E.TS)
            jk = z["jk"]
        if baseline_pb is None:
            baseline_pb = pm
        emitted[name] = geometry(pm-control, baseline_pb-control, mask)
        mse, budget, rho, denominator = B.mse_bootstrap(
            counts, pm, truth, control, jk, r["jk"][:-1], r["jk"][-1], own)
        mse_error = abs(mse[0]-metrics["raw"]["mse"])
        expected_budget = _across_pert_budget(
            pm, dict(enumerate(map(int, own))), np.full(n, 1/(len(genes)-1)))
        budget_error = abs(budget[0]-expected_budget)
        if mse_error > 1e-10 or budget_error > 1e-10 or np.any(denominator <= 0):
            raise AssertionError(f"{name}: MSE identity or denominator failure")
        raw = np.empty((len(counts), len(members)))
        table = pl.read_csv(A.WORK / f"per-target-{name}.csv")
        for j, member in enumerate(members):
            if member == "mse":
                raw[:, j] = mse
            else:
                rows = table.filter(pl.col("metric") == A.E.MEMBERS[member])
                mapping = dict(zip(rows["perturbation"], rows["value"]))
                values = np.array([mapping.get(str(t), np.nan) for t in targets])
                valid = np.isfinite(values)
                denominator_de = counts[:, valid].sum(1)
                if np.any(denominator_de <= 0):
                    raise ValueError("empty resampled DE cohort")
                raw[:, j] = counts[:, valid] @ values[valid] / denominator_de
            if abs(raw[0, j]-metrics["raw"][member]) > 1e-12:
                raise AssertionError(f"{name}/{member}: aggregate identity")
        if not np.isfinite(raw).all():
            raise ValueError("nonfinite bootstrap")
        scaled = B.scaled(raw, base, replicate)
        raw_by_arm[name], total_by_arm[name] = raw, scaled.mean(1)
        records.append({"arm": name, "raw": metrics["raw"],
                        "local_h1_overall": float(scaled[0].mean()),
                        "local_scaled": dict(zip(members, map(float, scaled[0])))})
        identity[name] = {
            "mse_identity_abs_difference": float(mse_error),
            "budget_identity_abs_difference": float(budget_error),
            "rho_range": [float(rho.min()), float(rho.max())]}
    comparisons = {}
    for other in ("cap15", "diagonal", "shuffled_fit", "shuffled_targets"):
        delta = total_by_arm["joint"]-total_by_arm[other]
        raw = raw_by_arm["joint"]-raw_by_arm[other]
        comparisons["joint-minus-"+other] = {
            "local_h1_overall_delta": float(delta[0]), "overall_ci": B.ci(delta[1:]),
            "raw_deltas": {m: {"point": float(raw[0, j]), **B.ci(raw[1:, j])}
                           for j, m in enumerate(members)}}
    primary = comparisons["joint-minus-cap15"]
    gates = {
        "overall_beats_cap15_interval": primary["overall_ci"]["lo"] > 0,
        "overall_beats_diagonal_interval": comparisons["joint-minus-diagonal"]["overall_ci"]["lo"] > 0,
        "overall_beats_shuffled_fit_interval": comparisons["joint-minus-shuffled_fit"]["overall_ci"]["lo"] > 0,
        "no_resolved_pds_harm": primary["raw_deltas"]["pds"]["hi"] >= 0,
        "no_resolved_raw_mse_harm": primary["raw_deltas"]["mse"]["lo"] <= 0,
        "no_resolved_nmae_harm": primary["raw_deltas"]["nmae"]["lo"] <= 0,
        "target_negative_control": comparisons["joint-minus-shuffled_targets"]["raw_deltas"]["pds"]["point"] > .10,
    }
    result = {
        "status": "EARNS_FURTHER_DEVELOPMENT_ONLY" if all(gates.values()) else "DOES_NOT_EARN_FURTHER_DEVELOPMENT",
        "scope": "one reused H1 context; no independent confirmation or A/B/C clearance",
        "decision_gates": gates,
        "bootstrap": {"draws": 10000, "seed": 20260909, "unit": "target",
                      "fixed": ["gallery", "anchors", "models"],
                      "mse": "panel ratio and global claim cap recomputed"},
        "arms": records, "comparisons": comparisons,
        "emitted_non_target_pseudobulk_geometry": emitted, "identity_gates": identity,
        "input_manifest_sha256": A.R.sha(A.HERE / "inputs.json"),
        "effects_sha256": A.R.sha(A.WORK / "effects.npz"),
        "anchors_sha256": A.R.sha(A.BASE / "work/h1-anchors.json"),
        "submission": False,
    }
    A.R.write_json(A.WORK / "results.json", result)
    pd.DataFrame([{"arm": r["arm"], **r["raw"], "local_h1_overall": r["local_h1_overall"]}
                  for r in records]).to_csv(A.WORK / "scores.csv", index=False)
    print(json.dumps({k: result[k] for k in ("status", "decision_gates", "comparisons")}, indent=2))


if __name__ == "__main__":
    main()
