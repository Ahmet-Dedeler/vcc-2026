#!/usr/bin/env python3
"""Paired target-bootstrap comparisons, with MSE recomputed as a panel ratio.

Condition on this reference gallery, frozen H1 anchors and constructed models.
Targets, not genes or cells, are resampled. This is not a new-context CI.
"""
import json

import numpy as np
import pandas as pd
import polars as pl

import evaluate323 as E
import run323 as R
from cell_eval2.metrics.delta import PRED_TRACE_CAP_K, _across_pert_budget
from cell_eval2.prep import bulk_lognorm_means


def mse_bootstrap(multiplicity, predicted, truth, control, pred_jk, real_jk,
                  control_jk, own_columns):
    """Recompute the global claim cap for a weighted panel without a genes x B tensor.

    B = sum_i a_i ||p_i,retained||² - sum_g S_g² / W_g.
    For non-target columns W_g is constant. Compute their sum using the small
    row Gram matrix, then patch the panel-target columns for own-gene exclusion.
    This is the weighted variance in cell-eval2's _across_pert_budget.
    """
    n, g = predicted.shape
    own = np.asarray(own_columns, dtype=int)
    if np.any(own < 0) or len(np.unique(own)) != n:
        raise ValueError("this bootstrap requires one distinct on-axis gene per target")
    a = np.asarray(multiplicity, dtype=float) / (g-1)
    total = a.sum(1)
    self_values = predicted[np.arange(n), own]
    norm2 = (predicted*predicted).sum(1) - self_values**2
    gram = predicted @ predicted.T
    squared_sum = np.einsum("bi,ij,bj->b", a, gram, a, optimize=True)
    panel_sum = a @ predicted[:, own]
    patch_weight = total[:, None]-a
    patch_sum = panel_sum-a*self_values
    patched = np.divide(patch_sum**2, patch_weight,
                        out=np.zeros_like(patch_sum), where=patch_weight > 0).sum(1)
    mean_term = (squared_sum-(panel_sum**2).sum(1))/total + patched
    budget = np.maximum(0., a @ norm2-mean_term)
    cap = np.minimum(pred_jk, PRED_TRACE_CAP_K*real_jk)
    claim = a @ cap
    rho = np.minimum(1., np.divide(budget, claim, out=np.ones_like(budget), where=claim > 0))
    diff = predicted-truth
    sq = (diff*diff).sum(1)-diff[np.arange(n), own]**2
    delta = truth-control
    dsq = (delta*delta).sum(1)-delta[np.arange(n), own]**2
    numerator = a @ (sq-real_jk)-rho*claim
    denominator = a @ (dsq-real_jk-control_jk)
    return numerator/denominator, budget, rho, denominator


def ci(values):
    v = np.asarray(values)
    finite = v[np.isfinite(v)]
    return {"lo": float(np.percentile(finite, 2.5)),
            "hi": float(np.percentile(finite, 97.5)),
            "nonfinite_resamples": int(v.size-finite.size)}


def scaled(raw, base, replicate):
    out = (raw-base)/(replicate-base)
    out[:, 1] = np.clip(out[:, 1], 0., 1.)  # member order below: PDS, MSE, ...
    out[:, 2] = np.maximum(out[:, 2], -6.)
    return out


def main():
    E.verify()
    anchors = json.loads((R.WORK / "h1-anchors.json").read_text())
    r = np.load(R.WORK / "reference.npz", allow_pickle=False)
    targets, genes = r["targets"], r["genes"]
    truth_all = bulk_lognorm_means(r["sums"].astype(float), E.TS)
    truth, control = truth_all[:-1], truth_all[-1]
    real_jk, control_jk = r["jk"][:-1], r["jk"][-1]
    own = np.array([{g: i for i, g in enumerate(genes)}[t] for t in targets])
    members = list(E.MEMBERS)
    n = len(targets)
    # Same 10,000 paired target draws for every arm and every metric.
    rng = np.random.default_rng(20260909)
    draws = rng.integers(0, n, (10000, n))
    counts = np.zeros((10000, n), dtype=np.int64)
    np.add.at(counts, (np.repeat(np.arange(10000), n), draws.ravel()), 1)
    counts = np.vstack([np.ones(n, dtype=np.int64), counts])
    base = np.array([anchors["baseline"][m] for m in members])
    anchor = np.array([anchors["replicate"][m] for m in members])
    all_boot, all_raw, records, gates = {}, {}, [], {}
    for name in [*R.ARMS, "shuffled_cap15"]:
        metrics = json.loads((R.WORK / f"metrics-{name}.json").read_text())
        z = np.load(R.WORK / f"stats-{name}.npz", allow_pickle=False)
        if R.sha(R.WORK / f"stats-{name}.npz") != metrics["statistics_sha256"]:
            raise ValueError(f"statistics changed: {name}")
        pm = bulk_lognorm_means(z["sums"].astype(float), E.TS)
        raw = np.empty((len(counts), len(members)))
        mse, budget, rho, denominator = mse_bootstrap(
            counts, pm, truth, control, z["jk"], real_jk, control_jk, own)
        if abs(mse[0]-metrics["raw"]["mse"]) > 1e-10:
            raise AssertionError(f"bootstrap MSE identity failed: {name}")
        columns = dict(enumerate(map(int, own)))
        expected_budget = _across_pert_budget(pm, columns, np.full(n, 1/(len(genes)-1)))
        if abs(budget[0]-expected_budget) > 1e-10:
            raise AssertionError(f"bootstrap budget identity failed: {name}")
        if np.any(denominator <= 0):
            raise ValueError(f"nonpositive MSE denominator in bootstrap: {name}")
        pt = pl.read_csv(R.WORK / f"per-target-{name}.csv")
        for j, member in enumerate(members):
            if member == "mse":
                raw[:, j] = mse
                continue
            rows = pt.filter(pl.col("metric") == E.MEMBERS[member])
            mapping = dict(zip(rows["perturbation"], rows["value"]))
            values = np.array([mapping.get(str(t), np.nan) for t in targets])
            mask = np.isfinite(values)
            raw[:, j] = (counts[:, mask] @ values[mask])/counts[:, mask].sum(1)
            if abs(raw[0, j]-metrics["raw"][member]) > 1e-12:
                raise AssertionError(f"per-target mean identity failed: {name}/{member}")
        ss = scaled(raw, base, anchor)
        all_boot[name] = ss.mean(1)
        all_raw[name] = raw
        records.append({"arm": name, **metrics["raw"],
                        "local_h1_overall": float(ss[0].mean())})
        gates[name] = {"mse_identity_abs_diff": float(abs(mse[0]-metrics["raw"]["mse"])),
                       "budget_identity_abs_diff": float(abs(budget[0]-expected_budget)),
                       "rho_point": float(rho[0]),
                       "rho_bootstrap_range": [float(rho[1:].min()), float(rho[1:].max())]}
    comparisons = {}
    for parent, child, sign in (("uniform", "plain15", 1), ("plain15", "cap15", 1),
                                 ("plain15", "plain30", -1), ("plain15", "parity15", -1),
                                 ("shuffled_cap15", "cap15", 1)):
        key = f"{child}-minus-{parent}"
        delta = all_boot[child]-all_boot[parent]
        point, interval = float(delta[0]), ci(delta[1:])
        point_agrees = bool(point*sign > 0)
        resolved = interval["lo"] > 0 if sign > 0 else interval["hi"] < 0
        metrics = {}
        for j, member in enumerate(members):
            d = all_raw[child][:, j]-all_raw[parent][:, j]
            metrics[member] = {"point": float(d[0]), **ci(d[1:])}
        comparisons[key] = {
            "local_h1_overall_delta": point, "overall_ci": interval,
            "known_board_expected_sign": sign if parent != "shuffled_cap15" else None,
            "point_sign_agrees": point_agrees,
            "resolved_in_expected_direction": bool(resolved),
            "raw_deltas": metrics,
        }
    nc = comparisons["cap15-minus-shuffled_cap15"]["raw_deltas"]["pds"]
    direction_pass = nc["lo"] > 0
    historical = [v for v in comparisons.values() if v["known_board_expected_sign"] is not None]
    qualified = direction_pass and all(v["resolved_in_expected_direction"] for v in historical)
    report = {
        "status": "QUALIFIED_ON_THIS_PANEL" if qualified else "NOT_QUALIFIED_AS_STANDALONE_SELECTOR",
        "scope": "one reused H1 context, union OFF; not an independent-context performance claim",
        "negative_control_pds_pass": bool(direction_pass),
        "known_board_point_signs_correct": sum(v["point_sign_agrees"] for v in historical),
        "known_board_contrasts_resolved": sum(v["resolved_in_expected_direction"] for v in historical),
        "bootstrap": {"n": 10000, "seed": 20260909, "unit": "target",
                      "fixed": ["reference PDS gallery", "H1 anchors", "constructed model"],
                      "mse": "ratio of resampled sums; global across-perturbation claim cap recomputed"},
        "anchors_sha256": R.sha(R.WORK / "h1-anchors.json"),
        "arms": records, "comparisons": comparisons, "identity_gates": gates,
    }
    pd.DataFrame(records).to_csv(R.WORK / "scores.csv", index=False)
    R.write_json(R.WORK / "results.json", report)
    print(json.dumps({k: report[k] for k in ("status", "negative_control_pds_pass",
          "known_board_point_signs_correct", "known_board_contrasts_resolved")}, indent=2))


if __name__ == "__main__":
    main()
