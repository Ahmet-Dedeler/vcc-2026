#!/usr/bin/env python3
"""346 step 6 — paired 10,000-target bootstrap of the H1-local overall and all six members.

Same construction as report342.py / report323.py: one fixed set of 10,000 paired target
resamples (seed 20260909), fixed reference PDS gallery and fixed H1-local anchors, with the
whole-panel MSE ratio and its jackknife claim cap recomputed inside every resample. Reads only
this experiment's arm outputs plus 323's reference and anchors; writes only into work/.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BASE = ROOT / "experiments/323-current-family-validation"
for _p in (str(BASE), str(ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import run323 as R            # noqa: E402
import evaluate323 as E       # noqa: E402
import report323 as B         # noqa: E402
from cell_eval2.metrics.delta import _across_pert_budget  # noqa: E402
from cell_eval2.prep import bulk_lognorm_means            # noqa: E402

ORDER = ("cap15", "diag_cross_b1", "same_ridge_b025", "same_ridge_b05", "same_ridge_b1",
         "shuffled_b1")
CONTRASTS = (("same_ridge_b1", "cap15"), ("same_ridge_b05", "cap15"),
             ("same_ridge_b025", "cap15"),
             ("diag_cross_b1", "cap15"), ("shuffled_b1", "cap15"),
             ("same_ridge_b1", "shuffled_b1"), ("same_ridge_b1", "diag_cross_b1"))
# same_ridge_b025 is POST-HOC (chosen from step 9's design-stage curve, after the two
# registered blend weights were already priced). Its interval is not selection-corrected.
POST_HOC = ("same_ridge_b025",)
LOWER_BETTER = ("mse", "nmae")
DRAWS, SEED = 10000, 20260909
OVERALL_FLOOR = 0.002
SIGN_FLOOR = 0.015


def resample_counts(n: int) -> np.ndarray:
    rng = np.random.default_rng(SEED)
    draws = rng.integers(0, n, (DRAWS, n))
    counts = np.zeros((DRAWS, n), dtype=np.int64)
    np.add.at(counts, (np.repeat(np.arange(DRAWS), n), draws.ravel()), 1)
    return np.vstack([np.ones(n, dtype=np.int64), counts])


def arm_raw(name, r, truth, control, own, counts, members):
    metrics = json.loads((WORK / f"metrics-{name}.json").read_text())
    path = WORK / f"stats-{name}.npz"
    if R.sha(path) != metrics["statistics_sha256"]:
        raise ValueError(f"statistics changed: {name}")
    with np.load(path, allow_pickle=False) as z:
        if not np.array_equal(z["targets"], r["targets"]):
            raise ValueError(f"target order drift: {name}")
        pm = bulk_lognorm_means(z["sums"].astype(float), E.TS)
        jk = z["jk"]
    n = len(r["targets"])
    mse, budget, rho, den = B.mse_bootstrap(counts, pm, truth, control, jk,
                                            r["jk"][:-1], r["jk"][-1], own)
    expected_budget = _across_pert_budget(pm, dict(enumerate(map(int, own))),
                                          np.full(n, 1 / (len(r["genes"]) - 1)))
    checks = {"mse_identity_abs": float(abs(mse[0] - metrics["raw"]["mse"])),
              "budget_identity_abs": float(abs(budget[0] - expected_budget)),
              "denominator_min": float(den.min()), "rho_point": float(rho[0]),
              "rho_bootstrap_range": [float(rho[1:].min()), float(rho[1:].max())]}
    if checks["mse_identity_abs"] > 1e-10 or checks["budget_identity_abs"] > 1e-10:
        raise AssertionError(f"bootstrap identity failed: {name}")
    if den.min() <= 0:
        raise ValueError(f"nonpositive MSE denominator in bootstrap: {name}")
    table = pl.read_csv(WORK / f"per-target-{name}.csv")
    raw = np.empty((len(counts), len(members)))
    for j, m in enumerate(members):
        if m == "mse":
            raw[:, j] = mse
            continue
        rows = table.filter(pl.col("metric") == E.MEMBERS[m])
        mapping = dict(zip(rows["perturbation"], rows["value"]))
        values = np.array([mapping.get(str(t), np.nan) for t in r["targets"]], dtype=float)
        valid = np.isfinite(values)
        denominator = counts[:, valid].sum(1)
        if not valid.any() or np.any(denominator <= 0):
            raise ValueError(f"empty cohort in a resample: {name}/{m}")
        raw[:, j] = (counts[:, valid] @ values[valid]) / denominator
        if abs(raw[0, j] - metrics["raw"][m]) > 1e-12:
            raise AssertionError(f"per-target mean identity failed: {name}/{m}")
        checks[f"eligible_targets_{m}"] = int(valid.sum())
    if not np.isfinite(raw).all():
        raise ValueError(f"non-finite resampled member: {name}")
    return metrics, raw, checks


def design_sign(signs, curve, name):
    """The knee arm is not in signs.json; step 9's curve holds its design-stage value."""
    if name in signs["arms"]:
        return signs["arms"][name]["design_w_all"]
    beta = {"same_ridge_b025": 0.25}[name]
    return next(r["design_weighted_sign"] for r in curve["curve"] if r["beta"] == beta)


def contrast(a, b, overall, raw, members):
    delta = overall[a] - overall[b]
    out = {"overall": float(delta[0]), "overall_ci": B.ci(delta[1:]), "raw": {}}
    for j, m in enumerate(members):
        d = raw[a][:, j] - raw[b][:, j]
        entry = {"point": float(d[0]), **B.ci(d[1:])}
        if m in LOWER_BETTER:
            entry["entire_interval_harmful"] = bool(entry["lo"] > 0)
            entry["entire_interval_beneficial"] = bool(entry["hi"] < 0)
        else:
            entry["entire_interval_harmful"] = bool(entry["hi"] < 0)
            entry["entire_interval_beneficial"] = bool(entry["lo"] > 0)
        out["raw"][m] = entry
    return out


def main() -> dict:
    ka = json.loads((WORK / "known-answer.json").read_text())
    if not ka["pass"]:
        raise ValueError("known-answer gate has not passed; modified arms are not interpreted")
    signs = json.loads((WORK / "signs.json").read_text())
    curve = json.loads((WORK / "beta-curve.json").read_text())
    if not signs["known_answer"]["pass"]:
        raise ValueError("343 sign known-answer has not passed")
    with np.load(BASE / "work/reference.npz", allow_pickle=False) as z:
        r = {k: z[k] for k in ("targets", "genes", "sums", "jk")}
    logmeans = bulk_lognorm_means(r["sums"].astype(float), E.TS)
    truth, control = logmeans[:-1], logmeans[-1]
    lookup = {g: i for i, g in enumerate(r["genes"])}
    own = np.array([lookup[t] for t in r["targets"]])
    counts = resample_counts(len(r["targets"]))
    members = list(E.MEMBERS)
    if members != ["pds", "mse", "nmae", "fid", "reach", "jac"]:
        raise AssertionError("member order drift")
    anchors = json.loads((BASE / "work/h1-anchors.json").read_text())
    base = np.array([anchors["baseline"][m] for m in members])
    replicate = np.array([anchors["replicate"][m] for m in members])

    raw, overall, rows, checks = {}, {}, [], {}
    for name in ORDER:
        metrics, raw[name], checks[name] = arm_raw(name, r, truth, control, own, counts, members)
        sc = B.scaled(raw[name], base, replicate)
        overall[name] = sc.mean(1)
        rows.append({"arm": name, **metrics["raw"],
                     "local_h1_overall": float(overall[name][0]),
                     "local_scaled": {m: float(sc[0, j]) for j, m in enumerate(members)},
                     "design_weighted_sign": design_sign(signs, curve, name),
                     "design_weighted_sign_on_covered_axis":
                         signs["arms"].get(name, {}).get("design_w_covered"),
                     "post_hoc": name in POST_HOC,
                     "panel_total_digest_sha256": metrics["panel_total_digest_sha256"],
                     "emitted_count_digest_sha256": metrics["emitted_count_digest_sha256"]})
    contrasts = {f"{a}-minus-{b}": contrast(a, b, overall, raw, members) for a, b in CONTRASTS}

    primary = contrasts["same_ridge_b1-minus-cap15"]
    shuffled = contrasts["shuffled_b1-minus-cap15"]
    sign_gain = (signs["arms"]["same_ridge_b1"]["design_w_all"]
                 - signs["arms"]["cap15"]["design_w_all"])
    share = (shuffled["overall"] / primary["overall"]) if primary["overall"] > 0 else (
        float("inf") if shuffled["overall"] > 0 else 0.0)
    kills = {
        "overall_lower_bound_at_or_below_plus_0_002":
            bool(primary["overall_ci"]["lo"] <= OVERALL_FLOOR),
        "shuffled_control_reproduces_at_least_half_the_gain":
            bool(primary["overall"] <= 0 or shuffled["overall"] >= 0.5 * primary["overall"]),
        "design_sign_gain_below_plus_0_015": bool(sign_gain < SIGN_FLOOR),
    }
    fired = [k for k, v in kills.items() if v]
    nonfinite = sum(c["overall_ci"]["nonfinite_resamples"]
                    + sum(e["nonfinite_resamples"] for e in c["raw"].values())
                    for c in contrasts.values())
    knee = contrasts["same_ridge_b025-minus-cap15"]
    result = {
        "post_hoc_arms": list(POST_HOC),
        "post_hoc_note": "same_ridge_b025 was selected from the design-stage curve AFTER the "
                         "two registered blend weights were priced; its interval carries no "
                         "correction for that selection and it is not part of any kill "
                         "criterion, which are defined on the beta=1 arm only.",
        "knee_overall_gain": knee["overall"], "knee_overall_ci": knee["overall_ci"],
        "status": "REGIME_NOT_WORTH_BUILDING_FOR_A" if fired
                  else "LOCAL_H1_GAIN_REQUIRES_SEPARATE_ABC_REVIEW",
        "kill_criteria": kills, "kill_criteria_fired": fired,
        "primary_overall_gain": primary["overall"],
        "primary_overall_ci": primary["overall_ci"],
        "shuffled_share_of_primary_gain": share,
        "design_weighted_sign_gain": sign_gain,
        "arms": rows, "contrasts": contrasts,
        "identity_checks": checks,
        "nonfinite_resamples": int(nonfinite),
        "known_answer_pass": bool(ka["pass"]),
        "sign_known_answer_pass": bool(signs["known_answer"]["pass"]),
        "bootstrap": {"draws": DRAWS, "seed": SEED, "unit": "target", "paired": True,
                      "fixed": ["reference PDS gallery", "H1-local anchors",
                                "one emission/control seed"],
                      "panel_mse_and_claim_cap_recomputed": True},
        "scope": "85-target reused H1 development panel, cap15 current family, H1 union OFF; "
                 "known-line held-target cross-study calibration; NOT unseen-line "
                 "generalization, NOT an A/B/C score, no board conversion factor",
        "anchors_sha256": R.sha(BASE / "work/h1-anchors.json"),
        "reserve_read": False, "submission": False, "production_default_change": False,
    }
    R.write_json(WORK / "results.json", result)
    print(json.dumps({k: result[k] for k in
                      ("status", "kill_criteria", "kill_criteria_fired",
                       "primary_overall_gain", "primary_overall_ci",
                       "shuffled_share_of_primary_gain", "design_weighted_sign_gain")},
                     indent=2), flush=True)
    for key, value in contrasts.items():
        print(key, json.dumps({"overall": value["overall"], **value["overall_ci"]}), flush=True)
    return result


if __name__ == "__main__":
    main()
