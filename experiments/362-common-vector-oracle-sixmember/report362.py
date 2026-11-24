#!/usr/bin/env python3
"""362 step 3 — paired all-six-member analysis: five arms minus cap15.

Same paired 10,000-target resamples (seed 20260909) as 323/337/342, fixed reference
gallery, fixed H1-local anchors, whole-panel MSE ratio with the jackknife claim cap
recomputed in every resample. No member is assumed invariant and no board
conversion factor is applied.
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
import run323 as R  # noqa: E402
import evaluate323 as E  # noqa: E402
import report323 as B  # noqa: E402
from cell_eval2.metrics.delta import _across_pert_budget  # noqa: E402
from cell_eval2.prep import bulk_lognorm_means  # noqa: E402

ORDER = ("cap15", "oracle_a05", "oracle_a1", "oracle_a1_normpreserve",
         "generic_a1", "shuffled_a1")
CONTRASTS = (("oracle_a05", "cap15"), ("oracle_a1", "cap15"),
             ("oracle_a1_normpreserve", "cap15"), ("generic_a1", "cap15"),
             ("shuffled_a1", "cap15"))
LOWER_BETTER = ("mse", "nmae")
DRAWS, SEED = 10000, 20260909


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
    if checks["mse_identity_abs"] > 1e-10:
        raise AssertionError(f"bootstrap MSE identity failed: {name}")
    if checks["budget_identity_abs"] > 1e-10:
        raise AssertionError(f"bootstrap budget identity failed: {name}")
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
    if not np.isfinite(raw).all():
        raise ValueError(f"non-finite resampled member: {name}")
    return metrics, raw, checks


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
    if (WORK / "results.json").exists():
        raise FileExistsError("terminal report exists; preserved")
    ka = json.loads((WORK / "known-answer.json").read_text())
    if not ka["pass"]:
        raise ValueError("known-answer gate has not passed; changed arms are not interpreted")
    with np.load(BASE / "work/reference.npz", allow_pickle=False) as z:
        r = {k: z[k] for k in ("targets", "genes", "sums", "jk")}
    logmeans = bulk_lognorm_means(r["sums"].astype(float), E.TS)
    truth, control = logmeans[:-1], logmeans[-1]
    lookup = {g: i for i, g in enumerate(r["genes"])}
    own = np.array([lookup[t] for t in r["targets"]])
    n = len(r["targets"])
    counts = resample_counts(n)
    members = list(E.MEMBERS)
    if members != ["pds", "mse", "nmae", "fid", "reach", "jac"]:
        raise AssertionError("member order drift")
    anchors = json.loads((BASE / "work/h1-anchors.json").read_text())
    base = np.array([anchors["baseline"][m] for m in members])
    replicate = np.array([anchors["replicate"][m] for m in members])
    raw, overall, rows, checks = {}, {}, [], {}
    for name in ORDER:
        metrics, raw[name], checks[name] = arm_raw(name, r, truth, control, own,
                                                   counts, members)
        scaled = B.scaled(raw[name], base, replicate)
        overall[name] = scaled.mean(1)
        rows.append({"arm": name, **metrics["raw"],
                     "local_h1_overall": float(overall[name][0]),
                     "local_scaled": {m: float(scaled[0, j]) for j, m in enumerate(members)},
                     "panel_total_digest_sha256": metrics["panel_total_digest_sha256"],
                     "emitted_count_digest_sha256": metrics["emitted_count_digest_sha256"]})
    contrasts = {f"{a}-minus-{b}": contrast(a, b, overall, raw, members)
                 for a, b in CONTRASTS}
    nonfinite = sum(c["overall_ci"]["nonfinite_resamples"]
                    + sum(e["nonfinite_resamples"] for e in c["raw"].values())
                    for c in contrasts.values())
    # scaled-unit member deltas (point) for the kill-criterion accounting
    span = replicate - base
    scaled_deltas = {}
    for key, c in contrasts.items():
        a = key.split("-minus-")[0]
        b = key.split("-minus-")[1]
        ra = next(r_ for r_ in rows if r_["arm"] == a)
        rb = next(r_ for r_ in rows if r_["arm"] == b)
        sd = {}
        for m in members:
            d = (ra[m] - rb[m]) / span[members.index(m)]
            if m == "mse":
                d = float(np.clip(
                    (ra[m] - base[members.index(m)]) / span[members.index(m)], 0., 1.)
                    - np.clip(
                    (rb[m] - base[members.index(m)]) / span[members.index(m)], 0., 1.))
            sd[m] = float(d)
        scaled_deltas[key] = sd
    result = {
        "status": "COMPLETE",
        "scope": "85-target reused H1 development panel, cap15 current family, H1 union OFF; "
                 "not independent A/B/C validation; no board conversion factor",
        "arms": rows, "contrasts": contrasts, "scaled_member_deltas_point": scaled_deltas,
        "integrity": {"known_answer_exact": bool(ka["pass"]),
                      "bootstrap_identity_checks": True,
                      "nonfinite_resamples": int(nonfinite)},
        "bootstrap": {"draws": DRAWS, "seed": SEED, "unit": "target", "paired": True,
                      "fixed": ["reference PDS gallery", "H1-local anchors", "cap15 effects"],
                      "panel_mse_and_claim_cap_recomputed": True},
        "anchors_sha256": R.sha(BASE / "work/h1-anchors.json"),
        "reserve_read": False, "submission": False, "production_default_change": False,
    }
    R.write_json(WORK / "results.json", result)
    print(json.dumps(
        {"contrasts": {k: {"overall": v["overall"], "overall_ci": v["overall_ci"]}
                       for k, v in contrasts.items()},
         "scaled_deltas": scaled_deltas}, indent=2), flush=True)
    return result


if __name__ == "__main__":
    main()
