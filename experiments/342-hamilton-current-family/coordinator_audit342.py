"""Post-run coordinator audit of the intended rounding decisions from frozen profiles.

Does not import rounding342, seam342, report342 or the worker's audit. It reconstructs
the original intended counts and independent Hamilton/random selections, checks all
255 delivered total vectors, and records which numerical gains were actually earned.
This is an audit, not a new arm or a change to the frozen gates.
"""

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = ROOT/"experiments/323-current-family-validation"
WORK = HERE/"work"
sys.path.insert(0, str(BASE))
import run323 as R


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    return json.loads(path.read_text())


def main():
    for relative, digest in read(HERE/"inputs.json")["inputs_sha256"].items():
        assert sha(ROOT/relative) == digest, relative
    with np.load(BASE/"work/reference.npz") as z:
        r = {k: z[k] for k in ("targets", "genes", "libraries", "b_ctrl")}
    with np.load(BASE/"work/effects.npz") as z:
        np.testing.assert_array_equal(z["targets"], r["targets"])
        effects = R.map_axis(z["cap15"], z["genes"], r["genes"])
    builder = R.builder()
    assert builder.H1_UNION_G is None and builder.MASS_FIX is None
    pick = np.random.default_rng(0)
    expected, rows = {a: [] for a in ("floor", "hamilton", "random_up")}, []
    for i, target in enumerate(r["targets"]):
        profile = builder.profile_from_effect(effects[i], r["b_ctrl"], min_cpm=10.)
        libraries = r["libraries"][pick.choice(len(r["libraries"]), 400, replace=False)]
        depth = int(min(np.median(libraries), 1_000_000))
        intended = 400*depth*profile
        floor = np.floor(intended).astype(np.int64)
        fraction = intended-floor
        increments = int(round(float(intended.sum())))-int(floor.sum())
        assert 0 <= increments <= np.count_nonzero(fraction > 0)
        hamilton = floor.copy()
        # Explicit gene-index secondary key; unlike the implementation's stable argsort.
        chosen = np.lexsort((np.arange(len(fraction)), -fraction))[:increments]
        hamilton[chosen] += 1
        random = floor.copy()
        candidates = np.where(fraction > 0)[0]
        selected = np.random.default_rng(3420909+i).choice(
            candidates, increments, replace=False)
        random[selected] += 1
        assert not np.any(fraction[chosen] <= 0)
        assert hamilton.sum() == random.sum() == round(float(intended.sum()))
        for name, totals in (("floor", floor), ("hamilton", hamilton), ("random_up", random)):
            expected[name].append(totals)
        rows.append({"target": str(target), "picked_median_depth": depth,
                     "increments": increments,
                     "hamilton_l1_to_intended": float(np.abs(hamilton-intended).sum()),
                     "random_l1_to_intended": float(np.abs(random-intended).sum()),
                     "floor_l1_to_intended": float(np.abs(floor-intended).sum())})
        assert rows[-1]["hamilton_l1_to_intended"] <= rows[-1]["random_l1_to_intended"]+1e-10
    checks = {}
    for name in expected:
        path = WORK/f"stats-{name}.npz"
        metrics = read(WORK/f"metrics-{name}.json")
        assert sha(path) == metrics["statistics_sha256"]
        with np.load(path) as z:
            np.testing.assert_array_equal(z["targets"], r["targets"])
            np.testing.assert_array_equal(z["sums"], expected[name])
            digests = [hashlib.sha256(np.ascontiguousarray(row).tobytes()).hexdigest()
                       for row in z["sums"]]
        assert hashlib.sha256("".join(digests).encode()).hexdigest() == (
            metrics["panel_total_digest_sha256"])
        rounding = read(WORK/f"rounding-{name}.json")
        for actual, derived in zip(rounding["targets"], rows):
            assert actual["target"] == derived["target"]
            assert actual["residual_k"] == derived["increments"]
            assert actual["calls"] == 2
        checks[name] = {"all_85_total_vectors_exact": True, "panel_digest_exact": True,
                        "statistics_sha256": sha(path)}
    result = read(WORK/"results.json")
    hf = result["contrasts"]["hamilton-minus-floor"]
    hr = result["contrasts"]["hamilton-minus-random_up"]
    independently_derived_gates = {
        "integrity_checks_pass": read(WORK/"known-answer.json")["pass"] and
                                 read(WORK/"independent-audit.json")["status"] == "PASS",
        "hamilton_overall_point_ge_plus_0_001_vs_floor": hf["overall"] >= .001,
        "hamilton_overall_lower_bound_positive_vs_floor": hf["overall_ci"]["lo"] > 0,
        "hamilton_overall_lower_bound_positive_vs_random_up": hr["overall_ci"]["lo"] > 0,
        "no_resolved_raw_mse_or_nmae_harm":
            not any(c["raw"][m]["lo"] > 0 for c in (hf, hr) for m in ("mse", "nmae")),
        "no_resolved_fid_reach_jac_harm":
            not any(c["raw"][m]["hi"] < 0 for c in (hf, hr) for m in ("fid", "reach", "jac")),
    }
    assert independently_derived_gates == result["gates"]
    status = ("LOCAL_GAIN_REQUIRES_SEPARATE_ABC_REVIEW"
              if all(independently_derived_gates.values()) else "NO_LOCAL_ADVANCEMENT")
    assert result["status"] == status
    path = WORK/"coordinator-audit.json"
    if path.exists():
        raise FileExistsError(path)
    output = {"status": "PASS", "checks": checks, "targets": rows,
              "frozen_gate_status_reproduced": status,
              "result_sha256": sha(WORK/"results.json"),
              "worker_audit_sha256": sha(WORK/"independent-audit.json"),
              "audit_code_sha256": sha(Path(__file__)),
              "finished_utc": datetime.now(timezone.utc).isoformat(),
              "not_reproduced": ["full changed-arm count regeneration"],
              "scope": "independent allocation replay, not biological replication"}
    path.write_text(json.dumps(output, indent=2, sort_keys=True, allow_nan=False)+"\n")
    print(json.dumps({k: v for k, v in output.items() if k != "targets"}, indent=2))


if __name__ == "__main__":
    main()
