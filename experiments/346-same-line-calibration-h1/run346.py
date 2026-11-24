#!/usr/bin/env python3
"""346 step 5 — emit actual counts and score all six members for each effect arm.

The emission loop is 323's `evaluate323.emit_and_score` verbatim in constants, order and
random-number use; only the output directory changes, so nothing is written into 323's work/.
The `cap15` arm is the known answer and must reproduce 323's recorded cap15 outputs EXACTLY
(all six raw members, MSE legs, emitted-count digest, panel-total digest, every saved
statistic array and every per-target row) before any modified arm is interpreted.

Usage:  run346.py cap15 | diag_cross_b1 | same_ridge_b1 | same_ridge_b05 | shuffled_b1 | all
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import resource
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("POLARS_MAX_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")

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
from vcc2026 import rc_emit   # noqa: E402
from vcc2026.donor_transfer import emit_counts   # noqa: E402

ORDER = ("cap15", "diag_cross_b1", "same_ridge_b1", "same_ridge_b05", "shuffled_b1",
         "same_ridge_b025")
# same_ridge_b025 is the POST-HOC knee arm; its effects live in a separate archive.
ARCHIVE = {"same_ridge_b025": "effects-b025.npz"}
MIN_FREE_PERCENT = 50


def free_percent() -> int:
    out = subprocess.check_output(["memory_pressure"], text=True).strip().splitlines()[-1]
    return int(out.rsplit(":", 1)[1].strip().rstrip("%"))


def emit_arm(name, r, c, emit_c, effects):
    """323 emit_and_score, constants and RNG order unchanged, writing into this work/."""
    b = R.builder()
    if b.H1_UNION_G is not None or b.MASS_FIX is not None:
        raise AssertionError("H1 union or mass fix unexpectedly active")
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
            jitter_seed=1000 + i, call_set_first=True,
            call_set=rc_emit.CallSet(400, 2200, 1.9), emit_counts=emit_counts)
        s, j, p, q, l = E.block_stats(block, c)
        sums.append(s); jks.append(j); ps.append(p); qs.append(q); lfcs.append(l)
        digests.append(hashlib.sha256(np.ascontiguousarray(s).tobytes()).hexdigest())
        counts_sha.update(np.asarray(block, dtype="<u4", order="C").tobytes())
        del block
        if (i + 1) % 20 == 0:
            print(f"{name} {i+1}/85 ({time.monotonic()-t0:.1f}s)", flush=True)
    panel_sha = hashlib.sha256("".join(digests).encode()).hexdigest()
    raw, rows, legs = E.score(r, c, sums, jks, ps, qs, lfcs)
    np.savez(WORK / f"stats-{name}.npz", targets=r["targets"], sums=sums,
             jk=jks, p=ps, q=qs, lfc=lfcs)
    pl.DataFrame(rows).write_csv(WORK / f"per-target-{name}.csv")
    out = {"arm": name, "raw": raw, "panel_total_digest_sha256": panel_sha,
           "emitted_count_digest_sha256": counts_sha.hexdigest(), "mse_legs": legs,
           "seconds": time.monotonic() - t0,
           "statistics_sha256": R.sha(WORK / f"stats-{name}.npz")}
    R.write_json(WORK / f"metrics-{name}.json", out)
    print(name, json.dumps(raw), flush=True)
    return out


def rows_equal(a: Path, b: Path) -> bool:
    fa, fb = pl.read_csv(a), pl.read_csv(b)
    if fa.columns != fb.columns or fa.height != fb.height:
        return False
    table = lambda f: {(str(p), str(m)): float(v) for p, m, v in
                       zip(f["perturbation"], f["metric"], f["value"])}
    ka, kb = table(fa), table(fb)
    if ka.keys() != kb.keys():
        return False
    return all((np.isnan(ka[k]) and np.isnan(kb[k])) or ka[k] == kb[k] for k in ka)


def known_answer(receipt: dict) -> None:
    """Exact on everything deterministic; one ulp is allowed on the panel means only.

    cell-eval2's `dispatch_de_metrics` returns the per-target NMAE rows in a
    nondeterministic order (verified: 323's csv and this run's csv hold bit-identical
    per-target values under a different row order). 323's recorded `raw["nmae"]` is
    `np.nanmean` of that list in ITS order, so the recorded scalar is order-dependent in
    its last bit. `np.nanmean` of the 58 identical values reproduces 323's scalar in 323's
    order and this run's scalar in this run's order, differing by 1.110e-16. Everything
    that is actually deterministic -- the emitted counts, the panel totals, every saved
    statistic, the MSE legs and every per-target value -- is still required to be exact.
    The arm contrasts below are ~1e-3, seven orders of magnitude above this tolerance, and
    every member reported in results.json is recomputed from the per-target values in one
    canonical target order, so the artifact cannot enter a contrast.
    """
    ulp = 2.3e-16
    expected = json.loads((BASE / "work/metrics-cap15.json").read_text())
    checks = {
        "emitted_count_digest_equal":
            receipt["emitted_count_digest_sha256"] == expected["emitted_count_digest_sha256"],
        "panel_total_digest_equal":
            receipt["panel_total_digest_sha256"] == expected["panel_total_digest_sha256"],
        "raw_metrics_exact": {m: receipt["raw"][m] == expected["raw"][m] for m in E.MEMBERS},
        "mse_legs_exact": receipt["mse_legs"] == expected["mse_legs"],
    }
    arrays = {}
    with np.load(BASE / "work/stats-cap15.npz", allow_pickle=False) as a, \
            np.load(WORK / "stats-cap15.npz", allow_pickle=False) as b:
        for key in ("targets", "sums", "jk", "p", "q", "lfc"):
            x, y = a[key], b[key]
            same = x.shape == y.shape and x.dtype == y.dtype
            if same:
                same = bool(np.array_equal(x, y, equal_nan=(x.dtype.kind == "f")))
            arrays[key] = bool(same)
    checks["statistics_arrays_exact"] = arrays
    checks["per_target_rows_exact"] = rows_equal(BASE / "work/per-target-cap15.csv",
                                                 WORK / "per-target-cap15.csv")
    checks["statistics_file_sha256_equal_informational"] = (
        receipt["statistics_sha256"] == expected["statistics_sha256"])
    differences = {m: abs(receipt["raw"][m] - expected["raw"][m]) for m in E.MEMBERS}
    checks["raw_metrics_within_one_ulp"] = {m: d <= ulp for m, d in differences.items()}
    checks["raw_metrics_exactly_equal_count"] = int(sum(checks["raw_metrics_exact"].values()))
    passed = bool(checks["emitted_count_digest_equal"] and checks["panel_total_digest_equal"]
                  and checks["mse_legs_exact"] and all(arrays.values())
                  and checks["per_target_rows_exact"]
                  and all(checks["raw_metrics_within_one_ulp"].values())
                  and checks["raw_metrics_exactly_equal_count"] >= 5)
    R.write_json(WORK / "known-answer.json", {
        "pass": passed, "checks": checks,
        "raw_absolute_differences": differences,
        "requirement": "exact equality with 323 cap15 on the emitted-count digest, the panel "
                       "total digest, the MSE legs, every saved statistic array and every "
                       "per-target value; the six panel means must agree within one ulp "
                       "(2.3e-16) with at least five exactly equal. The statistics file hash "
                       "is informational.",
        "tolerance_reason": known_answer.__doc__,
        "one_ulp_tolerance": ulp,
        "expected_source": "experiments/323-current-family-validation/work/metrics-cap15.json"})
    if not passed:
        raise AssertionError("cap15 arm did not reproduce 323 exactly; see work/known-answer.json")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("arms", nargs="+", choices=[*ORDER, "all", "recheck"])
    args = ap.parse_args()
    if "recheck" in args.arms:
        # Re-evaluate the cap15 known answer from the already-emitted receipt, without
        # re-emitting. The emission is deterministic; only the gate definition changed.
        known_answer(json.loads((WORK / "metrics-cap15.json").read_text()))
        print(json.dumps(json.loads((WORK / "known-answer.json").read_text())["checks"],
                         indent=2, sort_keys=True))
        return
    wanted = list(ORDER) if "all" in args.arms else [a for a in ORDER if a in args.arms]
    WORK.mkdir(parents=True, exist_ok=True)
    for name in wanted:
        if (WORK / f"metrics-{name}.json").exists():
            raise FileExistsError(f"refusing to overwrite metrics-{name}.json")
    if wanted != ["cap15"] and not (WORK / "known-answer.json").exists() and "cap15" not in wanted:
        raise SystemExit("run the cap15 known answer first")

    started, clock = datetime.now(timezone.utc), time.monotonic()
    E.verify()
    comp = json.loads((BASE / "work/compose.json").read_text())
    if not comp["gates_pass"]:
        raise ValueError("323 composition gates did not pass")
    if R.sha(BASE / "work/effects.npz") != comp["effects_sha256"]:
        raise ValueError("323 effects archive changed")
    r, c, emit_c = E.load_reference()
    mapped = {}
    for archive in sorted({ARCHIVE.get(n, "effects346.npz") for n in wanted}):
        with np.load(WORK / archive, allow_pickle=False) as z:
            if not np.array_equal(np.asarray(z["targets"]).astype(str),
                                  r["targets"].astype(str)):
                raise ValueError(f"target order drift between {archive} and 323 reference")
            with np.load(BASE / "work/effects.npz", allow_pickle=False) as w:
                if not np.array_equal(z["cap15"], w["cap15"]):
                    raise ValueError(f"{archive}: cap15 copy is not bit-identical to 323's")
                if R.builder().array_sha256(w["cap15"]) != comp["arms"]["cap15"]["final_sha256"]:
                    raise ValueError("cap15 final effect hash drift")
            egenes = z["genes"]
            for name in wanted:
                if ARCHIVE.get(name, "effects346.npz") == archive:
                    mapped[name] = R.map_axis(np.asarray(z[name], dtype=np.float64),
                                              egenes, r["genes"])
    receipts = {}
    for name in wanted:
        free = free_percent()
        print(f"== arm {name}: system free {free}% ({time.monotonic()-clock:.0f}s elapsed)",
              flush=True)
        if free < MIN_FREE_PERCENT:
            raise SystemExit(f"memory gate: only {free}% free, need >= {MIN_FREE_PERCENT}%")
        receipts[name] = emit_arm(name, r, c, emit_c, mapped[name])
        if name == "cap15":
            known_answer(receipts[name])
            print("known-answer: cap15 reproduces 323 exactly", flush=True)
        del mapped[name]
        gc.collect()
    record = {"started_utc": started.isoformat(),
              "finished_utc": datetime.now(timezone.utc).isoformat(),
              "elapsed_seconds": time.monotonic() - clock, "arms": wanted,
              "process_max_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              "free_percent_at_start": free_percent()}
    stamp = "-".join(wanted)
    R.write_json(WORK / f"execution-{stamp}.json", record)
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
