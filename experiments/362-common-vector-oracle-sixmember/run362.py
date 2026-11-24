#!/usr/bin/env python3
"""362 step 2 — emit + score all six arms through the exact 323 evaluator.

Own-directory replica of evaluate323.emit_and_score (same builder, same constants,
same RNG streams), writing ONLY into this experiment's work/. The cap15 arm must
reproduce 323's recorded cap15 outputs EXACTLY before any changed arm is interpreted.
Sequential arms, one process, ~142 s / ~5 GiB per arm.
"""
from __future__ import annotations

import gc
import hashlib
import json
import os
import resource
import sys
import time
from pathlib import Path

os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("POLARS_MAX_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = ROOT / "experiments/323-current-family-validation"
WORK = HERE / "work"
for _p in (str(BASE), str(ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import run323 as R  # noqa: E402
import evaluate323 as E  # noqa: E402
from vcc2026 import rc_emit  # noqa: E402
from vcc2026.donor_transfer import emit_counts  # noqa: E402

ORDER = ("cap15", "oracle_a05", "oracle_a1", "oracle_a1_normpreserve",
         "generic_a1", "shuffled_a1")

EXPECTED_RAW = {
    "fid": 0.4975545182292535, "jac": 0.07475983298131939,
    "mse": 0.9749462418923823, "nmae": 0.9508364630325662,
    "pds": 0.8896358543417369, "reach": 0.25434726355678233,
}
EXPECTED_COUNT_DIGEST = "547cf0abdc984ce92cd778adf3c86e7ef212046ad434b5cf94f3f49fbc3a529f"
EXPECTED_PANEL_DIGEST = "f699b3d546c18d45ab4f142ec4cdf1bcb6973f1ab4bc16fadc6b4f4ac62e8566"


def _rows_equal(a: Path, b: Path) -> bool:
    fa, fb = pl.read_csv(a), pl.read_csv(b)
    if fa.columns != fb.columns or fa.height != fb.height:
        return False

    def table(frame):
        return {(str(p), str(m)): float(v) for p, m, v in
                zip(frame["perturbation"], frame["metric"], frame["value"])}
    ka, kb = table(fa), table(fb)
    if ka.keys() != kb.keys():
        return False
    return all((np.isnan(ka[k]) and np.isnan(kb[k])) or ka[k] == kb[k] for k in ka)


def known_answer(receipt: dict) -> None:
    expected = json.loads((BASE / "work/metrics-cap15.json").read_text())
    checks = {
        "emitted_count_digest_equal":
            receipt["emitted_count_digest_sha256"] == expected["emitted_count_digest_sha256"],
        "panel_total_digest_equal":
            receipt["panel_total_digest_sha256"] == expected["panel_total_digest_sha256"],
        "raw_metrics_exact":
            {m: receipt["raw"][m] == expected["raw"][m] for m in E.MEMBERS},
        "mse_legs_exact": receipt["mse_legs"] == expected["mse_legs"],
    }
    differences = {m: abs(receipt["raw"][m] - expected["raw"][m]) for m in E.MEMBERS}
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
    checks["per_gene_totals_exact"] = arrays["sums"]
    checks["per_target_rows_exact"] = _rows_equal(
        BASE / "work/per-target-cap15.csv", WORK / "per-target-cap15.csv")
    passed = bool(checks["emitted_count_digest_equal"] and checks["panel_total_digest_equal"]
                  and all(checks["raw_metrics_exact"].values()) and checks["mse_legs_exact"]
                  and all(arrays.values()) and checks["per_target_rows_exact"])
    R.write_json(WORK / "known-answer.json", {
        "pass": passed, "checks": checks, "raw_absolute_differences": differences,
        "expected_source": "experiments/323-current-family-validation/work/metrics-cap15.json"})
    if not passed:
        raise AssertionError("cap15 did not reproduce 323 cap15 exactly; "
                             "see work/known-answer.json (preserved, not tuned)")


def emit_arm(name: str, r, c, emit_c, effects, builder) -> dict:
    pick_rng, emit_rng = builder.make_rngs(0)
    call_set = rc_emit.CallSet(400, 2200, 1.9)
    sums, jks, ps, qs, lfcs, digests = [], [], [], [], [], []
    counts_sha = hashlib.sha256()
    t0 = time.monotonic()
    for i, target in enumerate(r["targets"]):
        profile = builder.profile_from_effect(effects[i], r["b_ctrl"], min_cpm=10.)
        picks = pick_rng.choice(len(r["libraries"]), 400, replace=False)
        spread = builder.emission_spread_vector(
            effects[i], r["b_ctrl"], min_cpm=10., tail_spread=1., keep_top=8000,
            keep_top_frac=None, keep_spread=.25, rank_score=None, spread_map="block",
            spread_lo=0., spread_hi=.60)
        gate = builder.emission_rank_axis(effects[i], r["b_ctrl"], min_cpm=10.)[0]
        block, _ = rc_emit.rc_hybrid_block(
            profile, r["libraries"][picks], .1, emit_rng, spread=spread,
            overdispersion=4., control=emit_c, key=np.where(gate, effects[i], 0.),
            jitter_seed=1000 + i, call_set_first=True,
            call_set=call_set, emit_counts=emit_counts)
        s, j, p, q, l = E.block_stats(block, c)
        sums.append(s)
        jks.append(j)
        ps.append(p)
        qs.append(q)
        lfcs.append(l)
        digests.append(hashlib.sha256(np.ascontiguousarray(s).tobytes()).hexdigest())
        counts_sha.update(np.asarray(block, dtype="<u4", order="C").tobytes())
        del block
        if (i + 1) % 10 == 0:
            print(f"{name} {i + 1}/{len(r['targets'])} ({time.monotonic() - t0:.1f}s)",
                  flush=True)
    panel_sha = hashlib.sha256("".join(digests).encode()).hexdigest()
    raw, rows, legs = E.score(r, c, sums, jks, ps, qs, lfcs)
    np.savez(WORK / f"stats-{name}.npz", targets=r["targets"], sums=sums,
             jk=jks, p=ps, q=qs, lfc=lfcs)
    pl.DataFrame(rows).write_csv(WORK / f"per-target-{name}.csv")
    out = {"arm": name, "raw": raw, "panel_total_digest_sha256": panel_sha,
           "emitted_count_digest_sha256": counts_sha.hexdigest(),
           "mse_legs": legs, "seconds": time.monotonic() - t0,
           "statistics_sha256": R.sha(WORK / f"stats-{name}.npz")}
    R.write_json(WORK / f"metrics-{name}.json", out)
    print(name, json.dumps(raw), flush=True)
    return out


def main() -> None:
    for name in ORDER:
        if (WORK / f"metrics-{name}.json").exists():
            raise FileExistsError(f"refusing to overwrite metrics-{name}.json")
    t0 = time.monotonic()
    builder = R.builder()
    if builder.H1_UNION_G is not None or builder.MASS_FIX is not None:
        raise AssertionError("H1 union or mass fix unexpectedly active")
    r, c, emit_c = E.load_reference()
    with np.load(WORK / "designs362.npz", allow_pickle=False) as z:
        genes_eff = np.asarray(z["genes"]).astype(str)
        assert (np.asarray(z["targets"]).astype(str) == r["targets"]).all()
        designs = {name: np.asarray(z[name], dtype=np.float64) for name in ORDER}
    receipts = {}
    for name in ORDER:
        effects = R.map_axis(designs[name], genes_eff, r["genes"])
        del designs[name]
        receipts[name] = emit_arm(name, r, c, emit_c, effects, builder)
        del effects
        if name == "cap15":
            known_answer(receipts[name])
            print("known-answer: cap15 reproduces 323 cap15 exactly", flush=True)
        gc.collect()
    R.write_json(WORK / "emission.json", {
        "arms": ORDER, "elapsed_seconds": time.monotonic() - t0,
        "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "known_answer_pass": True})
    print(f"done in {time.monotonic() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
