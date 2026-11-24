#!/usr/bin/env python3
"""416 known answers.
  uv run python ka416.py base395          re-emit 410's base395 (k3c5bq80, no ballast, seed 0) into work/h1 via 400's harness;
                                          digests + six members must equal 410's work/h1/metrics-base395.json
  uv run python ka416.py compare MINE THEIRS   compare work/h1/{metrics,stats}-MINE with 410's -THEIRS (bit-level arrays)
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work" / "h1"
T410 = ROOT / "experiments/410-pertscale-revisit/work/h1"
TOL = 1e-12   # members: the scorer's float reductions can differ by 1 ulp (410's own KA: |d| <= 1.1e-16); arrays/digests exact
os.environ.setdefault("NUMBA_CACHE_DIR", str(HERE / "work" / "numba-cache"))
sys.path.insert(0, str(ROOT / "experiments/400-uncalled-knn"))

import numpy as np  # noqa: E402


def compare(mine, theirs):
    a = json.loads((WORK / f"metrics-{mine}.json").read_text())
    t = json.loads((T410 / f"metrics-{theirs}.json").read_text())
    diffs = {m: abs(a["raw"][m] - t["raw"][m]) for m in t["raw"]}
    za, zt = np.load(WORK / f"stats-{mine}.npz"), np.load(T410 / f"stats-{theirs}.npz")
    arrays = {k: bool(np.array_equal(za[k], zt[k], equal_nan=(za[k].dtype.kind == "f"))) for k in ("sums", "jk", "p", "q", "lfc")}
    ok = all(v <= TOL for v in diffs.values()) and all(arrays.values())
    rec = dict(mine=mine, theirs=f"410/{theirs}", pass_=ok, raw_abs_diffs=diffs, arrays_equal=arrays)
    (WORK / f"known-answer-{mine}.json").write_text(json.dumps(rec, indent=1))
    print("KNOWN ANSWER", mine, "vs 410", theirs, "PASS" if ok else "FAIL", json.dumps(rec), flush=True)
    if not ok:
        raise SystemExit(1)


def base395():
    import run400 as H  # noqa: PLC0415
    H.WORK = WORK
    orig = H.block_for_arm
    H.block_for_arm = lambda arm, *a, **k: orig("k3c5bq80", *a, **k)
    if (WORK / "metrics-base395.json").exists():
        rec = json.loads((WORK / "metrics-base395.json").read_text())
    else:
        b = H.builder389()
        eff, genes, _ = H.native_cap15(b)
        rec = H.emit_arm("base395", eff, genes)
    exp = json.loads((T410 / "metrics-base395.json").read_text())
    ok = (rec["emitted_count_digest_sha256"] == exp["emitted_count_digest_sha256"]
          and rec["panel_total_digest_sha256"] == exp["panel_total_digest_sha256"]
          and rec["statistics_sha256"] == exp["statistics_sha256"]
          and all(abs(rec["raw"][m] - exp["raw"][m]) <= TOL for m in exp["raw"]))
    (WORK / "known-answer-base395.json").write_text(json.dumps(dict(pass_=ok, diffs={m: abs(rec["raw"][m] - exp["raw"][m])
                                                                                     for m in exp["raw"]}), indent=1))
    print("KNOWN ANSWER base395", "PASS" if ok else "FAIL", flush=True)
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    if sys.argv[1] == "base395":
        base395()
    else:
        compare(sys.argv[2], sys.argv[3])
