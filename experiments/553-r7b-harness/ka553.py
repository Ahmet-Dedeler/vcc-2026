#!/usr/bin/env python3
"""553 — KA-H1: the generalised harness pointed at h1-test must reproduce recorded runs BIT FOR BIT.

usage: uv run python experiments/553-r7b-harness/ka553.py <mode> <seed>
For each recorded pair (pre-zs and, where recorded, zs) it compares this directory's work/h1_<mode>/ stats with the
recorded stats (np.array_equal, equal_nan, on p, q, lfc, sums and jk) and reports |d| per member between this run's inline
evaluate323.score and the recorded metrics json (pass: arrays bit-equal and |d| <= 1e-12; the scorer's own nmae
non-determinism is ~1e-16). Prints `KA PASS ...` / `KA FAIL ...`, appends a row to work/ka553.json,
exits 0 only if every compared pair is bit-equal (a recorded pair that does not exist is an error, not a skip).

  cpm            -> 528's h1_cpm (== 435's production run)                        pre-zs + zs
  lcap-kept      -> 528's h1_lcap-kept                                           pre-zs + zs
  cpm+c4h        -> 517's h1_c4-c4t0.05 (pre-zs) and 537's h1_cpm+c4h (zs)
  cpm-n1600      -> 537's h1_cpm-n1600 (537 section 3.2's K 2,000 row)           pre-zs + zs
  lcap-kept+c4h  -> 537's h1_lcap-kept+c4h (CPH, 537 section 3.2 P1)             pre-zs + zs
  cpm-n1200      -> 528's h1_cpm-n1200                                           pre-zs + zs
  lcap-kept+c4h-s-> 537's h1_lcap-kept+c4h-s (P9's shuffled hook)                pre-zs + zs
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TAG = "base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse"
HOOKTAG = "base_s{seed}_bps_hookem_mde0.1cap0.1min0.3coarse"
E528 = ROOT / "experiments/528-callprop-block/work"
E517 = ROOT / "experiments/517-flex-programme/work"
E537 = ROOT / "experiments/537-arm-pipeline/work"
REF = {
    "cpm": [(E528 / "h1_cpm", TAG, False), (E528 / "h1_cpm", TAG, True)],
    "lcap-kept": [(E528 / "h1_lcap-kept", TAG, False), (E528 / "h1_lcap-kept", TAG, True)],
    "cpm+c4h": [(E517 / "h1_c4-c4t0.05", HOOKTAG, False), (E537 / "h1_cpm+c4h", TAG, True)],
    "cpm-n1600": [(E537 / "h1_cpm-n1600", TAG, False), (E537 / "h1_cpm-n1600", TAG, True)],
    "lcap-kept+c4h": [(E537 / "h1_lcap-kept+c4h", TAG, False), (E537 / "h1_lcap-kept+c4h", TAG, True)],
    "cpm-n1200": [(E528 / "h1_cpm-n1200", TAG, False), (E528 / "h1_cpm-n1200", TAG, True)],
    "lcap-kept+c4h-s": [(E537 / "h1_lcap-kept+c4h-s", TAG, False), (E537 / "h1_lcap-kept+c4h-s", TAG, True)],
}
MEMBERS = ("pds", "mse", "nmae", "fid", "reach", "jac")


def arrays_equal(a, b) -> dict:
    """p, q, lfc, sums, jk bit-equal (NaN == NaN), and the target lists equal."""
    eq = {k: bool(np.array_equal(np.asarray(a[k]), np.asarray(b[k]), equal_nan=True)) for k in ("p", "q", "lfc", "sums", "jk")}
    eq["targets"] = bool(np.array_equal(np.asarray(a["targets"]).astype(str), np.asarray(b["targets"]).astype(str)))
    return eq


def compare(dir_a: str, dir_b: str, seed: int) -> None:
    """--compare: two of this directory's run dirs (e.g. the live-bulk and cached-bulk runs of one arm) must hold bit-equal
    stats, pre-zs and zs (KA-B)."""
    ok_all = True
    for zs in (False, True):
        t = TAG.format(seed=seed) + ("_dnk250zs" if zs else "")
        a = np.load(HERE / "work" / dir_a / f"stats-{t}.npz", allow_pickle=True)
        b = np.load(HERE / "work" / dir_b / f"stats-{t}.npz", allow_pickle=True)
        eq = arrays_equal(a, b)
        ok = all(eq.values())
        ok_all &= ok
        print(f"KA {'PASS' if ok else 'FAIL'} compare {dir_a} == {dir_b} s{seed} {'zs' if zs else 'pre-zs'}: {eq}")
    sys.exit(0 if ok_all else 1)


def main():
    if sys.argv[1] == "--compare":
        compare(sys.argv[2], sys.argv[3], int(sys.argv[4]))
    mode, seed = sys.argv[1], int(sys.argv[2])
    surf = sys.argv[3] if len(sys.argv) > 3 else "h1"      # h1cb: R7b's code path on h1 (sealed stub: stats only)
    ok_all = True
    rows = []
    for d, tag, zs in REF[mode]:
        t = tag.format(seed=seed) + ("_dnk250zs" if zs else "")
        mt = TAG.format(seed=seed) + ("_dnk250zs" if zs else "")
        ref_s, ref_m = d / f"stats-{t}.npz", d / f"metrics-{t}.json"
        mine_s, mine_m = HERE / f"work/{surf}_{mode}" / f"stats-{mt}.npz", HERE / f"work/{surf}_{mode}" / f"metrics-{mt}.json"
        form = "zs" if zs else "pre-zs"
        if not ref_s.exists():
            print(f"KA FAIL {mode} s{seed} {form}: the recorded {ref_s.relative_to(ROOT)} does not exist")
            ok_all = False
            continue
        if not mine_s.exists():
            print(f"KA FAIL {mode} s{seed} {form}: missing {mine_s.relative_to(ROOT)}")
            ok_all = False
            continue
        a, b = np.load(ref_s, allow_pickle=True), np.load(mine_s, allow_pickle=True)
        eq = arrays_equal(a, b)
        ra = json.loads(ref_m.read_text())["raw"]
        rb = json.loads(mine_m.read_text())["raw"]
        if surf == "h1cb":      # the sealed stub scores nothing: arrays only, and the stub must have written nulls
            assert all(rb[m] is None for m in MEMBERS), "the sealed stub wrote a score"
            dmem = {m: 0.0 for m in MEMBERS}
        else:
            dmem = {m: abs(float(ra[m]) - float(rb[m])) for m in MEMBERS}
        # cell-eval2's nmae dispatch is itself non-deterministic at the last ulp (kascore553: 3 distinct values in 6
        # identical calls), so members pass at |d| <= 1e-12; the arrays must be bit-equal
        ok = all(eq.values()) and max(dmem.values()) <= 1e-12
        ok_all &= ok
        print(f"KA {'PASS' if ok else 'FAIL'} {mode} s{seed} {form}: {mine_s.relative_to(ROOT)} == "
              f"{ref_s.relative_to(ROOT)} {eq}; |d| per member {dmem}")
        rows.append(dict(surface=surf, mode=mode, seed=seed, form=form, recorded=str(ref_s.relative_to(ROOT)),
                         mine=str(mine_s.relative_to(ROOT)), arrays_equal=eq, abs_d=dmem, raw_recorded=ra, raw_mine=rb,
                         pass_=ok, checked=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())))
    out = HERE / "work/ka553.json"
    prev = json.loads(out.read_text()) if out.exists() else []
    prev = [r for r in prev if not (r.get("surface", "h1") == surf and r["mode"] == mode and r["seed"] == seed)] + rows
    out.write_text(json.dumps(prev, indent=1))
    sys.exit(0 if ok_all and rows else 1)


if __name__ == "__main__":
    main()
