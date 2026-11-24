#!/usr/bin/env python3
"""537 — known-answer checks for run537.py: a run's stats must equal a recorded run's stats bit for bit.

usage: uv run python experiments/537-arm-pipeline/ka537.py <mode> <seed>
Prints `KA PASS <mode> s<seed> ...` or `KA FAIL ...` and exits 0 / 1. Pairs (pre-zs and, where recorded, zs):
  cpm        -> 528's h1_cpm (== 435's production run; 528's own KA)              pre-zs + zs
  lcap-kept  -> 528's h1_lcap-kept                                                 pre-zs + zs
  cpm+c4h    -> 517's h1_c4-c4t0.05 (== 496's real hook run; 517's own KA)        pre-zs only (517 ran it without zs)
  cpm+mzg    -> 528's h1_cpm+mzg (== 517's c4-mzg0.05)                             pre-zs + zs
  cpm-n1200  -> 528's h1_cpm-n1200                                                 pre-zs + zs
  cpm+c4h0.05 -> 517's h1_c4-c4t0.05 (pre-zs) and this directory's h1_cpm+c4h (zs): the edited harness (-s, -bs)
The comparison is np.array_equal on p, q, lfc and sums (equal_nan).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TAG = "base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse"
E528 = ROOT / "experiments/528-callprop-block/work"
E517 = ROOT / "experiments/517-flex-programme/work"
REF = {
    "cpm": [(E528 / "h1_cpm", TAG), (E528 / "h1_cpm", TAG + "_dnk250zs")],
    "lcap-kept": [(E528 / "h1_lcap-kept", TAG), (E528 / "h1_lcap-kept", TAG + "_dnk250zs")],
    "cpm+c4h": [(E517 / "h1_c4-c4t0.05", "base_s{seed}_bps_hookem_mde0.1cap0.1min0.3coarse")],
    "cpm+mzg": [(E528 / "h1_cpm+mzg", TAG), (E528 / "h1_cpm+mzg", TAG + "_dnk250zs")],
    "cpm-n1200": [(E528 / "h1_cpm-n1200", TAG), (E528 / "h1_cpm-n1200", TAG + "_dnk250zs")],
    # the harness after its 15:50Z / 16:10Z edits (-s, -bs): the same construction as `cpm+c4h` under an explicit tau, so a
    # fresh directory; pre-zs == 517's recorded run, zs == this directory's own `cpm+c4h` s0 (the unedited harness)
    "cpm+c4h0.05": [(E517 / "h1_c4-c4t0.05", "base_s{seed}_bps_hookem_mde0.1cap0.1min0.3coarse"),
                    (HERE / "work/h1_cpm+c4h", TAG + "_dnk250zs")],
    # the harness after its 18:55Z edit (-s for +mzg): `cpm+mzg` under an explicit r 0 (no flip), a fresh directory
    "cpm+mzg-r0": [(E528 / "h1_cpm+mzg", TAG), (E528 / "h1_cpm+mzg", TAG + "_dnk250zs")],
}


def main():
    mode, seed = sys.argv[1], int(sys.argv[2])
    ok_all = True
    for i, (d, tag) in enumerate(REF[mode]):
        ref = d / f"stats-{tag.format(seed=seed)}.npz"
        mine = HERE / f"work/h1_{mode}" / f"stats-{TAG.format(seed=seed)}{'_dnk250zs' if i == 1 else ''}.npz"
        if not ref.exists():
            print(f"KA SKIP {mode} s{seed}: no recorded {ref.relative_to(ROOT)}")
            continue
        if not mine.exists():
            print(f"KA FAIL {mode} s{seed}: missing {mine.relative_to(ROOT)}")
            ok_all = False
            continue
        a, b = np.load(ref, allow_pickle=True), np.load(mine, allow_pickle=True)
        eq = {k: bool(np.array_equal(np.asarray(a[k]), np.asarray(b[k]), equal_nan=True)) for k in ("p", "q", "lfc", "sums")}
        ok = all(eq.values())
        ok_all &= ok
        print(f"KA {'PASS' if ok else 'FAIL'} {mode} s{seed} {'zs' if i == 1 else 'pre-zs'}: {mine.name} == "
              f"{ref.relative_to(ROOT)} {eq}")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
