"""241 Phase 0 (DE half) — reproduce 200's published DE legs on 200's OWN renderer.

Reads 200's on-disk `de200.*.npz` tables (which were emitted with `VCC_MASS_FIX=minnorm`
everywhere except `art_nomf`) and re-derives the four DE members with this experiment's own
arithmetic.  Every row here must land on the number `experiments/200-.../RESULTS.md` §"Phase 1"
publishes, or Phase 1 of this experiment is measuring something else and must not run.

    uv run python phase0_de241.py
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
E200 = os.path.abspath(os.path.join(HERE, "..", "200-production-ablation"))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, HERE)
from lib241 import paired_de, row_for, seed_names                       # noqa: E402

os.makedirs(WORK, exist_ok=True)

# (label, arm stem, base stem, n seeds) -- exactly six200.py's COMPARISONS for the ten rows
PAIRS = [
    ("a23 mass fix",             "art_mf",          "art_nomf",        5),
    ("tiered spread policy",     "art_mf",          "art_flat_mf",     5),
    ("--overdispersion 4.0",     "art_mf",          "art_phi0_mf",     5),
    ("whitening (gamma 0.65)",   "reb_ship",        "reb_noWhiten",    3),
    ("PerturbMatch blend 0.5",   "reb_ship",        "reb_noPM",        3),
    ("CD4 residual f=0.10",      "reb_ship",        "reb_noCD4",       3),
    ("equal donor weights",      "reb_ship",        "reb_w41",         3),
    ("CPM gate 10 (shipped) LOO", "reb_ship_gate10", "reb_ship_gate0", 3),
]

# 200's published DE legs, RESULTS.md Phase 1 table (and work/six200-tablev2.json)
EXPECT = {
    "a23 mass fix":              dict(reach=-0.000295, fid=+0.000107, jac=+0.000045, nmae=+0.000022),
    "tiered spread policy":      dict(reach=-0.002119, fid=-0.002743, jac=+0.004863, nmae=+0.000036),
    "--overdispersion 4.0":      dict(reach=-0.001931, fid=+0.000603, jac=+0.005824, nmae=+0.000043),
    "whitening (gamma 0.65)":    dict(reach=+0.000777, fid=-0.002407, jac=-0.000355, nmae=+0.001835),
    "PerturbMatch blend 0.5":    dict(reach=-0.001104, fid=+0.000075, jac=+0.000299, nmae=-0.002514),
    "equal donor weights":       dict(reach=-0.000408, fid=+0.000417, jac=+0.000401, nmae=+0.009895),
    "CPM gate 10 (shipped) LOO": dict(reach=+0.000562, fid=+0.000240, jac=-0.001523, nmae=-0.000571),
}
# 200 §"Phase 0": KA_a15_raw must read reach 0.18283 and median n_pred 6034.
KA = dict(arm="KA_a15_raw", reach=0.18283, n_pred=6034.0)


def main() -> None:
    t0 = time.time()
    stems = sorted({s for _, a, b, n in PAIRS for s in (a, b)})
    need = {KA["arm"]}
    for _, a, b, n in PAIRS:
        need |= set(seed_names(a, n)) | set(seed_names(b, n))

    cache = os.path.join(WORK, "phase0_rows.json")
    rows = json.loads(open(cache).read()) if os.path.exists(cache) else {}
    for arm in sorted(need):
        if arm in rows:
            continue
        p = os.path.join(E200, f"de200.{arm}.npz")
        if not os.path.exists(p):
            sys.exit(f"missing {p}")
        rows[arm] = row_for(p)
        json.dump(rows, open(cache, "w"), indent=1, default=float)
        print(f"  scored {arm:<26} [{time.time() - t0:.0f}s]", flush=True)

    print("\n" + "=" * 100)
    print("KNOWN-ANSWER ROW — 164/157's a15, byte-identical artifact")
    print("=" * 100)
    k = rows[KA["arm"]]
    ok_r = abs(k["reach"] - KA["reach"]) < 5e-5
    ok_n = k["n_pred"] == KA["n_pred"]
    print(f"  reach        got {k['reach']:.5f}  expected {KA['reach']:.5f}   "
          f"{'PASS' if ok_r else 'FAIL'}")
    print(f"  median n_pred got {k['n_pred']:.0f}      expected {KA['n_pred']:.0f}        "
          f"{'PASS' if ok_n else 'FAIL'}")

    print("\n" + "=" * 118)
    print("200's DE legs, RE-DERIVED  (positive = removing the component would cost us)")
    print("=" * 118)
    print(f"  {'row':<28}{'member':>8}{'here':>12}{'200':>12}{'|d|':>11}   ")
    out, worst = {}, 0.0
    for label, arm, base, n in PAIRS:
        d = paired_de(rows, seed_names(arm, n), seed_names(base, n))
        out[label] = d
        e = EXPECT.get(label)
        line = f"  {label:<28}"
        for m in ("reach", "fid", "jac", "nmae"):
            if e is None:
                print(f"{line}{m:>8}{d[m]:>+12.6f}{'—':>12}{'—':>11}")
            else:
                dd = abs(d[m] - e[m])
                worst = max(worst, dd)
                print(f"{line}{m:>8}{d[m]:>+12.6f}{e[m]:>+12.6f}{dd:>11.2e}"
                      f"   {'ok' if dd < 1e-5 else 'MISMATCH'}")
            line = " " * 30
        print(f"{' ' * 30}{'sum4':>8}{sum(d[m] for m in ('reach','fid','jac','nmae')):>+12.6f}"
              f"   raw dreach {d['d_reach_raw']:+.5f} sem {d['reach_sem_raw']:.5f} "
              f"t {d['reach_t']:+.2f}  n={d['n_seeds']}")
    print(f"\n  worst |delta| against 200's published DE legs: {worst:.3e}   "
          f"{'PASS' if worst < 1e-5 else 'FAIL'}")

    json.dump(dict(rows=rows, pairs=out, ka=k, worst=worst),
              open(os.path.join(WORK, "phase0_de241.json"), "w"), indent=1, default=float)
    print(f"\nwrote work/phase0_de241.json [{time.time() - t0:.0f}s]")


if __name__ == "__main__":
    main()
