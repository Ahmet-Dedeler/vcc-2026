"""191 gate 0 (scoring half) — read 185's row off the tables `ka191.py` wrote.

Every number here must reproduce 185 §3b/3c/3d/3e for the SAME arm, computed with 185's own
instruments (`debench.Surface`, `six185.walk`, `six185.members`).  A miss means the port is
a different emitter and the preregistration does not apply to it.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
from scipy.stats import spearmanr

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
sys.path.insert(0, f"{ROOT}/experiments/185-rank-controlled-emission")
from debench import Surface  # noqa: E402
import six185 as SIX  # noqa: E402

P = os.path.join(WORK, "de191.%s.npz")
ARMS = [a for a in ("a15", "rc_hybrid", "rc_hybrid_seq") if os.path.exists(P % a)]

# 185 §1 and §3b/3c/3d, for the two arms this experiment reproduces.
EXPECT = {
    ("a15", "reach"): (0.18283, 2e-5),
    ("a15", "n_pred"): (6034, 0.5),
    ("a15", "jac"): (0.06824, 5e-5),
    ("a15", "nmae"): (0.96659, 5e-5),
    ("a15", "P1"): (0.7317, 5e-4),
    ("a15", "H_call"): (-0.391, 2e-3),
    ("a15", "E_call"): (+0.486, 2e-3),
    ("rc_hybrid", "reach"): (0.21681, 2e-5),
    ("rc_hybrid", "n_pred"): (6043, 0.5),
    ("rc_hybrid", "jac"): (0.06827, 5e-5),
    ("rc_hybrid", "nmae"): (0.96670, 5e-5),
    ("rc_hybrid", "P1"): (0.8780, 5e-4),
    ("rc_hybrid", "H_call"): (+0.984, 2e-3),
    ("rc_hybrid", "E_call"): (+0.007, 5e-3),
}


def main() -> None:
    meta = np.load(os.path.join(
        ROOT, "experiments/185-rank-controlled-emission/work/ctrl_meta.npz"
    ), allow_pickle=True)
    expr = np.log1p(meta["cpm"][meta["keep"]])
    z15 = np.load(P % "a15", allow_pickle=True)
    absH = z15["absH"].astype(np.float64)

    rows: list[tuple[str, str, float, float, float, bool]] = []
    out: dict[str, dict] = {}
    print(f"{'arm':<16}{'n_pred':>7}{'reach':>10}{'fid':>9}{'jac':>9}{'nmae':>9}"
          f"{'recall':>8}{'P(1)':>8}{'P(3)':>8}{'P(10)':>8}{'H_call':>9}{'E_call':>9}")
    for arm in ARMS:
        surface = Surface(P % arm)
        members = SIX.members(surface)
        _reach, pk, _ks = SIX.walk(surface)
        table = np.load(P % arm, allow_pickle=True)
        spear = SIX.spearmans(table, surface, absH, expr)
        record = dict(members, P1=pk[1], P3=pk[3], P10=pk[10],
                      H_call=spear["H_call"], E_call=spear["E_call"],
                      H_univ=spear["H_univ"], E_univ=spear["E_univ"])
        out[arm] = record
        print(f"{arm:<16}{record['n_pred']:>7}{record['reach']:>10.5f}{record['fid']:>9.5f}"
              f"{record['jac']:>9.5f}{record['nmae']:>9.5f}{record['recall']:>8.4f}"
              f"{pk[1]:>8.4f}{pk[3]:>8.4f}{pk[10]:>8.4f}"
              f"{spear['H_call']:>+9.3f}{spear['E_call']:>+9.3f}")
        for key, (want, tol) in EXPECT.items():
            if key[0] != arm:
                continue
            got = float(record[key[1]])
            rows.append((arm, key[1], got, want, tol, abs(got - want) <= tol))

    print("\nKNOWN-ANSWER ROWS against 185 (the mandate's gate 0)")
    bad = []
    for arm, name, got, want, tol, ok in rows:
        print(f"  {arm:<16}{name:<8}{got:>11.5f}   185 {want:<10.5f} tol {tol:<8g} "
              f"{'PASS' if ok else 'FAIL'}")
        if not ok:
            bad.append(f"{arm}/{name}")

    if "rc_hybrid" in out and "rc_hybrid_seq" in out:
        print("\nSEEDING: 185's two-stream bookkeeping vs production's single stream")
        for name in ("reach", "jac", "nmae", "fid"):
            a, b = out["rc_hybrid"][name], out["rc_hybrid_seq"][name]
            print(f"  {name:<8}{a:>11.5f}{b:>11.5f}   |diff| {abs(a - b):.5f}"
                  f"   (185 §3e twin spread on this family: 0.00007)")
        drift = abs(out["rc_hybrid"]["reach"] - out["rc_hybrid_seq"]["reach"])
        print(f"  reach drift {drift:.5f} vs 185's a15 seed spread 0.00746 "
              f"and its rc twin spread 0.00007")

    with open(os.path.join(WORK, "score191.json"), "w") as handle:
        json.dump(dict(members=out, known=[list(r) for r in rows]), handle, default=float)
    if bad:
        sys.exit(f"ABORT (rule 10 known-answer): {bad}")
    print("\nall known-answer rows PASS")


if __name__ == "__main__":
    main()
