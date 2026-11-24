#!/usr/bin/env python3
"""631f — held-out DOSE check for the called-DOWN magnitude lever (the coordinator's gate 3; 633's design).

proxy631 found, on H1's REAL emitter, that scaling every called DOWN gene's delivered lfc by m = 2 is worth −0.0227 (f 1)
/ −0.0251 (f 0.35) raw nmae, with m = 2 chosen on both cross-fit halves. Is that dose an H1 artefact, or does it
transfer to clean cancer lines? Model: 612's call-set-aware fold machinery, as in 633:
  base   = the plain-hook form (BETA * cons0 + s * add);
  called = production's head 400 + CPM block 2,200 (ds612c.callset); DOWN = base < 0;
  lever  = P -> m * P on called DOWN cells; m in 0.5, 1 (base), 1.5, 2, 2.5, 3, 4;
  score  = per-target nmae at 521's pools f 1 / 0.5 / 0.35; Δ vs base (m = 1), raw, negative = better.
Surfaces: h1 (retrodiction), k562_gwps (BINDING clean cancer line), nadig_jurkat (633j: NOT clean by its preregistered
block check, so descriptive), orion_hct116 / orion_hek293t (descriptive only: 609's depth skews on high-expression genes,
which the CPM block contains).

RETRODICTION GATE (fixed before the run): the h1 row must reproduce the real-emitter result in shape: Δ(m = 2) < 0 at
both f 1 and f 0.35, and the argmin m at f 0.35 in {1.5, 2, 2.5, 3}. FAIL -> every fold row is descriptive only.
DECISION (fixed before the run), on k562_gwps at f 0.5 and f 0.35 (both must agree):
  argmin m >= 2.5 and beats m = 2 by > 0.002 raw  -> a HIGHER dose transfers (reopen the dose upward);
  argmin m <= 1.0 or Δ(m = 2) >= 0                -> the x2 does NOT transfer (October keeps m = 1 for this lever);
  otherwise                                       -> m = 2 transfers (the H1 dose stands).
Analysis only, < 4 GB.   uv run python experiments/631-nmae-headroom/dose631f.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/612-key-gated-hook"))
import ds612 as D  # noqa: E402
import ds612b as B  # noqa: E402
import ds612c as C  # noqa: E402

MS = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]
REG = (1.0, 0.5, 0.35)
J633 = ROOT / "experiments/633-updose-folds"
SURFS = ["h1", "k562_gwps", "nadig_jurkat", "orion_hct116", "orion_hek293t"]


def main() -> int:
    gate_j = np.load(J633 / "work/gate-nadig_jurkat.npz", allow_pickle=True)
    _pv = B.pvals
    B.pvals = lambda name, T: (np.asarray(gate_j["p"], float) if name == "nadig_jurkat" else _pv(name, T))
    D.BANKS["nadig_jurkat"] = ["k562_gwps", "k562_essential", "nadig_hepg2"]
    D.PROFKEY["nadig_jurkat"] = "nadig_jurkat"
    BETA, prof = D.setup()
    prof["nadig_jurkat"] = dict(zip(gate_j["genes"].astype(str), gate_j["ctrl_cpm"].astype(float)))
    F_HERE0 = D.F.HERE
    out, lines = {}, []

    def say(s=""):
        print(s, flush=True)
        lines.append(s)

    for name in SURFS:
        D.F.HERE = J633 if name == "nadig_jurkat" else F_HERE0
        X = D.surface_forms(name, BETA, prof)
        T, S = X["T"], X["S"]
        base = X["P"]["plain"]
        called = C.callset(X["cons0"], S["cnt"], 400, 2200)
        sel = called & (base < 0)
        p = B.pvals(name, T)
        res = {}
        for f in REG:
            Tf = B.pool_T(T, p, f)
            r0 = np.nanmean(D.nmae_cells(base, Tf)[0])
            res[f] = {m: float(np.nanmean(D.nmae_cells(np.where(sel, m * base, base), Tf)[0]) - r0) for m in MS}
        out[name] = res
        say(f"\n== {name}: called DOWN cells/target {sel.sum(1).mean():.0f}; Δnmae raw vs m = 1 (− = better)")
        for f in REG:
            say(f"  f {f:<4g} " + "  ".join(f"m{m:g} {res[f][m]:+.5f}" for m in MS)
                + f"  | argmin m {min(MS, key=lambda m: res[f][m]):g}")
    h = out["h1"]
    gate = (h[1.0][2.0] < 0) and (h[0.35][2.0] < 0) and (min(MS, key=lambda m: h[0.35][m]) in (1.5, 2.0, 2.5, 3.0))
    say(f"\nRETRODICTION (h1 vs the real emitter: m = 2 best, Δ < 0 at f 1 and f 0.35): Δ(m2) f1 {h[1.0][2.0]:+.5f} "
        f"f0.35 {h[0.35][2.0]:+.5f}, argmin f0.35 {min(MS, key=lambda m: h[0.35][m]):g} -> {'PASS' if gate else 'FAIL'}")
    k = out["k562_gwps"]
    verdicts = []
    for f in (0.5, 0.35):
        am = min(MS, key=lambda m: k[f][m])
        if am >= 2.5 and k[f][am] < k[f][2.0] - 0.002:
            verdicts.append("higher")
        elif am <= 1.0 or k[f][2.0] >= 0:
            verdicts.append("no")
        else:
            verdicts.append("x2")
    v = verdicts[0] if verdicts[0] == verdicts[1] else "mixed"
    say(f"DECISION (fixed before the run; {'binding' if gate else 'descriptive: gate FAILED'}): K562 f0.5 / f0.35 -> "
        f"{verdicts} -> {v}")
    (HERE / "work/dose631f.json").write_text(json.dumps({"gate": gate, "verdict": v, "rows": {
        n: {str(f): {str(m): x for m, x in r.items()} for f, r in rr.items()} for n, rr in out.items()}}, indent=1))
    (HERE / "work/dose631f.out").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
