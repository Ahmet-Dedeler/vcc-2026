#!/usr/bin/env python3
"""633 — held-out DOSE reading for 630's rank-neutral UP pass: is the L1-optimal K on cancer lines above H1's?

Why: the L1-optimal magnitude rises with sign accuracy (m* = F^-1((2a-1)/2a)), and 613e measured the uncalled-UP sign
accuracy at 0.70-0.89 on the cancer folds vs 0.54-0.58 on H1. If the folds' optimum sits above K 1.0, H1 under-states
the dose for D/E/F; if it sits at or below, K 1.0 is the October dose. This is 11.1c.2's held-out reading for the dose,
which the board alone cannot give.

Model (612's call-set-aware fold machinery, ds612 / ds612b / ds612c, read-only):
  base   = the "plain" hook form: BETA * cons0 + s * add, s = CLEAN4 sign (tau 0.05) else the key sign; add = 480's
           coarse MDE add (c 0.1, cap 0.1, min count 0.3)
  called = ds612c.callset(cons0, cnt, 400, 2200) (production's head + CPM block); called genes keep the base
  UP pass at dose K: on UNCALLED genes with base > 0, P = max(base, K * 0.24 / sqrt(cnt_K)) (630's target);
           cnt_K = the surface's own control count per cell ("own"), or rescaled to a 20,000-count library ("d20",
           A/B/C's median depth; folds only)
  score  = per-target nmae (ds612.nmae_cells) at 521's pools f 1 / 0.5 / 0.35 (ds612b.pool_T); Δ = mean(K) - mean(base),
           negative = better; target-bootstrap 95 % CI (common480.boot_delta)
Doses K: 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0. Full delivery (an oracle): only the dose's SHAPE transfers, not its size.

RETRODICTION GATE (fixed here, before any fold row is read; rule 3). The h1 row must reproduce the real-emitter oracle
(630/work/oracle630.log, raw Δnmae, 3 seeds) at K 0.5 / 0.75 / 1.0 / 1.5 x f 1 / f 0.35:
  real: f 1   -0.00579 / -0.01085 / -0.01067 / +0.03666
        f 0.35 -0.00708 / -0.01415 / -0.02198 / -0.02800
  PASS iff (a) the sign agrees in all 8 cells, (b) at f 0.35 the model's best K among {0.5, 0.75, 1.0, 1.5} is 1.5 or 1.0
  (real: 1.5), (c) at f 1 the model's K 1.5 is worse than its K 1.0 (real: overshoot).
  FAIL -> the fold rows are reported as descriptive only, and October keeps H1's K 1.0.
Decision rule (fixed here): a cancer fold "prefers a higher dose" if, at both f 0.5 and f 0.35 (d20 rows), its
argmin K is > 1.0 AND that K's Δ beats K 1.0's by more than 0.002 raw.
AMENDMENT (13:49Z, before the first run): a dose optimum is an ABSOLUTE magnitude, and the Orion truths carry a
depth-linked technical blob (AGENTS 6.5: relative direction only, never absolutes; 609: Orion's block gate genes are
72-74 % DOWN). So the Orion rows (HCT116, HEK293T) are DESCRIPTIVE ONLY. The binding cancer surface is k562_gwps. If
the gate PASSES and K562 prefers a higher dose, the October dose question is reopened with a preregistered pair (and
needs a second clean cancer surface before any switch); otherwise K 1.0 stands.
Analysis only; < 4 GB; seconds to minutes.   uv run python experiments/633-updose-folds/dose633.py
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
from common480 import boot_delta  # noqa: E402

KS = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0]
REG = (1.0, 0.5, 0.35)
REAL = {1.0: {0.5: -0.00579, 0.75: -0.01085, 1.0: -0.01067, 1.5: +0.03666},
        0.35: {0.5: -0.00708, 0.75: -0.01415, 1.0: -0.02198, 1.5: -0.02800}}


def main() -> int:
    BETA, prof = D.setup()
    out, lines = {}, []

    def say(s=""):
        print(s, flush=True)
        lines.append(s)

    for name in B.SURFS:
        X = D.surface_forms(name, BETA, prof)
        T, S = X["T"], X["S"]
        base = X["P"]["plain"]
        called = C.callset(X["cons0"], S["cnt"], 400, 2200)
        p = B.pvals(name, T)
        cnt = np.asarray(S["cnt"], float)
        depths = {"own": cnt}
        if name != "h1":
            depths["d20"] = cnt * 20000.0 / float(S["lib_med"])
        sel = (~called) & (base > 0)
        res = {}
        for dname, cK in depths.items():
            mde = 0.24 / np.sqrt(np.maximum(cK, 1e-12))
            rows0 = {f: D.nmae_cells(base, B.pool_T(T, p, f))[0] for f in REG}
            rec = {}
            for K in KS:
                P = np.where(sel, np.maximum(base, K * mde[None, :]), base)
                rec[K] = {}
                for f in REG:
                    r1 = D.nmae_cells(P, B.pool_T(T, p, f))[0]
                    d, lo, hi, _ = boot_delta(rows0[f], r1)       # mean(r1 - rows0): raw Δnmae, negative = better
                    rec[K][f] = (d, lo, hi)
            res[dname] = rec
            say(f"\n== {name} [{dname}] uncalled UP-base cells/target {sel.sum(1).mean():.0f}; Δnmae raw vs base (− = better)")
            for f in REG:
                say(f"  f {f:<4g} " + "  ".join(f"K{K:g} {rec[K][f][0]:+.5f}" for K in KS)
                    + f"  | argmin K {min(KS, key=lambda k: rec[k][f][0]):g}")
        out[name] = res
    # retrodiction gate on h1
    h = out["h1"]["own"]
    sign_ok = all(np.sign(h[K][f][0]) == np.sign(REAL[f][K]) for f in REAL for K in REAL[f])
    best35 = min((0.5, 0.75, 1.0, 1.5), key=lambda k: h[k][0.35][0])
    over1 = h[1.5][1.0][0] > h[1.0][1.0][0]
    gate = sign_ok and best35 in (1.0, 1.5) and over1
    say(f"\nRETRODICTION (h1 vs the real-emitter oracle): signs {sign_ok}, best K at f 0.35 {best35:g} (real 1.5), "
        f"K1.5 worse than K1.0 at f 1 {over1} -> {'PASS' if gate else 'FAIL'}")
    for f in REAL:
        say(f"  f {f:<4g} model " + " ".join(f"K{K:g} {h[K][f][0]:+.5f}" for K in REAL[f])
            + "   real " + " ".join(f"K{K:g} {REAL[f][K]:+.5f}" for K in REAL[f]))
    prefer = {}
    for name in ("orion_hct116", "orion_hek293t", "k562_gwps"):
        r = out[name]["d20"]
        pf = True
        for f in (0.5, 0.35):
            kb = min(KS, key=lambda k: r[k][f][0])
            pf &= (kb > 1.0) and (r[kb][f][0] < r[1.0][f][0] - 0.002)
        prefer[name] = bool(pf)
    n_pref = int(prefer["k562_gwps"])                       # amendment: Orion rows are descriptive only
    say(f"\nDECISION (rule fixed in the docstring; {'binding' if gate else 'descriptive only: gate FAILED'}): folds preferring "
        f"a dose above 1.0 (d20): {prefer} (binding: k562_gwps only; Orion descriptive) -> "
        f"{'reopen the October dose with a preregistered pair' if (gate and n_pref >= 1) else 'K 1.0 stands'}")
    (HERE / "work/dose633.json").write_text(json.dumps({"gate": gate, "prefer": prefer, "rows": {
        n: {d: {str(K): {str(f): v for f, v in fr.items()} for K, fr in r.items()} for d, r in v.items()}
        for n, v in out.items()}}, indent=1, default=float))
    (HERE / "work/dose633.out").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
