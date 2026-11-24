"""214 — PLOLO `pds` and transfer cosine for the PerturbMatch x CD4-dose grid and the
rounding rule, priced against the CURRENT incumbent `cd4f034`.

INSTRUMENT.  162's panel-restricted leave-one-library-out validator, exactly as 200 and 208
used it: four donor folds over the real 300-target panel, rendered through
`profile_from_effect` into the challenge control frames A/B/C, scored against the held
library's own pseudobulk.  It is the only `pds` instrument validated against a board point
(6/6 on sign, LOO sigma 0.00149 raw).  `debench.Surface.pds()` is never used (188 Phase 3).

WHAT IS NEW.  200 ablated PerturbMatch at ONE CD4 dose (`f = 0.10`).  207 then tripled the
dose and `cd4f034` shipped `f = 0.34`.  064 composed PerturbMatch WITH CD4 at `f = 0.05` and
found them additive.  Nobody has measured whether the two are SUBSTITUTES -- i.e. whether
PerturbMatch's value falls as the CD4 dose rises.  This sweeps the 2-D grid.

THE INCUMBENT'S RENDERER.  `submissions/210-cd4f034.h5ad.meta.json` records
`mass_fix: null` -- a23's mass fix was NOT carried forward past a23.  200's `SHIP` row is
a23 (mass = "minnorm"); the incumbent is mass = None.  Both are carried here, because the
renderer is one of the three candidate explanations for PerturbMatch's sign.

KNOWN-ANSWER ROWS (rule 10; every one aborts):
  KA1  200's `SHIP` (= a23) absolute chal-mean per fold:
       orion_hct116 0.791823 / orion_hek293t 0.777375 / k562_gwps 0.789121 / kolf2 0.610128
  KA2  200's mass-fix LOO, `noMassFix` - `SHIP` = -0.004014 scaled
  KA3  200's PerturbMatch LOO, `noPM` - `SHIP` = -0.011187 scaled (= +0.001885 overall)
  KA4  200's whitening LOO, `noWhiten` - `SHIP`  (+0.025593 overall) and the CPM-gate LOO
  KA5  196/208's rounding rows on chalA at the a23 renderer, `f = 0.10`, pm 0.5:
       floor 0.739016, ham 0.742224, inf 0.742137
  KA6  169's cached chalA rows at m = inf: 0.790712 / 0.790896 / 0.771621 / 0.615317

    uv run python pds214.py
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
P162 = os.path.join(ROOT, "experiments/162-retrodicting-validator")
P200 = os.path.join(ROOT, "experiments/200-production-ablation")
P196 = os.path.join(ROOT, "experiments/196-emitted-depth")
sys.path.insert(0, P162)
sys.path.insert(0, P200)
sys.path.insert(0, P196)
os.makedirs(WORK, exist_ok=True)

import plolo                                                            # noqa: E402
from plolo import (FOLDS, PROD_OPS, Operator, cd4_for, deliver,         # noqa: E402
                   endpoints, mass_consistent_effect, solve_fraction)
from folds200 import consensus, frames, parts                           # noqa: E402
from ham196 import hamilton                                             # noqa: E402

TS = 50_000.0
PDS_SPAN = 0.95228 - 0.50020        # debench ANCHOR, board-fitted
K_RULE13 = 1.011                    # 162's LOO calibration (rule 13)
SIGMA_SCALED = 0.001494             # 162's leave-one-arm-out RMS on scaled `pds`
# per-context emitted library size, 196 `recon196.main` known-answer row
NL = {"chalA": 400 * 20_109, "chalB": 400 * 19_946, "chalC": 400 * 20_034}

BASE = dict(alpha=0.982, gamma=0.65, lam=0.0, cd4_f=0.10, cd4_frac=None, pm=0.50,
            weights=None, mass=None, min_cpm=10.0, op="union_pruned")


def arm(**kw):
    return dict(BASE, **kw)


DOSES = (None, 0.05, 0.10, 0.20, 0.34, 0.50, 0.80)
PMS = (0.0, 0.50)

def _tag(pm, f):
    return f"pm{int(round(pm * 100)):03d}_f{0 if f is None else int(round(f * 100)):03d}"


ARMS: dict[str, dict] = {}
# --- 1. the substitution grid: PerturbMatch x CD4 dose, incumbent renderer ---------------
for f in DOSES:
    for pm in PMS:
        ARMS[f"g_{_tag(pm, f)}"] = arm(pm=pm, cd4_f=f)
# --- 2. the same two doses on a23's renderer (200's own convention) ----------------------
for f in (0.10, 0.34):
    for pm in PMS:
        ARMS[f"mf_{_tag(pm, f)}"] = arm(pm=pm, cd4_f=f, mass="minnorm")
# --- 3. PerturbMatch's shape in the blend fraction, at both doses -------------------------
for f in (0.10, 0.34):
    for pm in (0.25, 0.75, 1.00):
        ARMS[f"g_{_tag(pm, f)}"] = arm(pm=pm, cd4_f=f)
# --- 4. the operator: is PerturbMatch's value a function of the whitener? ----------------
for f in (None, 0.10, 0.34):
    for pm in PMS:
        ARMS[f"nw_{_tag(pm, f)}"] = arm(pm=pm, cd4_f=f, gamma=0.0)
# --- 5. 200's own leave-one-out rows, for the known-answer gates -------------------------
ARMS["SHIP"] = arm(mass="minnorm")                       # 200's SHIP == a23
ARMS["noMassFix"] = arm(mass=None)                       # == a15 == g_pm050_f010
ARMS["noPM_mf"] = arm(mass="minnorm", pm=0.0)
ARMS["noWhiten_mf"] = arm(mass="minnorm", gamma=0.0)
ARMS["gate20_mf"] = arm(mass="minnorm", min_cpm=20.0)
ARMS["noCD4_mf"] = arm(mass="minnorm", cd4_f=None)

# Which (arm, frame) pairs also get the two integer rounding rules.  196/208 priced the
# rounding rule on chalA only; B and C are added here with each context's own emitted L.
ROUND_ARMS = ("g_pm050_f010", "g_pm000_f010", "g_pm050_f034", "g_pm000_f034",
              "mf_pm050_f010", "mf_pm000_f010", "mf_pm050_f034", "mf_pm000_f034")
assert set(ROUND_ARMS) <= set(ARMS), sorted(set(ROUND_ARMS) - set(ARMS))

if os.environ.get("ONLY"):
    keep = set(os.environ["ONLY"].split(","))
    ARMS = {k: v for k, v in ARMS.items() if k in keep}


def effects_for(p, cfg, op):
    """One arm's effect matrix out of the fold's per-donor parts (200's `effects_for`)."""
    pre = consensus(p, pm_blend=cfg["pm"], weights=cfg["weights"]) * cfg["alpha"]
    if cfg["gamma"] is None:
        whitened = np.asarray(pre, dtype=np.float64).copy()
    else:
        whitened = op.apply(pre, p["wcols"], cfg["gamma"], cfg["lam"])
    aud = dict(pre_median_norm=float(np.median(np.linalg.norm(pre, axis=1))),
               whitened_median_norm=float(np.median(np.linalg.norm(whitened, axis=1))))
    if cfg["cd4_f"] is not None:
        frac, per_unit = solve_fraction(whitened, p["cd4"], cfg["cd4_f"])
        final, _ = plolo.E062.add_post_whiten_residual(whitened, p["cd4"], frac, p["wcols"])
        aud["cd4_fraction"] = float(frac)
        aud["cd4_f_per_unit"] = float(per_unit)
    elif cfg["cd4_frac"] is not None:
        frac = float(cfg["cd4_frac"])
        final, _ = plolo.E062.add_post_whiten_residual(whitened, p["cd4"], frac, p["wcols"])
        aud["cd4_fraction"] = frac
    else:
        final = whitened
        aud["cd4_fraction"] = 0.0
    final = np.array(final, dtype=np.float64, copy=True)
    final[:, p["panel_positions"]] = 0.0
    aud["final_median_norm"] = float(np.median(np.linalg.norm(final, axis=1)))
    return final, aud


def deliver_rounded(effects, b, gate, mass, nL, rule):
    """Render, then quantise the delivered per-gene panel totals with `rule`.

    `floor` is `emit_counts`' shipped `totals = floor(n L p)`; `ham` is largest-remainder
    over the same `n L p` (196's `hamilton`).  Identical to `pdsham196.deliver_r` except
    that the mass-fix flag is a parameter here rather than hardcoded to "minnorm".
    """
    out = np.empty_like(effects)
    rel = np.zeros(effects.shape[0])
    for i in range(effects.shape[0]):
        g = np.where(gate, effects[i], 0.0)
        if mass is not None:
            g = mass_consistent_effect(g, b, gate, mode=mass)
        p = np.clip(np.expm1(b + g), 0.0, None) / TS
        p = p / float(p.sum())
        x = nL * p
        T = np.floor(x) if rule == "floor" else hamilton(x).astype(np.float64)
        nz = x > 0
        rel[i] = float(np.abs((T[nz] - x[nz]) / x[nz]).mean())
        out[i] = np.log1p(TS * (T / T.sum())) - b
    return out, float(rel.mean())


def main() -> None:
    t0 = time.time()
    out_path = os.path.join(WORK, "pds214.json")
    rows = json.loads(open(out_path).read()) if os.path.exists(out_path) else []
    done = {(r["fold"], r["arm"], r["frame"], r["rule"]) for r in rows}

    for held in FOLDS:
        p = parts(held)
        p["cd4"] = cd4_for(p["targets"])
        truth = np.asarray(p["truth"], dtype=np.float64)[:, p["scored"]]
        todo = [a for a in ARMS
                if any((held, a, f, r) not in done
                       for f in ("chalA", "chalB", "chalC")
                       for r in (("inf", "floor", "ham") if a in ROUND_ARMS else ("inf",)))]
        if not todo:
            print(f"[{held}] all arms done", flush=True)
            del p
            continue
        op = Operator(p["op_path"])
        print(f"[{held}] {len(p['targets'])} targets, {len(p['scored'])} scored, "
              f"operator {os.path.basename(str(p['op_path']))} maxeig {op.evals.max():.4f} "
              f"n {op.evals.size} [{time.time()-t0:.0f}s]", flush=True)
        for a in todo:
            cfg = ARMS[a]
            final, aud = effects_for(p, cfg, op)
            for fname, b, gate in frames(p, cfg["min_cpm"]):
                # effect-level ("intended") -- 058/064's instrument, no renderer
                if (held, a, fname, "intended") not in done:
                    ei = endpoints(np.where(gate[None, :], final, 0.0)[:, p["scored"]], truth)
                    rows.append(dict(fold=held, arm=a, frame=fname, rule="intended",
                                     **ei, **aud, gated=int(gate.sum())))
                if (held, a, fname, "inf") not in done:
                    mat = deliver(final, b, gate, cfg["mass"])
                    e = endpoints(mat[:, p["scored"]], truth)
                    rows.append(dict(fold=held, arm=a, frame=fname, rule="inf",
                                     **e, **aud, gated=int(gate.sum())))
                    print(f"  {held:<14} {a:<18} {fname} inf   pds {e['pds']:.6f} "
                          f"cos {e['cos']:.6f}", flush=True)
                if a in ROUND_ARMS:
                    for rule in ("floor", "ham"):
                        if (held, a, fname, rule) in done:
                            continue
                        mat, rel = deliver_rounded(final, b, gate, cfg["mass"],
                                                   NL[fname], rule)
                        e = endpoints(mat[:, p["scored"]], truth)
                        rows.append(dict(fold=held, arm=a, frame=fname, rule=rule,
                                         **e, **aud, gated=int(gate.sum()),
                                         mean_rel_round_err=rel))
                        print(f"  {held:<14} {a:<18} {fname} {rule:<5} pds "
                              f"{e['pds']:.6f} cos {e['cos']:.6f}", flush=True)
            open(out_path, "w").write(json.dumps(rows, indent=1))
        del op, p
    open(out_path, "w").write(json.dumps(rows, indent=1))
    print(f"wrote {out_path}  ({len(rows)} rows) [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    sys.exit(main())
