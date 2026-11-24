"""241 — PLOLO `pds` for the ten components under THREE renderer bases.

  base `MF`     a15 geometry + `VCC_MASS_FIX=minnorm`      == 200's SHIP.  Phase-0 only.
  base `NOMF`   the same geometry, mass fix UNSET          == the renderer we ship.
  base `SHIPPED` NOMF at the shipped CD4 dose f = 0.34     == `cd4f034`/`h1union`'s dose.

Instrument, verbatim from 200/162: the panel-restricted leave-one-library-out validator, four
donor folds over the real 300-target panel, rendered through production's
`profile_from_effect` into the CHALLENGE control frames A/B/C, scored against the held
library's own pseudobulk.  Nothing here is new physics — `folds200` and `plolo` are imported,
not copied — and every output goes to `experiments/241-ablation-reprice/work`.

FALSIFIER 3.  `kolf2` is pluripotent and gains 3-5x on this instrument; it is reported
SEPARATELY throughout and the headline is the 3-cancer-fold mean.  200 published the 4-fold
mean, so both are carried.

    ONLY=SHIP,noPM uv run python pds241.py            # subset
    uv run python pds241.py
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
E200 = os.path.join(ROOT, "experiments/200-production-ablation")
sys.path.insert(0, os.path.join(ROOT, "experiments/162-retrodicting-validator"))
sys.path.insert(0, E200)
sys.path.insert(0, HERE)

import plolo                                                              # noqa: E402
from plolo import FOLDS, PROD_OPS, Operator, cd4_for, deliver, endpoints, solve_fraction  # noqa: E402
import folds200                                                           # noqa: E402
from folds200 import consensus, frames, parts                             # noqa: E402

os.makedirs(WORK, exist_ok=True)

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
ALL_FOLDS = CANCER + ("kolf2",)

# 200's SHIP, field for field (pds200.py:58-61).
MF = dict(alpha=0.982, gamma=0.65, lam=0.0, cd4_f=0.10, cd4_frac=None, pm=0.50,
          weights=None, mass="minnorm", min_cpm=10.0, op="union_pruned")
NOMF = dict(MF, mass=None)                       # the shipped renderer
SHIPPED = dict(NOMF, cd4_f=0.34)                 # ...at the shipped CD4 dose


def loo_arms(base: dict, tag: str) -> dict:
    """One leave-one-out arm per component of the stack, plus the two controls."""
    A = {
        f"{tag}:BASE":       dict(base),
        f"{tag}:BASE-dup":   dict(base),                       # state-carry control
        f"{tag}:noWhiten":   dict(base, gamma=0.0),
        f"{tag}:noPM":       dict(base, pm=0.0),
        f"{tag}:noCD4":      dict(base, cd4_f=None),
        f"{tag}:alpha_v11":  dict(base, alpha=1.0796),
        f"{tag}:weights_4_1_1": dict(base, weights={"k562_gwps": 4.0, "orion_hct116": 1.0,
                                                    "orion_hek293t": 1.0}),
        f"{tag}:noGate":     dict(base, min_cpm=0.0),
        f"{tag}:noPrune":    dict(base, op="union"),           # kolf2 fold only
    }
    # the mass-fix row, expressed in whichever direction the base makes meaningful
    if base["mass"] is None:
        A[f"{tag}:addMassFix"] = dict(base, mass="minnorm")
    else:
        A[f"{tag}:noMassFix"] = dict(base, mass=None)
    return A


ARMS: dict[str, dict] = {}
ARMS.update(loo_arms(MF, "MF"))
ARMS.update(loo_arms(NOMF, "NOMF"))
ARMS.update(loo_arms(SHIPPED, "SHIPPED"))
OP_SWAP = {a for a in ARMS if a.endswith(":noPrune")}

if os.environ.get("ONLY"):
    _only = set(os.environ["ONLY"].split(","))
    ARMS = {k: v for k, v in ARMS.items() if k in _only}


def effects_for(p, cfg, op):
    """pds200.effects_for verbatim — compose one arm out of the fold's per-donor parts."""
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


def main() -> None:
    t0 = time.time()
    out_path = os.path.join(WORK, "pds241.json")
    rows = json.loads(open(out_path).read()) if os.path.exists(out_path) else []
    done = {(r["fold"], r["arm"], r["space"]) for r in rows}

    for held in ALL_FOLDS:
        production_fold = FOLDS[held]["op_kind"] == "production"
        todo = [a for a in ARMS if (a not in OP_SWAP or production_fold)]
        todo = [a for a in todo
                if any((held, a, f"chal{c}/delivered") not in done for c in "ABC")]
        if not todo:
            print(f"[{held}] all arms already done", flush=True)
            continue
        p = parts(held)
        p["cd4"] = cd4_for(p["targets"])
        truth = np.asarray(p["truth"], dtype=np.float64)[:, p["scored"]]
        by_op: dict[str, list[str]] = {}
        for a in todo:
            by_op.setdefault(ARMS[a]["op"] if production_fold else "held", []).append(a)
        for opkey, names in by_op.items():
            path = (PROD_OPS / f"whitener_{opkey}.npz") if opkey != "held" else p["op_path"]
            op = Operator(path)
            print(f"[{held}] operator {os.path.basename(str(path))} "
                  f"maxeig {op.evals.max():.4f} n {op.evals.size} [{time.time()-t0:.0f}s]",
                  flush=True)
            for a in names:
                cfg = ARMS[a]
                final, aud = effects_for(p, cfg, op)
                for fname, b, gate in frames(p, cfg["min_cpm"]):
                    space = f"{fname}/delivered"
                    if (held, a, space) in done:
                        continue
                    mat = deliver(final, b, gate, cfg["mass"])
                    e = endpoints(mat[:, p["scored"]], truth)
                    rows.append(dict(fold=held, arm=a, space=space, **e, **aud,
                                     gated=int(gate.sum())))
                print(f"  {held:<14} {a:<26} pds "
                      f"{np.mean([r['pds'] for r in rows if r['fold'] == held and r['arm'] == a]):.6f}"
                      f"  [{time.time()-t0:.0f}s]", flush=True)
                open(out_path, "w").write(json.dumps(rows, indent=1))
            del op
        del p
    print(f"wrote {out_path} [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    sys.exit(main())
