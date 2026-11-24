"""200 — leave-one-out and add-one-in on the SHIPPED stack, priced on PLOLO `pds`.

INSTRUMENT.  162's panel-restricted leave-one-library-out validator: four donor folds over the
real 300-target panel, rendered through production's `profile_from_effect` into the CHALLENGE
control frames A/B/C, scored against the held library's own pseudobulk.  It retrodicts board
`pds` 6/6 on sign with a leave-one-arm-out RMS of 0.001494 scaled, and it is the ONLY `pds`
instrument in this project validated against a board point.  `debench.Surface.pds()` is a
different statistic (log2 fold change, CPM-5 gate, 85 decoys) and 188 Phase 3's standing rule
forbids sizing a `pds` arm with it.

WHAT IS NEW HERE.  162 ran ARMS -- historical configurations, each differing from its own base.
This runs ABLATIONS: one component removed from the full shipped stack at a time, and the
mirror-image add-one-in from a bare consensus.  Sequential addition can only ever measure a
component against the stack that existed below it; the two directions together separate
marginal value from interaction.

KNOWN-ANSWER ROWS (rule 10's corollary -- a screen with no row whose answer is already known
cannot detect a whole class of silent bug).  All three abort the run:

  1. `SHIP` (a15 geometry + a23's mass fix) must reproduce 162's own `a23` rows fold by fold
     and frame by frame, from `experiments/162-retrodicting-validator/work/arms.json`.  This
     is the row that certifies the per-donor decomposition in `folds200.py`: if the
     recombination of the cached per-source effects were not the cached `pre`, this is where
     it shows.
  2. `noMassFix` (= a15) minus `SHIP` (= a23) must reproduce 162's measured mass-fix price,
     +0.004014 scaled, which 169 independently reproduced to six decimals.
  3. 169's q = 1.00 identity: 16/16 bit-identical rows, worst |d pds| 0.00e+00.  Reproduced
     here as the |d| between this script's SHIP rows and 162's, on the same 16 (fold, frame)
     cells.

No cell-eval2, no build, no submission, nothing written outside this directory.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
P162 = os.path.join(ROOT, "experiments/162-retrodicting-validator")
sys.path.insert(0, P162)
sys.path.insert(0, HERE)

import plolo                                                        # noqa: E402
from plolo import FOLDS, PROD_OPS, Operator, cd4_for, deliver, endpoints, solve_fraction  # noqa: E402
import folds200                                                     # noqa: E402
from folds200 import consensus, frames, parts                       # noqa: E402

PDS_SPAN = 0.95228 - 0.50020          # debench ANCHOR, the board-fitted pds anchors
K_RULE13 = 1.011                      # rule 13's calibration (162's LOO scale without a23)
K_ALLARM = 0.899                      # 162's all-arm fitted scale, reported as a sensitivity

# The shipped stack, field by field, exactly as `submissions/a15-cd4f010.vcc.meta.json` and
# `submissions/a23-massfix.h5ad.meta.json` record it (see work/inputs200.log).
SHIP = dict(alpha=0.982, gamma=0.65, lam=0.0, cd4_f=0.10, cd4_frac=None, pm=0.50,
            weights=None, mass="minnorm", min_cpm=10.0, op="union_pruned")
# A bare consensus: every optional component off, amplitude and gate held at the shipped
# value so the add-one-in ladder moves ONE thing per rung.
BARE = dict(SHIP, gamma=0.0, cd4_f=None, cd4_frac=None, pm=0.0, mass=None)


def arm(**kw) -> dict:
    return dict(SHIP, **kw)


def bare(**kw) -> dict:
    return dict(BARE, **kw)


ARMS: dict[str, dict] = {
    # ---- the incumbent and its two known-answer neighbours -------------------------------
    "SHIP":            arm(),
    "SHIP-dup":        arm(),                              # field-identical: state-carry check
    "noMassFix":       arm(mass=None),                     # == a15;  KA: -0.004014 scaled
    # ---- leave one out -------------------------------------------------------------------
    "noWhiten":        arm(gamma=0.0),
    "noPM":            arm(pm=0.0),
    "noCD4":           arm(cd4_f=None),
    "cd4_a14":         arm(cd4_f=None, cd4_frac=0.05),     # == a14's CD4;  KA: a14 -> a15
    "alpha_v11":       arm(alpha=1.0796),                  # undo 187/a14's calibration
    "weights_4_1_1":   arm(weights={"k562_gwps": 4.0, "orion_hct116": 1.0,
                                    "orion_hek293t": 1.0}),
    "weights_2_1_1":   arm(weights={"k562_gwps": 2.0, "orion_hct116": 1.0,
                                    "orion_hek293t": 1.0}),   # 175's K562-heavy mix
    "noGate":          arm(min_cpm=0.0),
    "gate5":           arm(min_cpm=5.0),                   # the pre-088 value
    "gate20":          arm(min_cpm=20.0),
    "noPrune":         arm(op="union"),                    # kolf2 fold only
    # ---- the gamma neighbourhood (170 says 0.65 is a six-member local optimum) ------------
    "gamma_050":       arm(gamma=0.50),
    "gamma_080":       arm(gamma=0.80),
    "a19ridge":        arm(gamma=0.80, lam=0.007),         # the a19 retrodiction
    # ---- add one in, from the bare consensus ---------------------------------------------
    "BARE":            bare(),
    "BARE+whiten":     bare(gamma=0.65),
    "BARE+pm":         bare(pm=0.50),
    "BARE+cd4":        bare(cd4_f=0.10),
    "BARE+massfix":    bare(mass="minnorm"),
    "BARE+whiten_unpruned": bare(gamma=0.65, op="union"),  # kolf2 only: AOI prune comparison
    # ---- the four board-measured pairs, reconstructed EXACTLY as the board ran them -------
    # (audit 181 rule 1: an arm family with a board point must retrodict it first.)
    "a14":             arm(mass=None, cd4_f=None, cd4_frac=0.05),   # a14 -> a15 is the CD4 arm
    "a19_exact":       arm(mass=None, gamma=0.80, lam=0.007),       # a19 vs a15
    "a11":             arm(mass=None, cd4_f=None, cd4_frac=0.05, alpha=0.838),
    "a12":             arm(mass=None, cd4_f=None, cd4_frac=0.05, alpha=1.0108),
    # control: is `gamma = 0` through the operator really the identity?  The eigenbasis is
    # complete on the operator columns, so it should be -- measured, not assumed.
    "noOperator":      arm(gamma=None),
}
# Arms that need the un-pruned operator can only be expressed on the fold whose operator is
# the production one (162: `kolf2` is the leak-free fold and the only one where a
# `union -> union_pruned` swap is expressible at all).
OP_SWAP = {"noPrune", "BARE+prune"}

if os.environ.get("ONLY"):
    _only = set(os.environ["ONLY"].split(","))
    ARMS = {k: v for k, v in ARMS.items() if k in _only}


def effects_for(p, cfg, op) -> tuple[np.ndarray, dict]:
    """Compose one arm's effect matrix out of the fold's per-donor parts."""
    pre = consensus(p, pm_blend=cfg["pm"], weights=cfg["weights"]) * cfg["alpha"]
    if cfg["gamma"] is None:                 # the operator is not applied at all
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
    out_path = os.path.join(WORK, "pds200.json")
    rows = json.loads(open(out_path).read()) if os.path.exists(out_path) else []
    done = {(r["fold"], r["arm"], r["space"]) for r in rows}

    for held in FOLDS:
        p = parts(held)
        p["cd4"] = cd4_for(p["targets"])
        truth = np.asarray(p["truth"], dtype=np.float64)[:, p["scored"]]
        production_fold = FOLDS[held]["op_kind"] == "production"
        todo = [a for a in ARMS if (a not in OP_SWAP or production_fold)]
        todo = [a for a in todo
                if any((held, a, f"chal{c}/delivered") not in done for c in "ABC")]
        if not todo:
            print(f"[{held}] all arms already done", flush=True)
            continue
        # group by operator so each 2.5 GB eigenbasis is loaded at most once per fold
        by_op: dict[str, list[str]] = {}
        for a in todo:
            key = ARMS[a]["op"] if production_fold else "held"
            by_op.setdefault(key, []).append(a)
        for opkey, names in by_op.items():
            path = (PROD_OPS / f"whitener_{opkey}.npz") if opkey != "held" else p["op_path"]
            op = Operator(path)
            print(f"[{held}] operator {os.path.basename(str(path))} "
                  f"maxeig {op.evals.max():.4f} n {op.evals.size}", flush=True)
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
                    print(f"  {held:<14} {a:<16} {space:<18} pds {e['pds']:.6f} "
                          f"cos {e['cos']:.5f} sse {e['sse']:.6f} sgnw {e['sgn_w']:.5f}",
                          flush=True)
                open(out_path, "w").write(json.dumps(rows, indent=1))
            del op
        del p
    print(f"wrote {out_path}")


if __name__ == "__main__":
    sys.exit(main())
