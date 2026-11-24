"""470b — for the STRONG two thirds, Arc's programme from its non-essential vs essential perturbations (H1; 467's gate)

A copy of run470.py (which chain470 executes; left untouched) with two more priors. The 2026 panel is essential-free
(0/300 K562-essential) and Arc's H1 bank is 78/300 essential, so is a programme over the non-essential perturbations a
better prior for an essential-free panel's strong targets?
  arcne   P_down over the non-h1-test, non-K562-essential perturbations
  arcess  P_down over the non-h1-test K562-essential perturbations
(modes s<prior>: <prior> on the strong two thirds, CLEAN4 on the weak third; report the non-essential subset.)

Original docstring follows.
470 — which block prior do design-WEAK targets want? (H1, production emitter; 467's q33 strength gate)

Why: 467/467b on H1, with Arc's programme on the design-strong two thirds, the weak third's reach ranks CLEAN4 (+0.0015
raw) > the design's own signs (−0.0099) > Arc's all-perturbation programme (−0.0388). reach lives in weak targets (406).
CLEAN4 is much more down-skewed than Arc's programme (66–68 % vs 51–52 % down on the A/B/C blocks), and 361 found all-down
buys fid on H1. So: is the weak third's best prior "down", a consensus, or a programme estimated from the bank's own weak
perturbations (the lab signature lives in weak responses, 466)?

modes (q33 and the strength measure are 467's, read from 467's strength record):
  w<prior>   Arc's all-perturbation programme on the strong two thirds, <prior> on the weak third
  s<prior>   <prior> on the strong two thirds, CLEAN4 on the weak third
priors:
  c4    CLEAN4 (the incumbent; `wc4` == 467's gated33, the known answer)
  down  every gene down
  arc   Arc's programme over all of h1_2025_full's non-h1-test perturbations (467's)
  arcw  the same over the WEAKER half of those perturbations (response size below the median)
  arcs  the same over the STRONGER half
Response size of a bank perturbation = count of |delta| > 0.2 over the 4,000 genes most expressed in h1_2025's control.

usage: RUN426_EMITTER=prod uv run python run470.py <mode> <run435 args>   -> work/h1_<mode>/ (metrics as 435/467)
"""
from __future__ import annotations

import glob
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
os.environ.setdefault("RUN426_EMITTER", "prod")
MODE = sys.argv.pop(1)
SIDE, PNAME = MODE[0], MODE[1:]
assert SIDE in "ws" and PNAME in ("c4", "down", "arc", "arcw", "arcs", "arcne", "arcess"), MODE
sys.path.insert(0, str(ROOT / "experiments/435-zs-on-coarse"))
import run435 as M  # noqa: E402

WORK = HERE / "work" / f"h1_{MODE}"
WORK.mkdir(parents=True, exist_ok=True)
M.WORK = WORK
M.H.WORK = WORK
_clean4 = M.H.prior_signs_full
REF = np.load(glob.glob(str(ROOT / "experiments/323-*/work/reference.npz"))[0], allow_pickle=True)
H1_TARGETS = set(REF["targets"].astype(str))
STRENGTH_JSON = ROOT / "experiments/467-strength-gated-programme/work/strength467.json"
Q = 0.33
PRIORS: dict[str, np.ndarray] = {}
GATE_ON: set[int] = set()
SPLIT: dict[str, object] = {}


def arc_pdowns(bank="h1_2025_full"):
    """P_down over the bank's non-h1-test perturbations: all, the weaker half and the stronger half by response size."""
    d = ROOT / "data/donors" / bank
    delta = np.load(d / "delta.npy", mmap_mode="r"); perts = [l.strip() for l in open(d / "perts.txt")]
    genes = np.asarray([l.strip() for l in open(d / "genes.txt")])
    rows = [i for i, p in enumerate(perts) if p not in H1_TARGETS]
    X = np.asarray(delta[rows], dtype=np.float64)
    b = np.load(ROOT / "data/donors/h1_2025/b_ctrl.npy")
    assert [l.strip() for l in open(ROOT / "data/donors/h1_2025/genes.txt")] == list(genes), "bank gene axes differ"
    expressed = np.argsort(-b)[:4000]
    size = (np.abs(X[:, expressed]) > 0.2).sum(1)
    med = float(np.median(size))
    strong = size > med
    SPLIT.update({"n_perts": int(len(rows)), "median_size": med, "n_strong": int(strong.sum()), "n_weak": int((~strong).sum()),
                  "size_q": [float(x) for x in np.quantile(size, [0, .25, .5, .75, 1])]})

    def pdown(Y):
        nz = (Y != 0).sum(0); dn = (Y < 0).sum(0)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)
    ess = {l.strip() for l in open(ROOT / "data/donors/k562_essential/perts.txt") if l.strip()}
    is_ess = np.array([perts[i] in ess for i in rows])
    SPLIT.update({"n_essential": int(is_ess.sum()), "n_nonessential": int((~is_ess).sum())})
    return genes, {"arc": pdown(X), "arcw": pdown(X[~strong]), "arcs": pdown(X[strong]),
                   "arcne": pdown(X[~is_ess]), "arcess": pdown(X[is_ess])}


def prior_signs_full(genes_full):
    base, cov = _clean4(genes_full)
    vg, pds = arc_pdowns()
    pos = {g: j for j, g in enumerate(vg)}
    col = np.array([pos.get(str(g), -1) for g in genes_full])
    PRIORS["c4"] = base
    PRIORS["down"] = -np.ones_like(base)
    for k, pd in pds.items():
        val = np.where(col >= 0, pd[np.maximum(col, 0)], np.nan)
        PRIORS[k] = np.where(np.isfinite(val), np.where(val > 0.5, -1.0, 1.0), -1.0)
    SPLIT["down_frac"] = {k: float((v < 0).mean()) for k, v in PRIORS.items()}
    SPLIT["agree_with_c4"] = {k: float((v == base).mean()) for k, v in PRIORS.items()}
    print(f"[470] split {json.dumps(SPLIT)}", flush=True)
    return base, cov


_carrier = M.carrier_block


def carrier_block(profile, libs, emit_rng, spread, key, i, emit_c):
    strong = i in GATE_ON
    if SIDE == "w":
        M.H.PRIOR_FULL = PRIORS["arc"] if strong else PRIORS[PNAME]
    else:
        M.H.PRIOR_FULL = PRIORS[PNAME] if strong else PRIORS["c4"]
    return _carrier(profile, libs, emit_rng, spread, key, i, emit_c)


M.H.prior_signs_full = prior_signs_full
M.carrier_block = carrier_block

st = json.loads(STRENGTH_JSON.read_text())
meas = st["chosen_measure"]
vals = np.array([st["per_target"][str(i)][meas] for i in range(len(st["per_target"]))], dtype=float)
thr = float(np.quantile(vals, Q))
GATE_ON.update(int(i) for i in np.flatnonzero(vals >= thr))
print(f"[470] {MODE}: measure {meas}, threshold {thr:.4g}, strong {len(GATE_ON)} of {len(vals)} targets", flush=True)

if __name__ == "__main__":
    M.main()
    (WORK / "split470.json").write_text(json.dumps(SPLIT, indent=1))
