"""467 — Arc's programme applied only to targets the design predicts to be STRONG (a per-target gate on the block prior).

Why (462l + 406): on H1 Arc's clean programme buys fid +0.0198 and nmae −0.0148 raw but costs reach −0.039, and the reach
loss is four small-gate targets (n_real 2–17). reach lives in weak targets (406: reach 0.67 weak, 0.06 mid, 0.02 strong),
while nmae ignores targets with fewer than 10 gate genes. So a prior applied to strong targets only should keep most of
the fid/nmae gain and avoid most of the reach cost.

Strength = a build-time function of each target's design row (the emitter's `key` = the DE-channel effect row inside the
emission gate): the number of genes with |key| > tau (tau 0.05 / 0.1 / 0.2) and the row's L2 norm. `record` saves them per
target (seed 0, clean4 prior); `gated` reads them, picks the measure with the best Spearman vs h1-test's n_real (reported,
and fixed before any gated number is seen), and gives a target Arc's programme when its strength is >= the q-quantile,
CLEAN4 otherwise.

usage: RUN426_EMITTER=prod uv run python run467.py record <run435 args>
       RUN426_EMITTER=prod uv run python run467.py gated<q> <run435 args>      (q in 33, 50, 67 = the share gated OFF)
       -> work/h1_<mode>/ (metrics as 435/450/462)
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
assert MODE == "record" or MODE.startswith(("gated", "c4off", "arcoff")), MODE
# gated<q>: Arc on strong targets, CLEAN4 on weak | c4off<q>: CLEAN4 on strong, NO prior (design signs) on weak |
# arcoff<q>: Arc on strong, no prior on weak. q = the share of targets (by strength) treated as weak.
sys.path.insert(0, str(ROOT / "experiments/435-zs-on-coarse"))
import run435 as M  # noqa: E402

WORK = HERE / "work" / f"h1_{MODE}"
WORK.mkdir(parents=True, exist_ok=True)
M.WORK = WORK
M.H.WORK = WORK
_clean4 = M.H.prior_signs_full
REF = np.load(glob.glob(str(ROOT / "experiments/323-*/work/reference.npz"))[0], allow_pickle=True)
H1_TARGETS = set(REF["targets"].astype(str))
STRENGTH_JSON = HERE / "work/strength467.json"
TAUS = (0.05, 0.1, 0.2)
PRIORS: dict[str, np.ndarray] = {}
GATE_ON: set[int] = set()
STRENGTH: dict[int, dict] = {}


def arc_pdown(bank="h1_2025_full"):
    d = ROOT / "data/donors" / bank
    delta = np.load(d / "delta.npy", mmap_mode="r"); perts = [l.strip() for l in open(d / "perts.txt")]
    genes = np.asarray([l.strip() for l in open(d / "genes.txt")])
    rows = [i for i, p in enumerate(perts) if p not in H1_TARGETS]
    X = np.asarray(delta[rows], dtype=np.float64); nz = (X != 0).sum(0); dn = (X < 0).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        pd = np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)
    return genes, pd


def prior_signs_full(genes_full):
    base, cov = _clean4(genes_full)
    vg, pd = arc_pdown()
    pos = {g: j for j, g in enumerate(vg)}
    col = np.array([pos.get(str(g), -1) for g in genes_full])
    val = np.where(col >= 0, pd[np.maximum(col, 0)], np.nan)
    PRIORS["clean4"] = base
    PRIORS["arc"] = np.where(np.isfinite(val), np.where(val > 0.5, -1.0, 1.0), -1.0)
    return base, cov


_carrier = M.carrier_block


def carrier_block(profile, libs, emit_rng, spread, key, i, emit_c):
    k = np.abs(np.asarray(key, dtype=np.float64))
    STRENGTH[i] = {**{f"n_gt_{t:g}": int((k > t).sum()) for t in TAUS}, "l2": float(np.sqrt((k * k).sum()))}
    strong = i in GATE_ON
    if MODE.startswith("gated"):
        M.H.PRIOR_FULL = PRIORS["arc"] if strong else PRIORS["clean4"]
    elif MODE.startswith("c4off"):
        M.H.PRIOR_FULL = PRIORS["clean4"]
    elif MODE.startswith("arcoff"):
        M.H.PRIOR_FULL = PRIORS["arc"]
    else:
        M.H.PRIOR_FULL = PRIORS["clean4"]
    off = MODE.startswith(("c4off", "arcoff")) and not strong
    saved = M.PROTECT_Q
    if off:
        # the weak target keeps its design signs: protect the WHOLE block (protect_q 0 -> only the block's minimum-|E0|
        # gene is flippable), with the rest of the spec (coarse MDE, shrink floor) unchanged. PRIOR_MODE "off" would also
        # drop the MDE, a second variable.
        M.PROTECT_Q = 0.0
    try:
        return _carrier(profile, libs, emit_rng, spread, key, i, emit_c)
    finally:
        M.PROTECT_Q = saved


M.H.prior_signs_full = prior_signs_full
M.carrier_block = carrier_block

if MODE.startswith(("gated", "c4off", "arcoff")):
    q = int(MODE.lstrip("gatedc4ofarc")) / 100.0
    st = json.loads(STRENGTH_JSON.read_text())
    meas = st["chosen_measure"]
    vals = np.array([st["per_target"][str(i)][meas] for i in range(len(st["per_target"]))], dtype=float)
    thr = float(np.quantile(vals, q))
    GATE_ON.update(int(i) for i in np.flatnonzero(vals >= thr))
    print(f"[467] {MODE}: measure {meas}, threshold {thr:.4g}, Arc's programme on {len(GATE_ON)} of "
          f"{len(vals)} targets", flush=True)

if __name__ == "__main__":
    M.main()
    if MODE == "record":
        from scipy.stats import spearmanr
        targets = list(REF["targets"].astype(str)); q_ = REF["q"]
        nreal = {t: int(n) for t, n in zip(targets, (q_ < 0.05).sum(1))}
        stats = sorted(WORK.glob("stats-*.npz"))
        assert stats, "no stats file written by run435"
        order = [str(t) for t in np.load(stats[0], allow_pickle=True)["targets"]]   # the harness's loop order
        per = {str(i): STRENGTH[i] for i in sorted(STRENGTH)}
        nr = np.array([nreal.get(order[i], np.nan) if i < len(order) else np.nan for i in sorted(STRENGTH)], dtype=float)
        rho = {}
        for m in list(next(iter(STRENGTH.values()))):
            v = np.array([STRENGTH[i][m] for i in sorted(STRENGTH)], dtype=float)
            ok = np.isfinite(nr)
            rho[m] = float(spearmanr(v[ok], nr[ok])[0])
        chosen = max(rho, key=rho.get)
        STRENGTH_JSON.write_text(json.dumps({"per_target": per, "targets_in_order": order, "spearman_vs_nreal": rho,
                                             "chosen_measure": chosen, "target_order_source": str(stats[0].name)}, indent=1))
        print(f"[467] strengths recorded for {len(per)} targets; Spearman vs n_real {json.dumps(rho)}; chosen {chosen}",
              flush=True)
