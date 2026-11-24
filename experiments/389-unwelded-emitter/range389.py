#!/usr/bin/env python3
"""389 achievable-range precomputation (kill-check input, before the S2 panel).

Reads 323's recorded cap15 per-gene totals (the S0 gate proves our totals are
bit-identical) + reference keep/ref_lin/b_ctrl/genes + cap15 effects for call-set
masks. Per target: class levels from S, optimistic caps (full class mass; sub-gate
mass only tightens), feasible (Amin, Amax) per gated gene, log2-range distribution
over call-set genes, and predicted 1.25x/1.48x/2.0x clip rates on the prescribed
call set. Writes work/range389.json. Lightweight numpy only.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("NUMBA_NUM_THREADS", "4")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = ROOT / "experiments/323-current-family-validation"
for _p in (str(HERE), str(BASE), str(ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import run323 as R  # noqa: E402
import evaluate323 as E  # noqa: E402
import k_emit389 as K  # noqa: E402
from vcc2026 import rc_emit  # noqa: E402  (CallSet only)

CLASSES3 = ((1 / 3, 16078.75), (1 / 3, 20190.5), (1 / 3, 25330.5))


def q(x):
    x = np.asarray(x, dtype=float)
    return {"n": int(x.size), "min": float(x.min()),
            "p10": float(np.percentile(x, 10)), "p25": float(np.percentile(x, 25)),
            "median": float(np.median(x)), "p75": float(np.percentile(x, 75)),
            "p90": float(np.percentile(x, 90)), "max": float(x.max()),
            "mean": float(x.mean())}


def main():
    WORK = HERE / "work"
    WORK.mkdir(parents=True, exist_ok=True)
    r, _, emit_c = E.load_reference()
    with np.load(BASE / "work/stats-cap15.npz", allow_pickle=False) as z:
        sums = z["sums"]
        assert np.array_equal(z["targets"], r["targets"])
    with np.load(BASE / "work/effects.npz", allow_pickle=False) as z:
        eff_native = z["cap15"]
        egenes = z["genes"]
    with R.model_env("uniform"):
        builder = R.load_module("builder389range",
                                ROOT / "scripts/make_submission_pm_cd4.py")
    effects = R.map_axis(eff_native, egenes, r["genes"])
    gidx = np.flatnonzero(emit_c.keep)
    cs = rc_emit.CallSet(400, 2200, 1.9)
    all_range, call_range, per_s = [], [], {}
    for s in (1.25, 1.48, 2.0):
        per_s[str(s)] = {"n_call": 0, "clip_hi": 0, "clip_lo": 0}
    Lk_first = None
    for i, target in enumerate(r["targets"]):
        Ms = sums[i][gidx].astype(np.int64)
        S = int(sums[i].sum())
        _, Lk, _, nk = K.class_cells(400, list(CLASSES3), S)
        if Lk_first is None:
            Lk_first = Lk.tolist()
        caps = np.broadcast_to((nk * Lk).astype(np.int64), (Ms.shape[0], 3)).copy()
        Amin, Amax = K.arith_range(Ms, Lk, nk, caps)
        pos = Ms > 0
        with np.errstate(divide="ignore", invalid="ignore"):
            lr = np.log2(Amax[pos] / np.maximum(Amin[pos], 1e-12))
        lr = lr[np.isfinite(lr)]
        all_range.append(lr)
        gate = builder.emission_rank_axis(effects[i], r["b_ctrl"], min_cpm=10.)[0]
        mask = K.prescribed_mask(np.where(gate, effects[i], 0.), emit_c, cs)
        call_range.append(lr[mask[pos]])
        for s in (1.25, 1.48, 2.0):
            A_t = K.scale_targets(Ms.astype(float),
                                  np.asarray(emit_c.ref_lin, dtype=float),
                                  S, s, mask)
            steered = np.isfinite(A_t)
            _, _, c_hi, c_lo = K.split_for_A(
                Ms[steered], A_t[steered], Lk, nk, caps[steered])
            per_s[str(s)]["n_call"] += int(steered.sum())
            per_s[str(s)]["clip_hi"] += int(c_hi.sum())
            per_s[str(s)]["clip_lo"] += int(c_lo.sum())
    all_range = np.concatenate(all_range)
    call_range = np.concatenate(call_range)
    out = {
        "classes_raw": [list(c) for c in CLASSES3],
        "levels_first_target": Lk_first,
        "log2_range_all_gated_Mpos": q(all_range),
        "log2_range_callset_Mpos": q(call_range),
        "clip_predicted_optimistic_caps": per_s,
        "note": "caps are full class mass (optimistic: sub-gate mass tightens); "
                "Ms are 323's recorded cap15 totals (S0 gate proves identity).",
    }
    R.write_json(WORK / "range389.json", out)
    print(__import__("json").dumps(out, indent=2), flush=True)


if __name__ == "__main__":
    main()
