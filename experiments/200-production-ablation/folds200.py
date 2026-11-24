"""200 — PLOLO folds, decomposed into their PER-DONOR parts so a component can be removed.

WHY THIS EXISTS
---------------
162's `build_fold` bakes three of the components this experiment has to ablate straight into
its cached `pre`:  the PerturbMatch K562 replacement, the equal donor weights, and the
consensus itself.  You cannot turn any of them off from the cache.  So this module rebuilds
the same folds from the same E062 code path and caches the pieces BEFORE they are combined:

    per-source effects (raw K562 / HCT116 / HEK293T), the PerturbMatch effects and the
    shared-column index, the held library's truth, its b_ctrl, the CPM gate, the scored
    columns, the operator columns, and the CD4 residual.

`compose()` then assembles an arm out of those pieces.  With every component at its shipped
value it must reproduce 162's own `arms.json` rows -- that is the known-answer gate in
`pds200.py`, and it is the only thing that certifies the decomposition.

NOTHING in 162 is modified: its `plolo.py` is imported, its fold cache is read where it is
already correct (truth / b_ctrl / gate / scored / wcols), and every artifact this module
writes lands in `experiments/200-production-ablation/work/`.
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

import plolo  # noqa: E402
from plolo import (CD4_ROOT, DONOR_PATHS, E062, E064, FOLDS, OFFICIAL, PANEL,  # noqa: E402
                   PANEL_SET, PM_BLEND, PM_PATH, Operator, build_fold, cd4_for)

os.makedirs(WORK, exist_ok=True)

# Every source a fold can carry, in a fixed order so a weight vector is unambiguous.
SOURCE_ORDER = ("k562_gwps", "orion_hct116", "orion_hek293t")


def _panel_positions() -> np.ndarray:
    pos = {g: i for i, g in enumerate(OFFICIAL)}
    return np.asarray([pos[g] for g in PANEL if g in pos], dtype=np.int64)


def parts(held: str, limit_gib: float = 24.0) -> dict:
    """Per-donor pieces for one PLOLO fold, cached.  Same E062 calls as 162's `build_fold`."""
    cache = os.path.join(WORK, f"parts-{held}.npz")
    fold = build_fold(held)                      # truth / b_ctrl / gate / scored / wcols
    targets = tuple(fold["targets"])
    spec = FOLDS[held]
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        if tuple(str(t) for t in z["targets"]) == targets:
            out = dict(fold)
            out["src"] = {n: z[f"src_{n}"].astype(np.float64) for n in spec["sources"]}
            out["avail"] = {n: z[f"avail_{n}"] for n in spec["sources"]}
            out["mean_on_axis"] = {n: z[f"mean_{n}"] for n in spec["sources"]}
            out["pm_effects"] = z["pm_effects"].astype(np.float64)
            out["pm_avail"] = z["pm_avail"]
            out["shared"] = z["shared"]
            out["audit"] = json.loads(str(z["audit"]))
            print(f"[{held}] parts cache hit", flush=True)
            return out

    pp = _panel_positions()
    banks = {n: E062.Bank.load(n, DONOR_PATHS[n]) for n in spec["sources"]}
    prepared = {n: E062.prepare_main_source(b, targets, OFFICIAL, pp,
                                            chunk_rows=256, max_rss_gib=limit_gib)
                for n, b in banks.items()}
    pm_bank = E064.load_perturbmatch(PM_PATH)
    pm_prep = E062.prepare_main_source(pm_bank, targets, OFFICIAL, pp,
                                       chunk_rows=256, max_rss_gib=limit_gib)
    shared = (E064.shared_k562_pm_columns(banks["k562_gwps"], pm_bank, OFFICIAL, PANEL_SET)
              if "k562_gwps" in banks else np.zeros(0, dtype=np.int64))

    audit = dict(held=held, sources=list(spec["sources"]), targets=len(targets),
                 shared_columns=int(np.asarray(shared).size),
                 per_source_available={n: int(prepared[n].target_available.sum())
                                       for n in spec["sources"]},
                 pm_available=int(pm_prep.target_available.sum()),
                 raw_median_norm={n: float(prepared[n].raw_median_norm)
                                  for n in spec["sources"]})
    np.savez_compressed(
        cache, targets=np.asarray(targets), audit=json.dumps(audit),
        pm_effects=pm_prep.effects.astype(np.float32),
        pm_avail=pm_prep.target_available, shared=np.asarray(shared),
        **{f"src_{n}": prepared[n].effects.astype(np.float32) for n in spec["sources"]},
        **{f"avail_{n}": prepared[n].target_available for n in spec["sources"]},
        **{f"mean_{n}": prepared[n].mean_on_axis for n in spec["sources"]},
    )
    print(f"[{held}] parts built: {json.dumps(audit)}", flush=True)
    out = dict(fold)
    out["src"] = {n: prepared[n].effects for n in spec["sources"]}
    out["avail"] = {n: prepared[n].target_available for n in spec["sources"]}
    out["mean_on_axis"] = {n: prepared[n].mean_on_axis for n in spec["sources"]}
    out["pm_effects"] = pm_prep.effects
    out["pm_avail"] = pm_prep.target_available
    out["shared"] = np.asarray(shared)
    out["audit"] = audit
    return out


def consensus(p: dict, *, pm_blend: float, weights: dict[str, float] | None) -> np.ndarray:
    """The pre-whitening effect matrix, with PerturbMatch and the donor weights as knobs.

    `weights=None` is production's equal weighting, which is `E062.build_baseline` verbatim.
    A weight dict takes the same weighted mean, renormalised over the donors that actually
    carry each target -- the rule `MultiDonorModel.effect` uses.
    """
    names = [n for n in SOURCE_ORDER if n in p["src"]]
    eff = {}
    for n in names:
        e = np.asarray(p["src"][n], dtype=np.float64)
        if n == "k562_gwps" and pm_blend > 0.0:
            e = e.copy()
            both = np.asarray(p["avail"][n], bool) & np.asarray(p["pm_avail"], bool)
            rows = np.flatnonzero(both)
            sh = p["shared"]
            if rows.size and sh.size:
                e[np.ix_(rows, sh)] = ((1.0 - pm_blend) * p["src"][n][np.ix_(rows, sh)]
                                       + pm_blend * p["pm_effects"][np.ix_(rows, sh)])
        eff[n] = e
    sources = [E062.PreparedSource(
        name=n, effects=eff[n],
        target_available=np.asarray(p["avail"][n], bool).copy(),
        mean_on_axis=np.asarray(p["mean_on_axis"][n], float).copy(),
        raw_median_norm=1.0, mapped_gene_count=0, duplicate_symbols=()) for n in names]
    if weights is None:
        pre, _ = E062.build_baseline(sources)
        return pre
    base = np.zeros_like(sources[0].effects)
    wsum = np.zeros(base.shape[0], dtype=np.float64)
    for s in sources:
        w = float(weights.get(s.name, 0.0))
        rows = s.target_available
        base[rows] += w * s.effects[rows]
        wsum += w * rows.astype(np.float64)
    covered = wsum > 0
    base[covered] /= wsum[covered, None]
    if (~covered).any():
        mean_response = sum(s.mean_on_axis for s in sources) / len(sources)
        mean_norm = float(np.linalg.norm(mean_response)) or 1.0
        base[~covered] = E062.EFFECT_SCALE * mean_response / mean_norm
    return base


def gate_for(b_ctrl: np.ndarray, min_cpm: float) -> np.ndarray:
    expr = np.expm1(b_ctrl)
    cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
    return (cpm > min_cpm) & (b_ctrl > 0)


def frames(fold, min_cpm: float):
    """Challenge frames A/B/C at a given CPM gate (the gate is one of the ablated flags)."""
    out = []
    for ctx in ("A", "B", "C"):
        b = plolo.challenge_bctrl(ctx)
        out.append((f"chal{ctx}", b, gate_for(b, min_cpm)))
    return out


if __name__ == "__main__":
    which = sys.argv[1:] or list(FOLDS)
    for held in which:
        p = parts(held)
        print(f"  {held}: {len(p['targets'])} targets, "
              f"{len(p['scored'])} scored, sources {list(p['src'])}", flush=True)
