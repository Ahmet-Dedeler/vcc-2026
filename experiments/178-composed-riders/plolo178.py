"""178 — the composed riders through PLOLO, the one board-validated `pds` instrument.

162 measured PLOLO (panel-restricted leave-one-library-out, scored after the renderer, in the
challenge frame) at **6/6 signs, fitted scale 0.899, leave-one-arm-out RMS 0.001494 scaled
`pds`, permutation p = 0.0028** across six historical geometry arms — and REJECTED it on every
other member. So it is used here for `pds` and for nothing else.

WHAT HAD TO BE ADDED
--------------------
162 caches one combined pre-whitening matrix per fold (`work/fold-<held>.npz:pre`), which is
`E062.build_baseline`'s UNWEIGHTED mean over the fold's sources. A donor-weighting arm is not
expressible on that object, so this script rebuilds each fold keeping the sources separate and
mixes them with the arm's weights. `weighted_baseline` at equal weights must reproduce
`build_baseline` bit-for-bit, and that is asserted per fold.

Note what the fold structure does to the weighting arm: `k562_gwps` is the up-weighted donor,
so the fold that HOLDS OUT K562 has no K562 to up-weight and reads exactly 0.0 by construction.
Three of the four folds carry the arm. That is reported, not averaged away.

KNOWN-ANSWER ROWS (rule 10's corollary)
---------------------------------------
`a15` (equal weights, gamma 0.65, no repair, no mass fix) and `a23` (the same plus
`mass_fix=minnorm`, which is the incumbent every arm here sits on) must reproduce 162's own
`chal{A,B,C}/delivered` endpoints from `work/arms.json`, e.g.

    a23 / orion_hct116 / chalA   pds 0.790896      a23 / kolf2 / chalA   pds 0.615317
    a15 / orion_hct116 / chalA   pds 0.786660      a15 / kolf2 / chalA   pds 0.615896

A mismatch means the fold rebuild is not 162's fold, and everything below is void.

No cell-eval2, no build, no submission.
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402
from plolo import (DONOR_PATHS, FOLDS, OFFICIAL, PANEL, PANEL_SET, PM_BLEND, PM_PATH,
                   PROD_OPS, Operator, cd4_for, deliver, endpoints, solve_fraction)  # noqa: E402

SPAN_PDS = 0.452                      # 162/retrodict.py, docs/reference BRIEF §8
K_CAL = 1.011                         # rule 13's calibration (162's a23-held-out scale)
REACH_TOLL = -0.0146                  # rule 13, scaled, for an emitted-profile change
FOLD_ORDER = ("orion_hct116", "orion_hek293t", "k562_gwps", "kolf2")

# arm -> (weights over PRODUCTION donor names, gamma, sign repair, mass fix)
W111 = {"k562_gwps": 1.0, "orion_hct116": 1.0, "orion_hek293t": 1.0}
W211 = {"k562_gwps": 2.0, "orion_hct116": 1.0, "orion_hek293t": 1.0}
ARMS: dict[str, tuple[dict, float, bool, str | None]] = {
    "a15":        (W111, 0.65, False, None),        # KNOWN ANSWER vs 162's arms.json
    "a23":        (W111, 0.65, False, "minnorm"),   # KNOWN ANSWER, and the incumbent base
    "g045":       (W111, 0.45, False, "minnorm"),
    "g030":       (W111, 0.30, False, "minnorm"),
    "flip":       (W111, 0.65, True,  "minnorm"),
    "w211":       (W211, 0.65, False, "minnorm"),
    "flip_w211":  (W211, 0.65, True,  "minnorm"),
    "g045_flip":  (W111, 0.45, True,  "minnorm"),
    "g045_w211":  (W211, 0.45, False, "minnorm"),
    "all3_g045":  (W211, 0.45, True,  "minnorm"),
    "all3_g030":  (W211, 0.30, True,  "minnorm"),
}


def load_builder():
    spec = importlib.util.spec_from_file_location(
        "prod_builder_178_plolo", ROOT / "scripts/make_submission_pm_cd4.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


M = load_builder()


def weighted_baseline(effects, available, means, weights):
    """`E062.build_baseline` generalised to per-source weights.

    At equal weights this must be the unweighted mean over the sources that carry each target,
    i.e. exactly `build_baseline`; the caller asserts that. The uncovered-target fallback is
    kept verbatim (the unweighted mean response, scaled to `EFFECT_SCALE`), because no arm here
    changes it and a silent divergence there would be a second variable.
    """
    baseline = np.zeros_like(effects[0])
    wsum = np.zeros(baseline.shape[0], dtype=np.float64)
    for eff, avail, w in zip(effects, available, weights):
        rows = np.asarray(avail, dtype=bool)
        if w:
            baseline[rows] += w * eff[rows]
        wsum += w * rows
    covered = wsum > 0
    baseline[covered] /= wsum[covered, None]
    if (~covered).any():
        mean_response = sum(means) / len(means)
        mean_norm = float(np.linalg.norm(mean_response)) or 1.0
        baseline[~covered] = plolo.E062.EFFECT_SCALE * mean_response / mean_norm
    return baseline


def build_fold_sources(held: str):
    """162's `build_fold`, but keeping the sources separate. Cached per fold."""
    cache = WORK / f"fold178-{held}.npz"
    spec = FOLDS[held]
    if cache.exists():
        z = np.load(cache, allow_pickle=True)
        return dict(held=held, targets=tuple(str(t) for t in z["targets"]),
                    sources=tuple(str(s) for s in z["sources"]),
                    eff=[z[f"eff{i}"] for i in range(len(z["sources"]))],
                    avail=[z[f"avail{i}"] for i in range(len(z["sources"]))],
                    mean=[z[f"mean{i}"] for i in range(len(z["sources"]))],
                    pre_ref=z["pre_ref"], truth=z["truth"], b_ctrl=z["b_ctrl"],
                    gate=z["gate"], scored=z["scored"], wcols=z["wcols"],
                    panel_positions=z["panel_positions"], op=spec["op"],
                    op_kind=spec["op_kind"])

    banks = {n: plolo.E062.Bank.load(n, DONOR_PATHS[n]) for n in spec["sources"]}
    held_bank = plolo.E062.Bank.load(held, DONOR_PATHS[held])
    targets = tuple(g for g in PANEL if g in held_bank.pert_pos)
    official_pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([official_pos[g] for g in PANEL if g in official_pos],
                                 dtype=np.int64)

    prepared = {n: plolo.E062.prepare_main_source(b, targets, OFFICIAL, panel_positions,
                                                  chunk_rows=256, max_rss_gib=12.0)
                for n, b in banks.items()}
    pm_bank = plolo.E064.load_perturbmatch(PM_PATH)
    pm_prepared = plolo.E062.prepare_main_source(pm_bank, targets, OFFICIAL, panel_positions,
                                                 chunk_rows=256, max_rss_gib=12.0)
    shared = (plolo.E064.shared_k562_pm_columns(banks["k562_gwps"], pm_bank, OFFICIAL,
                                                PANEL_SET)
              if "k562_gwps" in banks else None)

    source_prepared = []
    for n in spec["sources"]:
        if n == "k562_gwps":
            eff = np.asarray(prepared[n].effects, dtype=np.float64).copy()
            both = prepared[n].target_available & pm_prepared.target_available
            rows = np.flatnonzero(both)
            if len(rows):
                eff[np.ix_(rows, shared)] = (
                    (1 - PM_BLEND) * prepared[n].effects[np.ix_(rows, shared)]
                    + PM_BLEND * pm_prepared.effects[np.ix_(rows, shared)])
            source_prepared.append(plolo.E062.PreparedSource(
                name="k562_gwps_pm50", effects=eff,
                target_available=np.asarray(prepared[n].target_available, bool).copy(),
                mean_on_axis=np.asarray(prepared[n].mean_on_axis, float).copy(),
                raw_median_norm=float(prepared[n].raw_median_norm),
                mapped_gene_count=int(prepared[n].mapped_gene_count),
                duplicate_symbols=tuple(prepared[n].duplicate_symbols)))
        else:
            source_prepared.append(prepared[n])
    pre_ref, _ = plolo.E062.build_baseline(source_prepared)

    held_prep = plolo.E062.prepare_main_source(held_bank, targets, OFFICIAL, panel_positions,
                                               chunk_rows=256, max_rss_gib=12.0)
    truth = plolo.E062.raw_truth(held_prep)
    b_ctrl, held_cols = plolo._b_ctrl_on_axis(held_bank)
    expression = np.expm1(b_ctrl)
    cpm = 1e6 * expression / max(float(expression.sum()), 1e-12)
    gate = (cpm > plolo.MIN_CPM) & (b_ctrl > 0)
    scored_mask = np.zeros(len(OFFICIAL), dtype=bool)
    scored_mask[held_cols] = True
    scored_mask[panel_positions] = False
    scored = np.flatnonzero(scored_mask)

    op_path = Path(spec["op"])
    if spec["op_kind"] == "production":
        wcols = plolo.E064.operator_columns([banks[n] for n in spec["sources"]],
                                            OFFICIAL, PANEL_SET)
    else:
        _, wcols = plolo.E064.load_frozen_operator(
            op_path, [banks[n] for n in spec["sources"]], OFFICIAL, PANEL_SET,
            expected_sha256=plolo.E062.sha256_file(op_path))

    payload = dict(targets=np.asarray(targets), sources=np.asarray(spec["sources"]),
                   pre_ref=pre_ref.astype(np.float64), truth=truth.astype(np.float64),
                   b_ctrl=b_ctrl, gate=gate, scored=scored, wcols=wcols,
                   panel_positions=panel_positions)
    for i, s in enumerate(source_prepared):
        payload[f"eff{i}"] = np.asarray(s.effects, dtype=np.float64)
        payload[f"avail{i}"] = np.asarray(s.target_available, dtype=bool)
        payload[f"mean{i}"] = np.asarray(s.mean_on_axis, dtype=np.float64)
    np.savez_compressed(cache, **payload)
    print(f"[{held}] built {len(targets)}/300 targets, sources {spec['sources']}, "
          f"op {op_path.name} sha {plolo.E062.sha256_file(op_path)[:8]}, "
          f"{len(wcols)} operator columns, {len(scored)} scored genes", flush=True)
    return build_fold_sources(held)


def main() -> None:
    t0 = time.time()
    rows: list[dict] = []
    out_path = WORK / "plolo178.json"
    for held in FOLD_ORDER:
        fold = build_fold_sources(held)
        sources = fold["sources"]
        eff = [np.asarray(e, dtype=np.float64) for e in fold["eff"]]
        avail = [np.asarray(a, dtype=bool) for a in fold["avail"]]
        mean = fold["mean"]
        wcols = np.asarray(fold["wcols"], dtype=np.int64)
        panel_positions = np.asarray(fold["panel_positions"], dtype=np.int64)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]

        # gate: the equal-weight mix must be `build_baseline` bit-for-bit
        eq = weighted_baseline(eff, avail, mean, [1.0] * len(sources))
        ref = np.asarray(fold["pre_ref"], dtype=np.float64)
        drift = float(np.abs(eq.astype(np.float32) - ref.astype(np.float32)).max())
        print(f"[{held}] equal-weight rebuild vs 162's build_baseline: max |d| {drift:.3e}  "
              f"{'PASS' if drift == 0.0 else 'FAIL'}", flush=True)
        assert drift == 0.0

        cd4 = cd4_for(fold["targets"])
        op = Operator(Path(fold["op"]))
        frames = [(f"chal{c}", plolo.challenge_bctrl(c)) for c in "ABC"]
        frames = [(nm, b, (1e6 * np.expm1(b) / max(float(np.expm1(b).sum()), 1e-12)
                           > plolo.MIN_CPM) & (b > 0)) for nm, b in frames]

        pre_cache: dict[tuple, np.ndarray] = {}
        for arm, (weights, gamma, repair, mass) in ARMS.items():
            # `fold["sources"]` holds the ORIGINAL donor names (the k562 slot carries the
            # PerturbMatch blend inside its effects, not in its name), so the weights map
            # directly and a fold without k562 gets no up-weighting at all -- by construction.
            key = tuple(float(weights[s]) for s in sources)
            if key not in pre_cache:
                pre_cache[key] = weighted_baseline(eff, avail, mean, list(key))
            pre = pre_cache[key] * 0.982                      # a15/a23's VCC_ALPHA
            whitened = op.apply(pre, wcols, gamma, 0.0)
            n_flip = 0
            if repair:
                whitened, sr = M.sign_repair_to_consensus(whitened, pre, wcols)
                n_flip = int(sr["flipped_entries"])
            frac, per_unit = solve_fraction(whitened, cd4, 0.10)
            final, _ = plolo.E062.add_post_whiten_residual(whitened, cd4, frac, wcols)
            final[:, panel_positions] = 0.0
            for fname, b, g in frames:
                mat = deliver(final, b, g, mass)
                e = endpoints(mat[:, fold["scored"]], truth)
                rows.append(dict(fold=held, arm=arm, space=f"{fname}/delivered",
                                 weights=list(key), gamma=gamma, repair=repair,
                                 mass=mass, n_flip=n_flip, cd4_fraction=float(frac), **e))
            r = rows[-3:]
            print(f"  {held:<14} {arm:<11} pds(A/B/C) "
                  f"{r[0]['pds']:.6f} {r[1]['pds']:.6f} {r[2]['pds']:.6f}  "
                  f"cos {np.mean([x['cos'] for x in r]):.5f}  flips {n_flip}", flush=True)
            out_path.write_text(json.dumps(rows, indent=1))
        del op, cd4, eff, avail, pre_cache
    print(f"\nwrote {out_path} [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
