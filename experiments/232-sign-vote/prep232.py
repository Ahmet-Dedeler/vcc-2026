"""232 Phase 0a — rebuild the PER-SOURCE prepared effect matrices for every PLOLO fold.

162's cached folds hold only `pre`, the equal-weight MEAN over sources -- which is exactly
the object this experiment is about.  To ask "what would the donors' MAJORITY sign say?" the
sources have to be separated again, together with their COVERAGE, because on the official
18,533-gene axis a source that does not carry gene j writes a literal 0.0 into its effect
row.  That zero is "no vote", not "no response", and every statistic below has to know the
difference (falsifier 5, coverage honesty).

This is 203's `prep203.rebuild` re-derived inside 232's own directory so that its
known-answer row is this experiment's, not a borrowed one:

  KA-0  `E062.build_baseline(rebuilt sources)` must reproduce 162's cached `fold["pre"]`
        to float32 round-off on every fold.  If it does not, the separated sources are not
        the production input and nothing downstream may be believed.

Optionally (--extra kolf2) a NON-production source is prepared on the same target axis so
that the cancer folds, which only have TWO covering donors after the hold-out, can also be
read at production's real arity of three.  `kolf2` is never used on the `kolf2` fold.

    uv run python prep232.py
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "src"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
plolo.WORK.mkdir(parents=True, exist_ok=True)

from plolo import (DONOR_PATHS, FOLDS, OFFICIAL, PANEL, PANEL_SET,  # noqa: E402
                   PM_PATH, build_fold)

E062, E064 = plolo.E062, plolo.E064
LIMIT = 7.5          # GiB guard on PEAK RSS (062's `enforce_memory_limit` is cumulative);
                     # measured peak for the four folds in one process: see work/prep232.log


def gene_support(bank, panel_positions) -> np.ndarray:
    amap = E062.map_bank_to_axis(bank, OFFICIAL)
    m = np.zeros(len(OFFICIAL), dtype=bool)
    m[amap.unique_target_positions] = True
    m[panel_positions] = False
    return m


def rebuild(held: str, extra: tuple[str, ...] = ()) -> dict:
    cache = WORK / f"sources-{held}.npz"
    t0 = time.time()
    fold = build_fold(held, limit_gib=LIMIT)
    targets = fold["targets"]
    official_pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([official_pos[g] for g in PANEL if g in official_pos],
                                 dtype=np.int64)
    if cache.exists():
        z = np.load(cache, allow_pickle=True)
        if tuple(str(t) for t in z["targets"]) == tuple(targets):
            print(f"[{held}] source cache hit", flush=True)
            return dict(ka_rel=float(z["ka_rel"]), cached=True)

    spec = FOLDS[held]
    banks = {n: E062.Bank.load(n, DONOR_PATHS[n]) for n in spec["sources"]}
    prepared = {n: E062.prepare_main_source(b, targets, OFFICIAL, panel_positions,
                                            chunk_rows=256, max_rss_gib=LIMIT)
                for n, b in banks.items()}
    supp = {n: gene_support(b, panel_positions) for n, b in banks.items()}

    pm_bank = E064.load_perturbmatch(PM_PATH)
    pm_prepared = E062.prepare_main_source(pm_bank, targets, OFFICIAL, panel_positions,
                                           chunk_rows=256, max_rss_gib=LIMIT)
    shared = (E064.shared_k562_pm_columns(banks["k562_gwps"], pm_bank, OFFICIAL, PANEL_SET)
              if "k562_gwps" in banks else None)

    names, effs, avail, supps = [], [], [], []
    for n in spec["sources"]:
        if n == "k562_gwps":
            eff = np.asarray(prepared[n].effects, dtype=np.float64).copy()
            both = prepared[n].target_available & pm_prepared.target_available
            rows = np.flatnonzero(both)
            if len(rows):
                eff[np.ix_(rows, shared)] = (
                    (1 - plolo.PM_BLEND) * prepared[n].effects[np.ix_(rows, shared)]
                    + plolo.PM_BLEND * pm_prepared.effects[np.ix_(rows, shared)])
            names.append("k562_gwps_pm50")
        else:
            eff = np.asarray(prepared[n].effects, dtype=np.float64)
            names.append(n)
        effs.append(eff.astype(np.float32))
        avail.append(np.asarray(prepared[n].target_available, dtype=bool))
        supps.append(supp[n])

    # ---- KA-0: the equal-weight mean of these sources IS 162's cached `pre`
    src_objs = []
    for i, n in enumerate(spec["sources"]):
        src_objs.append(E062.PreparedSource(
            name=names[i], effects=np.asarray(effs[i], dtype=np.float64),
            target_available=avail[i].copy(),
            mean_on_axis=np.asarray(prepared[n].mean_on_axis, float).copy(),
            raw_median_norm=float(prepared[n].raw_median_norm),
            mapped_gene_count=int(prepared[n].mapped_gene_count),
            duplicate_symbols=tuple(prepared[n].duplicate_symbols)))
    pre_rebuilt, votes = E062.build_baseline(src_objs)
    pre_cached = np.asarray(fold["pre"], dtype=np.float64)
    err = float(np.abs(pre_rebuilt - pre_cached).max())
    rel = err / float(np.abs(pre_cached).max())
    print(f"[{held}] KA-0 build_baseline vs 162 cached pre: max abs {err:.3e} "
          f"rel {rel:.3e}", flush=True)
    assert rel < 1e-6, f"{held}: rebuilt sources are NOT the production input"

    # ---- optional non-production extra voters (never on their own fold)
    xnames, xeffs, xavail, xsupp = [], [], [], []
    for n in extra:
        if n == held:
            continue
        b = E062.Bank.load(n, DONOR_PATHS[n])
        p = E062.prepare_main_source(b, targets, OFFICIAL, panel_positions,
                                     chunk_rows=256, max_rss_gib=LIMIT)
        xnames.append(n)
        xeffs.append(np.asarray(p.effects, dtype=np.float32))
        xavail.append(np.asarray(p.target_available, dtype=bool))
        xsupp.append(gene_support(b, panel_positions))

    payload = dict(names=np.asarray(names), effects=np.stack(effs),
                   available=np.stack(avail), support=np.stack(supps),
                   votes=votes.astype(np.int32), targets=np.asarray(targets),
                   ka_max_abs=err, ka_rel=rel)
    if xnames:
        payload.update(extra_names=np.asarray(xnames), extra_effects=np.stack(xeffs),
                       extra_available=np.stack(xavail), extra_support=np.stack(xsupp))
    np.savez_compressed(cache, **payload)
    print(f"[{held}] wrote {cache.name} effects {np.stack(effs).shape} "
          f"extra {xnames} ({time.time() - t0:.0f}s)", flush=True)
    return dict(ka_rel=rel, ka_max_abs=err, names=names, extra=xnames,
                n_targets=len(targets), seconds=time.time() - t0)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    extra = ("kolf2",) if "--extra-kolf2" in sys.argv else ()
    folds = args or list(FOLDS)
    audit = {}
    for h in folds:
        audit[h] = rebuild(h, extra=extra)
    (WORK / "prep232.json").write_text(json.dumps(audit, indent=1, default=str))
    print("done")
