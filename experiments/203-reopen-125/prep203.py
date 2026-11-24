"""203 Phase 1 prep — rebuild the PER-SOURCE prepared effect matrices for each PLOLO fold.

162's cached folds hold only `pre`, the equal-weight MEAN over sources.  125's estimator
needs the sources separately, so this rebuilds them with 062's exact `prepare_main_source`
(including 064's PerturbMatch 50% blend on `k562_gwps`, which is what production ships)
and caches them under 203/work.

KNOWN-ANSWER ROW (rule 10): `build_baseline(rebuilt sources)` must reproduce 162's cached
`fold["pre"]` to float32 round-off.  If it does not, the rebuild is not the production
input and nothing downstream may be believed.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "src"))

import plolo  # noqa: E402
from plolo import DONOR_PATHS, FOLDS, OFFICIAL, PANEL, PANEL_SET, PM_PATH, build_fold  # noqa: E402

E062, E064 = plolo.E062, plolo.E064
OUT = HERE / "work"
OUT.mkdir(exist_ok=True)
LIMIT = 7.0


def rebuild(held: str) -> dict:
    cache = OUT / f"sources-{held}.npz"
    fold = build_fold(held, limit_gib=LIMIT)
    targets = fold["targets"]
    official_pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([official_pos[g] for g in PANEL if g in official_pos],
                                 dtype=np.int64)
    if cache.exists():
        print(f"[{held}] source cache hit", flush=True)
        return dict(cache=cache, fold=fold)

    spec = FOLDS[held]
    banks = {n: E062.Bank.load(n, DONOR_PATHS[n]) for n in spec["sources"]}
    prepared = {n: E062.prepare_main_source(b, targets, OFFICIAL, panel_positions,
                                            chunk_rows=256, max_rss_gib=LIMIT)
                for n, b in banks.items()}
    gene_support = {}
    for n, b in banks.items():
        amap = E062.map_bank_to_axis(b, OFFICIAL)
        m = np.zeros(len(OFFICIAL), dtype=bool)
        m[amap.unique_target_positions] = True
        m[panel_positions] = False
        gene_support[n] = m

    pm_prepared = E062.prepare_main_source(E064.load_perturbmatch(PM_PATH), targets,
                                           OFFICIAL, panel_positions,
                                           chunk_rows=256, max_rss_gib=LIMIT)
    pm_bank = E064.load_perturbmatch(PM_PATH)
    shared = (E064.shared_k562_pm_columns(banks["k562_gwps"], pm_bank, OFFICIAL, PANEL_SET)
              if "k562_gwps" in banks else None)

    names, effs, avail, supp = [], [], [], []
    for n in spec["sources"]:
        if n == "k562_gwps":
            eff = np.asarray(prepared[n].effects, dtype=np.float64).copy()
            both = prepared[n].target_available & pm_prepared.target_available
            rows = np.flatnonzero(both)
            if len(rows):
                eff[np.ix_(rows, shared)] = (
                    0.5 * prepared[n].effects[np.ix_(rows, shared)]
                    + 0.5 * pm_prepared.effects[np.ix_(rows, shared)])
            names.append("k562_gwps_pm50")
        else:
            eff = np.asarray(prepared[n].effects, dtype=np.float64)
            names.append(n)
        effs.append(eff.astype(np.float32))
        avail.append(np.asarray(prepared[n].target_available, dtype=bool))
        supp.append(gene_support[n])

    # ---- KNOWN-ANSWER: equal-mean of these sources must be 162's cached `pre`
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
    print(f"[{held}] KA build_baseline vs 162 cached pre: max abs {err:.3e} "
          f"rel {rel:.3e}  (float32 cache eps ~ {np.finfo(np.float32).eps:.1e})", flush=True)
    assert rel < 1e-6, f"{held}: rebuilt sources are NOT the production input"

    np.savez_compressed(
        cache, names=np.asarray(names), effects=np.stack(effs),
        available=np.stack(avail), support=np.stack(supp),
        votes=votes.astype(np.int32), targets=np.asarray(targets),
        ka_max_abs=err, ka_rel=rel)
    print(f"[{held}] wrote {cache.name}  effects {np.stack(effs).shape}", flush=True)
    return dict(cache=cache, fold=fold, ka_rel=rel)


if __name__ == "__main__":
    folds = sys.argv[1:] or ["k562_gwps", "orion_hct116", "orion_hek293t"]
    audit = {}
    for h in folds:
        r = rebuild(h)
        audit[h] = {k: v for k, v in r.items() if k == "ka_rel"}
    (OUT / "prep203.json").write_text(json.dumps(audit, indent=1, default=str))
    print("done")
