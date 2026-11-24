"""415 step 1 — the per-source prepared rows that `pre` averages, on the full 300-target panel, cached once.

`plolo.build_fold` prepares every source bank (`E062.prepare_main_source`: centre on the bank's column mean, scale by
EFFECT_SCALE / raw_median) and PM-blends k562 50/50 on the shared columns, then `build_baseline` averages them. The
fold cache keeps only that average, so a leave-one-donor-out prediction needs the individual rows back. Bank-level
statistics (column mean, raw median) do not depend on the target list, so preparing on the whole panel and indexing
by target is exact.

writes work/prep-{k562_gwps_pm50,k562_gwps,orion_hct116,orion_hek293t}.npz  (effects float64, available, raw_median)

    uv run python experiments/415-lodo-key/prep415.py
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo  # noqa: E402

plolo.WORK = WORK
from plolo import E062, E064, OFFICIAL, PANEL, PANEL_SET, PM_BLEND, PM_PATH, DONOR_PATHS  # noqa: E402

LIMIT = 8.0


def main() -> None:
    pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([pos[g] for g in PANEL if g in pos], dtype=np.int64)
    targets = list(PANEL)
    banks = {n: E062.Bank.load(n, DONOR_PATHS[n]) for n in ("k562_gwps", "orion_hct116", "orion_hek293t")}
    prepared = {}
    for n, b in banks.items():
        t0 = time.time()
        prepared[n] = E062.prepare_main_source(b, targets, OFFICIAL, panel_positions, chunk_rows=256, max_rss_gib=LIMIT)
        print(f"{n}: available {int(prepared[n].target_available.sum())}/300 raw_median "
              f"{prepared[n].raw_median_norm:.6f} ({time.time() - t0:.0f}s, peak {E062.peak_rss_gib():.2f} GiB)",
              flush=True)
    pm_bank = E064.load_perturbmatch(PM_PATH)
    pm = E062.prepare_main_source(pm_bank, targets, OFFICIAL, panel_positions, chunk_rows=256, max_rss_gib=LIMIT)
    shared = E064.shared_k562_pm_columns(banks["k562_gwps"], pm_bank, OFFICIAL, PANEL_SET)
    k = prepared["k562_gwps"]
    eff = np.asarray(k.effects, dtype=np.float64).copy()
    both = k.target_available & pm.target_available
    rows = np.flatnonzero(both)
    eff[np.ix_(rows, shared)] = ((1 - PM_BLEND) * k.effects[np.ix_(rows, shared)]
                                 + PM_BLEND * pm.effects[np.ix_(rows, shared)])
    print(f"pm blend: {int(both.sum())} rows x {len(shared)} shared columns", flush=True)
    out = {"k562_gwps_pm50": (eff, k.target_available, k.raw_median_norm, k.mean_on_axis),
           **{n: (p.effects, p.target_available, p.raw_median_norm, p.mean_on_axis) for n, p in prepared.items()}}
    for n, (e, a, rm, mo) in out.items():
        np.savez(WORK / f"prep-{n}.npz", effects=e, available=a, raw_median=rm, mean_on_axis=mo,
                 targets=np.asarray(targets))
    print(f"done, peak RSS {E062.peak_rss_gib():.2f} GiB", flush=True)


if __name__ == "__main__":
    main()
