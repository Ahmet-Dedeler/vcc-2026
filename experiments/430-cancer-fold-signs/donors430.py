"""430 step 0 — cache every donor table's PANEL rows on the official axis, raw and with the donor's all-perturbation
column mean (the centring production applies in `prepare_main_source`), plus the donor's control profile.

Per donor d in DON (+ kolf2 as the non-cancer truth line), writes work/donor430-<d>.npz:
  raw      (300, 18533) float32  delta row of each panel target (NaN where the gene is not on the donor's axis or the
                                  target is not covered); pert-minus-control, log1p CP50k (kolf2: z-score * 0.04)
  avail    (300,) bool           the donor covers the target
  colmean  (18533,) float64      mean over ALL the donor's perturbation rows (NaN off-axis)
  on_axis  (18533,) bool         the donor measures the gene
  b_ctrl   (18533,) float64      the donor's control log1p CP50k, mapped (NaN off-axis)
Mapping = E062.map_bank_to_axis (last source column wins for a duplicated symbol), identical to plolo/production.
Reads only; writes only experiments/430-cancer-fold-signs/work.

    uv run python experiments/430-cancer-fold-signs/donors430.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo  # noqa: E402

plolo.WORK = ROOT / "experiments/418-bulk-composition/work/plolo"     # read-only; nothing is built through plolo here
from plolo import E062, OFFICIAL, PANEL  # noqa: E402

DON = ("k562_gwps", "k562_essential", "rpe1", "nadig_jurkat", "orion_hct116", "orion_hek293t")
PATHS = {d: ROOT / "data/donors" / d for d in DON} | {"kolf2": ROOT / "data/external/kolf2"}


def build(name: str) -> dict:
    t0 = time.time()
    bank = E062.Bank.load(name, PATHS[name])
    amap = E062.map_bank_to_axis(bank, OFFICIAL)
    src, dst = amap.unique_source_columns, amap.unique_target_positions
    delta = np.load(bank.delta_path, mmap_mode="r")
    n = delta.shape[0]
    acc = np.zeros(src.size, dtype=np.float64)
    for i in range(0, n, 1000):
        acc += np.asarray(delta[i:i + 1000][:, src], dtype=np.float64).sum(0)
    colmean = np.full(len(OFFICIAL), np.nan)
    colmean[dst] = acc / n
    raw = np.full((len(PANEL), len(OFFICIAL)), np.nan, dtype=np.float32)
    avail = np.zeros(len(PANEL), bool)
    for i, t in enumerate(PANEL):
        r = bank.pert_pos.get(t)
        if r is None:
            continue
        raw[i, dst] = np.asarray(delta[r, src], dtype=np.float32)
        avail[i] = True
    on = np.zeros(len(OFFICIAL), bool)
    on[dst] = True
    b = np.full(len(OFFICIAL), np.nan)
    b[dst] = np.asarray(np.load(bank.b_ctrl_path, mmap_mode="r")[src], dtype=np.float64)
    np.savez_compressed(WORK / f"donor430-{name}.npz", raw=raw, avail=avail, colmean=colmean, on_axis=on, b_ctrl=b,
                        n_perts=n)
    info = dict(name=name, n_perts=int(n), n_on_axis=int(on.sum()), panel_covered=int(avail.sum()),
                duplicates=len(amap.duplicate_symbols), secs=round(time.time() - t0, 1))
    print(json.dumps(info), flush=True)
    return info


if __name__ == "__main__":
    out = [build(d) for d in (*DON, "kolf2")]
    (WORK / "donors430.json").write_text(json.dumps(out, indent=1))
