"""204 — the scale identity: de-attenuating a donor TABLE is an exact no-op in production.

`prepare_main_source` (062/run.py:239-297) computes, for donor i,

    effects = (delta_rows - mean_over_all_rows) * (EFFECT_SCALE / raw_median_norm)

and `raw_median_norm` is the median CENTRED row norm of the same matrix.  Every term is
degree-1 homogeneous in `delta`, so scaling a donor's whole table by any constant `c`
leaves `effects` **bit-identical**.  Production therefore already divides each donor by
its own amplitude before averaging: the median row norm IS its implicit attenuation
estimate.

Consequence for this experiment, and it is the reason the arms are what they are: the
mandate's "de-attenuate then average" cannot be applied at the table level -- it would
change nothing.  The only place a donor's attenuation can act is as a RELATIVE weight
between donors inside `build_baseline`, which is what every arm in `score204.py` moves.

Measured here rather than argued: `h1_2025` (255 x 18,080, the smallest bank on disk) is
copied with `delta` multiplied by 2.5, and both are pushed through the production
preparation.
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(HERE))

from prep204 import shared_targets  # noqa: E402
from plolo import E062, OFFICIAL, PANEL  # noqa: E402

C = 2.5
SRC = ROOT / "data/donors/h1_2025"
DST = WORK / "scaled-h1_2025"


def main() -> None:
    if DST.exists():
        shutil.rmtree(DST)
    DST.mkdir(parents=True)
    for f in ("genes.txt", "perts.txt", "b_ctrl.npy", "meta.json"):
        if (SRC / f).exists():
            shutil.copy(SRC / f, DST / f)
    np.save(DST / "delta.npy", (np.load(SRC / "delta.npy") * C))

    official_pos = {g: i for i, g in enumerate(OFFICIAL)}
    pp = np.asarray([official_pos[g] for g in PANEL if g in official_pos], dtype=np.int64)
    names = [p.strip() for p in (SRC / "perts.txt").read_text().splitlines() if p.strip()]
    tg = tuple(t for t in shared_targets() if t in set(names))[:120]

    a = E062.prepare_main_source(E062.Bank.load("h1_2025", SRC), tg, OFFICIAL, pp,
                                 chunk_rows=256, max_rss_gib=8.0)
    b = E062.prepare_main_source(E062.Bank.load("h1_scaled", DST), tg, OFFICIAL, pp,
                                 chunk_rows=256, max_rss_gib=8.0)
    A, B = np.asarray(a.effects), np.asarray(b.effects)
    out = dict(scale=C, n_targets=len(tg),
               raw_median_norm_original=float(a.raw_median_norm),
               raw_median_norm_scaled=float(b.raw_median_norm),
               raw_median_ratio=float(b.raw_median_norm / a.raw_median_norm),
               max_abs_diff_effects=float(np.max(np.abs(A - B))),
               max_rel_diff_effects=float(np.max(np.abs(A - B)) /
                                          max(float(np.max(np.abs(A))), 1e-300)),
               bit_identical=bool(np.array_equal(A, B)))
    (WORK / "homog204.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    shutil.rmtree(DST)


if __name__ == "__main__":
    main()
