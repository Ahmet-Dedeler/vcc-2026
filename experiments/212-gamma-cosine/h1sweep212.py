"""212 — 170's six-member h1-test sweep, unchanged, at the five gammas 170 did not run.

Everything is `experiments/170-gamma-mse/h1sweep.py`, imported not copied; the only changes
are that `WORK` points here and the grid is 0.30 / 0.40 / 0.525 / 0.575 / 0.625. The effect
files come from `gate212.py`, which reproduced 170's nine gammas BIT-IDENTICALLY (20/20), so
the five new rows are on the same surface as 170's nine and the two tables can be read as one.

The `mse` known-answer gate is 170's own: `msekit` on the frozen 095 h1-test artifact must
return **2.0551018299736388**, the value `cell-eval2`'s own `agg_results.csv` holds.

`nmae` here is `debench.Surface.nmae`, whose gate is **the reference's own BH<0.05
significant set** with finite reference lfc and at least 10 genes per row -- i.e. the
significance-selected pool 188 identified as the challenge's construction, not a
top-N-by-magnitude pool (which 188 measured INVERTS the sign). h1-test's median `n_real`
is ~90, against the mandate's stated ~89.

Geometry caveat, recorded because it is a real difference from the incumbent: `h1kit` forces
`VCC_MASS_FIX=minnorm` (a23's geometry) and RAISES otherwise, while `rc_callset` ships
without the mass fix. The PLOLO instrument carries both renders and they agree on every
gamma ordering, which is the control for it.

    uv run python h1sweep212.py 0.30,0.40,0.525,0.575,0.625
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "work"
E170 = REPO / "experiments/170-gamma-mse"
sys.path.insert(0, str(E170))

import h1kit  # noqa: E402

# clone 170's DE-reference cache so nothing can be written outside this experiment
for name in ("de_real_h1.npz", "pre.npz"):
    if not (WORK / name).exists():
        shutil.copyfile(E170 / "work" / name, WORK / name)
h1kit.WORK = WORK

import h1sweep  # noqa: E402

h1sweep.WORK = WORK

if __name__ == "__main__":
    grid = ([float(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1
            else [0.30, 0.40, 0.525, 0.575, 0.625])
    # `gate212.py` writes `h1_g0.525.npz` (three decimals); 170's loader asks for two.
    # Bridge with a clone rather than a code change, so `h1sweep.main` stays byte-unchanged.
    (WORK / "eff").mkdir(parents=True, exist_ok=True)
    for g in grid:
        src = WORK / "eff" / f"h1_g{g:.3f}.npz"
        dst = WORK / "eff" / f"h1_g{g:.2f}.npz"
        if not dst.exists():
            shutil.copyfile(src, dst)
    h1sweep.main(grid, 0, "h1sweep212.json")
