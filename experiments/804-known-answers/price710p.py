"""804 — row 710 on 801's ABCp folds (A/B/C, day-recipe held-out operators). Same arms and KA as price710.py.

price710.py looks up the design by fold id; ABCp shares A/B/C's panel, composition and targets with the 062 folds, so the
design file is ABC__<held>.npz. The KA (design base composition == the ABCp fold's fin, |d| < 1e-12) decides whether that
holds on the day operators; a failed KA stops the row (no re-cut here).
    KA804_OUT=/work/ka804p uv run python experiments/804-known-answers/price710p.py ABCp/<held>
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import price710 as P  # noqa: E402

_load = np.load


def _design_load(path, *a, **k):
    p = Path(path)
    # ABCp's own design file (801's re-cut with the day operators) when it exists; else 062's A/B/C design, which the KA
    # then rejects (10-16 09:26Z: |d| 1.2-4.0, the compositions live in the whitened space of 062's operators)
    if p.parent == P.DES and p.name.startswith("ABCp__") and not p.exists():
        p = p.with_name("ABC__" + p.name[len("ABCp__"):])
    return _load(p, *a, **k)


if __name__ == "__main__":
    assert sys.argv[1].startswith("ABCp/"), "ABCp folds only"
    P.np.load = _design_load          # price710 imported numpy as np; only its design lookup goes through here
    P.main(sys.argv[1])
