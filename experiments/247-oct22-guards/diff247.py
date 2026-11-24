"""247 — how far is the carrying-rule effects array from the shipped one at the matched dose?

Both staged with `stage247.py --save`: `on028` (flag off + H1 g=0.28 = the shipped h1union
geometry, sha abeb0749...) and `carry` (the same plus VCC_CD4_SIZING=carrying at g_cd4).  The
solved fraction differs by one ULP (5.6e-17), so the arrays differ at the last-bit level; this
quantifies it so "numeric, not byte, identity" is a number.

    uv run python experiments/247-oct22-guards/diff247.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

WORK = Path(__file__).resolve().parent / "work" / "stage"


def main() -> None:
    a = np.load(WORK / "effects-on028.npy")
    b = np.load(WORK / "effects-carry.npy")
    d = np.abs(a - b)
    denom = np.abs(a).max()
    out = {
        "shape": list(a.shape),
        "max_abs_diff": float(d.max()),
        "max_abs_diff_over_max_abs_value": float(d.max() / denom),
        "entries_differing": int(np.count_nonzero(d)),
        "entries_total": int(d.size),
        "frac_differing": float(np.count_nonzero(d) / d.size),
        "max_rel_diff_on_nonzero": float(np.max(d[a != 0] / np.abs(a[a != 0]))) if np.any(a != 0) else None,
        "rows_differing": int(np.count_nonzero(d.max(axis=1))),
        "on028_sha256_matches_shipped": None,
    }
    (WORK / "diff247.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
