"""551 -- compute ONE donor bank's panel-independent statistics (lib551.bank_stats: column mean over all mapped columns and
the median centred row norm) in its own process, then check the cache against 534's (lib534 == lib438 arithmetic):
the arrays must be equal and the npz bytes identical (np.savez is deterministic here). One bank per process keeps the
peak under 062's 9 GiB guard, which a single process computing three banks in a row exceeds (2026-10-01 00:09Z: the
memmapped pages of the first two banks plus the operator pushed orion_hek293t's centering to 9.013 GiB).

    uv run python experiments/551-r6-c10b/bankstats551.py <bank>
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard551  # noqa: E402,F401
import hashlib  # noqa: E402

import numpy as np  # noqa: E402

import lib551 as L  # noqa: E402

name = sys.argv[1]
mine = L.WORK / "bankstats" / f"{name}.npz"
ref = L.ROOT / "experiments/534-generalization-audit/work/partB/bankstats" / f"{name}.npz"
L.bank_stats(name)
a, b = np.load(mine), np.load(ref)
for k in ("all_source_columns", "mean_all", "raw_median"):
    assert np.array_equal(a[k], b[k]), f"{name} {k} differs from 534's cache"
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(mine) == sha(ref), f"{name}: npz bytes differ from 534's (arrays equal)"
print(f"{name}: bank stats == 534's cache (arrays and bytes, sha {sha(mine)[:16]}); peak RSS {L.E062.peak_rss_gib():.2f} GiB")
