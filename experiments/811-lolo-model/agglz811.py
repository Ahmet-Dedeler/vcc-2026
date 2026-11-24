"""811: agg811.py (unchanged, sha256 d34d1e26...) on the sealed reserve's six folds POOLED (LZ7/<held> + LZ8/<held>), for
LZ-CLAIM-811.md's one access. Same three lines as aggtag811.py, except the surface picks both LZ tags; no statistic is changed.
The run dir holds the six LZ folds only (copied from run-comp811kit-LZ7 / -LZ8; their KA folds stay out).
    uv run python experiments/811-lolo-model/agglz811.py <run-dir> --family 811-LZ --tries 1 --pairs C0:B0,C1:B1
Written 2026-10-22, before any LZ score existed."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import agg811 as G  # noqa: E402

G.SURF = {"LZ": lambda f: f.split("/")[0] in ("LZ7", "LZ8")}
if __name__ == "__main__":
    G.main()
