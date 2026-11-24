"""811: agg811.py (unchanged, sha256 d34d1e26...) restricted to ONE kit tag's folds (e.g. DEFo/<held>), so a run dir that also
holds the KA folds (L1/k562_gwps, ABCp/k562_gwps) aggregates only the tag. No statistic is changed.
    uv run python experiments/811-lolo-model/aggtag811.py <TAG> <run-dir> --family F --tries N --pairs cand:base"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import agg811 as G  # noqa: E402

TAG = sys.argv.pop(1)
G.SURF = {TAG: lambda f: f.split("/")[0] == TAG}
if __name__ == "__main__":
    G.main()
