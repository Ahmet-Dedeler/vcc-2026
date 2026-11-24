"""811: agg811.py (unchanged, sha256 d34d1e26...) with one more surface, "Aonly" = 801's A<k>/<held> folds. agg811's own
surfaces match fold ids starting "L" / "ABC", so A folds would be silently skipped. No statistic is changed.
    uv run python experiments/811-lolo-model/agga811.py <run-dir> --family F --tries N --pairs cand:base"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import agg811 as G  # noqa: E402

G.SURF = {"Aonly": lambda f: f[:1] == "A" and f[1:2].isdigit()}
if __name__ == "__main__":
    G.main()
