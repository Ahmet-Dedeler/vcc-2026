"""801 -- engine801's known answers on a dump that already exists on disk (522's local dump522-<panel>-<fold>.npz), so the
engine's arithmetic is checked on the Mac before any cloud dump exists. Light (~1-2 GB).

    uv run python experiments/801-lab-folds/test_local801.py efp k562_gwps      |  ... abc k562_gwps
"""
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import json  # noqa: E402
import time  # noqa: E402

import engine801 as E  # noqa: E402

panel, held = sys.argv[1], sys.argv[2]
t0 = time.time()
if panel == "abc":
    F = E.AbcFold(held)
else:
    F = E.DumpFold(f"{panel}/{held}", E.ROOT / f"experiments/522-october-bulk-key/work/dump522-{panel}-{held}.npz")
F.ka = E.known_answers(F)
print(json.dumps(F.ka, indent=1, default=float))
print(f"PASS {panel}/{held} ({time.time() - t0:.0f}s, peak {E.L.E062.peak_rss_gib():.2f} GiB)")
(HERE / "work" / "ka").mkdir(parents=True, exist_ok=True)
(HERE / "work" / "ka" / f"local-{panel}-{held}.json").write_text(json.dumps(F.ka, indent=1, default=float))
