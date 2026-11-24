"""233 — pre-build every (panel, fold) cache, serialised and memory-gated.

Independent of the operator fits (the fold's axis is a set operation on donor coverage and
the panel), so this runs while `whiten233.py` is still fitting.
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import panels233 as P
import fold233
from whiten233 import memory_ok

for panel in ("abc", "real", "pess"):
    for held in P.CANCER_FOLDS:
        for _ in range(120):
            ok, msg = memory_ok(min_swap_free_mb=600.0, min_free_gib=10.0)
            if ok:
                break
            print(f"[{panel}/{held}] HOLD ({msg})", flush=True)
            time.sleep(60)
        t0 = time.time()
        fold233.build(panel, held, limit_gib=6.0)
        print(f"[{panel}/{held}] done in {time.time()-t0:.0f}s", flush=True)
