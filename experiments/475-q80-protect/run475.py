"""475 — 484's q80 finding priced on H1: protected (|E0| > block q80) DOWN genes are still multiplied by the mass-balance kf
(median ~0.60 on real blocks), while protected UP genes keep E0. The fix leaves protected genes out of the shrink mask
(`rc_emit_k_475.py`, a copy of the production emitter edb5b8dc with 484's two-line patch).

modes: orig  the production emitter (known answer: stats == 435 seed 0 bit for bit)
       fix   the patched copy swapped in as the harness's emitter
usage: RUN426_EMITTER=prod uv run python run475.py <mode> <run435 args>   -> work/h1_<mode>/
"""
from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
os.environ["RUN426_EMITTER"] = "prod"
MODE = sys.argv.pop(1)
assert MODE in ("orig", "fix"), MODE
sys.path.insert(0, str(ROOT / "experiments/435-zs-on-coarse"))
import run435 as M  # noqa: E402

WORK = HERE / "work" / f"h1_{MODE}"
WORK.mkdir(parents=True, exist_ok=True)
M.WORK = WORK
M.H.WORK = WORK
if MODE == "fix":
    spec = importlib.util.spec_from_file_location("rc_emit_k_475", HERE / "rc_emit_k_475.py")
    K475 = importlib.util.module_from_spec(spec)
    sys.modules["rc_emit_k_475"] = K475
    spec.loader.exec_module(K475)
    M.K = K475
    print("[475] emitter: rc_emit_k_475 (protected genes leave the kf shrink mask)", flush=True)
else:
    print(f"[475] emitter: production {M.K.__file__}", flush=True)

if __name__ == "__main__":
    M.main()
