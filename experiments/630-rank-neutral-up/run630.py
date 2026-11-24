#!/usr/bin/env python3
"""630 — dump the REAL production emitter's final H1 blocks for 446bj2z's DE construction (537's `cpm` mode = 435's
production run: rc_emit_k 0f56209d, CallSet(400, 2200, 1.9), CLEAN4 q80 block prior, coarse MDE, ballast, zero repair,
+ the zs pass), so the rank-neutral UP pass (rnup630.py) can be developed and scored offline on identical emission draws.

usage: uv run python experiments/630-rank-neutral-up/run630.py <seed>
Runs run496.main exactly as 537's drive537 does for `cpm` (args `base <seed> --bulk ps --mde-c 0.1 --mde-cap 0.1
--mde-min-count 0.3 --mde-coarse --zs dn_k250_zs --dump 85`), with every output redirected into this directory:
work/h1/ (stats, metrics, diag) and work/blocks/ (the 85 pre-zs final blocks; rnup630 re-applies the zs pass itself and
must reproduce the zs stats bit for bit, a known answer). ~7.5 GB, ~5 min; run it through experiments/602-council-gate.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
os.environ["NUMBA_CACHE_DIR"] = str(HERE / "work" / "numba-cache")
os.environ["RUN426_EMITTER"] = "prod"
os.environ.pop("RUN496_EMITTER", None)
sys.dont_write_bytecode = True
seed = sys.argv.pop(1)
sys.path.insert(0, str(ROOT / "experiments/496-uncalled-sign-hook"))
import run496 as M  # noqa: E402
import hashlib  # noqa: E402
from vcc2026 import rc_emit_k as KPROD  # noqa: E402

sha = hashlib.sha256(Path(KPROD.__file__).read_bytes()).hexdigest()
assert sha.startswith("0f56209d57e8385a"), f"emitter is not the committed production 0f56209d: {sha}"
assert M.K is KPROD, "run496 did not import the production emitter"
WORK = HERE / "work" / "h1"
WORK.mkdir(parents=True, exist_ok=True)
M.WORK = WORK
M.H.WORK = WORK
M.HERE = HERE            # run496's block dump writes HERE / work/blocks/
sys.argv = [sys.argv[0], "base", seed, "--bulk", "ps", "--mde-c", "0.1", "--mde-cap", "0.1", "--mde-min-count", "0.3",
            "--mde-coarse", "--zs", "dn_k250_zs", "--dump", "85"]
if __name__ == "__main__":
    M.main()
