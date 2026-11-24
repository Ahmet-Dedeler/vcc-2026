#!/usr/bin/env python3
"""657 KA helper: d4cd5f's ORIGINAL rnup630 (630's dir, unmodified) with rnup618's swaps, output to 657's work/stacks_orig.
Used once: the same dumps through the original and through rnup630_657 (unsuffixed variant) must give identical stats."""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "experiments/630-rank-neutral-up"))
import rnup630 as U  # noqa: E402  (read-only)
U.R.WORK = ROOT / "experiments/618-h1-trainval-surface/work/ref_tv"
U.WORK = HERE / "work" / "stacks_orig"
U.WORK.mkdir(parents=True, exist_ok=True)
sys.argv = [str(ROOT / "experiments/630-rank-neutral-up/rnup630.py")] + sys.argv[1:]
U.main()
