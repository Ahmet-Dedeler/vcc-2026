#!/usr/bin/env python3
"""657: rnup630_657.py (d4cd5f's rnup630 + the `_o4` selection suffix) on the h1-trainval surface, rnup618's swaps exactly:
evaluate323's reference dir -> 618's work/ref_tv; output dir -> 657's work/stacks.
    uv run python rnup657.py <seed> <variant> [...] --src <run dir with blocks/ and stats> --name <prefix> --tag <block tag>"""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
import rnup630_657 as U  # noqa: E402
U.R.WORK = ROOT / "experiments/618-h1-trainval-surface/work/ref_tv"
U.WORK = HERE / "work" / "stacks"
U.WORK.mkdir(parents=True, exist_ok=True)
sys.argv = [str(HERE / "rnup630_657.py")] + sys.argv[1:]
U.main()
