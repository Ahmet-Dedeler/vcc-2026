#!/usr/bin/env python3
"""618 stacks: d4cd5f's 630/rnup630.py (UP pass, optionally + DN), run UNCHANGED as a module on the h1-trainval surface.
Swaps, all made before rnup630.main():
  * evaluate323's reference dir R.WORK -> work/ref_tv, so load_reference AND score (polars reads reference_de.parquet,
    which a numpy / io sandbox could not redirect) both read the trainval reference;
  * rnup630.WORK (its hard-coded output dir 630/work/h1) -> work/stacks: nothing is written into 630's dir;
  * block order: handled by rnup630 itself since 7ab0447e (numeric dump order, and each block's stored target asserted
    equal to the reference row; 21dd4c found that lexicographic order breaks at >= 100 targets).
    uv run python rnup618.py <seed> <variant> --src <run dir with blocks/ and stats> --name <prefix> --tag <block tag>
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "experiments/630-rank-neutral-up"))
import numpy as np  # noqa: E402
import rnup630 as U  # noqa: E402  (read-only)

U.R.WORK = HERE / "work" / "ref_tv"
U.WORK = HERE / "work" / "stacks"
U.WORK.mkdir(parents=True, exist_ok=True)
sys.argv = [str(ROOT / "experiments/630-rank-neutral-up/rnup630.py")] + sys.argv[1:]
U.main()
