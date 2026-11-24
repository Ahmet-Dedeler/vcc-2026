#!/usr/bin/env python3
"""618 stacks: d4cd5f's 630/price630.py run UNCHANGED as a module on the h1-trainval surface.
Swaps, before price630.main(): lib521's REF / EFF / ST435 (323's reference and key, 435's base stats) -> work/ref_tv,
work/tv/effects.npz and work/g2mde's stats (435's exact flags: --bulk ps, coarse MDE, no zs, no zero repair);
price630.W (its hard-coded arm-stats dir 630/work/h1) -> work/stacks (rnup618's outputs).
cmp keeps price630's H1-calibrated constants (frac_lt10 0.22, tail_f 0.425): f 1 / 0.5 / 0.35 are the surface's own.
    uv run python price618.py <variant> --base <zs-stats dir> --name <prefix> --tag <stats tag>
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
sys.path.insert(0, str(ROOT / "experiments/630-rank-neutral-up"))
import lib521 as L  # noqa: E402  (read-only; module constants patched before load_state's first call)

L.REF = HERE / "work/ref_tv/reference.npz"
L.EFF = HERE / "work/tv/effects.npz"
L.ST435 = HERE / "work/g2mde/stats-base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse.npz"
L.load_state.cache_clear()
import price630 as P  # noqa: E402  (read-only)

P.W = HERE / "work" / "stacks"
sys.argv = [str(ROOT / "experiments/630-rank-neutral-up/price630.py")] + sys.argv[1:]
P.main()
