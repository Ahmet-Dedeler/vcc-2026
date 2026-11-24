#!/usr/bin/env python3
"""552 (R7 owner [d80a10], 10-16): run one of 551's instrument scripts UNCHANGED, with 062's peak-RSS self-guard lifted to
>= 16 GiB in this process. An admission-only deviation, recorded in the CONFIRMATION-LEDGER's R7 section before use.

Why: 062's enforce_memory_limit (experiments/062-cd4-post-whiten/run.py) raises if this PROCESS's peak RSS (ru_maxrss) has
reached the cap its caller hard-codes (sens415:92 passes 8.0 GiB). R7a's T4 comp551 died there on 10-16 at 10:26Z ("cd4_rest
centering: peak RSS 8.430 GiB reached 8.000 GiB"), with ~25 GiB free+inactive: a guard on the process's history, not
memory pressure. R-NEW died on the same line the same morning, and readiness lifted it in scripts/def_lodo_key.py (690j).
sens415.py and 062/run.py are byte-pinned in R7a's instrument closure (and production code), so neither is edited.

What it does: wraps importlib's SourceFileLoader.exec_module so that, right after the module plolo registers as "e162_062"
(file 062-cd4-post-whiten/run.py; the one T4 died in, which sens415 and lib418 share) finishes executing, that ONE
module's enforce_memory_limit is replaced by readiness's 690j guard: the original function called with max(cap, 16 GiB).
Exactly one attribute of one module (coordinator [2d6f46], 10:33Z). 064/run.py's private copy of the same file
("exp062_composition_base", reached only from 064's own CLI main) is left alone. Nothing else is touched: no module is
imported early, nothing is written, and the script runs as `python <script>` would (sys.argv, sys.path[0], __main__).
The guard decides only whether the process may continue (062/run.py lines 48-52 are a pure check-and-raise, and
max_rss_gib reaches only that function), so every computed value is the pinned code's.
usage: python rss_lift552.py <path/to/551 script.py> [its arguments ...]
"""
from __future__ import annotations

import atexit
import importlib.machinery
import runpy
import sys
from pathlib import Path

CAP_GIB = 16.0
_LIFTED: list[str] = []
_orig_exec = importlib.machinery.SourceFileLoader.exec_module


def _exec_module(self, module):
    _orig_exec(self, module)
    fn = getattr(module, "enforce_memory_limit", None)
    if (getattr(module, "__name__", "") == "e162_062" and str(getattr(module, "__file__", "") or "").endswith(
            "062-cd4-post-whiten/run.py") and callable(fn) and not getattr(fn, "_lifted552", False) and not _LIFTED):
        def guard(limit_gib, stage, _orig=fn):
            return _orig(max(float(limit_gib), CAP_GIB) if float(limit_gib) > 0 else limit_gib, stage)
        guard._lifted552 = True
        module.enforce_memory_limit = guard
        _LIFTED.append(module.__name__)


importlib.machinery.SourceFileLoader.exec_module = _exec_module
atexit.register(lambda: print(f"[rss_lift552] 062's peak-RSS guard lifted to >= {CAP_GIB:g} GiB on {len(_LIFTED)} "
                              f"module(s): {_LIFTED}", file=sys.stderr, flush=True))


def main() -> None:
    script = Path(sys.argv[1]).resolve()
    sys.argv = [str(script)] + sys.argv[2:]
    sys.path[0] = str(script.parent)          # what `python <script>` puts first on sys.path
    runpy.run_path(str(script), run_name="__main__")


if __name__ == "__main__":
    main()
