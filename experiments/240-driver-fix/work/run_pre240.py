"""Run a pre-240 byte copy of a `scripts/` file under the path it had in `scripts/`.

Both scripts derive `ROOT` from `__file__`, so a copy kept in this experiment's `work/`
resolves the repo root two levels too deep. Nothing else is changed: the source executed is the
byte copy, compiled with `__file__` set to `<root>/scripts/<name>` -- 235's byte-copy trick.

    uv run python work/run_pre240.py work/def_postbuild.PRE240.py --meta ... --bundle ...
"""
import sys
from pathlib import Path

src = Path(sys.argv[1]).resolve()
root = src.parents[2].parent          # <root>/experiments/240-driver-fix/work/x.py -> <root>
asif = root / "scripts" / src.name
sys.argv = [str(asif), *sys.argv[2:]]
exec(compile(src.read_text(), str(asif), "exec"),
     {"__name__": "__main__", "__file__": str(asif)})
