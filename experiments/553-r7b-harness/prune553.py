#!/usr/bin/env python3
"""553 — delete regenerable intermediates (h5ad files of finished ablations) to save disk. Refuses any path outside this
directory's work/ and anything that is not an .h5ad; records each deletion (size, sha256 from the file's json) in
work/PRUNED.md.  usage: uv run python experiments/553-r7b-harness/prune553.py work/h1t7s3/h1t7s3.h5ad [...]"""
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = (HERE / "work").resolve()
for arg in sys.argv[1:]:
    p = (HERE / arg).resolve() if not Path(arg).is_absolute() else Path(arg).resolve()
    if WORK not in p.parents or p.suffix != ".h5ad":
        raise SystemExit(f"refused: {p} (only .h5ad files under {WORK})")
    if not p.exists():
        print(f"absent: {p}")
        continue
    size = p.stat().st_size
    p.unlink()
    with open(WORK / "PRUNED.md", "a") as f:
        f.write(f"- `{p.relative_to(WORK)}` ({size / 1e9:.2f} GB, {time.strftime('%H:%MZ', time.gmtime())}, prune553.py): "
                f"regenerable by its build step; sha256 in the json beside it.\n")
    print(f"pruned {p} ({size / 1e9:.2f} GB)")
