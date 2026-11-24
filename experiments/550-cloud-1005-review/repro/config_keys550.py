#!/usr/bin/env python3
"""Compare the config-key set of two submission metas. The values are printed; the exit code is the set.

  uv run python experiments/550-cloud-1005-review/repro/config_keys550.py built.meta.json ref.meta.json k1,k2
Exit 0 when the keys that differ are exactly that set. This is the check run_arm_generic.sh does
against one reference. It does not score pds or mse.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path


def config(path: Path) -> dict:
    return json.loads(path.read_text())["config"]


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 3:
        print("usage: config_keys550.py built.meta.json ref.meta.json k1,k2", file=sys.stderr)
        return 2
    built, ref = config(Path(args[0])), config(Path(args[1]))
    diff = {k: (ref.get(k, "<absent>"), built.get(k, "<absent>"))
            for k in sorted(set(ref) | set(built))
            if ref.get(k, "<absent>") != built.get(k, "<absent>")}
    want = {x for x in args[2].split(",") if x}
    for k, (a, b) in diff.items():
        print(f"    {k}")
    ok = set(diff) == want
    print(f"{'PASS' if ok else 'FAIL'}: got {sorted(diff)} want {sorted(want)}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
