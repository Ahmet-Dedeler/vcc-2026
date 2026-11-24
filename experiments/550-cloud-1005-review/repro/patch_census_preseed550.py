#!/usr/bin/env python3
"""NOT APPLIED. Shows the census fresh-dir check dying on recipe_dr524's pre-seeded whitener,
and the replacement that allows only that pair.

Refuses to edit a path under the live repo. --self-test runs the old and new checks on a temp dir.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

OLD = '''if [ -d "$W" ] && [ -n "$(ls -A "$W" 2>/dev/null)" ]; then
    die "DEF_DAY_CENSUS_MAX_ZERO=$DEF_DAY_CENSUS_MAX_ZERO is a one-shot override and needs a fresh work dir ($W is not empty)"
  fi'''

NEW = '''# 550: recipe_dr524.sh (and the October day) copy the whitener pair in before def_day.sh.
# Those two names are the pre-seed, not a reused build. Anything else is still a reused dir.
extra="$(find "$W" -mindepth 1 -maxdepth 1 ! -name whitener_union_pruned_def.npz ! -name whitener_union_pruned_def.json -print -quit 2>/dev/null)"
if [ -n "$extra" ]; then
    die "DEF_DAY_CENSUS_MAX_ZERO=$DEF_DAY_CENSUS_MAX_ZERO is a one-shot override and needs a fresh work dir ($W contains $extra)"
  fi'''

WRAP = '''
set -uo pipefail
W="$1"
DEF_DAY_CENSUS_MAX_ZERO=8
die(){ echo "STOP: $*"; exit 1; }
%s
echo OK
'''


def check(snippet: str, work: Path) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", "-c", WRAP % snippet, "census", str(work)], capture_output=True, text=True)


def self_test() -> int:
    root = Path(tempfile.mkdtemp(prefix="550census-"))
    pre = root / "pre"
    pre.mkdir()
    (pre / "whitener_union_pruned_def.npz").write_bytes(b"x")
    (pre / "whitener_union_pruned_def.json").write_text("{}\n")
    dirty = root / "dirty"
    dirty.mkdir()
    (dirty / "whitener_union_pruned_def.npz").write_bytes(b"x")
    (dirty / "def_day.log").write_text("leftover\n")
    old_pre = check(OLD, pre)
    new_pre = check(NEW, pre)
    old_dirty = check(OLD, dirty)
    new_dirty = check(NEW, dirty)
    print("old on whitener pair", old_pre.returncode, (old_pre.stdout or old_pre.stderr).strip())
    print("new on whitener pair", new_pre.returncode, (new_pre.stdout or new_pre.stderr).strip())
    print("old on leftover log", old_dirty.returncode, (old_dirty.stdout or old_dirty.stderr).strip())
    print("new on leftover log", new_dirty.returncode, (new_dirty.stdout or new_dirty.stderr).strip())
    ok = old_pre.returncode != 0 and new_pre.returncode == 0 and old_dirty.returncode != 0 and new_dirty.returncode != 0
    print("SELF-TEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--apply-to", type=Path, help="a COPY of patched def_day.sh or patch_extras538.py. Refused for the live repo.")
    a = ap.parse_args(argv)
    if a.self_test:
        return self_test()
    if a.apply_to is None:
        print("not applied. Pass --self-test, or --apply-to a copy outside the repo.")
        return 0
    live = Path(__file__).resolve().parents[3]
    target = a.apply_to.resolve()
    if live in target.parents or target == live:
        print(f"refusing to edit {target}: it is inside the live repo")
        return 2
    text = target.read_text()
    if OLD not in text:
        print("the fresh-dir check this patch replaces is not in that file")
        return 1
    target.write_text(text.replace(OLD, NEW, 1))
    print(f"patched copy {target}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
