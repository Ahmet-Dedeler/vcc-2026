#!/usr/bin/env python3
"""Apply the October pass-through in order: 530's bundle, 528's call-propensity block, 538's extras.

528 does not apply cleanly after 530 (its MDE-grep anchor is gone). That miss is reported and is not
a failure: 538 re-inserts the same snippets. A miss in 530 or in 538 is a failure and writes nothing.

  uv run python experiments/538-cloud-bundle-extras/apply_october538.py --root <repo> --check
  uv run python experiments/538-cloud-bundle-extras/apply_october538.py --root <repo>

--check prints the before/after sha256 of every file the sequence would change and does not replace
the destination. A second apply (and a --check of an already-patched tree) exits `already patched`
before writing. Destination writes are a copy next to the target, then os.replace.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT_DEFAULT = HERE.parents[1]
PATCH530 = ROOT_DEFAULT / "experiments/530-cloud-october-bundle/patch_bundle530.py"
PATCH528 = ROOT_DEFAULT / "experiments/528-callprop-block/patch_defday528.py"
PATCH538 = HERE / "patch_extras538.py"

# Every file 530, 528, or 538 replaces. 528 only edits def_day.sh, which is already in this list.
FILES = [
    Path("scripts/def_day.sh"),
    Path("scripts/def_day_helpers.py"),
    Path("experiments/448-def-rehearsal-q175/relevel_iter.py"),
    Path("scripts/def_sameline.py"),
    Path("scripts/october_recipe.py"),
    Path("docs/october-decisions.json"),
    Path("experiments/524-october-recipe-dr/gate524.py"),
    Path("scripts/def_postbuild.py"),
    Path("scripts/def_pipeline.py"),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _run(cmd: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True)


def already_patched(root: Path) -> bool:
    day = (root / "scripts/def_day.sh").read_text()
    post = (root / "scripts/def_postbuild.py").read_text()
    pipe = (root / "scripts/def_pipeline.py").read_text()
    helpers = (root / "scripts/def_day_helpers.py").read_text()
    recipe = (root / "scripts/october_recipe.py").read_text()
    decisions = (root / "docs/october-decisions.json").read_text()
    gate = (root / "experiments/524-october-recipe-dr/gate524.py").read_text()
    # DEF_DAY_LODO_HI is a prefix of DEF_DAY_LODO_HI_PCT, and "lodo_hi" of lodo_hi_pct.
    # Match the tokens this patch adds, not those prefixes.
    return (
        '--hi "$DEF_DAY_LODO_HI"' in day
        and 'REPO="$(CDPATH=' in day
        and "DEF_DAY_CENSUS_MAX_ZERO" in day
        and "DEF_DAY_CALLSET_BLOCK" in day
        and "--require-config-diff" in post
        and "--census-max-zero" in pipe
        and '{"lodo_hi", "callset"}' in helpers
        and '"callset_block": "DEF_DAY_CALLSET_BLOCK"' in recipe
        and '"lodo_hi": {' in decisions
        and '"DEF_DAY_LODO_HI":' in gate
    )


def partial_patched(root: Path) -> str:
    """A tree that has some 538 markers and not others. 530-only is not partial: 538 has not started."""
    flags = {
        "def_day LODO_HI": '--hi "$DEF_DAY_LODO_HI"' in (root / "scripts/def_day.sh").read_text(),
        "postbuild require": "--require-config-diff" in (root / "scripts/def_postbuild.py").read_text(),
        "pipeline census": "--census-max-zero" in (root / "scripts/def_pipeline.py").read_text(),
        "recipe lodo_hi": '"lodo_hi": "DEF_DAY_LODO_HI"' in (root / "scripts/october_recipe.py").read_text(),
        "gate LODO_HI": '"DEF_DAY_LODO_HI":' in (root / "experiments/524-october-recipe-dr/gate524.py").read_text(),
    }
    if any(flags.values()) and not all(flags.values()):
        have = ", ".join(k for k, v in flags.items() if v)
        missing = ", ".join(k for k, v in flags.items() if not v)
        return f"partial extras (have {have}; missing {missing})"
    return ""


def stage_tree(root: Path, stage: Path) -> None:
    for rel in FILES:
        src = root / rel
        if not src.is_file():
            raise SystemExit(f"missing {src}")
        dest = stage / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(src, dest)


def apply_sequence(stage: Path) -> None:
    r = _run([sys.executable, str(PATCH530), "--root", str(stage)])
    text = (r.stdout or "") + (r.stderr or "")
    if r.returncode != 0:
        if "already patched" in text:
            print("530 already applied on this tree; continuing")
        else:
            raise SystemExit(f"530 failed:\n{text}")
    else:
        print(r.stdout.strip())
    day = stage / "scripts/def_day.sh"
    before_528 = day.read_bytes()
    r = _run([sys.executable, str(PATCH528), str(day)])
    text = (r.stdout or "") + (r.stderr or "")
    if r.returncode != 0:
        if day.read_bytes() != before_528:
            raise SystemExit(f"528 failed and rewrote def_day.sh:\n{text}")
        if "anchor found" not in text:
            raise SystemExit(f"528 failed:\n{text}")
        print("528 did not apply (anchor miss after 530); extras re-anchors the call-propensity block")
    else:
        print(r.stdout.strip() or "528 applied cleanly")
    r = _run([sys.executable, str(PATCH538), "--root", str(stage)])
    text = (r.stdout or "") + (r.stderr or "")
    if r.returncode != 0:
        raise SystemExit(f"538 failed:\n{text}")
    print(r.stdout.strip())


def print_shas(before: dict[str, str], after: dict[str, str]) -> None:
    for rel in FILES:
        flag = "CHANGE" if before[str(rel)] != after[str(rel)] else "same"
        print(f"{rel}")
        print(f"  before {before[str(rel)]}")
        print(f"  after  {after[str(rel)]}  {flag}")


def install(stage: Path, root: Path) -> None:
    """Copy each staged file beside its destination, then os.replace. All temps exist before any replace."""
    temps = []
    for rel in FILES:
        dest = root / rel
        tmp = dest.with_name(dest.name + ".538apply")
        shutil.copyfile(stage / rel, tmp)
        temps.append((tmp, dest))
    for tmp, dest in temps:
        os.replace(tmp, dest)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=ROOT_DEFAULT)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args(argv)
    root = a.root.resolve()
    if already_patched(root):
        raise SystemExit("already patched")
    partial = partial_patched(root)
    if partial:
        raise SystemExit(partial + "; restore and re-apply")
    before = {str(rel): sha256(root / rel) for rel in FILES}
    stage = Path(tempfile.mkdtemp(prefix="october538-"))
    try:
        stage_tree(root, stage)
        apply_sequence(stage)
        after = {str(rel): sha256(stage / rel) for rel in FILES}
        print_shas(before, after)
        if a.check:
            now = {str(rel): sha256(root / rel) for rel in FILES}
            if now != before:
                raise SystemExit("check wrote the destination")
            print("patches apply cleanly (not written)")
            return 0
        install(stage, root)
        print("patched " + ", ".join(str(root / rel) for rel in FILES))
        return 0
    finally:
        shutil.rmtree(stage, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
