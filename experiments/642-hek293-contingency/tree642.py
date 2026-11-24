"""642: build a POST-BUNDLE tree (copies; the live repo is never written): 538's applier (530 + 528 + 538 extras), then
554, 590 and 641, exactly the order of the 10-05 checklist; data/, submissions/ and experiments/ are read-only symlinks.

    from tree642 import make_tree; make_tree(dest)            # -> dest (post-bundle, 642 NOT applied)
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PY = ROOT / ".venv/bin/python"
sys.path.insert(0, str(ROOT / "experiments/554-dose-guard"))
import test_guard554 as G  # noqa: E402

EXTRA642 = [Path("scripts/def_identity.py")]


def make_tree(dest: Path) -> Path:
    for rel in G.REL538 + G.EXTRA + EXTRA642:
        t = dest / rel
        t.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / rel, t)
    steps = [[str(PY), str(G.APPLIER538), "--root", str(dest)],
             [str(PY), str(G.PATCH), "--root", str(dest)],
             [str(PY), str(ROOT / "experiments/590-identity-validation/patch_identity590.py"), "--root", str(dest), "--apply"],
             [str(PY), str(ROOT / "experiments/641-depth-row/patch_depth641.py"), "--root", str(dest)]]
    for cmd in steps:
        # B2 ([4d0207] for the coordinator): after the 10-05 apply the live files already carry 538 / 554 / 590 / 641, so a
        # step whose --check says "already" is skipped (a PARTIALLY patched tree still fails at the apply below)
        chk = [x for x in cmd if x != "--apply"] + ["--check"]
        c = subprocess.run(chk, capture_output=True, text=True)
        co = (c.stdout + c.stderr).lower()
        if c.returncode != 0 and "already" in co and "partially" not in co:
            continue
        r = subprocess.run(cmd, capture_output=True, text=True)
        assert r.returncode == 0, (cmd, r.stdout[-2000:], r.stderr[-2000:])
    for d in ("data", "submissions"):
        os.symlink(ROOT / d, dest / d)
    # experiments/: some bundle files live there (copied above), so link every other entry one level down
    (dest / "experiments").mkdir(exist_ok=True)
    for e in sorted((ROOT / "experiments").iterdir()):
        if not (dest / "experiments" / e.name).exists():
            os.symlink(e, dest / "experiments" / e.name)
    return dest


if __name__ == "__main__":
    print(make_tree(Path(sys.argv[1])))
