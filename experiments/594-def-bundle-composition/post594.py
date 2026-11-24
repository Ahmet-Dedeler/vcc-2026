"""594 — post-apply check for 10-05: after 538 → 554 → 590 are applied to the LIVE tree, confirm the result equals the
composition test's. Read-only on the live tree: it writes only temp files.

1. The 11 touched files equal a canonical application of 538 -> 554 -> 590 to git HEAD's copies of the same files,
   rebuilt in a temp tree (HEAD is still the pre-bundle commit at checklist step 6, before step 9 commits). This is
   robust to legitimate HEAD changes since 594 (e.g. october-decisions.json gains decisions daily); 594's frozen table
   is printed beside it for reference.
2. `scripts/def_identity.py --verify` passes and prints the lowest top-1 correlation for A/B/C (0.772 / 0.719 / 0.733).
3. The 590 gate on the live (applied) `def_sameline.py`, against the pre-590 file rebuilt in a temp tree from git HEAD's
   `def_sameline.py` + 538's applier (its sha must be 03cd591a…, the post-538 file 594 saw): a registered unanimous call at
   0.554 is blocked; at 0.70, without the key, and with `--min-top1 0` it switches on with the pre-590 file's exact output;
   a non-unanimous call stays off.
Usage: uv run python experiments/594-def-bundle-composition/post594.py   (exit 0 = all PASS)
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

R = Path(__file__).resolve().parents[2]
PY = str(R / ".venv/bin/python")
FINAL = {
    "scripts/def_day.sh": "5851f84b9f238c932816d22fbe6cd810c6db68c863db6bdb41d8746d10932f6c",
    "scripts/def_day_helpers.py": "e3ff45bcefc1a6d3af651126903726e51972c61343386752f36be19bfd8215f3",
    "experiments/448-def-rehearsal-q175/relevel_iter.py": "e360a1d9907801a1b5f67832a1e0922a10cf12f605b765eaabbd38739bf71b35",
    "scripts/def_sameline.py": "22c6f0f08c1f02883887a34690bc5642f172026433d90b19275a578f24d86259",
    "scripts/october_recipe.py": "503c45a00bf74164a206f13cfcd878c468d6d5d8783a1a04068c32f42bb2c793",
    "docs/october-decisions.json": "8ee9e74dee8ec0c84de2e8d263f4bc2cc05288ca7c112feda57a10c11e82e3db",
    "experiments/524-october-recipe-dr/gate524.py": "7a43d2f5ae0f03c017a6903b9d524e6e7f9352de0f32f4641b0d3a4dba0cf5ae",
    "scripts/def_postbuild.py": "557c2b522b565636261ea297cc4dcc2e5e1af8f1ea0bbc8c788790388f0b68f6",
    "scripts/def_pipeline.py": "9a718b9eaa24f10cfbf255ce830ea4a748f61252ec503bbb6db69d3ddf5ca9fd",
    "scripts/def_preflight.py": "2d1811945ebeb506cf2bde739c03c57f5350e84311a4476d11337a2355b1622d",
    "scripts/def_identity.py": "c98a608fe2838d1ca6b68dbf090874548f286c911f472becf8ef3dcfcc157f59",
}
PRE590_SAMELINE = "03cd591a82cc"
REL538 = ["scripts/def_day.sh", "scripts/def_day_helpers.py", "experiments/448-def-rehearsal-q175/relevel_iter.py",
          "scripts/def_sameline.py", "scripts/october_recipe.py", "docs/october-decisions.json",
          "experiments/524-october-recipe-dr/gate524.py", "scripts/def_postbuild.py", "scripts/def_pipeline.py"]


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> None:
    fails = []
    # 1: expected = HEAD's files + 538 -> 554 -> 590 in a temp tree
    extra = ["scripts/frozen_inputs.json", "scripts/make_submission_pm_cd4.py"]
    with tempfile.TemporaryDirectory() as td:
        canon = Path(td) / "canon"
        for rel in list(FINAL) + extra:
            (canon / rel).parent.mkdir(parents=True, exist_ok=True)
            (canon / rel).write_bytes(subprocess.run(["git", "show", f"HEAD:{rel}"], cwd=R, capture_output=True, check=True).stdout)
        steps = [[PY, str(R / "experiments/538-cloud-bundle-extras/apply_october538.py"), "--root", str(canon)],
                 [PY, str(R / "experiments/554-dose-guard/patch_guard554.py"), "--root", str(canon)],
                 [PY, str(R / "experiments/590-identity-validation/patch_identity590.py"), "--root", str(canon), "--apply"]]
        for cmd in steps:
            q = subprocess.run(cmd, capture_output=True, text=True)
            if q.returncode != 0:
                print(f"[1] canonical application failed at {Path(cmd[1]).name}: {q.stdout[-300:]} {q.stderr[-300:]}")
                fails.append("canonical application")
                break
        for rel, frozen in FINAL.items():
            want = sha(canon / rel) if (canon / rel).exists() else "?"
            got = sha(R / rel)
            ok = got == want
            note = "" if want == frozen else f"   (594's frozen table: {frozen[:16]}; HEAD moved since 594)"
            print(f"[1] {'OK  ' if ok else 'DIFF'} {rel} live {got[:16]} canonical {want[:16]}{note}")
            if not ok:
                fails.append(f"sha {rel}")
    # 2
    p = subprocess.run([PY, "scripts/def_identity.py", "--verify"], cwd=R, capture_output=True, text=True)
    lines = [ln.strip() for ln in p.stdout.splitlines() if "lowest top-1" in ln]
    ok = p.returncode == 0 and "VERIFY PASS" in p.stdout and len(lines) == 3
    print(f"[2] def_identity --verify {'PASS' if ok else 'FAIL'}: " + " | ".join(lines))
    if not ok:
        fails.append("def_identity --verify")
    # 3
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        rebuild = tmp / "rebuild"
        for rel in REL538:
            (rebuild / rel).parent.mkdir(parents=True, exist_ok=True)
            blob = subprocess.run(["git", "show", f"HEAD:{rel}"], cwd=R, capture_output=True).stdout
            (rebuild / rel).write_bytes(blob)
        r = subprocess.run([PY, str(R / "experiments/538-cloud-bundle-extras/apply_october538.py"), "--root", str(rebuild)],
                           capture_output=True, text=True)
        pre_sha = sha(rebuild / "scripts/def_sameline.py")
        ok = r.returncode == 0
        print(f"[3] pre-590 def_sameline rebuilt from HEAD + 538: {pre_sha[:12]} "
              f"({'as in 594' if pre_sha.startswith(PRE590_SAMELINE) else 'HEAD moved since 594'}; applier rc {r.returncode})")
        if not ok:
            fails.append("pre-590 rebuild")
        pre = tmp / "pre"
        (pre / "scripts").mkdir(parents=True)
        for name in ("data", "experiments"):
            os.symlink(R / name, pre / name)
        shutil.copy(rebuild / "scripts/def_sameline.py", pre / "scripts/def_sameline.py")
        for name, unan, s1, xargs, want_on in [("a", True, 0.554, [], False), ("b", True, 0.70, [], True),
                                               ("c", True, None, [], True), ("d", True, 0.554, ["--min-top1", "0"], True),
                                               ("f", False, 0.70, [], False)]:
            d = {"call": {"A": "HAP1_ENGINEERED"}, "call_unanimous_across_methods": {"A": unan}}
            if s1 is not None:
                d["min_top1_score"] = {"A": s1}
            idf = tmp / f"{name}.json"
            idf.write_text(json.dumps(d))
            res = {}
            for side, root, extra in (("live", R, xargs), ("pre590", pre, [a for a in xargs if a not in ("--min-top1", "0")])):
                o = tmp / f"{name}_{side}.npz"
                q = subprocess.run([PY, "scripts/def_sameline.py", "--bundle", str(R / "data/challenge"), "--identity", str(idf),
                                    "--out", str(o), *extra], cwd=root, capture_output=True, text=True)
                res[side] = (q.returncode, o, q.stdout)
            rc, o, out = res["live"]
            on = rc == 0 and "A" in list(np.load(o)["contexts"])
            same = None
            if want_on and rc == 0 and res["pre590"][0] == 0:
                za, zb = np.load(o), np.load(res["pre590"][1])
                same = sorted(za.files) == sorted(zb.files) and all(
                    np.array_equal(za[k], zb[k], equal_nan=za[k].dtype.kind == "f") for k in za.files)
            ok = on == want_on and same in (None, True) and (name != "c" or "WARNING" in out) and (name != "a" or "590" in out)
            print(f"[3{name}] on={on} want={want_on} same_as_pre590={same} -> {'PASS' if ok else 'FAIL'}")
            if not ok:
                fails.append(f"gate {name}")
    print("ALL PASS" if not fails else f"FAILED: {fails}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
