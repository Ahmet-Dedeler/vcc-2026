#!/usr/bin/env python3
"""NOT APPLIED. On a 538-patched october_recipe.py, refuse two October combinations.

  * DEF_DAY_MDE_SIGN on while DEF_DAY_CALLSET_BLOCK is not prop (the hook without the block).
    Override: the mde_sign switch's evidence contains the token hook-alone-explicit.
  * DEF_DAY_MDE_SIGN whose source is not clean4 (programme, hookbank, file: Arc or a file hook).
    Override: the evidence contains the token arc-source-explicit.

clean4:0.05 together with callset_block=prop needs neither token. `decide` does not add a token
on its own. Refuses to edit a path inside the live repo.

  uv run python experiments/550-cloud-1005-review/repro/patch_hook_conditional550.py --self-test
  uv run python experiments/550-cloud-1005-review/repro/patch_hook_conditional550.py --apply-to <copy>
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

CALL = "    _problems_538(env, problems)\n"
CALL_NEW = "    _problems_538(env, problems, s.get(\"mde_sign\") or {})\n"

DEF = '''def _problems_538(env: dict, problems: list) -> None:
    """Range and dependency checks for 538's knobs. PENDING stays out of env."""
'''

DEF_NEW = '''def _problems_538(env: dict, problems: list, mde_sign: dict | None = None) -> None:
    """Range and dependency checks for 538's knobs. PENDING stays out of env.

    550: the hook belongs with the block, and CLEAN4 is the default source (537 §5.4).
    A non-clean4 source, or a hook with the block off, needs an explicit token in the
    switch's evidence. `decide` does not write either token.
    """
    sign = env.get("DEF_DAY_MDE_SIGN", "")
    if sign:
        block = env.get("DEF_DAY_CALLSET_BLOCK", "")
        evidence = str((mde_sign or {}).get("evidence") or "")
        on_block = block == "prop" or block.startswith("prop:")
        if not on_block and "hook-alone-explicit" not in evidence:
            problems.append(
                "mde_sign is on without callset_block=prop "
                "(the hook belongs with the block; evidence token hook-alone-explicit overrides)"
            )
        source = sign.split(":", 1)[0].lower()
        if source != "clean4" and "arc-source-explicit" not in evidence:
            problems.append(
                f"mde_sign source {source!r} is not clean4 "
                "(CLEAN4 is the default hook source; evidence token arc-source-explicit overrides)"
            )
'''


def live_root() -> Path:
    return Path(__file__).resolve().parents[3]


def patch_text(src: str) -> str:
    if CALL not in src or src.count(CALL) != 1:
        raise SystemExit("call anchor missing or not unique")
    if DEF not in src or src.count(DEF) != 1:
        raise SystemExit("def anchor missing or not unique")
    if "hook-alone-explicit" in src:
        raise SystemExit("already patched")
    return src.replace(CALL, CALL_NEW, 1).replace(DEF, DEF_NEW, 1)


def apply_to(path: Path) -> None:
    target = path.resolve()
    root = live_root()
    if target == root or root in target.parents:
        raise SystemExit(f"refusing to edit {target}: it is inside the live repo")
    target.write_text(patch_text(target.read_text()))
    print(f"patched copy {target}")


def self_test() -> int:
    import json
    import os
    import shutil
    import subprocess
    import tempfile

    root = live_root()
    files = [
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
    before = {str(p): (root / p).read_bytes() for p in files}
    dest = Path(tempfile.mkdtemp(prefix="550hook-"))
    for rel in files:
        t = dest / rel
        t.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(root / rel, t)
    py = sys.executable
    r = subprocess.run([py, str(root / "experiments/538-cloud-bundle-extras/apply_october538.py"),
                        "--root", str(dest)], capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout, r.stderr)
        return 1
    recipe = dest / "scripts/october_recipe.py"
    decisions = dest / "docs/october-decisions.json"

    def check(label: str) -> subprocess.CompletedProcess:
        env = os.environ.copy()
        env.update({
            "OCTOBER_DECISIONS": str(decisions),
            "OCTOBER_DEF_DAY": str(dest / "scripts/def_day.sh"),
            "OCTOBER_BUILDER": str(root / "scripts/make_submission_pm_cd4.py"),
            "OCTOBER_EMITTER": str(root / "src/vcc2026/rc_emit_k.py"),
            "OCTOBER_CALLBLOCK": str(root / "experiments/528-callprop-block/callblock528.py"),
            "OCTOBER_SAMELINE": str(dest / "scripts/def_sameline.py"),
        })
        proc = subprocess.run([py, str(recipe), "check", "--allow-pending"], capture_output=True, text=True, env=env)
        print(f"[{label}] rc {proc.returncode}")
        print((proc.stdout or "").strip())
        err = (proc.stderr or "").strip()
        if err:
            print(err)
        return proc

    def set_switches(mde_value, mde_evidence, block_value):
        d = json.loads(decisions.read_text())
        d["switches"]["mde_sign"].update({
            "status": "DECIDED", "value": mde_value, "evidence": mde_evidence,
        })
        if block_value is None:
            d["switches"]["callset_block"].update({"status": "PENDING", "value": None, "evidence": ""})
        else:
            d["switches"]["callset_block"].update({
                "status": "DECIDED", "value": block_value, "evidence": "550 block on",
            })
        decisions.write_text(json.dumps(d, indent=1) + "\n")

    set_switches("clean4:0.05", "496a transfers, block still pending", None)
    bare = check("unpatched hook without block")
    env_line = bare.stdout.split("pending")[0]
    if bare.returncode != 0 or "DEF_DAY_MDE_SIGN" not in env_line or "DEF_DAY_CALLSET_BLOCK" in env_line:
        print("FAIL: unpatched recipe should accept the hook with the block off")
        return 1
    apply_to(recipe)
    refused = check("patched hook without block")
    if refused.returncode == 0 or "hook belongs with the block" not in refused.stdout:
        print("FAIL: patched recipe did not refuse the hook without the block")
        return 1
    set_switches("clean4:0.05", "block and CLEAN4", "prop")
    ok = check("patched clean4 with prop")
    if ok.returncode != 0 or "DEF_DAY_MDE_SIGN" not in ok.stdout or "DEF_DAY_CALLSET_BLOCK" not in ok.stdout:
        print("FAIL: clean4 + prop should pass")
        return 1
    set_switches("hookbank:h1_2025_full:0.05:0.05:le1", "496a and the block", "prop")
    arc = check("patched arc source, no token")
    if arc.returncode == 0 or "not clean4" not in arc.stdout:
        print("FAIL: arc source without the token should be refused")
        return 1
    set_switches("hookbank:h1_2025_full:0.05:0.05:le1", "arc-source-explicit D/E/F read", "prop")
    arc_ok = check("patched arc source, explicit token")
    if arc_ok.returncode != 0:
        print("FAIL: explicit arc source with the block on should pass")
        return 1
    set_switches("clean4:0.05", "hook-alone-explicit 496a and a D/E/F entry", None)
    alone = check("patched hook alone, explicit token")
    if alone.returncode != 0:
        print("FAIL: explicit hook-alone should pass")
        return 1
    for rel, blob in before.items():
        if (root / rel).read_bytes() != blob:
            print(f"FAIL: live {rel} changed")
            return 1
    print("HOOK CONDITIONAL SELF-TEST PASS")
    shutil.rmtree(dest, ignore_errors=True)
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--apply-to", type=Path)
    a = ap.parse_args(argv)
    if a.self_test:
        return self_test()
    if a.apply_to is None:
        print("not applied. Pass --self-test, or --apply-to a copy of the 538-patched october_recipe.py outside the repo.")
        return 0
    apply_to(a.apply_to)
    return 0


if __name__ == "__main__":
    sys.exit(main())
