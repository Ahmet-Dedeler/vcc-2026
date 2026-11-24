"""594 — composition test of the 10-05 sequence: 538's applier → 554's patch → 590's patch, end to end, on scratch trees.

Never touches the live files. The main run uses a git worktree of HEAD (`git worktree add --detach`) as the scratch root S,
with the gitignored inputs the suites read (control bundles, CCLE cache, donor arrays, rehearsal metas, audit dumps, .venv)
added as symlinks into S. Patches only write the 11 files listed in TOUCHED, all of them real (tracked) files in S.

Steps (S is checked after each):
  0  S == HEAD == live on the touched files; run every suite (538 test_extras538, 554 test_guard554 + identity554,
     590 test_identity590); bash -n def_day.sh; py_compile; `october_recipe.py check --allow-pending`.
  1  538: apply_october538.py --root S --check, then apply. Checks + 590's suite (now against the post-538 tree) + the
     post-538 shas equal a fresh reference application of 538 to copies of the live files.
  2  554: patch_guard554.py --root S --check, then apply. Checks + 590's suite (post-538+554) + def_preflight --verify +
     shas equal test_guard554.make_tree(with554=True).
  3  590: its suite has run against this exact tree in step 2; patch_identity590.py --root S --check, then apply. Checks +
     def_identity --verify on the applied tree + the def_sameline gate cases on the applied files, each compared with the
     step-2 (pre-590) def_sameline.
  4  Ordering hazards: every permutation of {538, 554, 590} on fresh temp trees (copies of the live files); which steps
     apply or are refused, and whether the final shas equal the canonical order's.
The 538 and 554 suites copy the LIVE-side files and apply their own patch(es) themselves, so they are pre-application
suites: they run at step 0 (on S == live) and are recorded n/a after their own patch is in.

Writes work/compose594.json and work/compose594.log. Usage: uv run python experiments/594-def-bundle-composition/compose594.py
"""

from __future__ import annotations

import hashlib
import itertools
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

R = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PY = str(R / ".venv/bin/python")
SCRATCH = Path(os.environ.get("SCRATCH594", "/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/"
                                          "e132c515-fa77-477b-acf1-db04d2a4f831/scratchpad"))
S = SCRATCH / "wt594"
LOG = HERE / "work/compose594.log"
TOUCHED = [
    "scripts/def_day.sh", "scripts/def_day_helpers.py", "experiments/448-def-rehearsal-q175/relevel_iter.py",
    "scripts/def_sameline.py", "scripts/october_recipe.py", "docs/october-decisions.json",
    "experiments/524-october-recipe-dr/gate524.py", "scripts/def_postbuild.py", "scripts/def_pipeline.py",
    "scripts/def_preflight.py", "scripts/def_identity.py",
]
EXTRA = ["scripts/frozen_inputs.json", "scripts/make_submission_pm_cd4.py"]  # 554's make_tree copies these too
OVERLAY = [
    "data/challenge", "data/external", "data/donors", "data/harness",
    "experiments/180-def-readiness/work/def-bundle", "experiments/238-essential-free-panel/work/ef-bundle",
    "experiments/238-essential-free-panel/work/ef-bundle-r", "experiments/534-generalization-audit/work/partB/S1/bundle",
    "experiments/552-r7-reserve/work/T1/bundle", "experiments/552-r7-reserve/work/r7b/bundle",
    "experiments/526-cloud-october-review2/work/efp_inert",
    "experiments/497-october-switches/work/efr497", "experiments/510-october-recipe-rehearsal/work/efr510",
    "experiments/452-def-day-rehearsal/work/efr", "experiments/474-october-hardening/work/efr_retrace",
    "experiments/433-coarse-mde-arm/work/rc_audit_433a.jsonl", "experiments/429-zs-production/work/rc_audit_429a.jsonl",
    "experiments/225-source-sweep/work/tables", "experiments/223-silent-zero-audit/work/tables",
    "experiments/166-donor-acquisition/donors", "experiments/120-myllia-a549/work/compact-donor", ".venv",
]
APPLY = {
    "538": lambda root, check: [PY, str(root_of(root) / "experiments/538-cloud-bundle-extras/apply_october538.py"),
                                "--root", str(root), *(["--check"] if check else [])],
    "554": lambda root, check: [PY, str(root_of(root) / "experiments/554-dose-guard/patch_guard554.py"),
                                "--root", str(root), *(["--check"] if check else [])],  # applies by default (no --apply)
    "590": lambda root, check: [PY, str(root_of(root) / "experiments/590-identity-validation/patch_identity590.py"),
                                "--root", str(root), "--check" if check else "--apply"],
}


def root_of(root: Path) -> Path:
    """The patch scripts run from S when S is the target (isolation), from the live repo for temp trees."""
    return S if root == S else R


def log(msg: str) -> None:
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {msg}"
    print(line, flush=True)
    with LOG.open("a") as fh:
        fh.write(line + "\n")


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def shas(root: Path) -> dict[str, str]:
    return {rel: sha(root / rel) for rel in TOUCHED}


def run(cmd: list[str], cwd: Path, timeout: int = 2400, env: dict | None = None) -> tuple[int, str]:
    e = {k: v for k, v in os.environ.items() if not k.startswith(("DEF_DAY_", "VCC_"))}
    e.update(env or {})
    try:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout, env=e)
        return p.returncode, (p.stdout + p.stderr)
    except subprocess.TimeoutExpired as ex:
        return 124, f"TIMEOUT after {timeout}s: {ex}"


def tail(text: str, n: int = 6) -> str:
    return "\n".join(text.strip().splitlines()[-n:])


def fill(rel: str) -> int:
    """Symlink every live path under `rel` that S lacks (recursively into dirs S already has). Returns links made."""
    src, dst = R / rel, S / rel
    if not src.exists() and not src.is_symlink():
        return 0
    if not dst.exists() and not dst.is_symlink():
        dst.parent.mkdir(parents=True, exist_ok=True)
        os.symlink(src, dst)
        return 1
    if dst.is_dir() and not dst.is_symlink() and src.is_dir():
        return sum(fill(f"{rel}/{c.name}") for c in sorted(src.iterdir()))
    return 0


def common_checks(root: Path, label: str) -> dict:
    out = {}
    rc, t = run(["/bin/bash", "-n", str(root / "scripts/def_day.sh")], root)
    out["bash -n def_day.sh"] = "PASS" if rc == 0 else f"FAIL: {tail(t)}"
    bad = []
    for rel in TOUCHED:
        if rel.endswith(".py"):
            rc, t = run([PY, "-m", "py_compile", str(root / rel)], root)
            if rc != 0:
                bad.append(f"{rel}: {tail(t, 2)}")
    out["py_compile (touched .py)"] = "PASS" if not bad else "FAIL: " + "; ".join(bad)
    rc, t = run([PY, str(root / "scripts/october_recipe.py"), "check", "--allow-pending"], root)
    out["october_recipe check --allow-pending"] = ("PASS" if rc == 0 else f"FAIL rc {rc}") + f" | {tail(t, 3)}"
    for k, v in out.items():
        log(f"[{label}] {k}: {v}")
    return out


def suite(name: str, script: str, label: str) -> str:
    t0 = time.time()
    rc, t = run([PY, str(S / script)], S)
    verdict = ("PASS" if rc == 0 else f"FAIL rc {rc}") + f" ({time.time() - t0:.0f} s) | {tail(t, 4)}"
    log(f"[{label}] suite {name}: {verdict}")
    return verdict


def gate_cases(label: str, pre590_sameline: Path) -> dict:
    """590's T2 cases on the APPLIED def_sameline in S, each against the pre-590 (post-538+554) def_sameline."""
    res = {}
    with tempfile.TemporaryDirectory(dir=HERE / "work") as td:
        tmp = Path(td)
        pre = tmp / "pre"
        (pre / "scripts").mkdir(parents=True)
        for name in ("data", "experiments"):
            os.symlink(S / name, pre / name)
        shutil.copy(pre590_sameline, pre / "scripts/def_sameline.py")
        cases = [("a", True, 0.554, [], False), ("b", True, 0.70, [], True), ("c", True, None, [], True),
                 ("d", True, 0.554, ["--min-top1", "0"], True), ("f", False, 0.70, [], False)]
        for name, unan, s1, xargs, want_on in cases:
            d = {"call": {"A": "HAP1_ENGINEERED"}, "call_unanimous_across_methods": {"A": unan}}
            if s1 is not None:
                d["min_top1_score"] = {"A": s1}
            idf = tmp / f"{name}.json"
            idf.write_text(json.dumps(d))
            outs = {}
            for side, root, extra in (("applied", S, xargs), ("pre590", pre, [a for a in xargs if a not in ("--min-top1", "0")])):
                o = tmp / f"{name}_{side}.npz"
                rc, t = run([PY, "scripts/def_sameline.py", "--bundle", str(S / "data/challenge"), "--identity", str(idf),
                             "--out", str(o), *extra], root)
                outs[side] = (rc, o, t)
            import numpy as np
            rc, o, t = outs["applied"]
            on = rc == 0 and "A" in list(np.load(o)["contexts"])
            same = None
            if want_on and rc == 0 and outs["pre590"][0] == 0:
                za, zb = np.load(o), np.load(outs["pre590"][1])
                same = sorted(za.files) == sorted(zb.files) and all(
                    np.array_equal(za[k], zb[k], equal_nan=za[k].dtype.kind == "f") for k in za.files)
            ok = (on == want_on) and (same in (None, True)) and (name != "c" or "WARNING" in t) and (name != "a" or "590" in t)
            res[name] = {"on": on, "want_on": want_on, "same_as_pre590": same, "PASS": ok}
            log(f"[{label}] gate case {name}: on={on} want={want_on} same_as_pre590={same} -> {'PASS' if ok else 'FAIL'}")
    return res


def permutations() -> dict:
    res = {}
    canon = None
    for order in itertools.permutations(["538", "554", "590"]):
        with tempfile.TemporaryDirectory(prefix="perm594-") as td:
            root = Path(td)
            for rel in TOUCHED + EXTRA:
                (root / rel).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy(R / rel, root / rel)
            steps = []
            for p in order:
                rc, t = run(APPLY[p](root, False), R)
                steps.append({"patch": p, "rc": rc, "msg": tail(t, 2)})
                if rc != 0:
                    break
            final = shas(root) if all(s["rc"] == 0 for s in steps) else None
            key = "→".join(order)
            res[key] = {"steps": steps, "final": final}
            if order == ("538", "554", "590"):
                canon = final
            log(f"[perm] {key}: " + ", ".join(f"{s['patch']} rc {s['rc']}" for s in steps))
    for key, v in res.items():
        v["final_equals_canonical"] = (v["final"] == canon) if v["final"] is not None else None
    return res


def main() -> None:
    (HERE / "work").mkdir(parents=True, exist_ok=True)
    LOG.write_text("")
    out: dict = {"started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "head": None, "steps": {}}
    live_before = shas(R)
    out["head"] = subprocess.run(["git", "rev-parse", "HEAD"], cwd=R, capture_output=True, text=True).stdout.strip()
    # the scratch worktree
    if S.exists():
        subprocess.run(["git", "worktree", "remove", "--force", str(S)], cwd=R, capture_output=True, text=True)
        shutil.rmtree(S, ignore_errors=True)
    rc = subprocess.run(["git", "worktree", "add", "--detach", str(S), out["head"]], cwd=R, capture_output=True, text=True)
    if rc.returncode != 0:
        raise SystemExit(f"worktree add failed: {rc.stderr}")
    links = sum(fill(rel) for rel in OVERLAY)
    log(f"worktree {S} at {out['head'][:12]}; {links} overlay symlinks to live gitignored inputs")
    try:
        s0 = shas(S)
        assert s0 == live_before, "live differs from HEAD on a touched file: the composition would not describe 10-05"
        st = {"shas": s0, "checks": common_checks(S, "step0")}
        st["suites"] = {n: suite(n, sc, "step0") for n, sc in (
            ("538 test_extras538", "experiments/538-cloud-bundle-extras/test_extras538.py"),
            ("554 test_guard554", "experiments/554-dose-guard/test_guard554.py"),
            ("554 identity554", "experiments/554-dose-guard/identity554.py"),
            ("590 test_identity590", "experiments/590-identity-validation/test_identity590.py"))}
        out["steps"]["0 (HEAD)"] = st

        # step 1: 538
        rcc, tc = run(APPLY["538"](S, True), S)
        rca, ta = run(APPLY["538"](S, False), S)
        st = {"check": f"rc {rcc} | {tail(tc, 2)}", "apply": f"rc {rca} | {tail(ta, 2)}", "shas": shas(S)}
        with tempfile.TemporaryDirectory(prefix="ref538-") as td:
            ref = Path(td)
            for rel in TOUCHED:
                (ref / rel).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy(R / rel, ref / rel)
            rcr, _ = run(APPLY["538"](ref, False), R)
            st["equals_reference_538_application"] = rcr == 0 and shas(ref) == st["shas"]
        st["checks"] = common_checks(S, "step1")
        st["suites"] = {"538 test_extras538": "n/a (pre-application suite; ran at step 0)",
                        "554 test_guard554": "n/a (applies 538 to copies itself; ran at step 0)",
                        "554 identity554": "n/a (ran at step 0)",
                        "590 test_identity590": suite("590 test_identity590", "experiments/590-identity-validation/test_identity590.py", "step1")}
        log(f"[step1] check rc {rcc}, apply rc {rca}, equals reference: {st['equals_reference_538_application']}")
        out["steps"]["1 (538)"] = st

        # step 2: 554
        rcc, tc = run(APPLY["554"](S, True), S)
        rca, ta = run(APPLY["554"](S, False), S)
        st = {"check": f"rc {rcc} | {tail(tc, 2)}", "apply": f"rc {rca} | {tail(ta, 2)}", "shas": shas(S)}
        sys.path.insert(0, str(R / "experiments/554-dose-guard"))
        import test_guard554 as G  # live code: make_tree copies the LIVE files and applies 538 then 554
        with tempfile.TemporaryDirectory(prefix="ref554-") as td:
            ref = Path(td) / "t"
            G.make_tree(ref, True)
            st["equals_test_guard554_make_tree"] = {rel: sha(ref / rel) == st["shas"][rel]
                                                    for rel in TOUCHED if (ref / rel).exists()}
        st["checks"] = common_checks(S, "step2")
        rcp, tp = run([PY, str(S / "scripts/def_preflight.py"), "--verify"], S)
        st["checks"]["def_preflight --verify"] = ("PASS" if rcp == 0 else f"FAIL rc {rcp}") + f" | {tail(tp, 2)}"
        log(f"[step2] def_preflight --verify: {st['checks']['def_preflight --verify']}")
        st["suites"] = {"538 test_extras538": "n/a", "554 test_guard554": "n/a (pre-application suite; ran at step 0)",
                        "554 identity554": "n/a (ran at step 0)",
                        "590 test_identity590": suite("590 test_identity590", "experiments/590-identity-validation/test_identity590.py", "step2")}
        log(f"[step2] check rc {rcc}, apply rc {rca}")
        out["steps"]["2 (554)"] = st
        pre590 = HERE / "work/def_sameline.pre590.py"
        shutil.copy(S / "scripts/def_sameline.py", pre590)

        # step 3: 590
        rcc, tc = run(APPLY["590"](S, True), S)
        rca, ta = run(APPLY["590"](S, False), S)
        st = {"check": f"rc {rcc} | {tail(tc, 3)}", "apply": f"rc {rca} | {tail(ta, 2)}", "shas": shas(S)}
        st["checks"] = common_checks(S, "step3")
        rcv, tv = run([PY, str(S / "scripts/def_identity.py"), "--verify"], S)
        lines = [ln.strip() for ln in tv.splitlines() if "lowest top-1" in ln]
        st["checks"]["def_identity --verify"] = ("PASS" if rcv == 0 and "VERIFY PASS" in tv else f"FAIL rc {rcv}") + \
            " | " + "; ".join(lines)
        log(f"[step3] def_identity --verify: {st['checks']['def_identity --verify']}")
        st["gate_cases"] = gate_cases("step3", pre590)
        rc2, t2 = run(APPLY["590"](S, True), S)
        st["second_application"] = "refused" if rc2 != 0 and "already" in t2 else f"NOT refused (rc {rc2})"
        log(f"[step3] check rc {rcc}, apply rc {rca}; second application {st['second_application']}")
        out["steps"]["3 (590)"] = st

        # step 4: order permutations
        out["permutations"] = permutations()
        out["final_shas"] = out["steps"]["3 (590)"]["shas"]
    finally:
        live_after = shas(R)
        out["live_untouched"] = live_after == live_before
        log(f"live files untouched: {out['live_untouched']}")
        subprocess.run(["git", "worktree", "remove", "--force", str(S)], cwd=R, capture_output=True, text=True)
        out["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        (HERE / "work/compose594.json").write_text(json.dumps(out, indent=1) + "\n")
        log("wrote work/compose594.json; worktree removed")


if __name__ == "__main__":
    main()
