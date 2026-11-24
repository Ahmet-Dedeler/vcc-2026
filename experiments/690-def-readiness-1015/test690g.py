#!/usr/bin/env python3
"""690g tests. Light (a real 89 MB stack dump is copied a few times; no build). Scratch tree = the live files with 690a applied.
  G1 --check/apply/second apply refused; refused without 690a     G2 bash -n
  G3 KA: `command --allow-pending --json` env and args unchanged with dump_hook PENDING
  G4 def_dumphook_check on a real stack dump: identical, moved, env changed, perts changed, NaN, per-context good/bad
  G5 the inserted hook block (cut from the patched def_day.sh) in a bash harness: good hook; good + bulk_ctx with and
     without 672 in the builder; failing hook; invalid output; unset knob; the validation block (git-committed hook only)
  G6 recipe: a decided existing path reaches env; a missing path is a problem
  G7 scripts/test_october_recipe.py green on the patched tree, pending and all-decided
usage: uv run python experiments/690-def-readiness-1015/test690g.py
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PY = sys.executable
DUMP = ROOT / "experiments/658-prereg-1007-builds/work/dr2/run/dump_h1x2_cd4x1.npz"
BUNDLE = ROOT / "experiments/180-def-readiness/work/def-bundle"   # contexts D, E, F
FILES = ["scripts/def_day.sh", "scripts/october_recipe.py", "docs/october-decisions.json",
         "scripts/make_submission_pm_cd4.py", "scripts/def_sameline.py", "src/vcc2026/rc_emit_k.py",
         "experiments/528-callprop-block/callblock528.py", "scripts/test_october_recipe.py",
         "experiments/497-october-switches/work/def_day.497.sh"]
RESULTS: list[tuple[bool, str]] = []


def rec(ok: bool, name: str) -> None:
    RESULTS.append((bool(ok), name))
    print(f"{'PASS' if ok else 'FAIL'} {name}", flush=True)


def tree(td: Path, name: str, with690a: bool = True) -> Path:
    t = td / name
    for f in FILES:
        (t / f).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / f, t / f)
    if with690a:
        r = subprocess.run([PY, str(HERE / "patch690a.py"), "--root", str(t)], capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
    return t


def run(*a, **kw) -> subprocess.CompletedProcess:
    return subprocess.run([str(x) for x in a], capture_output=True, text=True, **kw)


def write_dump(p: Path, base: dict, **over) -> Path:
    d = {k: base[k] for k in ("bulk", "carriers", "perts", "genes", "env")}
    d.update(over)
    np.savez(p, **d)
    return p


def main() -> int:
    td = Path(tempfile.mkdtemp(prefix="t690g_"))
    try:
        # ---- G1, G2
        t = tree(td, "a")
        rec(run(PY, HERE / "patch690g.py", "--root", t, "--check").returncode == 0, "G1 --check on a post-690a tree")
        before = (t / "scripts/def_day.sh").read_text()
        u = tree(td, "u")   # post-690a, unpatched by 690g: the KA reference
        rec(run(PY, HERE / "patch690g.py", "--root", t).returncode == 0, "G1 apply")
        rec(run(PY, HERE / "patch690g.py", "--root", t).returncode == 1, "G1 second apply refused")
        n = tree(td, "n", with690a=False)
        r = run(PY, HERE / "patch690g.py", "--root", n, "--check")
        rec(r.returncode == 1 and "690a" in r.stderr, "G1 refused without 690a")
        rec(run("bash", "-n", t / "scripts/def_day.sh").returncode == 0, "G2 bash -n")
        # ---- G3
        args = ["command", "--bundle", "B", "--work", "W", "--name", "N", "--allow-pending", "--json"]
        j0 = json.loads(run(PY, u / "scripts/october_recipe.py", *args, cwd=u).stdout)
        j1 = json.loads(run(PY, t / "scripts/october_recipe.py", *args, cwd=t).stdout)
        rec(j0["env"] == j1["env"] and j0["args"] == j1["args"]
            and set(j1["pending_treated_off"]) - set(j0["pending_treated_off"]) == {"dump_hook"},
            "G3 KA: env and args unchanged; dump_hook only extends PENDING")
        # ---- G4
        base = dict(np.load(DUMP, allow_pickle=False))
        chk = t / "scripts/def_dumphook_check.py"
        def check(hook: Path, ctx: Path | None = None) -> subprocess.CompletedProcess:
            return run(PY, chk, DUMP, hook, BUNDLE, *([ctx] if ctx else []))
        same = write_dump(td / "same.npz", base)
        rec(check(same).returncode == 0, "G4 identical dump PASS")
        rng = np.random.default_rng(0)
        b = base["bulk"].copy()
        v = rng.normal(size=b.shape)
        nb = np.linalg.norm(b, axis=1, keepdims=True)
        mv = b / np.where(nb > 0, nb, 1) * 0.7 + v / np.linalg.norm(v, axis=1, keepdims=True) * 0.3
        mv = mv / np.linalg.norm(mv, axis=1, keepdims=True) * nb
        moved = write_dump(td / "moved.npz", base, bulk=mv)
        r = check(moved)
        rec(r.returncode == 0 and "mean cosine" in r.stdout, f"G4 moved dump PASS ({[l for l in r.stdout.splitlines() if 'cosine' in l][:1]})")
        rec(check(write_dump(td / "env.npz", base, env=np.array("other env"))).returncode == 1, "G4 env changed FAILS")
        p2 = base["perts"].copy(); p2[0] = "XXX"
        rec(check(write_dump(td / "perts.npz", base, perts=p2)).returncode == 1, "G4 perts changed FAILS")
        b2 = base["bulk"].copy(); b2[3, 5] = np.nan
        rec(check(write_dump(td / "nan.npz", base, bulk=b2)).returncode == 1, "G4 NaN FAILS")
        ctx_ok = td / "ctx_ok.npz"
        np.savez(ctx_ok, perts=base["perts"], genes=base["genes"], contexts=np.array(["D", "E", "F"]),
                 bulk_D=mv, bulk_E=base["bulk"], bulk_F=mv)
        rec(check(same, ctx_ok).returncode == 0, "G4 per-context npz PASS")
        ctx_bad = td / "ctx_bad.npz"
        np.savez(ctx_bad, perts=base["perts"], genes=base["genes"], contexts=np.array(["A", "B", "C"]),
                 bulk_A=mv, bulk_B=mv, bulk_C=mv)
        rec(check(same, ctx_bad).returncode == 1, "G4 per-context npz with the wrong contexts FAILS")
        # ---- G5: the hook block in a harness (cwd = the scratch tree, made a git repo so the hook can be committed)
        dd = (t / "scripts/def_day.sh").read_text()
        hblk = dd[dd.index("# 690g: dump hook -- the hook's dump replaces"):dd.index('[ -n "$BLC" ] && { export VCC_BULK_LOAD_CTX')]
        hblk += dd[dd.index('[ -n "$BLC" ] && { export VCC_BULK_LOAD_CTX'):].split("\n", 1)[0] + "\n"
        vblk = dd[dd.index("# 690g: dump hook -- DEF_DAY_DUMP_HOOK=<script> (validated"):dd.index("# 502: the DE-channel riders")]
        run("git", "init", "-q", cwd=t)
        hooks = {
            "good.sh": 'cp "$3" "$5/$(basename "$3")"\n',
            "ctx.sh": f'cp "$3" "$5/$(basename "$3")"; cp "{ctx_ok}" "$5/bulk_ctx.npz"\n',
            "fail.sh": "exit 3\n",
            "badenv.sh": f'cp "{td / "env.npz"}" "$5/$(basename "$3")"\n',
        }
        (t / "hooks").mkdir()
        for k, s in hooks.items():
            (t / "hooks" / k).write_text("#!/bin/bash\n" + s)
        run("git", "add", "-A", cwd=t)
        run("git", "-c", "user.email=t@t", "-c", "user.name=t", "commit", "-qm", "t", cwd=t)
        w = td / "w"
        def harness(hook: str | None, builder672: bool = False) -> subprocess.CompletedProcess:
            if w.exists():
                shutil.rmtree(w)
            w.mkdir()
            shutil.copy2(DUMP, w / "dump_h1x2_cd4x1.npz")
            (w / "identity.json").write_text("{}")
            bpath = t / "scripts/make_submission_pm_cd4.py"
            btxt = (ROOT / "scripts/make_submission_pm_cd4.py").read_text()
            bpath.write_text(btxt + ("\n# 672: per-context bulk load (test marker)\n" if builder672 else ""))
            env = {k: x for k, x in os.environ.items() if not k.startswith(("DEF_DAY_", "VCC_"))}
            env.update(W=str(w), B=str(BUNDLE), LOG=str(w / "def_day.log"), SQ="1.75")
            if hook:
                env["DEF_DAY_DUMP_HOOK"] = f"hooks/{hook}"
            env["PYBIN"] = PY
            sh = ('log(){ echo "LOG: $*" >> "$LOG"; }\ndie(){ echo "DIE: $*"; exit 7; }\n'
                  'uv(){ if [ "$1" = run ] && [ "$2" = python ]; then shift 2; "$PYBIN" "$@"; else command uv "$@"; fi; }\n'
                  f'DUMP="$W/dump_h1x2_cd4x1.npz"\nDSHA="$(shasum -a 256 "$DUMP" | cut -d\' \' -f1)"\n'
                  + vblk + hblk +
                  'echo "LOAD=$VCC_BULK_LOAD"; echo "CTX=${VCC_BULK_LOAD_CTX:-unset}"; echo "DUMP=$DUMP"\n')
            return run("bash", "-c", sh, cwd=t, env=env)
        r = harness(None)
        rec(r.returncode == 0 and "CTX=unset" in r.stdout and r.stdout.count("dump_h1x2_cd4x1.npz") >= 2
            and "dumphook" not in r.stdout, "G5 unset: the stack dump is exported unchanged")
        r = harness("good.sh")
        rec(r.returncode == 0 and "dumphook/dump_h1x2_cd4x1.npz" in r.stdout and "CTX=unset" in r.stdout,
            "G5 good hook: the hook's dump (same name) is exported")
        r = harness("ctx.sh", builder672=True)
        rec(r.returncode == 0 and "CTX=" in r.stdout and "bulk_ctx.npz:" in r.stdout, "G5 hook + bulk_ctx with 672: VCC_BULK_LOAD_CTX exported")
        r = harness("ctx.sh", builder672=False)
        rec(r.returncode == 7 and "no VCC_BULK_LOAD_CTX (672)" in r.stdout, "G5 bulk_ctx without 672 in the builder: STOP")
        r = harness("fail.sh")
        rec(r.returncode == 7 and "the dump hook failed" in r.stdout, "G5 failing hook: STOP")
        r = harness("badenv.sh")
        rec(r.returncode == 7 and "failed validation" in r.stdout, "G5 invalid hook output: STOP")
        (t / "hooks/good.sh").write_text("#!/bin/bash\n# edited\n" + hooks["good.sh"])
        r = harness("good.sh")
        rec(r.returncode == 7 and "uncommitted changes" in r.stdout, "G5 an uncommitted hook edit: STOP")
        (t / "hooks/new.sh").write_text("#!/bin/bash\ntrue\n")
        r = harness("new.sh")
        rec(r.returncode == 7 and "must be a committed script" in r.stdout, "G5 an uncommitted hook: STOP")
        # ---- G6
        d = json.loads((t / "docs/october-decisions.json").read_text())
        for hp, want_ok in (("hooks/fail.sh", True), ("hooks/nope.sh", False)):
            d6 = json.loads(json.dumps(d))
            d6["switches"]["dump_hook"].update(status="DECIDED", value=hp)
            f = td / "d6.json"
            f.write_text(json.dumps(d6))
            r = run(PY, t / "scripts/october_recipe.py", *args, cwd=t, env=dict(os.environ, OCTOBER_DECISIONS=str(f)))
            got = r.returncode == 0 and json.loads(r.stdout)["env"].get("DEF_DAY_DUMP_HOOK") == hp
            rec(got == want_ok, f"G6 dump_hook {hp}: {'reaches env' if want_ok else 'is a problem'}")
        # ---- G7
        r = run(PY, t / "scripts/test_october_recipe.py", cwd=t)
        rec(r.returncode == 0, "G7 test_october_recipe green (pending state)")
        dfile = t / "docs/october-decisions.json"
        d7 = json.loads(dfile.read_text())
        for v in d7["switches"].values():
            if v.get("status") == "PENDING":
                v.update(status="DECIDED", value="off", evidence="test690g")
        dfile.write_text(json.dumps(d7, indent=1) + "\n")
        r = run(PY, t / "scripts/test_october_recipe.py", cwd=t)
        rec(r.returncode == 0, "G7 test_october_recipe green (all-decided state)")
    finally:
        shutil.rmtree(td, ignore_errors=True)
    k = sum(ok for ok, _ in RESULTS)
    print(f"\n{k}/{len(RESULTS)} PASS")
    return 0 if k == len(RESULTS) else 1


if __name__ == "__main__":
    sys.exit(main())
