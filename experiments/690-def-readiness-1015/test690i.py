#!/usr/bin/env python3
"""690i tests on a scratch tree (live files copied; the repo is never written).
  I1 patch --check / apply OK; a second apply is refused; decisions json differs only by + h1_stack (DECIDED x2)
  I2 bash -n def_day.sh; py_compile helpers + recipe
  I3 the DUMP SUBSHELL, run for real from the extracted def_day text with a stub driver that prints its VCC_* env:
     every combination the unpatched script supports (H1 unset/x2 x CD4 unset/x1/x1.5 x strength unset/0.25) gives the
     SAME env (VCC_H1_UNION_G, VCC_CD4_TARGET_G, VCC_BULK_DUMP, VCC_PERT_SCALE) and the same log lines; every new grid point
     gives recommend438's values and the dump_h1x<d>_cd4x<d> name
  I4 validation: the 5 x 6 grid passes, x4 / x1.1 / x0.5 die
  I5 the riders spec: unset -> no October spec, as before; CD4 x1.5 -> the same spec as before; H1 x2.5 -> spec + h1=x2.5
  I6 helpers: the old _dump_name530 answers unchanged; new names; spec with h1 parses; riders_match CLI on fake metas
  I7 recipe: live decisions -> `command` env unchanged by the patch (plain and --allow-pending); decide h1_stack x2.5 /
     cd4_stack x1.25 -> exported and installed on the patched def_day, refused (rc 3) on the unpatched one; h1_stack x4 -> rc 3
  I8 scripts/test_october_recipe.py passes on the patched scratch tree
usage: uv run python experiments/690-def-readiness-1015/test690i.py
"""
from __future__ import annotations

import importlib.util
import json
import os
import py_compile
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PY = sys.executable
FILES = ["scripts/def_day.sh", "scripts/def_day_helpers.py", "scripts/october_recipe.py", "scripts/test_october_recipe.py",
         "docs/october-decisions.json", "experiments/497-october-switches/work/def_day.497.sh",
         "scripts/make_submission_pm_cd4.py", "scripts/def_sameline.py", "src/vcc2026/rc_emit_k.py",
         "experiments/528-callprop-block/callblock528.py"]
RES: list[tuple[bool, str]] = []
H1V = {"x1": "0.28", "x1.5": "0.42", "x2": "0.56", "x2.5": "0.7", "x3": "0.84"}
C4V = {"x0.75": "0.2582741458987231", "x1": None, "x1.25": "0.4304569098312051", "x1.5": "0.5165482917974462",
       "x1.75": "0.6026396737636872", "x2": "0.6887310557299282"}
TOK = lambda d: d.replace(".", "")


def rec(ok: bool, name: str, detail: str = "") -> None:
    RES.append((bool(ok), name))
    print(f"{'PASS' if ok else 'FAIL'} {name}" + (f"  [{detail}]" if detail and not ok else ""), flush=True)


def tree(td: Path, name: str) -> Path:
    t = td / name
    for f in FILES:
        (t / f).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / f, t / f)
    return t


def cut(text: str, start: str, end: str) -> str:
    i = text.index(start)
    j = text.index(end, i)
    return text[i:j + len(end)]


def dump_env(dd: str, td: Path, env: dict) -> tuple[str, str]:
    """run def_day's own dump block (CD4/H1 case -> the subshell -> rc=$?) with a stub driver; return (env lines, log)"""
    start = 'case "${DEF_DAY_H1_STACK:-x2}" in\n  x1) H1T' if "690i: the doses" in dd else 'case "${DEF_DAY_CD4_STACK:-x1}" in\n  x1.5) DUMP'
    blk = cut(dd, start, "rc=$?\n")
    w = td / f"w{abs(hash(json.dumps(env, sort_keys=True))) % 10**9}"
    if w.exists():
        shutil.rmtree(w)
    w.mkdir()
    sh = w / "h.sh"
    sh.write_text('set -uo pipefail\nW="$1"; N=x; LODO=""\nlog(){ printf "%s\\n" "$*" >> "$W/log"; }\n'
                  "DRV=(bash -c 'env | grep -E \"^VCC_(H1_UNION_G|CD4_TARGET_G|BULK_DUMP|PERT_SCALE)=\" | sort')\n" + blk)
    e = {k: v for k, v in os.environ.items() if not k.startswith(("DEF_DAY_", "VCC_"))}
    e.update(env)
    subprocess.run(["bash", str(sh), str(w)], env=e, capture_output=True, text=True)
    out = (w / "driver-dump.log").read_text() if (w / "driver-dump.log").exists() else "<no driver log>"
    lg = (w / "log").read_text() if (w / "log").exists() else ""
    return out.replace(str(w), "$W"), lg.replace(str(w), "$W")


def run_block(script: str, env: dict) -> str:
    e = {k: v for k, v in os.environ.items() if not k.startswith(("DEF_DAY_", "VCC_"))}
    e.update(env)
    r = subprocess.run(["bash", "-c", script], env=e, capture_output=True, text=True)
    return (r.stdout + r.stderr).strip() + f" rc={r.returncode}"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def recipe(t: Path, dec: Path, *args, def_day: Path | None = None) -> tuple[int, str, str]:
    env = dict(os.environ, OCTOBER_DECISIONS=str(dec))
    if def_day:
        env["OCTOBER_DEF_DAY"] = str(def_day)
    r = subprocess.run([PY, str(t / "scripts/october_recipe.py"), *args], env=env, capture_output=True, text=True, cwd=t)
    return r.returncode, r.stdout, r.stderr


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="t690i_") as tds:
        td = Path(tds)
        o, p = tree(td, "orig"), tree(td, "patched")
        pat = lambda *a: subprocess.run([PY, str(HERE / "patch690i.py"), "--root", str(p), *a], capture_output=True, text=True)
        r1, r2, r3 = pat("--check"), pat(), pat()
        rec(r1.returncode == 0 and r2.returncode == 0 and r3.returncode == 1 and "REFUSED" in r3.stderr,
            "I1 --check / apply OK, second apply refused", r1.stderr + r2.stderr)
        jo = json.loads((o / "docs/october-decisions.json").read_text())["switches"]
        jp = json.loads((p / "docs/october-decisions.json").read_text())["switches"]
        rec(set(jp) - set(jo) == {"h1_stack"} and all(jp[k] == jo[k] for k in jo)
            and jp["h1_stack"]["status"] == "DECIDED" and jp["h1_stack"]["value"] == "x2",
            "I1 decisions: only + h1_stack (DECIDED x2)")
        ok2 = subprocess.run(["bash", "-n", str(p / "scripts/def_day.sh")]).returncode == 0
        try:
            for f in ("scripts/def_day_helpers.py", "scripts/october_recipe.py"):
                py_compile.compile(str(p / f), doraise=True)
        except py_compile.PyCompileError:
            ok2 = False
        rec(ok2, "I2 bash -n def_day.sh; py_compile helpers + recipe")

        # I3 the dump subshell
        do, dp = (o / "scripts/def_day.sh").read_text(), (p / "scripts/def_day.sh").read_text()
        same, bad = 0, []
        for h1 in (None, "x2"):
            for c4 in (None, "x1", "x1.5"):
                for st in (None, "0.25"):
                    env = {k: v for k, v in (("DEF_DAY_H1_STACK", h1), ("DEF_DAY_CD4_STACK", c4), ("DEF_DAY_LODO_STRENGTH", st)) if v}
                    a, b = dump_env(do, td, env), dump_env(dp, td, env)
                    if a == b and "VCC_H1_UNION_G=0.56" in a[0]:
                        same += 1
                    else:
                        bad.append((env, a, b))
        rec(same == 12 and not bad, f"I3 dump subshell identical to the unpatched script on all 12 shared combinations", str(bad[:1]))
        newbad = []
        for h1, hv in H1V.items():
            for c4, cv in C4V.items():
                env = {"DEF_DAY_H1_STACK": h1, "DEF_DAY_CD4_STACK": c4}
                out, lg = dump_env(dp, td, env)
                want = [f"VCC_BULK_DUMP=$W/dump_h1{TOK(h1)}_cd4{TOK(c4)}.npz"] + ([f"VCC_CD4_TARGET_G={cv}"] if cv else []) \
                    + [f"VCC_H1_UNION_G={hv}"]
                got = [ln for ln in out.splitlines() if ln]
                logs_ok = (("H1 stack" in lg) == (h1 != "x2")) and (("CD4 stack" in lg) == (cv is not None))
                if sorted(got) != sorted(want) or not logs_ok:
                    newbad.append((env, got, want, lg))
        rec(not newbad, "I3 all 30 grid points: recommend438's values, dump_h1x<d>_cd4x<d> names, log lines", str(newbad[:1]))

        # I4 validation
        vblk = cut(dp, 'case "${DEF_DAY_CD4_STACK:-x1}" in\n  x0.75|x1', "esac\ncase \"${DEF_DAY_H1_STACK:-x2}\" in   # 690i: 438's H1 grid")
        vblk = vblk.rsplit("case", 1)[0] + cut(dp, 'case "${DEF_DAY_H1_STACK:-x2}" in   # 690i: 438', "esac\n")
        harness = 'die(){ echo "DIE $*"; exit 9; }\n' + vblk + "echo VALID\n"
        okv = all(run_block(harness, {"DEF_DAY_H1_STACK": h, "DEF_DAY_CD4_STACK": c}).startswith("VALID") for h in H1V for c in C4V)
        okv &= run_block(harness, {}).startswith("VALID")
        okv &= all("DIE" in run_block(harness, e) for e in ({"DEF_DAY_H1_STACK": "x4"}, {"DEF_DAY_CD4_STACK": "x1.1"},
                                                            {"DEF_DAY_H1_STACK": "x0.5"}, {"DEF_DAY_CD4_STACK": "x3"}))
        rec(okv, "I4 validation: 30 grid points + unset pass; x4 / x1.1 / x0.5 / CD4 x3 die")

        # I5 the October riders spec
        def spec(dd: str, env: dict) -> str:
            blk = cut(dd, '_r530=""\n', '  [ -n "${DEF_DAY_CALLSET_BLOCK:-}" ] && _spec="${_spec},callset=${DEF_DAY_CALLSET_BLOCK}"\n')
            if "690i" in dd:
                blk += '  [ "${DEF_DAY_H1_STACK:-x2}" != "x2" ] && _spec="${_spec},h1=${DEF_DAY_H1_STACK}"   # 690i\n'
            return run_block(blk + '  echo "SPEC $_spec"\nelse\n  echo NOSPEC\nfi\n', env)
        s_ok = (spec(do, {}) == spec(dp, {}) == "NOSPEC rc=0"
                and spec(do, {"DEF_DAY_CD4_STACK": "x1.5"}) == spec(dp, {"DEF_DAY_CD4_STACK": "x1.5"})
                and spec(dp, {"DEF_DAY_H1_STACK": "x2.5"}) == "SPEC rc=-,stack=x1,rho=-,strength=-,h1=x2.5 rc=0"
                and spec(dp, {"DEF_DAY_H1_STACK": "x2"}) == "NOSPEC rc=0")
        rec(s_ok, "I5 riders spec: unset / CD4 x1.5 unchanged; H1 x2.5 adds h1=x2.5; H1 x2 stays unset",
            f"{spec(dp, {})} | {spec(dp, {'DEF_DAY_H1_STACK': 'x2.5'})}")

        # I6 helpers
        ho, hp = load(o / "scripts/def_day_helpers.py", "h_orig"), load(p / "scripts/def_day_helpers.py", "h_pat")
        old_cases = [("x1", "-"), ("x1.5", "-"), ("x1.5", "0.25"), ("x1", "-0.25"), ("x1", "0.25")]
        h_ok = all(ho._dump_name530(a, b) == hp._dump_name530(a, b) for a, b in old_cases)
        h_ok &= hp._dump_name530("x0.75", "-", "x2.5") == "dump_h1x25_cd4x075.npz"
        h_ok &= hp._dump_name530("x1.5", "0.25", "x3") == "dump_h1x3_cd4x15_s0.25.npz"
        sp = hp._october530_spec("rc=-,stack=x1,rho=-,strength=-,h1=x2.5")
        h_ok &= sp is not None and sp["h1"] == "x2.5" and ho._october530_spec("rc=-,stack=x1,rho=-,strength=-,h1=x2.5") is None
        cfg = {"bulk_load": "dump_h1x25_cd4x1.npz", "cd4_target_g": hp.X1_G, "rc_callset_cpm": 2200}
        h_ok &= hp._october530_mismatch(cfg, sp) is None
        h_ok &= "bulk_load" in (hp._october530_mismatch(dict(cfg, bulk_load="dump_h1x2_cd4x1.npz"), sp) or "")
        h_ok &= "h1" in (hp._october530_mismatch(cfg, dict(sp, h1="x4")) or "")
        old_spec = ho._october530_spec("rc=-,stack=x1.5,rho=-,strength=0.25")
        c15 = {"bulk_load": "dump_h1x2_cd4x15_s0.25.npz", "cd4_target_g": ho.X1_G, "rc_callset_cpm": 2200}
        h_ok &= ho._october530_mismatch(c15, old_spec) == hp._october530_mismatch(c15, hp._october530_spec("rc=-,stack=x1.5,rho=-,strength=0.25")) is None
        rec(h_ok, "I6 helpers: old dump names unchanged, new names, spec h1 parses and is checked")
        mt = td / "meta.json"
        mt.write_text(json.dumps({"config": cfg}))
        cli = lambda s: subprocess.run([PY, str(p / "scripts/def_day_helpers.py"), "riders_match", str(mt), "-", "-", s],
                                       capture_output=True, text=True).returncode
        rec(cli("rc=-,stack=x1,rho=-,strength=-,h1=x2.5") == 0 and cli("rc=-,stack=x1,rho=-,strength=-,h1=x3") != 0
            and cli("rc=-,stack=x1,rho=-,strength=-") != 0, "I6 riders_match CLI: right h1 rc 0, wrong / missing h1 refused")

        # I7 recipe
        live_o, live_p = o / "docs/october-decisions.json", p / "docs/october-decisions.json"
        envs = []
        for extra in ([], ["--allow-pending"]):
            a = recipe(o, live_o, "command", "--bundle", "B", "--work", "W", "--name", "N", "--json", *extra)
            b = recipe(p, live_p, "command", "--bundle", "B", "--work", "W", "--name", "N", "--json", *extra)
            envs.append((a[0], b[0], json.loads(a[1])["env"] if a[0] == 0 else a[2][-200:],
                         json.loads(b[1])["env"] if b[0] == 0 else b[2][-200:]))
        rec(all(x[0] == x[1] and x[2] == x[3] for x in envs), "I7 recipe: live decisions -> identical env (plain, --allow-pending)",
            str(envs))
        dec2 = td / "d2.json"
        shutil.copy(live_p, dec2)
        d1 = recipe(p, dec2, "decide", "h1_stack", "x2.5", "--evidence", "test690i")
        d2 = recipe(p, dec2, "decide", "cd4_stack", "x1.25", "--evidence", "test690i")
        a = recipe(p, dec2, "command", "--bundle", "B", "--work", "W", "--name", "N", "--json", "--allow-pending")
        b = recipe(p, dec2, "command", "--bundle", "B", "--work", "W", "--name", "N", "--json", "--allow-pending",
                   def_day=o / "scripts/def_day.sh")
        e7 = json.loads(a[1])["env"] if a[0] == 0 else {}
        rec(d1[0] == 0 and d2[0] == 0 and a[0] == 0 and e7.get("DEF_DAY_H1_STACK") == "x2.5" and e7.get("DEF_DAY_CD4_STACK") == "x1.25"
            and b[0] == 3 and "DEF_DAY_H1_STACK" in b[2] and "predates 690i" in b[2],
            "I7 decide h1 x2.5 + cd4 x1.25: exported on the patched def_day; rc 3 on the unpatched one", f"{a[2][-300:]} | {b[2][-300:]}")
        recipe(p, dec2, "decide", "h1_stack", "x4", "--evidence", "test690i")
        c = recipe(p, dec2, "command", "--bundle", "B", "--work", "W", "--name", "N", "--json", "--allow-pending")
        rec(c[0] == 3 and "h1_stack 'x4'" in c[2], "I7 h1_stack x4 -> rc 3", c[2][-200:])

        # I8 the recipe's own test suite on the patched tree
        r8 = subprocess.run([PY, "-m", "pytest", "-q", str(p / "scripts/test_october_recipe.py")], capture_output=True, text=True, cwd=p)
        if "No module named pytest" in r8.stderr:
            r8 = subprocess.run([PY, str(p / "scripts/test_october_recipe.py")], capture_output=True, text=True, cwd=p)
        rec(r8.returncode == 0, "I8 scripts/test_october_recipe.py on the patched tree", (r8.stdout + r8.stderr)[-500:])
    n = sum(ok for ok, _ in RES)
    print(f"\n{n}/{len(RES)} PASS")
    return 0 if n == len(RES) else 1


if __name__ == "__main__":
    sys.exit(main())
