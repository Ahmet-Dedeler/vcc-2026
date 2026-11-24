#!/usr/bin/env python3
"""690a tests. Light: no data, no build. Each tree is a scratch copy of the live files 690a and october_recipe touch.

  T1  --check, apply, a second apply refused; refusals with a pre-631 builder and with a def_day.sh lacking 630dd
  T2  bash -n on the patched def_day.sh
  T3  KA "today's command unchanged": `command --allow-pending --json` env and args identical before and after; the new
      switches only extend the PENDING list; plain `command` refuses with them listed
  T4  the inserted def_day.sh text, cut from the patched file and run in a bash harness with stubbed die/log:
      validation accepts the builder grammar (m in [1.25, 4]) and refuses everything else, and refuses without
      DEF_DAY_DE_MAG; unset knob -> no export, no final.env line, no ALLOW change; set -> all three
  T5  the final-meta check passes only with the requested variant AND the reviewed module sha
  T6  recipe: decided values reach env; dependency and grammar problems fire; installed_problems flags a builder
      without VCC_DE_DNMAG
  T7  scripts/test_october_recipe.py passes on the patched tree
usage: uv run python experiments/690-def-readiness-1015/test690a.py
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PY = sys.executable
PATCH = HERE / "patch690a.py"
sys.path.insert(0, str(HERE))
import patch690a as P  # noqa: E402

FILES = ["scripts/def_day.sh", "scripts/october_recipe.py", "docs/october-decisions.json",
         "scripts/make_submission_pm_cd4.py", "scripts/def_sameline.py", "src/vcc2026/rc_emit_k.py",
         "experiments/528-callprop-block/callblock528.py", "scripts/test_october_recipe.py",
         "experiments/497-october-switches/work/def_day.497.sh"]
PRE631 = ROOT / "experiments/671-dnmag-arm/work/builder.pre631.py"
RESULTS: list[tuple[bool, str]] = []


def rec(ok: bool, name: str) -> None:
    RESULTS.append((bool(ok), name))
    print(f"{'PASS' if ok else 'FAIL'} {name}", flush=True)


def tree(td: Path, name: str) -> Path:
    t = td / name
    for f in FILES:
        (t / f).parent.mkdir(parents=True, exist_ok=True)
        if (ROOT / f).exists():
            shutil.copy2(ROOT / f, t / f)
    return t


def patch(t: Path, *extra: str) -> subprocess.CompletedProcess:
    return subprocess.run([PY, str(PATCH), "--root", str(t), *extra], capture_output=True, text=True)


def recipe(t: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run([PY, str(t / "scripts/october_recipe.py"), *args], capture_output=True, text=True, cwd=t)


def cut(text: str, start: str, end: str) -> str:
    i = text.index(start)
    j = text.index(end, i) + len(end)
    return text[i:j]


def bash(snippet: str, env: dict) -> subprocess.CompletedProcess:
    harness = 'die(){ echo "DIE: $*"; exit 7; }\nlog(){ echo "LOG: $*"; }\n' + snippet
    e = {k: v for k, v in os.environ.items() if not k.startswith(("DEF_DAY_", "VCC_"))}
    e.update(env)
    return subprocess.run(["bash", "-c", harness], capture_output=True, text=True, env=e)


def main() -> int:
    td = Path(tempfile.mkdtemp(prefix="t690a_"))
    try:
        # ---- T1
        t = tree(td, "a")
        r = patch(t, "--check")
        rec(r.returncode == 0 and "CHECK OK" in r.stdout, "T1 --check on the live files")
        before_dd = (t / "scripts/def_day.sh").read_text()
        r = patch(t)
        rec(r.returncode == 0 and "APPLIED" in r.stdout, "T1 apply")
        r2 = patch(t)
        rec(r2.returncode == 1 and "already applied" in r2.stderr, "T1 a second apply is refused")
        b = tree(td, "b")
        shutil.copy2(PRE631, b / "scripts/make_submission_pm_cd4.py")
        r = patch(b, "--check")
        rec(r.returncode == 1 and "631" in r.stderr, "T1 refused on a pre-631 builder")
        c = tree(td, "c")
        p = c / "scripts/def_day.sh"
        p.write_text(p.read_text().replace("630dd: rank-neutral UP pass-through", "630dd-removed"))
        r = patch(c, "--check")
        rec(r.returncode == 1 and "630dd" in r.stderr, "T1 refused on a def_day.sh without 630dd")

        # ---- T2
        r = subprocess.run(["bash", "-n", str(t / "scripts/def_day.sh")], capture_output=True, text=True)
        rec(r.returncode == 0, "T2 bash -n on the patched def_day.sh")
        dd = (t / "scripts/def_day.sh").read_text()
        import difflib  # noqa: PLC0415
        diff = list(difflib.unified_diff(before_dd.splitlines(), dd.splitlines(), lineterm="", n=0))
        added = [ln[1:] for ln in diff if ln.startswith("+") and not ln.startswith("+++")]
        removed = [ln for ln in diff if ln.startswith("-") and not ln.startswith("---")]
        rec(not removed and len(added) == 19 and all(("DNMAG" in ln or "dnmag" in ln or "690a" in ln or ln.strip() in
                 {"fi", "PY", "import json, sys"} or "json.load(open(sys.argv[1]))" in ln or "sys.exit(0 if u.get" in ln
                 or ln.startswith("#    pass after zs")) for ln in added),
            f"T2 the patch only adds DNMAG lines ({len(added)} added, {len(removed)} removed)")

        # ---- T3
        u = tree(td, "u")
        args = ["command", "--bundle", "data/challenge", "--work", "W", "--name", "N", "--allow-pending", "--json"]
        r0 = recipe(u, *args)
        r1 = recipe(t, *args)
        j0, j1 = json.loads(r0.stdout), json.loads(r1.stdout)
        rec(r0.returncode == 0 and r1.returncode == 0 and j0["env"] == j1["env"] and j0["args"] == j1["args"],
            "T3 KA: the resolved env and args are unchanged with the new switches PENDING")
        rec(set(j1["pending_treated_off"]) - set(j0["pending_treated_off"]) == {"de_uprn", "mde_sign_keygate", "de_dnmag"},
            "T3 the new switches only extend the PENDING list")
        r = recipe(t, "command", "--bundle", "data/challenge", "--work", "W", "--name", "N")
        rec(r.returncode != 0 and "de_dnmag" in (r.stdout + r.stderr), "T3 plain command refuses with them listed")
        r = recipe(t, "show")
        rec(r.returncode == 0 and all(k in r.stdout for k in ("de_uprn", "mde_sign_keygate", "de_dnmag")),
            "T3 show lists the three switches")

        # ---- T4
        v = cut(dd, "# 690a: DOWN magnitude pass-through -- DEF_DAY_DE_DNMAG", "(631)\"\nfi\n")
        os.chdir(t)  # the snippet greps scripts/make_submission_pm_cd4.py relative to the repo
        ok_vals = ["d_m125", "d_m200", "d_m400", "d_m200_q5", "d_m150_q12"]
        bad_vals = ["d_m124", "d_m401", "d_m020", "d_m2000", "d_m200_q0", "d_m200_q01", "d_m200_", "dm200", "d_m200q5",
                    "D_M200", "d_m200;id", " d_m200"]
        okv = [bash(v, {"DEF_DAY_DE_MAG": "dn_k250_zs", "DEF_DAY_DE_DNMAG": x}).returncode == 0 for x in ok_vals]
        badv = [bash(v, {"DEF_DAY_DE_MAG": "dn_k250_zs", "DEF_DAY_DE_DNMAG": x}).returncode == 7 for x in bad_vals]
        rec(all(okv), f"T4 validation accepts {ok_vals}")
        rec(all(badv), f"T4 validation refuses {bad_vals}")
        r = bash(v, {"DEF_DAY_DE_DNMAG": "d_m200"})
        rec(r.returncode == 7 and "needs DEF_DAY_DE_MAG" in r.stdout, "T4 refused without DEF_DAY_DE_MAG")
        r = bash(v, {})
        rec(r.returncode == 0 and r.stdout == "", "T4 unset: validation is silent")
        e_line = next(ln for ln in dd.splitlines() if ln.startswith('[ -n "${DEF_DAY_DE_DNMAG:-}" ] && { export'))
        f_line = next(ln for ln in dd.splitlines() if ln.startswith('[ -n "${VCC_DE_DNMAG:-}" ] && printf'))
        a_line = next(ln for ln in dd.splitlines() if ln.startswith('[ -n "${VCC_DE_DNMAG:-}" ] && ALLOW'))
        efa = (f'W="$TMPW"; ALLOW=base\n{e_line}\n{f_line}\n{a_line}\n'
               'echo "ENV=${VCC_DE_DNMAG:-unset}"; echo "ALLOW=$ALLOW"; cat "$W/def_day.final.env" 2>/dev/null || echo NOENV\n')
        wdir = td / "w_unset"
        wdir.mkdir()
        r = bash(efa, {"TMPW": str(wdir)})
        rec(r.returncode == 0 and "ENV=unset" in r.stdout and "ALLOW=base\n" in r.stdout and "NOENV" in r.stdout,
            "T4 unset: no export, no ALLOW change, no final.env line")
        wdir = td / "w_set"
        wdir.mkdir()
        r = bash(efa, {"TMPW": str(wdir), "DEF_DAY_DE_DNMAG": "d_m200"})
        rec(r.returncode == 0 and "ENV=d_m200" in r.stdout and "ALLOW=base,de_dnmag" in r.stdout
            and "VCC_DE_DNMAG=d_m200" in r.stdout, "T4 set: export, ALLOW += de_dnmag, final.env line")

        # ---- T5
        cblk = cut(dd, 'if [ -n "${VCC_DE_DNMAG:-}" ]; then   # 690a', "PY\nfi\n")
        def meta(cfg: dict) -> Path:
            m = td / f"meta_{len(list(td.glob('meta_*')))}.json"
            m.write_text(json.dumps({"config": cfg}))
            return m
        good = meta({"de_dnmag": {"variant": "d_m200", "module_sha256": P.DNMAG_MODULE_SHA256}})
        wrongsha = meta({"de_dnmag": {"variant": "d_m200", "module_sha256": "0" * 64}})
        wrongvar = meta({"de_dnmag": {"variant": "d_m150", "module_sha256": P.DNMAG_MODULE_SHA256}})
        missing = meta({})
        res = {}
        for name, m in (("good", good), ("wrongsha", wrongsha), ("wrongvar", wrongvar), ("missing", missing)):
            fin = str(m)[:-len(".h5ad.meta.json")] if str(m).endswith(".h5ad.meta.json") else None
            # the block reads "$FINAL.h5ad.meta.json": give it a FINAL whose meta is m
            link = td / f"F_{name}"
            os.link(m, str(link) + ".h5ad.meta.json")
            res[name] = bash(cblk, {"FINAL": str(link), "VCC_DE_DNMAG": "d_m200"}).returncode
        rec(res == {"good": 0, "wrongsha": 7, "wrongvar": 7, "missing": 7}, f"T5 meta check {res}")
        r = bash(cblk, {"FINAL": str(td / "nonexistent")})
        rec(r.returncode == 0, "T5 unset: no meta check")

        # ---- T6
        dj = json.loads((t / "docs/october-decisions.json").read_text())
        def resolve_with(updates: dict, builder: Path | None = None) -> dict:
            d6 = json.loads(json.dumps(dj))
            for k, v in updates.items():
                d6["switches"][k].update({"status": "DECIDED", "value": v})
            f = td / f"d6_{len(list(td.glob('d6_*')))}.json"
            f.write_text(json.dumps(d6))
            env = dict(os.environ, OCTOBER_DECISIONS=str(f))
            if builder:
                env["OCTOBER_BUILDER"] = str(builder)
            r = subprocess.run([PY, str(t / "scripts/october_recipe.py"), "command", "--bundle", "B", "--work", "W",
                                "--name", "N", "--allow-pending", "--json"], capture_output=True, text=True, env=env, cwd=t)
            try:
                return {"rc": r.returncode, **json.loads(r.stdout)}
            except json.JSONDecodeError:
                return {"rc": r.returncode, "out": r.stdout + r.stderr}
        o = resolve_with({"de_dnmag": "d_m200", "de_uprn": "b_k100_q20r3p10", "mde_sign_keygate": "0.5"})
        rec(o["rc"] == 0 and o.get("env", {}).get("DEF_DAY_DE_DNMAG") == "d_m200"
            and o["env"].get("DEF_DAY_DE_UPRN") == "b_k100_q20r3p10" and o["env"].get("DEF_DAY_MDE_SIGN_KEYGATE") == "0.5",
            "T6 decided values reach the env (de_mag and mde_sign already DECIDED)")
        o = resolve_with({"de_dnmag": "d_m200", "de_mag": "off"})
        rec(o["rc"] != 0 and "de_dnmag runs after the zs pass" in json.dumps(o), "T6 de_dnmag without de_mag is a problem")
        o = resolve_with({"de_uprn": "b_k100", "de_mag": "off"})
        rec(o["rc"] != 0 and "de_uprn runs after the zs pass" in json.dumps(o), "T6 de_uprn without de_mag is a problem")
        o = resolve_with({"mde_sign_keygate": "0.5", "mde_sign": "off"})
        rec(o["rc"] != 0 and "needs mde_sign" in json.dumps(o), "T6 keygate without mde_sign is a problem")
        o = resolve_with({"de_dnmag": "d_m500"})
        rec(o["rc"] != 0 and "must be d_m<m*100>" in json.dumps(o), "T6 de_dnmag grammar")
        o = resolve_with({"mde_sign_keygate": "1.5"})
        rec(o["rc"] != 0 and "quantile" in json.dumps(o), "T6 keygate range")
        o = resolve_with({"de_dnmag": "d_m200"}, builder=PRE631)
        rec(o["rc"] != 0 and "VCC_DE_DNMAG" in json.dumps(o), "T6 installed_problems: a builder without VCC_DE_DNMAG")

        # ---- T7
        r = subprocess.run([PY, str(t / "scripts/test_october_recipe.py")], capture_output=True, text=True, cwd=t)
        tail = (r.stdout + r.stderr).strip().splitlines()[-3:]
        rec(r.returncode == 0, f"T7 test_october_recipe on the patched tree: {' | '.join(tail)}")
        # the all-decided state (Oct 21): decide every PENDING switch off in the scratch tree; the suite stays green
        dfile = t / "docs/october-decisions.json"
        keep = dfile.read_text()
        d7 = json.loads(keep)
        for k, v in d7["switches"].items():
            if v.get("status") == "PENDING":
                v.update({"status": "DECIDED", "value": "off", "evidence": "test690a all-decided state"})
        dfile.write_text(json.dumps(d7, indent=1) + "\n")
        r = subprocess.run([PY, str(t / "scripts/test_october_recipe.py")], capture_output=True, text=True, cwd=t)
        tail = (r.stdout + r.stderr).strip().splitlines()[-3:]
        rec(r.returncode == 0, f"T7 test_october_recipe in the all-decided state: {' | '.join(tail)}")
        dfile.write_text(keep)
    finally:
        os.chdir(ROOT)
        shutil.rmtree(td, ignore_errors=True)
    n_ok = sum(ok for ok, _ in RESULTS)
    print(f"\n{n_ok}/{len(RESULTS)} PASS")
    return 0 if n_ok == len(RESULTS) else 1


if __name__ == "__main__":
    sys.exit(main())
