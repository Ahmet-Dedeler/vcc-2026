#!/usr/bin/env python3
"""673 tests: gate673 on 695's REAL finished R-ABC2 run (the reference W1 is judged against). The mirrored run holds a copy of
R-ABC2's stack dump renamed to 690i's x2.5 name, its env's VCC_H1_UNION_G set per case; def_day.log gets the H1 stack line.
  B1 the recipe adds DEF_DAY_H1_STACK=x2.5, the meta names the x2.5 dump (moved sha), G 0.7 -> ALL GATES PASS, adopted == [H1]
  B2 the export is there but the meta still names R-ABC2's dump unchanged -> FAIL
  B3 the meta names the x2.5 dump but its env records G 0.56 (the dose did not act) -> FAIL
  B4 no export (R-ABC2's recipe) but the meta moved -> FAIL
  B5 the export + an unadopted key (de_mag) moved -> FAIL
usage: uv run python experiments/673-h1x25-abc/test673.py
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
R2_RUN = ROOT / "experiments/695-rabc2/work/run"
R2_META = ROOT / "submissions/695rabc2-october-recipe-abc.h5ad.meta.json"
R2_RECIPE = ROOT / "experiments/695-rabc2/work/recipe695.json"
D2, D25 = "dump_h1x2_cd4x1_s0.25.npz", "dump_h1x25_cd4x1_s0.25.npz"
RES: list[tuple[bool, str]] = []


def rec(ok, name, detail=""):
    RES.append((bool(ok), name))
    print(f"{'PASS' if ok else 'FAIL'} {name}" + (f"  [{detail}]" if detail and not ok else ""), flush=True)


def case(td: Path, tag: str, h1: str | None, g: str | None, meta_edit) -> str:
    w = td / tag
    run = w / "run"
    run.mkdir(parents=True)
    for p in R2_RUN.iterdir():
        if p.name == "def_day.log":
            t = p.read_text()
            if h1:
                t += f"2026-10-20T00:00:00Z H1 stack {h1}: dump {D25} VCC_H1_UNION_G={g} (level builds keep DEF_ENV)\n"
            (run / p.name).write_text(t)
        else:
            (run / p.name).symlink_to(p)
    if g is not None:   # the renamed dump: R-ABC2's arrays, env G per case, its own run dir in the env string
        z = dict(np.load(R2_RUN / D2, allow_pickle=False))
        e = json.loads(str(z["env"])); e["VCC_H1_UNION_G"] = g
        e["VCC_BULK_DUMP"] = str(run / D25)
        z["env"] = np.array(json.dumps(e))
        z["bulk"] = z["bulk"].copy(); z["bulk"].flat[0] += 1.0
        np.savez(run / D25, **z)
    base = json.loads(R2_RECIPE.read_text())
    env = dict(base["env"], **({"DEF_DAY_H1_STACK": h1} if h1 else {}))
    (w / "recipe673.json").write_text(json.dumps({"env": env, "args": base["args"], "pending_treated_off": [], "problems": []}))
    m = json.loads(R2_META.read_text())
    meta_edit(m["config"], run)
    stem = td / f"{tag}-stem"
    Path(str(stem) + ".h5ad.meta.json").write_text(json.dumps(m))
    Path(str(stem) + ".vcc").write_bytes(b"pkg")
    cl = w / "w1_673.log"
    cl.write_text(f"2026-10-20T00:00:00Z package {stem}.vcc sha256 {hashlib.sha256(b'pkg').hexdigest()}\n2026-10-20T00:00:01Z novelty rc=0\n")
    p = subprocess.run([sys.executable, str(HERE / "gate673.py"), "--work", str(w), "--stem", str(stem), "--chainlog", str(cl)],
                       capture_output=True, text=True, cwd=ROOT)
    return p.stdout + p.stderr


def renamed(c, run):
    c["bulk_load"] = D25
    c["bulk_load_sha256"] = hashlib.sha256((run / D25).read_bytes()).hexdigest()


def main() -> int:
    if not (R2_RUN.is_dir() and R2_META.is_file() and R2_RECIPE.is_file()):
        print("695's run / meta / recipe not on disk"); return 2
    with tempfile.TemporaryDirectory(prefix="t673_") as t:
        td = Path(t)
        o = case(td, "b1", "x2.5", "0.7", renamed)
        line = next((ln for ln in o.splitlines() if ln.startswith("[w1] PASS: config diff vs R-ABC2")), "")
        rec(bool(line) and "adopted ['DEF_DAY_H1_STACK']" in line and "ALL GATES PASS" in o, "B1 H1 x2.5, renamed dump, G 0.7 -> ALL GATES PASS", o[-1400:])
        o = case(td, "b2", "x2.5", None, lambda c, r: None)
        rec("ALL GATES FAIL" in o, "B2 export but the meta unchanged -> FAIL", o[-600:])
        o = case(td, "b3", "x2.5", "0.56", renamed)
        rec("FAIL H1 stack dose (690i)" in o and "ALL GATES FAIL" in o, "B3 renamed dump but G 0.56 (dose did not act) -> FAIL", o[-600:])
        o = case(td, "b4", None, "0.7", renamed)
        rec("[w1] FAIL: config diff vs R-ABC2 (adopted [])" in o and "ALL GATES FAIL" in o, "B4 no export but the meta moved -> FAIL", o[-600:])
        def mag(c, r):
            renamed(c, r)
            c["de_mag"] = dict(c["de_mag"], variant="other")
        o = case(td, "b5", "x2.5", "0.7", mag)
        rec("[w1] FAIL: config diff vs R-ABC2" in o and "ALL GATES FAIL" in o, "B5 export + de_mag moved -> FAIL", o[-600:])
    n = sum(ok for ok, _ in RES)
    print(f"\n{n}/{len(RES)} PASS")
    return 0 if n == len(RES) else 1


if __name__ == "__main__":
    sys.exit(main())
