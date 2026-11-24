#!/usr/bin/env python3
"""674 tests: gate674 against a W1 STAND-IN made from 695's real R-ABC2 run (W1 itself is not built yet): reference recipe =
R-ABC2's + DEF_DAY_H1_STACK=x2.5, reference run = 695's mirrored with a renamed x2.5 dump (G 0.7). The W2 case adds
DEF_DAY_CD4_STACK=x1.5 and its renamed dump dump_h1x25_cd4x15_s0.25.npz (G 0.7, CD4 G per case). Overrides GATE674_REF_*.
  B1 + CD4 x1.5, the meta names the x2.5/x1.5 dump, CD4 G 0.5165... -> ALL GATES PASS, adopted == [CD4]
  B2 the CD4 export but the meta still names W1's dump -> FAIL
  B3 the renamed dump records CD4 G 0.3443... (the CD4 dose did not act) -> FAIL
  B4 no CD4 export (W1's recipe) but the meta moved -> FAIL
  B5 the CD4 export + an unadopted key (de_mag) moved -> FAIL
usage: uv run python experiments/674-h1x25cd4x15-abc/test674.py
"""
from __future__ import annotations

import hashlib
import json
import os
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
D2, D25, D2515 = "dump_h1x2_cd4x1_s0.25.npz", "dump_h1x25_cd4x1_s0.25.npz", "dump_h1x25_cd4x15_s0.25.npz"
C1, C15 = "0.3443655278649641", "0.5165482917974462"
RES: list[tuple[bool, str]] = []


def rec(ok, name, detail=""):
    RES.append((bool(ok), name))
    print(f"{'PASS' if ok else 'FAIL'} {name}" + (f"  [{detail}]" if detail and not ok else ""), flush=True)


def mirror(run: Path, extra_log: str):
    run.mkdir(parents=True)
    for p in R2_RUN.iterdir():
        if p.name == "def_day.log":
            (run / p.name).write_text(p.read_text() + extra_log)
        else:
            (run / p.name).symlink_to(p)


def put_dump(run: Path, name: str, g: str, cg: str, bump: float) -> str:
    z = dict(np.load(R2_RUN / D2, allow_pickle=False))
    e = json.loads(str(z["env"])); e["VCC_H1_UNION_G"] = g; e["VCC_CD4_TARGET_G"] = cg; e["VCC_BULK_DUMP"] = str(run / name)
    z["env"] = np.array(json.dumps(e)); z["bulk"] = z["bulk"].copy(); z["bulk"].flat[0] += bump
    np.savez(run / name, **z)
    return hashlib.sha256((run / name).read_bytes()).hexdigest()


def w1_standin(td: Path):
    run = td / "w1" / "run"
    mirror(run, f"2026-10-20T00:00:00Z H1 stack x2.5: dump {D25} VCC_H1_UNION_G=0.7 (level builds keep DEF_ENV)\n")
    sha = put_dump(run, D25, "0.7", C1, 1.0)
    base = json.loads(R2_RECIPE.read_text())
    (td / "w1" / "recipe673.json").write_text(json.dumps({"env": dict(base["env"], DEF_DAY_H1_STACK="x2.5"), "args": base["args"]}))
    m = json.loads(R2_META.read_text()); m["config"]["bulk_load"] = D25; m["config"]["bulk_load_sha256"] = sha
    (td / "w1" / "meta.json").write_text(json.dumps(m))
    return run, td / "w1" / "meta.json", td / "w1" / "recipe673.json", m


def case(td: Path, tag: str, ref, cd4: bool, cg: str | None, meta_edit) -> str:
    ref_run, ref_meta, ref_recipe, w1m = ref
    w = td / tag
    run = w / "run"
    mirror(run, f"2026-10-20T00:00:00Z H1 stack x2.5: dump {D2515} VCC_H1_UNION_G=0.7 (level builds keep DEF_ENV)\n"
                + (f"2026-10-20T00:00:00Z CD4 stack x1.5: dump {D2515} VCC_CD4_TARGET_G={C15} (level builds keep DEF_ENV)\n" if cd4 else ""))
    (run / D25).symlink_to(ref_run / D25)   # W1's dump is reachable by name from this run too
    sha = put_dump(run, D2515, "0.7", cg, 2.0) if cg else None
    env = json.loads(ref_recipe.read_text())["env"]
    env = dict(env, **({"DEF_DAY_CD4_STACK": "x1.5"} if cd4 else {}))
    base = json.loads(R2_RECIPE.read_text())
    (w / "recipe674.json").write_text(json.dumps({"env": env, "args": base["args"], "pending_treated_off": [], "problems": []}))
    m = json.loads(json.dumps(w1m))
    meta_edit(m["config"], sha)
    stem = td / f"{tag}-stem"
    Path(str(stem) + ".h5ad.meta.json").write_text(json.dumps(m))
    Path(str(stem) + ".vcc").write_bytes(b"pkg")
    cl = w / "w2_674.log"
    cl.write_text(f"2026-10-20T00:00:00Z package {stem}.vcc sha256 {hashlib.sha256(b'pkg').hexdigest()}\n2026-10-20T00:00:01Z novelty rc=0\n")
    e = dict(os.environ, GATE674_REF_META=str(ref_meta), GATE674_REF_RUN=str(ref_run), GATE674_REF_RECIPE=str(ref_recipe))
    p = subprocess.run([sys.executable, str(HERE / "gate674.py"), "--work", str(w), "--stem", str(stem), "--chainlog", str(cl)],
                       capture_output=True, text=True, cwd=ROOT, env=e)
    return p.stdout + p.stderr


def renamed(c, sha):
    c["bulk_load"] = D2515; c["bulk_load_sha256"] = sha


def main() -> int:
    if not (R2_RUN.is_dir() and R2_META.is_file() and R2_RECIPE.is_file()):
        print("695's run / meta / recipe not on disk"); return 2
    with tempfile.TemporaryDirectory(prefix="t674_") as t:
        td = Path(t); ref = w1_standin(td)
        o = case(td, "b1", ref, True, C15, renamed)
        line = next((ln for ln in o.splitlines() if ln.startswith("[w2] PASS: config diff vs W1")), "")
        rec(bool(line) and "adopted ['DEF_DAY_CD4_STACK']" in line and "ALL GATES PASS" in o, "B1 + CD4 x1.5, renamed dump, CD4 G x1.5 -> ALL GATES PASS", o[-1500:])
        o = case(td, "b2", ref, True, None, lambda c, s: None)
        rec("ALL GATES FAIL" in o, "B2 CD4 export but the meta still names W1's dump -> FAIL", o[-600:])
        o = case(td, "b3", ref, True, C1, renamed)
        rec("FAIL CD4 stack dose (690i)" in o and "ALL GATES FAIL" in o, "B3 renamed dump but CD4 G x1 (dose did not act) -> FAIL", o[-600:])
        o = case(td, "b4", ref, False, C15, renamed)
        rec("[w2] FAIL: config diff vs W1 (adopted [])" in o and "ALL GATES FAIL" in o, "B4 no CD4 export but the meta moved -> FAIL", o[-600:])
        def mag(c, s):
            renamed(c, s); c["de_mag"] = dict(c["de_mag"], variant="other")
        o = case(td, "b5", ref, True, C15, mag)
        rec("[w2] FAIL: config diff vs W1" in o and "ALL GATES FAIL" in o, "B5 CD4 export + de_mag moved -> FAIL", o[-600:])
    n = sum(ok for ok, _ in RES)
    print(f"\n{n}/{len(RES)} PASS")
    return 0 if n == len(RES) else 1


if __name__ == "__main__":
    sys.exit(main())
