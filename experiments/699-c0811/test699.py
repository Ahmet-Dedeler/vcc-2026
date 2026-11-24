#!/usr/bin/env python3
"""699 tests: gate699's reference logic on 695's REAL finished R-ABC2 run (the reference C0 is judged against).
A test run dir mirrors 695's run (symlinks), except def_day.log, which gets 690g's two hook lines appended, as C0's will have.
  B1 the recipe adds only the hook (dumphook811c0.sh) and the meta moves bulk_load_sha256 -> ALL GATES PASS, adopted == [hook]
  B2 the hook is added but the meta's dump is unmoved -> FAIL (a hook that changes nothing is not C0)
  B3 the recipe is R-ABC2's exactly (no hook) but the meta moves -> FAIL (adopted is not [hook])
  B4 the hook is added and an unadopted key (de_mag) moves too -> FAIL
usage: uv run python experiments/699-c0811/test699.py
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
R2_RUN = ROOT / "experiments/695-rabc2/work/run"
R2_META = ROOT / "submissions/695rabc2-october-recipe-abc.h5ad.meta.json"
R2_RECIPE = ROOT / "experiments/695-rabc2/work/recipe695.json"
HC0 = "experiments/811-lolo-model/dumphook811c0.sh"
RES: list[tuple[bool, str]] = []


def rec(ok, name, detail=""):
    RES.append((bool(ok), name))
    print(f"{'PASS' if ok else 'FAIL'} {name}" + (f"  [{detail}]" if detail and not ok else ""), flush=True)


def case(td: Path, tag: str, hook: str | None, meta_edit) -> str:
    w = td / tag
    run = w / "run"
    run.mkdir(parents=True)
    for p in R2_RUN.iterdir():   # mirror 695's run; def_day.log gets the hook lines when the case adds a hook
        if p.name == "def_day.log":
            t = p.read_text()
            if hook:
                t += (f"2026-10-20T00:00:00Z 690g: dump hook {hook} (0000000000000000) on dump_h1x2_cd4x1_s0.25.npz\n"
                      "2026-10-20T00:00:01Z 690g: dump replaced by the hook's 0000000000000000 (the stack dump stays at run)\n")
            (run / p.name).write_text(t)
        else:
            (run / p.name).symlink_to(p)
    base = json.loads(R2_RECIPE.read_text())
    env = dict(base["env"], **({"DEF_DAY_DUMP_HOOK": hook} if hook else {}))
    (w / "recipe699.json").write_text(json.dumps({"env": env, "args": base["args"], "pending_treated_off": [], "problems": []}))
    m = json.loads(R2_META.read_text())
    meta_edit(m["config"])
    stem = td / f"{tag}-stem"
    Path(str(stem) + ".h5ad.meta.json").write_text(json.dumps(m))
    Path(str(stem) + ".vcc").write_bytes(b"pkg")
    cl = w / "c0811_699.log"
    cl.write_text(f"2026-10-20T00:00:00Z package {stem}.vcc sha256 {hashlib.sha256(b'pkg').hexdigest()}\n2026-10-20T00:00:01Z novelty rc=0\n")
    p = subprocess.run([sys.executable, str(HERE / "gate699.py"), "--work", str(w), "--stem", str(stem), "--chainlog", str(cl)],
                       capture_output=True, text=True, cwd=ROOT)
    return p.stdout + p.stderr


def moved(c):
    c["bulk_load_sha256"] = "9" * 64


def main() -> int:
    if not (R2_RUN.is_dir() and R2_META.is_file() and R2_RECIPE.is_file()):
        print("695's run / meta / recipe not on disk"); return 2
    with tempfile.TemporaryDirectory(prefix="t699_") as t:
        td = Path(t)
        o = case(td, "b1", HC0, moved)
        line = next((ln for ln in o.splitlines() if ln.startswith("[c0811] PASS: config diff vs R-ABC2")), "")
        rec(bool(line) and "adopted ['DEF_DAY_DUMP_HOOK']" in line and "ALL GATES PASS" in o,
            "B1 only the hook is added and the dump moves -> ALL GATES PASS", o[-1200:])
        o = case(td, "b2", HC0, lambda c: None)
        rec("[c0811] FAIL: config diff vs R-ABC2" in o and "ALL GATES FAIL" in o, "B2 hook added, dump unmoved -> FAIL", o[-500:])
        o = case(td, "b3", None, moved)
        rec("[c0811] FAIL: config diff vs R-ABC2 (adopted [])" in o and "ALL GATES FAIL" in o,
            "B3 no hook (R-ABC2's recipe) but the dump moved -> FAIL", o[-500:])
        def mag(c):
            moved(c)
            c["de_mag"] = dict(c["de_mag"], variant="other")
        o = case(td, "b4", HC0, mag)
        rec("[c0811] FAIL: config diff vs R-ABC2" in o and "ALL GATES FAIL" in o, "B4 hook added + de_mag moved -> FAIL", o[-500:])
    n = sum(ok for ok, _ in RES)
    print(f"\n{n}/{len(RES)} PASS")
    return 0 if n == len(RES) else 1


if __name__ == "__main__":
    sys.exit(main())
