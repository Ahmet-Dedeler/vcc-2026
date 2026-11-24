#!/usr/bin/env python3
"""698 tests: gate698's reference logic on 697's REAL finished R-ABC2H run (the reference C1 is judged against).
A test run dir mirrors 697's run (symlinks), except def_day.log, whose 690g hook lines name dumphook811.sh, as C1's will.
  B1 the recipe changes only the hook path (707 -> 811) and the meta moves bulk_load_sha256 -> ALL GATES PASS, adopted == [hook]
  B2 the hook path changes but the meta's dump is unmoved -> FAIL (a hook that changes nothing is not C1)
  B3 the recipe is R-ABC2H's exactly (hook unchanged) but the meta moves -> FAIL (adopted is not [hook])
  B4 the hook changes and an unadopted key (de_mag) moves too -> FAIL
usage: uv run python experiments/698-c1811/test698.py
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
R2H_RUN = ROOT / "experiments/697-rabc2h/work/run"
R2H_META = ROOT / "submissions/697rabc2h-october-recipe-abc.h5ad.meta.json"
R2H_RECIPE = ROOT / "experiments/697-rabc2h/work/recipe697.json"
H707 = "experiments/707-strong-tier-direction/dumphook707.sh"
H811 = "experiments/811-lolo-model/dumphook811.sh"
RES: list[tuple[bool, str]] = []


def rec(ok, name, detail=""):
    RES.append((bool(ok), name))
    print(f"{'PASS' if ok else 'FAIL'} {name}" + (f"  [{detail}]" if detail and not ok else ""), flush=True)


def case(td: Path, tag: str, hook: str, meta_edit) -> str:
    w = td / tag
    run = w / "run"
    run.mkdir(parents=True)
    for p in R2H_RUN.iterdir():   # mirror 697's run; def_day.log names the hook this case uses
        if p.name == "def_day.log":
            (run / p.name).write_text(p.read_text().replace(H707, hook))
        else:
            (run / p.name).symlink_to(p)
    base = json.loads(R2H_RECIPE.read_text())
    new = {"env": dict(base["env"], DEF_DAY_DUMP_HOOK=hook), "args": base["args"], "pending_treated_off": [], "problems": []}
    (w / "recipe698.json").write_text(json.dumps(new))
    m = json.loads(R2H_META.read_text())
    meta_edit(m["config"])
    stem = td / f"{tag}-stem"
    Path(str(stem) + ".h5ad.meta.json").write_text(json.dumps(m))
    Path(str(stem) + ".vcc").write_bytes(b"pkg")
    cl = w / "c1811_698.log"
    cl.write_text(f"2026-10-19T00:00:00Z package {stem}.vcc sha256 {hashlib.sha256(b'pkg').hexdigest()}\n2026-10-19T00:00:01Z novelty rc=0\n")
    p = subprocess.run([sys.executable, str(HERE / "gate698.py"), "--work", str(w), "--stem", str(stem), "--chainlog", str(cl)],
                       capture_output=True, text=True, cwd=ROOT)
    return p.stdout + p.stderr


def moved(c):
    c["bulk_load_sha256"] = "8" * 64


def main() -> int:
    if not (R2H_RUN.is_dir() and R2H_META.is_file() and R2H_RECIPE.is_file()):
        print("697's run / meta / recipe not on disk"); return 2
    with tempfile.TemporaryDirectory(prefix="t698_") as t:
        td = Path(t)
        o = case(td, "b1", H811, moved)
        line = next((ln for ln in o.splitlines() if ln.startswith("[c1811] PASS: config diff vs R-ABC2H")), "")
        rec(bool(line) and "adopted ['DEF_DAY_DUMP_HOOK']" in line and "ALL GATES PASS" in o,
            "B1 only the hook path changes and the dump moves -> ALL GATES PASS", o[-1200:])
        o = case(td, "b2", H811, lambda c: None)
        rec("[c1811] FAIL: config diff vs R-ABC2H" in o and "ALL GATES FAIL" in o, "B2 hook changed, dump unmoved -> FAIL", o[-500:])
        o = case(td, "b3", H707, moved)
        rec("[c1811] FAIL: config diff vs R-ABC2H (adopted [])" in o and "ALL GATES FAIL" in o,
            "B3 hook unchanged (R-ABC2H's recipe) but the dump moved -> FAIL", o[-500:])
        def mag(c):
            moved(c)
            c["de_mag"] = dict(c["de_mag"], variant="other")
        o = case(td, "b4", H811, mag)
        rec("[c1811] FAIL: config diff vs R-ABC2H" in o and "ALL GATES FAIL" in o, "B4 hook changed + de_mag moved -> FAIL", o[-500:])
    n = sum(ok for ok, _ in RES)
    print(f"\n{n}/{len(RES)} PASS")
    return 0 if n == len(RES) else 1


if __name__ == "__main__":
    sys.exit(main())
