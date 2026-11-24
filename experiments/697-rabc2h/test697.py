#!/usr/bin/env python3
"""697 tests: gate697's reference logic. R-ABC2 does not exist yet, so R-ABC (694) stands in as the reference run (overrides
GATE697_REF_META / _RUN / _RECIPE). Its finished run supplies def_day.log, postbuild, ntclevel, audit428 and the npz.
  A1 reference recipe = R-ABC + θ (as if R-ABC2 adopted θ); new recipe = that + K 3000; the meta moves rc_callset_cpm only
     -> ALL PASS: θ is the reference, not a change, so it is not required to move
  A2 the same, but the meta does not move rc_callset_cpm -> FAIL
  A3 the same, plus an unadopted key (de_mag) moved -> FAIL
usage: uv run python experiments/697-rabc2h/test697.py
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RABC_RUN = ROOT / "experiments/694-rabc/work/run"
RABC_META = ROOT / "submissions/694rabc-october-recipe-abc.h5ad.meta.json"
RABC_RECIPE = ROOT / "experiments/694-rabc/work/recipe694.json"
RES: list[tuple[bool, str]] = []


def rec(ok, name, detail=""):
    RES.append((bool(ok), name))
    print(f"{'PASS' if ok else 'FAIL'} {name}" + (f"  [{detail}]" if detail and not ok else ""), flush=True)


def case(td: Path, tag: str, meta_edit) -> str:
    w = td / tag
    w.mkdir()
    (w / "run").symlink_to(RABC_RUN)
    base = json.loads(RABC_RECIPE.read_text())
    ref = {"env": dict(base["env"], DEF_DAY_LODO_STRENGTH="0.25"), "args": base["args"]}   # a stand-in "R-ABC2" with θ
    (w / "ref_recipe.json").write_text(json.dumps(ref))
    new = {"env": dict(ref["env"], DEF_DAY_DUMP_HOOK="experiments/707-strong-tier-direction/dumphook707.sh"), "args": base["args"], "pending_treated_off": [], "problems": []}
    (w / "recipe697.json").write_text(json.dumps(new))
    m = json.loads(RABC_META.read_text())
    meta_edit(m["config"])
    stem = td / f"{tag}-stem"
    Path(str(stem) + ".h5ad.meta.json").write_text(json.dumps(m))
    Path(str(stem) + ".vcc").write_bytes(b"pkg")
    cl = w / "rabc2h_697.log"
    cl.write_text(f"2026-10-18T00:00:00Z package {stem}.vcc sha256 {hashlib.sha256(b'pkg').hexdigest()}\n2026-10-18T00:00:01Z novelty rc=0\n")
    env = dict(os.environ, GATE697_REF_META=str(RABC_META), GATE697_REF_RUN=str(RABC_RUN), GATE697_REF_RECIPE=str(w / "ref_recipe.json"))
    p = subprocess.run([sys.executable, str(HERE / "gate697.py"), "--work", str(w), "--stem", str(stem), "--chainlog", str(cl)],
                       capture_output=True, text=True, cwd=ROOT, env=env)
    return p.stdout + p.stderr


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="t697_") as t:
        td = Path(t)
        o = case(td, "a1", lambda c: c.__setitem__("bulk_load_sha256", "3" * 64))
        line = next((ln for ln in o.splitlines() if ln.startswith("[rabc2h] PASS: config diff vs R-ABC2")), "")
        rec(bool(line) and "adopted ['DEF_DAY_DUMP_HOOK']" in line and re.search(r"^\[rabc2h\] CLIFF OK", o, re.M) is not None,
            "A1 θ in the reference, the hook new: the diff row PASSes with adopted == [hook] only", o[-900:])
        o = case(td, "a2", lambda c: None)
        rec("[rabc2h] FAIL: config diff vs R-ABC2" in o, "A2 hook new but the dump unmoved -> FAIL", o[-500:])
        def k_mag(c):
            c["bulk_load_sha256"] = "3" * 64
            c["de_mag"] = dict(c["de_mag"], variant="other")
        o = case(td, "a3", k_mag)
        rec("[rabc2h] FAIL: config diff vs R-ABC2" in o, "A3 an unadopted key moved -> FAIL", o[-500:])
    n = sum(ok for ok, _ in RES)
    print(f"\n{n}/{len(RES)} PASS")
    return 0 if n == len(RES) else 1


if __name__ == "__main__":
    sys.exit(main())
