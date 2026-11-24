#!/usr/bin/env python3
"""695 tests: gate695.py's new logic on R-ABC's finished run (694 work/run: def_day.log, postbuild, ntclevel, audit428, npz) with
synthetic R-ABC2 metas. Scratch dirs only.
  K1 K adopted (DEF_DAY_RC_CALLSET_CPM=3800) and the meta moves rc_callset_cpm only -> ALL PASS, "[rabc2] PASS" anchored
  K2 K adopted but rc_callset_cpm unchanged -> FAIL (its key must move)
  L1 θ adopted but nothing moved -> FAIL
  L2 θ adopted, the dump renamed _s0.25 (+sha) and the tercilec key sha moved -> the config-diff row PASSes
  H1 the hook adopted, def_day's log lacks the hook lines -> FAIL (the hook row)
  X1 K adopted and an unadopted key (de_mag) also changed -> FAIL
usage: uv run python experiments/695-rabc2/test695.py
"""
from __future__ import annotations

import copy
import hashlib
import json
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


def case(td: Path, tag: str, extra_env: dict, meta_edit) -> str:
    w = td / tag
    w.mkdir()
    (w / "run").symlink_to(RABC_RUN)          # same run files: the npz equal R-ABC's, so only config keys can differ
    r = json.loads(RABC_RECIPE.read_text())
    r = {"env": dict(r["env"], **extra_env), "args": r["args"], "pending_treated_off": [], "problems": []}
    (w / "recipe695.json").write_text(json.dumps(r))
    m = json.loads(RABC_META.read_text())
    meta_edit(m["config"])
    stem = td / f"{tag}-stem"
    Path(str(stem) + ".h5ad.meta.json").write_text(json.dumps(m))
    Path(str(stem) + ".vcc").write_bytes(b"pkg")
    cl = w / "rabc2_695.log"
    cl.write_text(f"2026-10-18T00:00:00Z package {stem}.vcc sha256 {hashlib.sha256(b'pkg').hexdigest()}\n2026-10-18T00:00:01Z novelty rc=0\n")
    p = subprocess.run([sys.executable, str(HERE / "gate695.py"), "--work", str(w), "--stem", str(stem), "--chainlog", str(cl)],
                       capture_output=True, text=True, cwd=ROOT)
    return p.stdout + p.stderr


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="t695_") as t:
        td = Path(t)
        o = case(td, "k1", {"DEF_DAY_RC_CALLSET_CPM": "3800"}, lambda c: c.__setitem__("rc_callset_cpm", 3800))
        anchored = any(re.match(r"^\[rabc2\] PASS", ln) for ln in o.splitlines()) and any(re.match(r"^\[rabc2\] CLIFF OK", ln) for ln in o.splitlines())
        rec("ALL GATES PASS" in o and "[rabc2] PASS: config diff vs R-ABC (adopted ['DEF_DAY_RC_CALLSET_CPM'])" in o and anchored,
            "K1 K adopted, rc_callset_cpm moved -> ALL PASS, anchored lines", o[-900:])
        o = case(td, "k2", {"DEF_DAY_RC_CALLSET_CPM": "3800"}, lambda c: None)
        rec("ALL GATES FAIL" in o and "[rabc2] FAIL: config diff vs R-ABC" in o, "K2 K adopted, key unmoved -> FAIL", o[-600:])
        o = case(td, "l1", {"DEF_DAY_LODO_STRENGTH": "0.25"}, lambda c: None)
        rec("ALL GATES FAIL" in o and "[rabc2] FAIL: config diff vs R-ABC" in o, "L1 θ adopted, nothing moved -> FAIL", o[-600:])
        def lstr(c):   # what def_day's θ does to the final meta: the dump renamed _s0.25 (+sha), the tercilec key sha
            c["bulk_load"] = "dump_h1x2_cd4x1_s0.25.npz"
            c["bulk_load_sha256"] = "1" * 64
            c["bulk_sharpen_pt"] = dict(c["bulk_sharpen_pt"], carried_key_sha256="2" * 64)
        o = case(td, "l2", {"DEF_DAY_LODO_STRENGTH": "0.25"}, lstr)
        rec("[rabc2] PASS: config diff vs R-ABC (adopted ['DEF_DAY_LODO_STRENGTH'])" in o,
            "L2 θ adopted, dump renamed + tercilec key sha -> the config-diff row PASSes (pert_scale not required)", o[-700:])
        o = case(td, "h1", {"DEF_DAY_DUMP_HOOK": "experiments/707-strong-tier-direction/dumphook707.sh"},
                 lambda c: c.__setitem__("bulk_load_sha256", "0" * 64))
        rec("ALL GATES FAIL" in o and "FAIL dump hook (690g / 707)" in o, "H1 hook adopted without its def_day lines -> FAIL", o[-600:])
        def k_and_mag(c):
            c["rc_callset_cpm"] = 3800
            c["de_mag"] = dict(c["de_mag"], variant="other")
        o = case(td, "x1", {"DEF_DAY_RC_CALLSET_CPM": "3800"}, k_and_mag)
        rec("ALL GATES FAIL" in o and "[rabc2] FAIL: config diff vs R-ABC" in o, "X1 an unadopted key changed -> FAIL", o[-600:])
    n = sum(ok for ok, _ in RES)
    print(f"\n{n}/{len(RES)} PASS")
    return 0 if n == len(RES) else 1


if __name__ == "__main__":
    sys.exit(main())
