#!/usr/bin/env python3
"""694 tests: gate694.py and rabc694.sh's recipe check, on 524dr2's finished run (658 work/dr2/run: def_day.log, postbuild,
ntclevel, audit428) and synthetic metas. Scratch dirs only; nothing in submissions/ is touched.
  G1 on 524dr2 itself (no riders): the strict diff vs 524dr2 FAILs (diff []), everything else PASSes -> ALL GATES FAIL
  G2 synthetic R-ABC meta (524dr2 + de_uprn + de_dnmag at the pinned module shas, a different sameline.path) -> ALL PASS,
     the evidence lines start with "[rabc] PASS" / "[rabc] CLIFF OK" (arm1015's anchored greps), --write writes gates694.log
  G3 dn_mag v1's module sha in the meta -> FAIL;  G4 a moved re-level factor -> FAIL (only the 524dr2 diff catches it);
  G5 a package sha that does not match the chain log -> FAIL
  R1 the recipe check: today's resolved recipe OK; a dropped rider, a changed old value, an extra export -> STOP
usage: uv run python experiments/694-rabc/test694.py
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
PY = sys.executable
DR2 = ROOT / "experiments/658-prereg-1007-builds/work/dr2"
META2 = ROOT / "submissions/524dr2-october-recipe-abc.h5ad.meta.json"
UPRN = "7e6b3e52e2ad374b7b6e21e92327fcf65771433693eb3846c09c0189f3f2959f"
DN_V1 = "0a8cc560231a16f0fff7fbbf316e3a5e190cbdf20490239ef354b49b3ff60dac"
DN_V2 = "b9d6cfcf17ad783c1ea2607f95366173a0fc1f7920c4104604f9ee925ac55c87"
RES: list[tuple[bool, str]] = []


def rec(ok: bool, name: str, detail: str = "") -> None:
    RES.append((bool(ok), name))
    print(f"{'PASS' if ok else 'FAIL'} {name}" + (f"  [{detail}]" if detail and not ok else ""), flush=True)


def setup(td: Path, tag: str, recipe: dict, meta: dict, pkg_bytes: bytes = b"pkg", log_sha: str | None = None):
    w = td / tag
    w.mkdir()
    (w / "run").symlink_to(DR2 / "run")
    (w / "recipe694.json").write_text(json.dumps(recipe))
    stem = td / f"{tag}-stem"
    Path(str(stem) + ".h5ad.meta.json").write_text(json.dumps(meta))
    Path(str(stem) + ".vcc").write_bytes(pkg_bytes)
    sha = log_sha or hashlib.sha256(pkg_bytes).hexdigest()
    cl = w / "rabc694.log"
    cl.write_text(f"2026-10-15T00:00:00Z package {stem}.vcc sha256 {sha}\nOK: novel\n2026-10-15T00:00:01Z novelty rc=0\n")
    return w, stem, cl


def gate(w: Path, stem: Path, cl: Path, write: bool = False):
    a = [PY, str(HERE / "gate694.py"), "--work", str(w), "--stem", str(stem), "--chainlog", str(cl)] + (["--write"] if write else [])
    return subprocess.run(a, capture_output=True, text=True, cwd=ROOT)


def main() -> int:
    r2 = json.loads((DR2 / "recipe524.json").read_text())
    m2 = json.loads(META2.read_text())
    new_env = dict(r2["env"], DEF_DAY_DE_UPRN="b_k100_q20r3p10", DEF_DAY_DE_DNMAG="d_m200")
    rabc_recipe = {"env": new_env, "args": r2["args"], "pending_treated_off": [], "problems": []}

    def rabc_meta(dn_sha=DN_V2, ctx=None):
        m = copy.deepcopy(m2)
        c = m["config"]
        c["de_uprn"] = {"module_sha256": UPRN, "variant": "b_k100_q20r3p10", "where": "x"}
        c["de_dnmag"] = {"module_sha256": dn_sha, "variant": "d_m200", "where": "x"}
        c["sameline"]["path"] = "/elsewhere/experiments/694-rabc/work/run/sameline.npz"
        if ctx:
            c["bulk_ctx_scale"] = ctx
        return m

    with tempfile.TemporaryDirectory(prefix="t694_") as t:
        td = Path(t)
        # G1
        w, s, cl = setup(td, "g1", r2, m2)
        r = gate(w, s, cl)
        o = r.stdout
        rec(r.returncode == 0 and "ALL GATES FAIL" in o and "[rabc] FAIL: config diff vs 524dr2 == []" in o
            and "[rabc] PASS: config diff vs 453dr" in o and "[rabc] CLIFF OK" in o and o.count("  FAIL ") == 0,
            "G1 524dr2 itself: only the strict 524dr2 diff fails", o[-600:])
        # G2
        w, s, cl = setup(td, "g2", rabc_recipe, rabc_meta())
        r = gate(w, s, cl, write=True)
        o = r.stdout
        lines = o.splitlines()
        anchored = (any(re.match(r"^\[rabc\] PASS", ln) for ln in lines) and any(re.match(r"^\[rabc\] CLIFF OK", ln) for ln in lines)
                    and "[rabc] PASS: config diff vs 524dr2 == ['de_dnmag', 'de_uprn']" in o)
        rec(r.returncode == 0 and "ALL GATES PASS" in o and anchored and (w / "gates694.log").is_file()
            and "ok   rank-neutral UP rider (630)" in o and "ok   DOWN magnitude rider (631, dn_mag2)" in o,
            "G2 synthetic R-ABC meta: ALL PASS, anchored evidence lines, gates694.log written", o[-900:])
        # BSD grep, exactly as arm1015 runs it, on a log that holds the gate's raw output
        rl = td / "g2-run.log"
        rl.write_text(o + "2026-10-15T00:00:02Z novelty rc=0\n")
        g = lambda pat: subprocess.run(["/usr/bin/grep", "-q", pat, str(rl)]).returncode == 0
        rec(g("^\\[rabc\\] PASS") and g("^\\[rabc\\] CLIFF OK") and g("novelty rc=0"), "G2 arm1015's three BSD greps match")
        # G3
        w, s, cl = setup(td, "g3", rabc_recipe, rabc_meta(dn_sha=DN_V1))
        r = gate(w, s, cl, write=True)
        rec(r.returncode == 1 and "ALL GATES FAIL" in r.stdout and "FAIL DOWN magnitude rider" in r.stdout
            and not (w / "gates694.log").exists(), "G3 dn_mag v1 module sha -> FAIL, no gates694.log", r.stdout[-500:])
        # G4
        moved = dict(m2["config"]["bulk_ctx_scale"], E=1.1114)
        w, s, cl = setup(td, "g4", rabc_recipe, rabc_meta(ctx=moved))
        r = gate(w, s, cl)
        rec("ALL GATES FAIL" in r.stdout and "[rabc] FAIL: config diff vs 524dr2 == ['bulk_ctx_scale', 'de_dnmag', 'de_uprn']" in r.stdout
            and "[rabc] PASS: config diff vs 453dr" in r.stdout, "G4 a moved re-level factor -> FAIL (524dr2 diff only)", r.stdout[-500:])
        # G5
        w, s, cl = setup(td, "g5", rabc_recipe, rabc_meta(), log_sha="0" * 64)
        r = gate(w, s, cl)
        rec("ALL GATES FAIL" in r.stdout and "MISMATCH" in r.stdout, "G5 package sha mismatch -> FAIL", r.stdout[-300:])

        # G6-G9 (gate v2): a rebuilt run in a NEW work dir. Its dump `env` and same-line `provenance` carry that dir, so the
        # npz file shas differ while every array is 524dr2's.
        import numpy as np

        def provrun(tag: str, bump_array: bool = False, bump_text: bool = False, w_change: bool = False):
            w = td / tag
            (w / "run").mkdir(parents=True)
            src = DR2 / "run"
            for f in src.iterdir():
                if f.name not in ("dump_h1x2_cd4x1.npz", "sameline.npz"):
                    (w / "run" / f.name).symlink_to(f)
            newrun = str((w / "run").resolve())
            shas = {}
            for name, key in (("dump_h1x2_cd4x1.npz", "env"), ("sameline.npz", "provenance")):
                z = dict(np.load(src / name, allow_pickle=False))
                old = str(z[key]).replace(str(src.resolve().relative_to(ROOT)), newrun)
                if bump_text and name == "sameline.npz":
                    old = old.replace("same-line", "programme")
                z[key] = np.array(old)
                if bump_array and name == "dump_h1x2_cd4x1.npz":
                    k = next(k for k in z if z[k].dtype.kind == "f" and z[k].size)
                    z[k] = z[k].copy(); z[k].flat[0] += 1.0
                np.savez(w / "run" / name, **z)
                shas[name] = hashlib.sha256((w / "run" / name).read_bytes()).hexdigest()
            (w / "recipe694.json").write_text(json.dumps(rabc_recipe))
            m = rabc_meta()
            m["config"]["bulk_load_sha256"] = shas["dump_h1x2_cd4x1.npz"]
            m["config"]["sameline"]["sha256"] = shas["sameline.npz"]
            m["config"]["sameline"]["path"] = newrun + "/sameline.npz"
            if w_change:
                m["config"]["sameline"]["w"] = 0.6
            stem = td / f"{tag}-stem"
            Path(str(stem) + ".h5ad.meta.json").write_text(json.dumps(m))
            Path(str(stem) + ".vcc").write_bytes(b"pkg")
            cl = w / "rabc694.log"
            cl.write_text(f"2026-10-16T00:00:00Z package {stem}.vcc sha256 {hashlib.sha256(b'pkg').hexdigest()}\n"
                          "2026-10-16T00:00:01Z novelty rc=0\n")
            return gate(w, stem, cl).stdout

        o6 = provrun("g6")
        rec("ALL GATES PASS" in o6 and "provenance-only ['bulk_load_sha256', 'sameline.sha256']" in o6,
            "G6 new work dir, arrays identical: provenance-only shas accepted -> ALL PASS", o6[-700:])
        o7 = provrun("g7", bump_array=True)
        rec("ALL GATES FAIL" in o7 and "FAIL bulk_load_sha256 vs 524dr2 is provenance only: array" in o7,
            "G7 one dump array changed -> FAIL", o7[-500:])
        o8 = provrun("g8", bump_text=True)
        rec("ALL GATES FAIL" in o8 and "provenance differs beyond the work dir" in o8, "G8 provenance text changed beyond the path -> FAIL",
            o8[-500:])
        o9 = provrun("g9", w_change=True)
        rec("ALL GATES FAIL" in o9 and "the same-line record differs beyond its sha256" in o9, "G9 same-line w changed -> FAIL",
            o9[-500:])

        # R1: rabc694.sh's recipe check, extracted verbatim from the script
        sh = (HERE / "rabc694.sh").read_text()
        blk = sh.split('CMP="$(uv run python - "$W/${PFX}recipe694.json" "$PREV_RECIPE" "$BUNDLE" "$RW" "$NAME" <<\'PY\'\n', 1)[1].split("\nPY\n", 1)[0]
        chk = td / "cmp.py"
        chk.write_text(blk)
        bundle, rw = r2["args"][0], "experiments/694-rabc/work/run"

        def cmp(recipe: dict) -> str:
            p = td / "r.json"
            rr = copy.deepcopy(recipe)
            rr["args"] = [bundle, rw, "rabc"] + r2["args"][3:]
            p.write_text(json.dumps(rr))
            return subprocess.run([PY, str(chk), str(p), str(DR2 / "recipe524.json"), bundle, rw, "rabc"],
                                  capture_output=True, text=True).stdout.strip()
        ok_case = cmp(rabc_recipe)
        drop = copy.deepcopy(rabc_recipe); del drop["env"]["DEF_DAY_DE_DNMAG"]
        chg = copy.deepcopy(rabc_recipe); chg["env"]["DEF_DAY_MDE_SIGN"] = "clean4:0.10"
        extra = copy.deepcopy(rabc_recipe); extra["env"]["DEF_DAY_RHO_SCALE"] = "0.99"
        live = json.loads((HERE / "work/dry-recipe694.json").read_text()) if (HERE / "work/dry-recipe694.json").is_file() else None
        rec(ok_case.startswith("OK ") and cmp(drop).startswith("STOP") and cmp(chg).startswith("STOP") and cmp(extra).startswith("STOP")
            and (live is None or cmp(live).startswith("OK ")),
            "R1 recipe check: today's recipe OK; dropped rider / changed value / extra export STOP",
            f"{ok_case} | {cmp(drop)} | {cmp(chg)} | {cmp(extra)}")
    n = sum(ok for ok, _ in RES)
    print(f"\n{n}/{len(RES)} PASS")
    return 0 if n == len(RES) else 1


if __name__ == "__main__":
    sys.exit(main())
