#!/usr/bin/env python3
"""693p tests on a scratch tree (live files copied; synthetic KA records in a scratch dir until the real KA ends).
  P1 refused without the six KA records; refused with a record for another v2 sha
  P2 --check / apply with six matching records; second apply refused
  P3 the patched builder compiles (py_compile), keeps v1 on disk, imports dn_mag2 twice, pins v2's sha; def_day.sh bash -n
     and its meta check names v2
  P4 the installed src/vcc2026/dn_mag2.py is byte-identical to 693's dn_mag2.py
usage: uv run python experiments/693-dnmag-speed/test693p.py [--real-ka]   (--real-ka: use work/ka693_*.json)
"""
from __future__ import annotations

import hashlib
import json
import py_compile
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PY = sys.executable
V1 = "0a8cc560231a16f0fff7fbbf316e3a5e190cbdf20490239ef354b49b3ff60dac"
RESULTS: list[tuple[bool, str]] = []


def rec(ok: bool, name: str) -> None:
    RESULTS.append((bool(ok), name))
    print(f"{'PASS' if ok else 'FAIL'} {name}", flush=True)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    real = "--real-ka" in sys.argv
    td = Path(tempfile.mkdtemp(prefix="t693p_"))
    try:
        t = td / "repo"
        for f in ("scripts/make_submission_pm_cd4.py", "scripts/def_day.sh", "src/vcc2026/dn_mag.py"):
            (t / f).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / f, t / f)
        v2 = sha(HERE / "dn_mag2.py")
        kd = HERE / "work" if real else td / "ka"
        run = lambda *a: subprocess.run([PY, str(HERE / "patch693p.py"), "--root", str(t), "--ka-dir", str(kd), *a],
                                        capture_output=True, text=True)
        if not real:
            kd.mkdir()
            r = run("--check")
            rec(r.returncode == 1 and "missing KA record" in r.stderr, "P1 refused without KA records")
            for s in (0, 1, 2):
                for m in ("dn", "updn"):
                    (kd / f"ka693_s{s}_{m}.json").write_text(json.dumps({"blocks": 85, "v1_sha": V1, "v2_sha": "0" * 64}))
            r = run("--check")
            rec(r.returncode == 1 and "does not cover" in r.stderr, "P1 refused with records for another v2")
            for s in (0, 1, 2):
                for m in ("dn", "updn"):
                    (kd / f"ka693_s{s}_{m}.json").write_text(json.dumps({"blocks": 85, "v1_sha": V1, "v2_sha": v2}))
        rec(run("--check").returncode == 0, f"P2 --check OK ({'real' if real else 'synthetic'} KA records)")
        rec(run().returncode == 0, "P2 apply")
        rec(run().returncode == 1, "P2 second apply refused")
        b = (t / "scripts/make_submission_pm_cd4.py")
        try:
            py_compile.compile(str(b), doraise=True)
            ok = True
        except py_compile.PyCompileError:
            ok = False
        bt = b.read_text()
        rec(ok and bt.count("from vcc2026 import dn_mag2 as _dn_mag") == 2 and f'DE_DNMAG_MODULE_SHA256 = "{v2}"' in bt
            and 'ROOT / "src/vcc2026/dn_mag2.py"' in bt and sha(t / "src/vcc2026/dn_mag.py") == V1,
            "P3 builder compiles, imports dn_mag2 twice, pins v2; v1 kept")
        dd = (t / "scripts/def_day.sh").read_text()
        rec(subprocess.run(["bash", "-n", str(t / "scripts/def_day.sh")]).returncode == 0 and v2 in dd and V1 not in dd,
            "P3 def_day.sh bash -n; its meta check names v2 only")
        rec(sha(t / "src/vcc2026/dn_mag2.py") == v2, "P4 installed dn_mag2.py == 693's v2")
    finally:
        shutil.rmtree(td, ignore_errors=True)
    n = sum(ok for ok, _ in RESULTS)
    print(f"\n{n}/{len(RESULTS)} PASS")
    return 0 if n == len(RESULTS) else 1


if __name__ == "__main__":
    sys.exit(main())
