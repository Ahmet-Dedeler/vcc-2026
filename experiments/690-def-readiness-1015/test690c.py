#!/usr/bin/env python3
"""690c tests. Light: the restore() trap cut from the patched def_day.sh runs in a bash harness against a scratch repo.
  C1  --check, apply, second apply refused; bash -n on the patched file
  C2  with a pin record: builder.pinned.py == the builder at trap time; the restore command still runs; rc preserved
  C3  without a pin record: nothing written; rc preserved
  C4  with DEF_DAY_KEEP_PINS=1: the copy is still saved, the builder is left pinned, rc preserved
usage: uv run python experiments/690-def-readiness-1015/test690c.py
"""
from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PY = sys.executable
RESULTS: list[tuple[bool, str]] = []


def rec(ok: bool, name: str) -> None:
    RESULTS.append((bool(ok), name))
    print(f"{'PASS' if ok else 'FAIL'} {name}", flush=True)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    td = Path(tempfile.mkdtemp(prefix="t690c_"))
    try:
        t = td / "repo"
        (t / "scripts").mkdir(parents=True)
        shutil.copy2(ROOT / "scripts/def_day.sh", t / "scripts/def_day.sh")
        r = subprocess.run([PY, str(HERE / "patch690c.py"), "--root", str(t), "--check"], capture_output=True, text=True)
        rec(r.returncode == 0, "C1 --check on the live def_day.sh")
        r = subprocess.run([PY, str(HERE / "patch690c.py"), "--root", str(t)], capture_output=True, text=True)
        rec(r.returncode == 0, "C1 apply")
        r = subprocess.run([PY, str(HERE / "patch690c.py"), "--root", str(t)], capture_output=True, text=True)
        rec(r.returncode == 1 and "already applied" in r.stderr, "C1 second apply refused")
        dd = (t / "scripts/def_day.sh").read_text()
        rec(subprocess.run(["bash", "-n", str(t / "scripts/def_day.sh")]).returncode == 0, "C1 bash -n")
        i = dd.index("restore(){")
        j = dd.index("\n}\n", i) + 3
        fn = dd[i:j]
        # harness: the trap function, with the restore command stubbed (it records that it ran and restores a file)
        harness = ('log(){ echo "LOG: $*" >> "$W/log.txt"; }\n'
                   'uv(){ echo restore-ran >> "$W/log.txt"; cp "$W/before.py" scripts/make_submission_pm_cd4.py; }\n'
                   + fn + '\ntrap restore EXIT\nexit "${RC:-0}"\n')
        def case(name: str, pin: bool, keep: bool, rc: int) -> tuple[subprocess.CompletedProcess, Path]:
            w = td / name
            w.mkdir()
            (t / "scripts/make_submission_pm_cd4.py").write_text(f"# pinned builder for {name}\n")
            (w / "before.py").write_text("# pre-pin builder\n")
            if pin:
                (w / "pin-artifacts.json").write_text("{}")
            env = {k: v for k, v in os.environ.items() if not k.startswith("DEF_DAY_")}
            # the trap compares the post-restore builder with BSHA0: the pre-pin copy if a pin ran, else the untouched builder
            b0 = sha(w / "before.py") if pin else sha(t / "scripts/make_submission_pm_cd4.py")
            env.update(W=str(w), RC=str(rc), BSHA0=b0, LOG=str(w / "def_day.log"))
            if keep:
                env["DEF_DAY_KEEP_PINS"] = "1"
            r = subprocess.run(["bash", "-c", harness], capture_output=True, text=True, cwd=t, env=env)
            return r, w
        r, w = case("pin", True, False, 5)
        ok = (w / "builder.pinned.py").read_text() == "# pinned builder for pin\n" and "restore-ran" in (w / "log.txt").read_text() \
            and (t / "scripts/make_submission_pm_cd4.py").read_text() == "# pre-pin builder\n" and r.returncode == 5
        rec(ok, f"C2 pin record: copy == trap-time builder, restore ran, rc 5 kept (rc {r.returncode})")
        r, w = case("nopin", False, False, 0)
        rec(not (w / "builder.pinned.py").exists() and r.returncode == 0, "C3 no pin record: nothing written, rc 0")
        r, w = case("keep", True, True, 3)
        ok = (w / "builder.pinned.py").exists() and "restore-ran" not in (w / "log.txt").read_text() \
            and (t / "scripts/make_submission_pm_cd4.py").read_text() == "# pinned builder for keep\n" and r.returncode == 3
        rec(ok, f"C4 KEEP_PINS: copy saved, builder left pinned, rc 3 kept (rc {r.returncode})")
    finally:
        shutil.rmtree(td, ignore_errors=True)
    n = sum(ok for ok, _ in RESULTS)
    print(f"\n{n}/{len(RESULTS)} PASS")
    return 0 if n == len(RESULTS) else 1


if __name__ == "__main__":
    sys.exit(main())
