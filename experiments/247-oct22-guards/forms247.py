"""247 item B #8 — `def_pipeline.other_builders()` on all six launch forms, with a sleeping dummy.

The dummy is a file literally named `scripts/make_submission_pm_cd4.py` under `work/forms/`, whose
body is `time.sleep(90)`, so each launched process carries the builder's file name on its command
line exactly as a real build would while doing nothing.  Each form is launched, both predicates
(247 and the pre-247 byte copy) are run against `ps`, and every dummy is killed.

    uv run python experiments/247-oct22-guards/forms247.py
"""

from __future__ import annotations

import importlib.util
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work" / "forms"
PY = str(ROOT / ".venv/bin/python")


def load_script(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    module.__file__ = str(ROOT / "scripts/def_pipeline.py")
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    dummy_dir = WORK / "scripts"
    dummy_dir.mkdir(exist_ok=True)
    dummy = dummy_dir / "make_submission_pm_cd4.py"
    dummy.write_text("import time\ntime.sleep(90)\n")
    rel = str(dummy.relative_to(ROOT))
    forms = {
        "1 plain: python <script>": [PY, rel],
        "2 absolute script path": [PY, str(dummy)],
        "3 two interpreter flags: python -X dev -u <script>": [PY, "-X", "dev", "-u", rel],
        "4 runpy (how the test suite starts it): python -c \"import runpy; runpy.run_path('<script>')\"":
            [PY, "-c", f"import runpy; runpy.run_path({rel!r})"],
        "5 uv shim: uv run python <script> (the `uv` parent, before/after exec)": ["uv", "run", "python", rel],
        "6 wrapped: /usr/bin/time -l python <script>": ["/usr/bin/time", "-l", PY, rel],
    }
    new = load_script(ROOT / "scripts/def_pipeline.py", "def_pipeline_247")
    old = load_script(HERE / "work/pre247/def_pipeline.py", "def_pipeline_pre247")

    result = {}
    procs = []
    try:
        for name, cmd in forms.items():
            p = subprocess.Popen(cmd, cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                 start_new_session=True)
            procs.append(p)
            time.sleep(2.5 if cmd[0] == "uv" else 0.8)   # let uv resolve and spawn its child
            hits_new = new.other_builders()
            hits_old = old.other_builders()
            # only lines that belong to THIS form's process tree (the dummy path is unique to it)
            mine_new = [h for h in hits_new if rel in h or str(dummy) in h]
            mine_old = [h for h in hits_old if rel in h or str(dummy) in h]
            result[name] = {"cmd": cmd, "247": "DETECTED" if mine_new else "MISSED",
                            "pre247": "DETECTED" if mine_old else "MISSED",
                            "hits_247": mine_new, "hits_pre247": mine_old}
            print(f"{result[name]['247']:<9} 247   | {result[name]['pre247']:<9} pre-247 | {name}")
            for h in mine_new:
                print(f"            {h[:150]}")
            if cmd[0] in ("uv", "/usr/bin/time"):
                # the WRAPPER line alone (the `uv` shim before it execs, `/usr/bin/time`): 242's
                # third miss. The child interpreter is a separate process the old predicate sees.
                parent_new = [h for h in mine_new if h.split(None, 1)[1].startswith(cmd[0])]
                parent_old = [h for h in mine_old if h.split(None, 1)[1].startswith(cmd[0])]
                sub = f"{name.split(':')[0]}-parent: the `{cmd[0]}` line alone"
                result[sub] = {"cmd": cmd, "247": "DETECTED" if parent_new else "MISSED",
                               "pre247": "DETECTED" if parent_old else "MISSED",
                               "hits_247": parent_new, "hits_pre247": parent_old}
                print(f"{result[sub]['247']:<9} 247   | {result[sub]['pre247']:<9} pre-247 | {sub}")
            os.killpg(os.getpgid(p.pid), signal.SIGTERM)
            p.wait(timeout=10)
            time.sleep(0.3)
    finally:
        for p in procs:
            if p.poll() is None:
                try:
                    os.killpg(os.getpgid(p.pid), signal.SIGKILL)
                except ProcessLookupError:
                    pass
    # and with the dummies gone, neither predicate may still see one (no leftover false positive);
    # a REAL build running concurrently is a true positive and is reported, not counted against.
    quiet_new, quiet_old = new.other_builders(), old.other_builders()
    leftover = [h for h in quiet_new if "work/forms" in h]
    real_new = [h for h in quiet_new if "work/forms" not in h]
    real_old = [h for h in quiet_old if "work/forms" not in h]
    result["_quiet"] = {"247_leftover_dummies": leftover, "247_real_builds": real_new,
                        "pre247_real_builds": real_old}
    print(f"\ndummies gone: leftover false positives -> {leftover or 'none'}")
    if real_new or real_old:
        print(f"a REAL build is running: 247 sees {len(real_new)} process line(s), pre-247 sees {len(real_old)}:")
        for h in real_new:
            print(f"    {h[:110]}")
    # a false-positive control: the test suite's file name must NOT read as a build
    fp = subprocess.Popen([PY, "-c", "import time; time.sleep(5)  # scripts/test_make_submission_pm_cd4.py"],
                          cwd=ROOT, start_new_session=True)
    time.sleep(0.6)
    fp_hits = [h for h in new.other_builders() if "test_make_submission_pm_cd4" in h]
    os.killpg(os.getpgid(fp.pid), signal.SIGTERM)
    result["_false_positive_control"] = fp_hits
    print(f"false-positive control (test_make_submission_pm_cd4.py on a command line): "
          f"{'FLAGGED (bad)' if fp_hits else 'not flagged (good)'}")
    (WORK / "forms247.json").write_text(json.dumps(result, indent=1))
    rows = {k: r for k, r in result.items() if not k.startswith("_")}
    main_rows = {k: r for k, r in rows.items() if k in forms}
    ok = (all(r["247"] == "DETECTED" for r in rows.values()) and not leftover and not fp_hits)
    missed_old = [k for k, r in rows.items() if r["pre247"] == "MISSED"]
    print(f"\npre-247 missed {len(missed_old)} of {len(rows)} rows ({len(forms)} forms + 2 wrapper-parent sub-rows): {missed_old}")
    print(f"ITEM B#8 {'PASS' if ok else 'FAIL'}: 247 detects {sum(r['247'] == 'DETECTED' for r in main_rows.values())}/{len(forms)} forms "
          f"and {sum(r['247'] == 'DETECTED' for r in rows.values())}/{len(rows)} rows")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
