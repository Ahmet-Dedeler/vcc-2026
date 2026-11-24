#!/usr/bin/env python3
"""What `october_recipe.py command` does today, before 502's pass-through is applied.

Uses the live decisions file and the live def_day.sh, then a comment-only fake def_day via
OCTOBER_DEF_DAY. Does not modify either file.
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
SCRIPT = REPO / "scripts/october_recipe.py"
bad = 0


def expect(label, cond, detail=""):
    global bad
    print(f"{'OK' if cond else 'UNEXPECTED'}  {label}" + (f": {detail}" if detail != "" else ""))
    if not cond:
        bad += 1


def run(args, env):
    return subprocess.run([sys.executable, str(SCRIPT), *args], env=env, text=True, capture_output=True)


def main() -> int:
    base = dict(os.environ)
    base.pop("OCTOBER_DECISIONS", None)
    base.pop("OCTOBER_DEF_DAY", None)
    base.pop("OCTOBER_BUILDER", None)
    r = run(["command", "--bundle", "data/challenge", "--work", "W", "--name", "stem"], base)
    last = r.stdout.strip().splitlines()[-1]
    expect("live command exits 2 (PENDING wins over the missing knob)", r.returncode == 2, r.returncode)
    expect("stdout is a runnable line carrying DEF_DAY_DE_MAG",
           last.startswith("DEF_DAY_LODO=1 ") and "DEF_DAY_DE_MAG=dn_k250_zs" in last and last.endswith("noprog 4"),
           last)
    expect("stderr names both PENDING and the missing knob",
           "PENDING:" in r.stderr and "has no DEF_DAY_DE_MAG" in r.stderr, r.stderr.strip())

    r = run(["command", "--bundle", "data/challenge", "--work", "W", "--name", "stem", "--allow-pending"], base)
    last = r.stdout.strip().splitlines()[-1]
    expect("--allow-pending exits 3", r.returncode == 3, r.returncode)
    expect("--allow-pending still prints the rider", "DEF_DAY_DE_MAG=dn_k250_zs" in last, last)
    expect("--allow-pending stderr is only the missing knob",
           "has no DEF_DAY_DE_MAG" in r.stderr and "PENDING:" not in r.stderr, r.stderr.strip())

    # Pasting that line at current def_day.sh: the var is not VCC_*, the wipe does not export it.
    live = (REPO / "scripts/def_day.sh").read_text()
    expect("live def_day.sh has no DEF_DAY_DE_MAG", "DEF_DAY_DE_MAG" not in live, "")
    sim = r"""
        set -uo pipefail
        DEF_DAY_DE_MAG=dn_k250_zs
        for v in $(env | grep -o '^VCC_[A-Z0-9_]*'); do unset "$v"; done
        # current def_day.sh never reads DEF_DAY_DE_MAG and never exports VCC_DE_MAG
        echo "VCC_DE_MAG=${VCC_DE_MAG-<unset>}"
        echo "DEF_DAY_DE_MAG=${DEF_DAY_DE_MAG-<unset>}"
    """
    s = subprocess.run(["bash", "-c", sim], text=True, capture_output=True, check=True)
    expect("pasted command leaves VCC_DE_MAG unset", "VCC_DE_MAG=<unset>" in s.stdout, s.stdout.strip())
    expect("DEF_DAY_DE_MAG survives as an unread shell variable", "DEF_DAY_DE_MAG=dn_k250_zs" in s.stdout, "")

    # installed_problems is a substring test. A comment containing the token satisfies it.
    with tempfile.TemporaryDirectory() as td:
        fake = Path(td) / "def_day.sh"
        fake.write_text(
            live + "\n# note: DEF_DAY_DE_MAG is mentioned here and is not exported\n"
            'SLF="${6-pt:noprog}"\n'
        )
        env = dict(base, OCTOBER_DEF_DAY=str(fake))
        r = run(["command", "--bundle", "B", "--work", "W", "--name", "N", "--allow-pending"], env)
        expect("comment-only def_day makes installed_problems pass",
               r.returncode == 0 and "has no DEF_DAY_DE_MAG" not in r.stderr, f"rc={r.returncode} err={r.stderr.strip()}")
        expect("and stdout still tells you to export the rider",
               "DEF_DAY_DE_MAG=dn_k250_zs" in r.stdout, r.stdout.strip().splitlines()[-1])

    # The live-state unit test still expects --allow-pending rc 0 and no rider.
    t = subprocess.run([sys.executable, str(REPO / "scripts/test_october_recipe.py"),
                        "RecipeTests.test_live_state_refuses_pending_then_resolves_with_allow"],
                       text=True, capture_output=True)
    expect("test_live_state fails on the decided de_mag", t.returncode != 0, f"rc={t.returncode}")
    blob = t.stderr + t.stdout
    expect("the failure is the missing knob, not an import error",
           "has no DEF_DAY_DE_MAG" in blob or "AssertionError" in blob, blob[-400:])

    print("FAILS" if bad else "ALL DEMONSTRATED", bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
