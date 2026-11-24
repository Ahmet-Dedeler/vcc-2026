"""247 item B #10 — `def_census.py`'s coverage gate, FIRING on a collapse and PASSING on the draws.

    pass-1  A/B/C, --verify, production banks                    -> VERIFY PASS, GATE PASS
    pass-2  199's pessimistic def-bundle + its re-extracted banks -> GATE PASS (150/291/293, 174, 167, 0)
    fire-3  199's def-bundle counted against the A/B/C banks      -> GATE FAIL (CD4 union 0, complete 0)
    fire-4  a synthetic bundle whose 300 targets NO production donor carries
                                                                  -> GATE FAIL on every row
    BLIND   the pre-247 script on fire-4                          -> exit 0

    uv run python experiments/247-oct22-guards/census247.py
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work" / "census"
ABC = ROOT / "data/challenge"
DEF199 = ROOT / "experiments/199-def-runbook/work/def-bundle"
DEF199_CD4 = ROOT / "experiments/199-def-runbook/work/cd4_compact_def199"


def run(script: Path, args: list[str], pin_root: bool = False) -> tuple[int, str]:
    if pin_root:
        # the pre-247 byte copy computes ROOT from __file__; pin it to scripts/ (stage222's trick:
        # set module.__file__ BEFORE exec_module, since runpy.run_path overrides init_globals)
        code = ("import importlib.util, sys\n"
                f"spec = importlib.util.spec_from_file_location('def_census_pre247', {str(script)!r})\n"
                "m = importlib.util.module_from_spec(spec)\n"
                f"m.__file__ = {str(ROOT / 'scripts/def_census.py')!r}\n"
                "sys.modules['def_census_pre247'] = m\n"
                "spec.loader.exec_module(m)\n"
                "sys.argv = ['def_census.py'] + sys.argv[1:]\n"
                "m.main()\n")
        cmd = ["uv", "run", "python", "-c", code, *args]
    else:
        cmd = ["uv", "run", "python", str(script), *args]
    proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    return proc.returncode, proc.stdout + proc.stderr


def collapsed_bundle() -> Path:
    """300 targets drawn from the 953 axis genes NO production donor carries (247 measured)."""
    out = WORK / "collapsed-bundle"
    out.mkdir(parents=True, exist_ok=True)
    genes = [str(g) for g in pd.read_csv(ABC / "gene_names.csv")["gene_name"]]
    covered: set[str] = set()
    for d in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        covered |= {ln.strip() for ln in (ROOT / "data/donors" / d / "perts.txt").read_text().splitlines() if ln.strip()}
    pool = [g for g in genes if g not in covered]
    assert len(pool) >= 300, len(pool)
    panel = pool[:300]
    shutil.copyfile(ABC / "gene_names.csv", out / "gene_names.csv")
    pd.DataFrame({"target_gene": panel, "n_cells": [400] * 300}).to_csv(out / "pert_counts.csv", index=False)
    manifest = json.loads((DEF199 / "manifest.json").read_text())
    manifest["panel_id"] = "vcc2026-synthetic-collapsed-247"
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return out


def gate_lines(text: str) -> list[str]:
    return [ln.rstrip() for ln in text.splitlines()
            if ln.strip().startswith(("PASS", "FAIL", "GATE", "VERIFY"))]


def main() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    new = ROOT / "scripts/def_census.py"
    old = HERE / "work/pre247/def_census.py"
    collapsed = collapsed_bundle()
    cases = {
        "pass-1 A/B/C --verify": (new, ["--verify", "--out", str(WORK / "abc.json")], False, 0),
        "pass-2 199 def-bundle + its banks": (new, ["--bundle", str(DEF199), "--cd4-root", str(DEF199_CD4),
                                                    "--out", str(WORK / "def199.json")], False, 0),
        "fire-3 199 def-bundle vs the A/B/C banks": (new, ["--bundle", str(DEF199), "--out", str(WORK / "def199-abcbanks.json")], False, 1),
        "fire-4 collapsed synthetic panel": (new, ["--bundle", str(collapsed), "--cd4-root", str(DEF199_CD4),
                                                    "--out", str(WORK / "collapsed.json")], False, 1),
        "BLIND pre-247 on fire-4": (old, ["--bundle", str(collapsed), "--cd4-root", str(DEF199_CD4)], True, 0),
    }
    result = {}
    for name, (script, args, pin, expect) in cases.items():
        rc, out = run(script, args, pin_root=pin)
        (WORK / (name.split()[0] + ".log")).write_text(out)
        lines = gate_lines(out)
        result[name] = {"rc": rc, "expect_rc": expect, "ok": rc == expect, "gate_lines": lines}
        print(f"\n=== {name}: rc={rc} (expected {expect})")
        for ln in lines:
            print(f"    {ln}")
        if rc not in (0, 1):
            print(out[-2000:])
    ok = all(r["ok"] for r in result.values())
    (WORK / "census247.json").write_text(json.dumps(result, indent=1))
    print(f"\nITEM B#10 {'PASS' if ok else 'FAIL'}")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
