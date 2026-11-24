"""248 -- import-time validation of the two new variables, every combination, in subprocesses.

Each case imports the builder (working tree by default, or `--builder <path>`) with one
environment and reports whether the import raised, and with what message.  The rule under test:
both variables or neither; names in the registry; no repeats; no overlap with an ACTIVE union 1.

    uv run python experiments/248-second-residual/envcheck248.py [--builder <path>]
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BUILDER = Path(sys.argv[sys.argv.index("--builder") + 1]) if "--builder" in sys.argv \
    else ROOT / "scripts/make_submission_pm_cd4.py"

BASE = {"VCC_ALPHA": "0.982", "VCC_CD4_TARGET_F": "0.34"}
CASES = [
    ("unset (default path)", {}, None),
    ("union2 without a dose", {"VCC_UNION2_SOURCES": "dld1_lowmoi"}, "VCC_UNION2_G is not"),
    ("dose without union2", {"VCC_UNION2_G": "0.35"}, "would do nothing"),
    ("empty name list", {"VCC_UNION2_SOURCES": " , ", "VCC_UNION2_G": "0.35"}, "names no source"),
    ("unknown name", {"VCC_UNION2_SOURCES": "not_a_bank", "VCC_UNION2_G": "0.35"}, "not in UNION_SOURCE_REGISTRY"),
    ("repeated name", {"VCC_UNION2_SOURCES": "dld1_lowmoi,dld1_lowmoi", "VCC_UNION2_G": "0.35"}, "repeats a source"),
    ("dose out of range", {"VCC_UNION2_SOURCES": "dld1_lowmoi", "VCC_UNION2_G": "5"}, "must be in (0, 5)"),
    ("overlap with ACTIVE union 1", {"VCC_H1_UNION_G": "0.28", "VCC_UNION2_SOURCES": "h1_2025,dld1_lowmoi", "VCC_UNION2_G": "0.35"}, "dosed twice"),
    ("H1 pair as union 2 with union 1 UNSET (row 2)", {"VCC_UNION2_SOURCES": "gse295214_h1,h1_2025", "VCC_UNION2_G": "0.28"}, None),
    ("the arm (CUM1)", {"VCC_H1_UNION_G": "0.28", "VCC_UNION2_SOURCES": "dld1_lowmoi", "VCC_UNION2_G": "0.35"}, None),
    ("CUM4", {"VCC_H1_UNION_G": "0.28", "VCC_UNION2_SOURCES": "dld1_lowmoi,jiang_5line_noK562,igvf_cardiomyocyte,myllia_a549", "VCC_UNION2_G": "0.35"}, None),
]
PROBE = (
    "import importlib.util, sys, json; "
    f"spec = importlib.util.spec_from_file_location('b', {str(BUILDER)!r}); m = importlib.util.module_from_spec(spec); "
    f"m.__file__ = {str(ROOT / 'scripts/make_submission_pm_cd4.py')!r}; sys.modules['b'] = m; spec.loader.exec_module(m); "
    "print(json.dumps({'UNION2_SOURCES': None if m.UNION2_SOURCES is None else [n for n, _ in m.UNION2_SOURCES], "
    "'UNION2_G': m.UNION2_G, 'H1_UNION_G': m.H1_UNION_G, 'union1': list(m.UNION_SOURCE_NAMES)}))"
)


def main() -> None:
    rows = []
    ok_all = True
    for label, env, want in CASES:
        e = {k: v for k, v in os.environ.items() if not k.startswith("VCC_")}
        e.update(BASE)
        e.update(env)
        p = subprocess.run([sys.executable, "-c", PROBE], env=e, capture_output=True, text=True, cwd=ROOT)
        raised = p.returncode != 0
        msg = p.stderr.strip().splitlines()[-1] if raised else p.stdout.strip().splitlines()[-1]
        if want is None:
            ok = not raised
        else:
            ok = raised and want in msg
        ok_all &= ok
        rows.append(dict(case=label, env=env, expect="import OK" if want is None else f"refused: {want}",
                         raised=raised, message=msg[:220], ok=ok))
        print(f"[{'PASS' if ok else 'FAIL'}] {label:48s} -> {'RAISED' if raised else 'ok'}: {msg[:150]}")
    (HERE / "work/envcheck248.json").write_text(json.dumps(dict(builder=str(BUILDER), rows=rows, all_ok=ok_all), indent=1))
    print(f"\nENVCHECK {'PASS' if ok_all else 'FAIL'} {sum(r['ok'] for r in rows)}/{len(rows)}  builder {BUILDER}")
    raise SystemExit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
