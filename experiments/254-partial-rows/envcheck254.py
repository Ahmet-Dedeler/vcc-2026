"""254 -- the import-time rules for `VCC_CD4_PARTIAL_G`, each in a fresh subprocess.

    uv run python experiments/254-partial-rows/envcheck254.py
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WT = ROOT / "scripts/make_submission_pm_cd4.py"
PRE = WORK / "make_submission_pm_cd4.PRE.py"

PROBE = r'''
import importlib.util, sys
spec = importlib.util.spec_from_file_location("envprobe_builder", sys.argv[1])
m = importlib.util.module_from_spec(spec); m.__file__ = sys.argv[2]; sys.modules["envprobe_builder"] = m
spec.loader.exec_module(m)
print("OK", repr(getattr(m, "CD4_PARTIAL_G", "<no-attribute>")), repr(getattr(m, "CD4_DEFLATE_K", None)), repr(getattr(m, "CD4_DEFLATE_BASIS", None)))
'''

CASES = [
    # (label, builder, env overrides, expect_ok, expected CD4_PARTIAL_G repr or None)
    ("unset", WT, {}, True, "None"),
    ("blank string", WT, {"VCC_CD4_PARTIAL_G": ""}, True, "None"),
    ("0.34", WT, {"VCC_CD4_PARTIAL_G": "0.34"}, True, "0.34"),
    ("0.34 with K=5 (the composed arm)", WT, {"VCC_CD4_PARTIAL_G": "0.34", "VCC_CD4_DEFLATE_K": "5"}, True, "0.34"),
    ("0.34 with the pinned basis (October composition)", WT,
     {"VCC_CD4_PARTIAL_G": "0.34", "VCC_CD4_DEFLATE_BASIS": "V5_abc"}, True, "0.34"),
    ("0 refused", WT, {"VCC_CD4_PARTIAL_G": "0"}, False, None),
    ("5 refused", WT, {"VCC_CD4_PARTIAL_G": "5"}, False, None),
    ("-0.3 refused", WT, {"VCC_CD4_PARTIAL_G": "-0.3"}, False, None),
    ("abc refused", WT, {"VCC_CD4_PARTIAL_G": "abc"}, False, None),
    ("K=5 alone still imports (249's path untouched)", WT, {"VCC_CD4_DEFLATE_K": "5"}, True, "None"),
    ("pre-254 byte copy with the flag set: imports, no attribute", PRE,
     {"VCC_CD4_PARTIAL_G": "0.34"}, True, "'<no-attribute>'"),
]


def main() -> None:
    base_env = {k: v for k, v in os.environ.items() if not k.startswith("VCC_")}
    base_env.update({"VCC_ALPHA": "0.982", "VCC_CD4_TARGET_F": "0.34", "VCC_H1_UNION_G": "0.28"})
    results = []
    n_ok = 0
    for label, builder, env, expect_ok, expect_g in CASES:
        run_env = {**base_env, **env}
        p = subprocess.run([sys.executable, "-c", PROBE, str(builder), str(WT)], env=run_env,
                           capture_output=True, text=True, cwd=str(ROOT))
        ok_import = p.returncode == 0 and p.stdout.startswith("OK")
        got_g = p.stdout.split()[1] if ok_import else None
        passed = (ok_import == expect_ok) and (not expect_ok or got_g == expect_g)
        err = p.stderr.strip().splitlines()[-1] if p.stderr.strip() else ""
        n_ok += passed
        results.append(dict(label=label, expect_ok=expect_ok, imported=ok_import, got=got_g,
                            expected=expect_g, passed=passed, stderr_last=err[:200]))
        print(f"[{'PASS' if passed else 'FAIL'}] {label:60s} imported={ok_import} got={got_g} "
              f"{'' if ok_import else err[:110]}")
    (WORK / "envcheck254.json").write_text(json.dumps(results, indent=1))
    print(f"\nENVCHECK {n_ok}/{len(CASES)}")
    raise SystemExit(0 if n_ok == len(CASES) else 1)


if __name__ == "__main__":
    main()
