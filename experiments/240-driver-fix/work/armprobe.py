"""What the BUILDER reads out of the environment it is handed.

Imports `scripts/make_submission_pm_cd4.py` as a module -- it parses `VCC_ALPHA`,
`VCC_CD4_TARGET_F`, `VCC_H1_UNION_G` and `VCC_UNION_SOURCES` at import time and nothing else
runs -- and prints its own parsed values plus the config keys those become.
"""
import importlib.util, json, os, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location("pm_cd4", ROOT / "scripts/make_submission_pm_cd4.py")
mod = importlib.util.module_from_spec(spec)
sys.argv = ["make_submission_pm_cd4.py"]
sys.modules["pm_cd4"] = mod
spec.loader.exec_module(mod)
print(json.dumps({
    "env_seen": {k: os.environ.get(k, "<absent>")
                 for k in ("VCC_ALPHA", "VCC_CD4_TARGET_F", "VCC_H1_UNION_G", "VCC_UNION_SOURCES")},
    "builder_parsed": {"alpha": mod.ALPHA, "cd4_target_f": mod.CD4_TARGET_F,
                       "h1_union_g": mod.H1_UNION_G,
                       "h1_union_sources": [n for n, _ in mod.H1_UNION_SOURCES]},
}, indent=1))
