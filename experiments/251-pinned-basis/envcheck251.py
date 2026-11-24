"""251 -- the import-time rules of `VCC_CD4_DEFLATE_BASIS`, each in its own subprocess (248's
`envcheck248.py` pattern).  Every row is a builder IMPORT under one environment; the expected
outcome is `ok` (imports; the flag's resolved value printed) or `refuse` (ValueError with the
given phrase).  The pre-251 byte copy is imported too, with the flag set, to show it is unknown
there (the 249/248 "flag unknown to the old builder" row).

    uv run python experiments/251-pinned-basis/envcheck251.py
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
NEW = ROOT / "scripts/make_submission_pm_cd4.py"
PRE = WORK / "make_submission_pm_cd4.PRE.py"

PROBE = r'''
import importlib.util, sys, json
from pathlib import Path
path = Path(sys.argv[1])
spec = importlib.util.spec_from_file_location("envcheck_builder", path)
m = importlib.util.module_from_spec(spec)
m.__file__ = str(Path(sys.argv[2]))
sys.modules["envcheck_builder"] = m
try:
    spec.loader.exec_module(m)
except ValueError as exc:
    print(json.dumps({"outcome": "refuse", "error": str(exc)}))
    raise SystemExit(0)
print(json.dumps({"outcome": "ok", "cd4_deflate_k": getattr(m, "CD4_DEFLATE_K", "<no attribute>"),
                  "cd4_deflate_basis": getattr(m, "CD4_DEFLATE_BASIS", "<no attribute>"),
                  "registry": sorted(getattr(m, "CD4_DEFLATION_BASIS_REGISTRY", {}))}))
'''

ROWS = [
    # (label, builder, env, expected outcome, phrase the refusal must contain)
    ("unset", NEW, {}, "ok", None),
    ("BASIS=V5_abc", NEW, {"VCC_CD4_DEFLATE_BASIS": "V5_abc"}, "ok", None),
    ("K=5", NEW, {"VCC_CD4_DEFLATE_K": "5"}, "ok", None),
    ("BASIS blank (' ') reads as unset", NEW, {"VCC_CD4_DEFLATE_BASIS": "  "}, "ok", None),
    ("BOTH K=5 and BASIS=V5_abc", NEW, {"VCC_CD4_DEFLATE_K": "5", "VCC_CD4_DEFLATE_BASIS": "V5_abc"},
     "refuse", "two variables"),
    ("BASIS=nope (unregistered)", NEW, {"VCC_CD4_DEFLATE_BASIS": "nope"}, "refuse",
     "not in CD4_DEFLATION_BASIS_REGISTRY"),
    ("BASIS=V5_ABC (case matters)", NEW, {"VCC_CD4_DEFLATE_BASIS": "V5_ABC"}, "refuse",
     "not in CD4_DEFLATION_BASIS_REGISTRY"),
    ("BOTH, K invalid (K=0)", NEW, {"VCC_CD4_DEFLATE_K": "0", "VCC_CD4_DEFLATE_BASIS": "V5_abc"},
     "refuse", "VCC_CD4_DEFLATE_K must be an integer >= 1"),
    ("pre-251 byte copy, BASIS=V5_abc (unknown to it)", PRE, {"VCC_CD4_DEFLATE_BASIS": "V5_abc"}, "ok", None),
    ("pre-251 byte copy, K=5", PRE, {"VCC_CD4_DEFLATE_K": "5"}, "ok", None),
]


def main() -> None:
    base_env = {k: v for k, v in os.environ.items()
                if not k.startswith("VCC_")}
    base_env["VCC_ALPHA"] = "0.982"
    results, fails = [], 0
    for label, builder, env, want, phrase in ROWS:
        proc = subprocess.run([sys.executable, "-c", PROBE, str(builder), str(NEW)],
                              env={**base_env, **env}, capture_output=True, text=True, cwd=ROOT)
        last = [ln for ln in proc.stdout.splitlines() if ln.startswith("{")]
        if proc.returncode != 0 or not last:
            got = {"outcome": "crash", "error": (proc.stderr or proc.stdout)[-600:]}
        else:
            got = json.loads(last[-1])
        ok = got["outcome"] == want and (phrase is None or phrase in got.get("error", ""))
        if builder == PRE and want == "ok":
            ok = ok and got.get("cd4_deflate_basis") == "<no attribute>"
        if label == "BASIS blank (' ') reads as unset":
            ok = ok and got.get("cd4_deflate_basis") is None
        if label == "BASIS=V5_abc":
            ok = ok and got.get("cd4_deflate_basis") == "V5_abc" and got.get("cd4_deflate_k") is None
        fails += 0 if ok else 1
        results.append({"label": label, "builder": str(builder.relative_to(ROOT)), "env": env,
                        "want": want, "phrase": phrase, "got": got, "ok": ok})
        detail = got.get("error", "")[:110] if got["outcome"] != "ok" else \
            f"K={got.get('cd4_deflate_k')!r} BASIS={got.get('cd4_deflate_basis')!r}"
        print(f"  {'PASS' if ok else 'FAIL'}  {label:52s} -> {got['outcome']:7s} {detail}", flush=True)
    (WORK / "envcheck251.json").write_text(json.dumps(results, indent=1))
    print(f"\nENVCHECK {len(results) - fails}/{len(results)} {'PASS' if not fails else 'FAIL'}")
    raise SystemExit(1 if fails else 0)


if __name__ == "__main__":
    main()
