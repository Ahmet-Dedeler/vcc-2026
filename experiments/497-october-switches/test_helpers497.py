"""497: meta_matches with the scalar and the per-target sharpening forms (unit test on synthetic metas).
usage: uv run python experiments/497-october-switches/test_helpers497.py <def_day_helpers.py>"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

H = sys.argv[1]
F = "D=1.1,E=1.2,F=1.3"
SC = {"D": 1.1, "E": 1.2, "F": 1.3}
cases = [
    ("scalar ok", {"bulk_sharpen_q": 1.75, "bulk_ctx_scale": SC}, "1.75", 0),
    ("scalar wrong q", {"bulk_sharpen_q": 1.5, "bulk_ctx_scale": SC}, "1.75", 1),
    ("scalar wrong scale", {"bulk_sharpen_q": 1.75, "bulk_ctx_scale": {"D": 1.1, "E": 1.2, "F": 1.4}}, "1.75", 1),
    ("pt ok", {"bulk_sharpen_q": None, "bulk_sharpen_pt": {"mode": "stercile", "q": [2.5, 2.0, 1.75], "counts": [100, 100, 100]}, "bulk_ctx_scale": SC}, "stercile:2.5,2.0,1.75", 0),
    ("pt wrong mode", {"bulk_sharpen_q": None, "bulk_sharpen_pt": {"mode": "stercile", "q": [2.5, 2.0, 1.75]}, "bulk_ctx_scale": SC}, "tercile:2.5,2.0,1.75", 1),
    ("pt wrong q", {"bulk_sharpen_q": None, "bulk_sharpen_pt": {"mode": "stercile", "q": [2.5, 2.0, 1.5]}, "bulk_ctx_scale": SC}, "stercile:2.5,2.0,1.75", 1),
    ("pt meta, scalar arg", {"bulk_sharpen_q": None, "bulk_sharpen_pt": {"mode": "stercile", "q": [2.5, 2.0, 1.75]}, "bulk_ctx_scale": SC}, "1.75", 1),
    ("scalar meta, pt arg", {"bulk_sharpen_q": 1.75, "bulk_ctx_scale": SC}, "stercile:2.5,2.0,1.75", 1),
    ("pt ok, no scale at unit", {"bulk_sharpen_q": None, "bulk_sharpen_pt": {"mode": "tercile", "q": [2.5, 2.0, 1.75]}}, "tercile:2.5,2.0,1.75", None),
]
bad = 0
with tempfile.TemporaryDirectory() as d:
    for name, cfg, qarg, want in cases:
        p = Path(d) / "m.json"; p.write_text(json.dumps({"config": cfg}))
        fact = "D=1,E=1,F=1" if want is None else F
        rc = subprocess.run([sys.executable, H, "meta_matches", str(p), fact, qarg], capture_output=True, text=True).returncode
        exp = 0 if want is None else want
        ok = rc == exp
        bad += not ok
        print(f"{'PASS' if ok else 'FAIL'}  {name}: rc {rc} (want {exp})")
print("ALL PASS" if not bad else f"{bad} FAILED")
sys.exit(1 if bad else 0)
