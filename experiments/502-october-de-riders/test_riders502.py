"""502: riders_match on synthetic metas. usage: uv run python experiments/502-october-de-riders/test_riders502.py <def_day_helpers.py>"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

H = sys.argv[1]
zs = {"variant": "dn_k250_zs", "module_sha256": "x"}
hk = {"source": "clean4", "tau": 0.05}
cases = [("neither requested, neither present", {}, "-", "-", 0),
         ("zs requested and present", {"de_mag": zs}, "dn_k250_zs", "-", 0),
         ("zs requested, absent", {}, "dn_k250_zs", "-", 1),
         ("zs present, not requested", {"de_mag": zs}, "-", "-", 1),
         ("hook requested and present", {"de_mde_sign": hk}, "-", "clean4:0.05", 0),
         ("hook tau differs", {"de_mde_sign": hk}, "-", "clean4:0.1", 1),
         ("hook present, not requested", {"de_mde_sign": hk}, "-", "-", 1),
         ("both requested and present", {"de_mag": zs, "de_mde_sign": hk}, "dn_k250_zs", "clean4:0.05", 0),
         ("wrong zs variant", {"de_mag": zs}, "up_k140_zs", "-", 1)]
bad = 0
with tempfile.TemporaryDirectory() as d:
    for name, cfg, mag, sign, want in cases:
        p = Path(d) / "m.json"; p.write_text(json.dumps({"config": cfg}))
        rc = subprocess.run([sys.executable, H, "riders_match", str(p), mag, sign], capture_output=True, text=True).returncode
        ok = rc == want; bad += not ok
        print(f"{'PASS' if ok else 'FAIL'}  {name}: rc {rc} (want {want})")
print("ALL PASS" if not bad else f"{bad} FAILED"); sys.exit(1 if bad else 0)
