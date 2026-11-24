"""437: the builder half of 432's trap T6. Under VCC_CD4_SIZING=carrying, CD4 is sized by VCC_CD4_TARGET_G and
VCC_CD4_TARGET_F is only RECORDED. So 417's stack-dump recipe (`VCC_CD4_TARGET_F=0.51` for CD4 x1.5) silently does
nothing on D/E/F, and every check passes. Refuse a TARGET_F other than the driver's shipped 0.34 under carrying sizing.
A/B/C (sizing unset) and the driver's own D/E/F default (F 0.34) are untouched, so every existing build is bit-identical.
usage: uv run python patch_t6.py [--builder PATH]"""
import argparse
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ANCHOR = '''if CD4_TARGET_G is not None and not 0.0 < CD4_TARGET_G < 5.0:
    raise ValueError(f"VCC_CD4_TARGET_G must be in (0, 5); got {CD4_TARGET_G}")
'''
ADD = '''if CD4_SIZING == "carrying" and CD4_TARGET_F is not None and CD4_TARGET_F != 0.34:
    # 432 T6 (2026-09-23): under carrying sizing CD4 is sized by VCC_CD4_TARGET_G and VCC_CD4_TARGET_F is only
    # recorded, so a dose change asked for through TARGET_F (417's `f=0.51` stack-dump recipe) would silently do
    # nothing on a new panel. 0.34 is the driver's shipped default and stays accepted (recorded, inert).
    raise ValueError(
        f"VCC_CD4_TARGET_F={CD4_TARGET_F} under VCC_CD4_SIZING=carrying would only be recorded; "
        "set VCC_CD4_TARGET_G (e.g. 1.5 x 0.3443655278649641 for CD4 x1.5) instead (432 T6)"
    )
'''

ap = argparse.ArgumentParser()
ap.add_argument("--builder", default=str(ROOT / "scripts/make_submission_pm_cd4.py"))
a = ap.parse_args()
p = Path(a.builder); s = p.read_text()
if "432 T6" in s:
    raise SystemExit(f"{p} already carries the T6 guard")
if s.count(ANCHOR) != 1:
    raise SystemExit(f"anchor found {s.count(ANCHOR)} times")
p.write_text(s.replace(ANCHOR, ANCHOR + ADD))
print(f"{p}: patched, sha {hashlib.sha256(p.read_bytes()).hexdigest()}")
