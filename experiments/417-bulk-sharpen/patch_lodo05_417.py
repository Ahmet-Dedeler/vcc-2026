"""417 window 3: register 420's LODO key at clip(0.5, 3) (8c82bef3; 415's V1_PP key re-clipped) as VCC_PERT_SCALE
`lodo415_lo0.5_hi3`. 420 priced it best among key variants on the sharpened bulk: +0.0031 (k 1) / +0.0038 (k 1.407) vs
417a. Apply only when no build has the builder pinned. Apply AFTER the window-2 patches.
usage: uv run python patch_lodo05_417.py [--write]      (PATCH_TARGET=<path> to dry-run on a copy)
"""
import difflib
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
P = Path(os.environ.get("PATCH_TARGET", ROOT / "scripts/make_submission_pm_cd4.py"))
src = P.read_text()
out = src
old = '''    "wshare_lo0.5_inf": (
        ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.5_hi1e+09.npz",
        "851a689de20dbd3794e05191efce356ea35baab4b58daaba75c14672c8363987",
    ),
}'''
new = '''    "wshare_lo0.5_inf": (
        ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.5_hi1e+09.npz",
        "851a689de20dbd3794e05191efce356ea35baab4b58daaba75c14672c8363987",
    ),
    # 420: 415's LODO V1_PP key re-clipped at (0.5, 3); best key variant on the sharpened bulk (+0.0031/+0.0038 vs 417a).
    "lodo415_lo0.5_hi3": (
        ROOT / "experiments/420-bulk-stack/work/pert_scale_415pp_lo0.5_hi3.npz",
        "8c82bef37d0a59b3ed1f4917ccad0d8fd24618a89fd9e2700c058b9aeff4339b",
    ),
}'''
if out.count(old) != 1:
    raise SystemExit(f"[lodo05] expected exactly one match, found {out.count(old)}")
out = out.replace(old, new)
sys.stdout.writelines(difflib.unified_diff(src.splitlines(True), out.splitlines(True), "a/" + P.name, "b/" + P.name, n=1))
if "--write" in sys.argv:
    P.write_text(out)
    print("\nWROTE", P)
