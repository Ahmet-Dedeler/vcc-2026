"""495: register the LODO key at cap 6 (`lodo415_lo0.5_hi6`, the third session's 540) in the builder's PERT_SCALE_ARTIFACTS.
A pure addition: with VCC_PERT_SCALE unset or naming another key, no code path changes. The sha is read from the file at
patch time and must equal the recorded dc374281...; the builder re-verifies it at build time.
usage: uv run python experiments/495-lodo-cap6/patch_builder495.py <builder.py> [--check]
"""
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "experiments/495-lodo-cap6/work/pert_scale_415pp_lo0.5_hi6.npz"
WANT = "dc374281452ed4f7ae3a25b392dbfa4ff0d5cc591bd7b5b71e2c0a2dca86a4c0"
got = hashlib.sha256(ART.read_bytes()).hexdigest()
assert got == WANT, f"artifact sha {got} != {WANT}"
path = Path(sys.argv[1]); src = path.read_text()
old = '''    "lodo415_lo0.5_hi3": (
        ROOT / "experiments/420-bulk-stack/work/pert_scale_415pp_lo0.5_hi3.npz",
        "8c82bef37d0a59b3ed1f4917ccad0d8fd24618a89fd9e2700c058b9aeff4339b",
    ),
}
'''
new = old[:-2] + f'''    # 495 (2026-09-27): the same LODO key at cap 6 (540: the cap of 3 binds on 49 of 294 carried targets and costs).
    # w = clip(key / median(key[carried]), 0.5, 6) on carried targets, 1 elsewhere; the rule at (0.5, 3) reproduces 8c82bef3.
    "lodo415_lo0.5_hi6": (
        ROOT / "experiments/495-lodo-cap6/work/pert_scale_415pp_lo0.5_hi6.npz",
        "{WANT}",
    ),
}}
'''
if src.count(old) != 1:
    raise SystemExit(f"anchor count {src.count(old)} != 1")
if '"lodo415_lo0.5_hi6"' in src:
    raise SystemExit("already registered")
src = src.replace(old, new)
compile(src, str(path), "exec")
if "--check" in sys.argv:
    print("patch applies cleanly (not written)")
else:
    path.write_text(src)
    print(f"patched {path}")
