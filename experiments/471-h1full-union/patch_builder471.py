"""471: register the denoised H1 bank (`data/donors/h1_2025_full`, 458d's delta + 471's b_ctrl) in the builder's
UNION_SOURCE_REGISTRY, so `VCC_UNION_SOURCES=gse295214_h1,h1_2025_full` can replace `h1_2025` in the H1 union. A pure
addition: with VCC_UNION_SOURCES unset (or naming other banks) no code path changes. The shas are read from the files
at patch time and verified by the builder at build time.
usage: uv run python experiments/471-h1full-union/patch_builder471.py <builder.py> [--check]
"""
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BANK = ROOT / "data/donors/h1_2025_full"
DELTA_SHA = "26662c527164ddd434ae19a13076e915d1f82de9683c75e27e0af1f5b28591c0"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 24), b""):
            h.update(c)
    return h.hexdigest()


shas = {n: sha(BANK / n) for n in ("delta.npy", "b_ctrl.npy", "genes.txt", "perts.txt")}
assert shas["delta.npy"] == DELTA_SHA, shas["delta.npy"]
path = Path(sys.argv[1]); src = path.read_text()
old = '''            "perts.txt": "977b808e723c73eea5dc4a60cec6e71512387dc520512bf133d3b18fce7287dd",
        },
    },
'''
new = old + f'''    # --- 471 (2026-09-26): the denoised H1 bank. 458d summed ALL public VCC 2025 H1 cells (300 perturbations, ~1,071
    #     cells each, 38,176 controls); 471 added its control profile b_ctrl. Select it in place of `h1_2025` with
    #     VCC_UNION_SOURCES=gse295214_h1,h1_2025_full (25 panel rows vs 20; 3 new union targets: HSBP1, SHPRH, ZFP62).
    "h1_2025_full": {{
        "path": ROOT / "data/donors/h1_2025_full",
        "sha256": {{
            "delta.npy": "{shas['delta.npy']}",
            "b_ctrl.npy": "{shas['b_ctrl.npy']}",
            "genes.txt": "{shas['genes.txt']}",
            "perts.txt": "{shas['perts.txt']}",
        }},
    }},
'''
if src.count(old) != 1:
    raise SystemExit(f"anchor count {src.count(old)} != 1")
if '"h1_2025_full": {' in src:
    raise SystemExit("already registered")
src = src.replace(old, new)
compile(src, str(path), "exec")
if "--check" in sys.argv:
    print("anchor found; compiles; not writing", shas)
else:
    path.write_text(src)
    print(f"patched {path}: h1_2025_full registered {shas}")
