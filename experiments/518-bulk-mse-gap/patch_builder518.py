"""518 (NOT APPLIED; a contingency for the coordinator): register the strength-key pert-scale artifacts in the builder's
PERT_SCALE_ARTIFACTS, exactly as 495 registered cap 6. A pure addition: with VCC_PERT_SCALE unset or naming another key no
code path changes. Anchored after 495's "lodo415_lo0.5_hi6" entry (the builder as of 2026-09-28). The shas are read from
the files at patch time and must equal the recorded ones; the builder re-verifies them at build time.
  strength518_S1_lo0.5_hi3 : clip(key/med, 0.5, 3) x exp(1.25 z)  (stand-alone; the alternative to 495's cap 6)
  strength518_S2_lo0.5_hi6 : clip(key/med, 0.5, 6) x exp(0.5 z)   (rider on cap 6; priced +0.0001 ... +0.0003 over it)
usage: uv run python experiments/518-bulk-mse-gap/patch_builder518.py <builder.py> [--check]
"""
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ARTS = {
    "strength518_S1_lo0.5_hi3": (ROOT / "experiments/518-bulk-mse-gap/work/pert_scale_518_S1_strengthmeta_th1.25_lo0.5_hi3.npz",
                                 "1443c1a41b285872769f20a99ee958d06266b0171e99328de22f13a4cdd00da6"),
    "strength518_S2_lo0.5_hi6": (ROOT / "experiments/518-bulk-mse-gap/work/pert_scale_518_S2_strengthmeta_th0.5_lo0.5_hi6.npz",
                                 "dc0856e07a4b826dcd7bd3bdc27a2fecef12133c539df0bdf999c0c542fff01a"),
}
for name, (p, want) in ARTS.items():
    got = hashlib.sha256(p.read_bytes()).hexdigest()
    assert got == want, f"{name}: artifact sha {got} != {want}"
path = Path(sys.argv[1])
src = path.read_text()
old = '''    "lodo415_lo0.5_hi6": (
        ROOT / "experiments/495-lodo-cap6/work/pert_scale_415pp_lo0.5_hi6.npz",
        "dc374281452ed4f7ae3a25b392dbfa4ff0d5cc591bd7b5b71e2c0a2dca86a4c0",
    ),
}
'''
add = ""
for name, (p, want) in ARTS.items():
    add += (f'    # 518 (2026-09-28): 8c82bef3\'s LODO key x exp(theta z), z = the robust response-strength rank (518 RESULTS).\n'
            f'    "{name}": (\n        ROOT / "{p.relative_to(ROOT)}",\n        "{want}",\n    ),\n')
new = old[:-2] + add + "}\n"
if src.count(old) != 1:
    raise SystemExit(f"anchor count {src.count(old)} != 1")
if '"strength518_S1_lo0.5_hi3"' in src:
    raise SystemExit("already registered")
src = src.replace(old, new)
compile(src, str(path), "exec")
if "--check" in sys.argv:
    print("patch applies cleanly (not written)")
else:
    path.write_text(src)
    print(f"patched {path}")
