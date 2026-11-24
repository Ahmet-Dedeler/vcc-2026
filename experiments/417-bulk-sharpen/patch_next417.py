"""417 builder window 2: (a) 415's production leave-one-donor-out key as drop-in VCC_PERT_SCALE artifacts; (b) the bulk-only
per-context scale VCC_BULK_CTX_SCALE ("A=0.97,B=1.01,C=1.03") for 414 section 6's cliff re-levelling.

(b) multiplies each context's BULK effect row (the ballast render only; the carriers are untouched) before rendering.
Sharpening restores each row's own norm, so the scale commutes with VCC_BULK_SHARPEN_Q. Under ballast it touches only
pds/mse, and pds is scale-invariant apart from the render, so the knob is an mse-only amplitude lever. Unset is a no-op.
Apply only when no build has the builder pinned. Apply AFTER patch_sharpen417.py.
usage: uv run python patch_next417.py [--write]      (PATCH_TARGET=<path> to dry-run on a copy)
"""
import difflib
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
P = Path(os.environ.get("PATCH_TARGET", ROOT / "scripts/make_submission_pm_cd4.py"))
src = P.read_text()
out = src


def sub(old, new, label):
    global out
    n = out.count(old)
    if n != 1:
        raise SystemExit(f"[{label}] expected exactly one match, found {n}")
    out = out.replace(old, new)


sub('''    "wshare_lo0.3_inf": (
        ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.3_hi1e+09.npz",
        "94c048e99aecd40004e70f559fe26bdd3b8b9aef9cb8c8cd0efbf8a72d2fe69e",
    ),
}''', '''    "wshare_lo0.3_inf": (
        ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.3_hi1e+09.npz",
        "94c048e99aecd40004e70f559fe26bdd3b8b9aef9cb8c8cd0efbf8a72d2fe69e",
    ),
    # 415: the production leave-one-donor-out key, key_p = mean_j <P^{-j}_p, d_j>/|P_p|^2 on the shipped composition
    # (base + CD4 + H1 residual); beats 410's key on 3/3 cancer folds and kolf2 (+1.4..+5.9 % rel cos at clip 3).
    "lodo415_lo0.3_hi3": (
        ROOT / "experiments/415-lodo-key/work/pert_scale_415pp_lo0.3_hi3.npz",
        "5f1ae21f530734735b23bc0a114bb7b84988f0b62221c6a11c7a64d8176547d9",
    ),
    "lodo415_lo0.3_inf": (
        ROOT / "experiments/415-lodo-key/work/pert_scale_415pp_lo0.3_hi1e+09.npz",
        "5a12079f83b607bb93e14ea30aa4a8443c56deb426d26a0c7f5854a019296dfd",
    ),
    # 414 section 6/7: 410's key with the floor raised to 0.5 (same key; w = max(w_lo0.3, 0.5)); removes most of the key's
    # board pds cost (board-emulated +0.0028 pds, +0.0012 overall).
    "wshare_lo0.5_hi3": (
        ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.5_hi3.npz",
        "387d6f1a76b866a314c633e6f879329036bee5a1f7d5bcadcad1e86cd04f1574",
    ),
    "wshare_lo0.5_inf": (
        ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.5_hi1e+09.npz",
        "851a689de20dbd3794e05191efce356ea35baab4b58daaba75c14672c8363987",
    ),
}''', "lodo_artifacts")

sub('''BULK_SHARPEN_Q = os.environ.get("VCC_BULK_SHARPEN_Q", "").strip() or None
''', '''BULK_SHARPEN_Q = os.environ.get("VCC_BULK_SHARPEN_Q", "").strip() or None
# `VCC_BULK_CTX_SCALE` (417, 414 section 6): one multiplier per context on the BULK effect row only (ballast render), e.g.
# "A=0.978,B=1.0136,C=1.0386" to put each context's bulk on its mse cliff. Requires VCC_PERT_SCALE_BALLAST. Unset is a no-op.
BULK_CTX_SCALE = os.environ.get("VCC_BULK_CTX_SCALE", "").strip() or None
if BULK_CTX_SCALE is not None:
    _pairs = [p.split("=") for p in BULK_CTX_SCALE.split(",")]
    if any(len(p) != 2 for p in _pairs):
        raise SystemExit(f"VCC_BULK_CTX_SCALE must look like 'A=0.98,B=1.01,C=1.04': {BULK_CTX_SCALE!r}")
    BULK_CTX_SCALE = {k.strip(): float(v) for k, v in _pairs}
    if not all(0.8 <= v <= 1.2 for v in BULK_CTX_SCALE.values()):
        raise SystemExit(f"VCC_BULK_CTX_SCALE factors must be in [0.8, 1.2]: {BULK_CTX_SCALE!r}")
''', "ctx_knob")

sub('''    if PERT_SCALE_BALLAST is None:
        raise SystemExit("VCC_BULK_SHARPEN_Q is a bulk-only arm and needs VCC_PERT_SCALE_BALLAST (welded it reaches the DE channel)")
''', '''    if PERT_SCALE_BALLAST is None:
        raise SystemExit("VCC_BULK_SHARPEN_Q is a bulk-only arm and needs VCC_PERT_SCALE_BALLAST (welded it reaches the DE channel)")
if BULK_CTX_SCALE is not None and PERT_SCALE_BALLAST is None:
    raise SystemExit("VCC_BULK_CTX_SCALE is a bulk-only arm and needs VCC_PERT_SCALE_BALLAST (welded it reaches the DE channel)")
''', "ctx_refuse")

sub('''    if BULK_SHARPEN_Q is not None:
        if bulk_effects is None:
            raise ValueError("VCC_BULK_SHARPEN_Q needs the ballast path (bulk_effects)")
''', '''    if BULK_SHARPEN_Q is not None:
        if bulk_effects is None:
            raise ValueError("VCC_BULK_SHARPEN_Q needs the ballast path (bulk_effects)")
    if BULK_CTX_SCALE is not None:
        if bulk_effects is None:
            raise ValueError("VCC_BULK_CTX_SCALE needs the ballast path (bulk_effects)")
        if set(BULK_CTX_SCALE) != set(map(str, spec.contexts)):
            raise ValueError(f"VCC_BULK_CTX_SCALE contexts {sorted(BULK_CTX_SCALE)} != bundle {sorted(map(str, spec.contexts))}")
        print(f"[bulk-ctx-scale] {BULK_CTX_SCALE} on the ballast bulk render only", flush=True)
''', "ctx_check")

sub('''                                bulk_effects[target_index],
''', '''                                bulk_effects[target_index] * (1.0 if BULK_CTX_SCALE is None else BULK_CTX_SCALE[str(context)]),
''', "ctx_apply")

sub('''            **({} if BULK_SHARPEN_Q is None else {"bulk_sharpen_q": BULK_SHARPEN_Q}),
''', '''            **({} if BULK_SHARPEN_Q is None else {"bulk_sharpen_q": BULK_SHARPEN_Q}),
            **({} if BULK_CTX_SCALE is None else {"bulk_ctx_scale": BULK_CTX_SCALE}),
''', "config")

sys.stdout.writelines(difflib.unified_diff(src.splitlines(True), out.splitlines(True), "a/" + P.name, "b/" + P.name, n=1))
if "--write" in sys.argv:
    P.write_text(out)
    print("\nWROTE", P)
