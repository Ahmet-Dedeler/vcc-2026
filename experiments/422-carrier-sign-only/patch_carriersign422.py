"""422 builder patch: VCC_CARRIER_SIGN_ONLY=1 (requires VCC_CARRIER_GAMMA). The carriers keep the SHIPPED design's
magnitudes |final| and take their signs from the second pass at VCC_CARRIER_GAMMA, falling back to the shipped sign where
that pass is exactly zero. Row norms are unchanged by construction. 417c showed the gamma-0 pass carries better SIGNS
(board fid +0.0012, reach +0.0024) and a worse MAGNITUDE profile for nmae (+0.0089); this keeps the first and drops the
second. Unset is a no-op. Apply only when no build has the builder pinned. Apply AFTER patch_carrier417.py.
usage: uv run python patch_carriersign422.py [--write]      (PATCH_TARGET=<path> to dry-run on a copy)
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


sub('''    if PERT_SCALE_BALLAST is None:
        raise SystemExit("VCC_CARRIER_GAMMA changes the carriers only and needs VCC_PERT_SCALE_BALLAST (welded it moves the bulk)")
''', '''    if PERT_SCALE_BALLAST is None:
        raise SystemExit("VCC_CARRIER_GAMMA changes the carriers only and needs VCC_PERT_SCALE_BALLAST (welded it moves the bulk)")
# `VCC_CARRIER_SIGN_ONLY` (422): with VCC_CARRIER_GAMMA, keep the shipped carrier MAGNITUDES and take only the SIGNS from the
# second pass (fallback: the shipped sign where the pass is exactly 0). Unset is a no-op.
CARRIER_SIGN_ONLY = os.environ.get("VCC_CARRIER_SIGN_ONLY", "").strip() or None
if CARRIER_SIGN_ONLY not in (None, "1"):
    raise SystemExit(f"VCC_CARRIER_SIGN_ONLY must be 1 when set: {CARRIER_SIGN_ONLY!r}")
if CARRIER_SIGN_ONLY is not None and CARRIER_GAMMA is None:
    raise SystemExit("VCC_CARRIER_SIGN_ONLY needs VCC_CARRIER_GAMMA (the pass whose signs it takes)")
''', "knob")

sub('''        n0 = np.linalg.norm(final, axis=1)
        n1 = np.linalg.norm(fc, axis=1)
        scale_rows = np.where(n1 > 0, n0 / np.maximum(n1, 1e-300), 0.0)
        carrier_final = fc * scale_rows[:, None]''', '''        n0 = np.linalg.norm(final, axis=1)
        n1 = np.linalg.norm(fc, axis=1)
        scale_rows = np.where(n1 > 0, n0 / np.maximum(n1, 1e-300), 0.0)
        carrier_final = fc * scale_rows[:, None]
        if CARRIER_SIGN_ONLY is not None:
            # 422: shipped magnitudes, second-pass signs (shipped sign where the pass is exactly 0).
            carrier_final = np.abs(final) * np.where(fc != 0, np.sign(fc), np.sign(final))''', "sign_only")

sub('''        carrier_audit = {"gamma": CARRIER_GAMMA, "whitening": wc_audit,''', '''        flips = (np.sign(carrier_final) != np.sign(final)) & (final != 0)
        carrier_audit = {"gamma": CARRIER_GAMMA, "whitening": wc_audit, "sign_only": CARRIER_SIGN_ONLY is not None,
                         "frac_entries_sign_changed": float(flips.mean()),
                         "abs_weighted_frac_sign_changed": float((np.abs(final) * flips).sum() / max(float(np.abs(final).sum()), 1e-300)),''', "audit")

sub('''            **({} if CARRIER_GAMMA is None else {"carrier_gamma": CARRIER_GAMMA}),
''', '''            **({} if CARRIER_GAMMA is None else {"carrier_gamma": CARRIER_GAMMA}),
            **({} if CARRIER_SIGN_ONLY is None else {"carrier_sign_only": True}),
''', "config")

sys.stdout.writelines(difflib.unified_diff(src.splitlines(True), out.splitlines(True), "a/" + P.name, "b/" + P.name, n=1))
if "--write" in sys.argv:
    P.write_text(out)
    print("\nWROTE", P)
