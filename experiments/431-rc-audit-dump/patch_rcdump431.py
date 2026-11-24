"""431: optional per-block rc diagnostics dump (`VCC_RC_AUDIT_DUMP=<path.jsonl>`). The meta keeps only median/min/max over
the 900 blocks; this writes one JSON line per block (context, target and every numeric rc diagnostic), so a dose question
like 426a's A/B/C `mde_lambda` median 0.84 can be answered per context and per block. Unset is a no-op (no code runs).
It does not touch the emission, so a knob-on build's X is bit-identical to the knob-off build's.
usage: uv run python patch_rcdump431.py [--builder PATH]"""
import argparse
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def rep(s, old, new, what):
    n = s.count(old)
    if n != 1:
        raise SystemExit(f"{what}: anchor found {n} times (need exactly 1)")
    return s.replace(old, new)


ap = argparse.ArgumentParser(); ap.add_argument("--builder", default=str(ROOT / "scripts/make_submission_pm_cd4.py"))
a = ap.parse_args()
b = Path(a.builder); s = b.read_text()
if "VCC_RC_AUDIT_DUMP" in s:
    raise SystemExit("already carries the 431 patch")
s = rep(s, "BALLAST_LB = 995_000\n", '''BALLAST_LB = 995_000
# `VCC_RC_AUDIT_DUMP=<path.jsonl>` (431): one JSON line of rc diagnostics per block (context, target, numeric values);
# diagnostics only, the emission is untouched. Unset is a no-op.
RC_AUDIT_DUMP = os.environ.get("VCC_RC_AUDIT_DUMP", "").strip() or None
''', "knob")
s = rep(s, "                    rc_audit.append(rc_diag)\n", '''                    rc_audit.append(rc_diag)
                    if RC_AUDIT_DUMP is not None:
                        import json as _json431  # noqa: PLC0415
                        with open(RC_AUDIT_DUMP, "a") as _fh431:
                            _fh431.write(_json431.dumps({"context": str(context), "target": str(target),
                                                         **{k: float(v) for k, v in rc_diag.items()
                                                            if isinstance(v, (int, float, np.integer, np.floating))}}) + "\\n")
''', "dump")
b.write_text(s)
print(f"builder {b}: patched, sha {hashlib.sha256(b.read_bytes()).hexdigest()}")
