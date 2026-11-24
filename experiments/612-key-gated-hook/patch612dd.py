#!/usr/bin/env python3
"""612dd -- the def_day.sh pass-through for 612's key gate: DEF_DAY_MDE_SIGN_KEYGATE=<q> -> VCC_DE_MDE_SIGN_KEYGATE.
UNAPPLIED. Order: 598 -> 612 -> 612dd (after the 10-05 bundle). Ownership agreed with [4d0207], who composes it into the
bundle -> 641 -> 642 -> 604 chain test.

ONE file, scripts/def_day.sh (the POST-598 file), six anchored insertions. Unset knob = no export, no check, no meta key,
so today's command is unchanged.
  DOC  the knob documented next to DEF_DAY_MDE_SIGN;
  V    validation after the DEF_DAY_MDE_SIGN case: needs DEF_DAY_MDE_SIGN, 0 < q < 1, and the builder carries
       VCC_DE_MDE_SIGN_KEYGATE (612);
  E    export VCC_DE_MDE_SIGN_KEYGATE for the level builds, next to VCC_DE_MDE_SIGN's export;
  C    the final meta must record config.de_mde_sign_keygate.q == the requested q (an inline check, like
       VCC_RC_CALLSET_BLOCK's; riders_match's interface is untouched);
  F    def_day.final.env records it;
  A    postbuild's allowed config diff gains de_mde_sign_keygate when set.
The identity hook itself needs no def_day edit: [4d0207]'s patch604 serves it as DEF_DAY_MDE_SIGN=hookbank:identity:0.

    uv run python experiments/612-key-gated-hook/patch612dd.py --root <repo> --check
    uv run python experiments/612-key-gated-hook/patch612dd.py --root <repo>
"""
from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path

MARK = "612dd: key gate pass-through"
MARK598 = "598: tercilec"
MARK612 = "612: key gate"
DAY = Path("scripts/def_day.sh")
BUILDER = Path("scripts/make_submission_pm_cd4.py")

DOC_OLD = '''#    DEF_DAY_HOOK_ESTIMATOR=meanz_ne (517, default pdown): the hook artifact's estimator (def_sameline.py >= 517's patch).
'''
DOC_NEW = DOC_OLD + '''#    DEF_DAY_MDE_SIGN_KEYGATE=<q> (612dd: key gate pass-through; default off; e.g. 0.5): 612's key gate on the hook's sign
#    (VCC_DE_MDE_SIGN_KEYGATE): the hook's sign only on each target's weak-key genes. Needs DEF_DAY_MDE_SIGN. Builder >= 612.
'''
V_OLD = '''    *) die "DEF_DAY_MDE_SIGN must be clean4:<tau>, programme:<tau>[:<tau_c4>] or hookbank:<bank>:<tau>[:<tau_c4>[:le<N>]] (got $DEF_DAY_MDE_SIGN)" ;;
  esac
fi
'''
V_NEW = V_OLD + '''# 612dd: key gate pass-through -- DEF_DAY_MDE_SIGN_KEYGATE=<q> (validated here, exported with the other DE riders)
if [ -n "${DEF_DAY_MDE_SIGN_KEYGATE:-}" ]; then
  [ -n "${DEF_DAY_MDE_SIGN:-}" ] || die "DEF_DAY_MDE_SIGN_KEYGATE gates the hook's sign: it needs DEF_DAY_MDE_SIGN"
  python3 -c 'import sys; q=float(sys.argv[1]); sys.exit(0 if 0.0 < q < 1.0 else 1)' "$DEF_DAY_MDE_SIGN_KEYGATE" 2>/dev/null \\
    || die "DEF_DAY_MDE_SIGN_KEYGATE must be a quantile 0 < q < 1 (got $DEF_DAY_MDE_SIGN_KEYGATE)"
  grep -q "VCC_DE_MDE_SIGN_KEYGATE" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_DE_MDE_SIGN_KEYGATE (612)"
fi
'''
E_OLD = '''[ -n "${DEF_DAY_MDE_SIGN:-}" ] && { export VCC_DE_MDE_SIGN="$MDS"; log "DE rider for the level builds: VCC_DE_MDE_SIGN=$VCC_DE_MDE_SIGN"; }
'''
E_NEW = E_OLD + '''[ -n "${DEF_DAY_MDE_SIGN_KEYGATE:-}" ] && { export VCC_DE_MDE_SIGN_KEYGATE="$DEF_DAY_MDE_SIGN_KEYGATE"; log "DE rider for the level builds: VCC_DE_MDE_SIGN_KEYGATE=$VCC_DE_MDE_SIGN_KEYGATE"; }   # 612dd: key gate pass-through
'''
C_OLD = '''else
  uv run python scripts/def_day_helpers.py riders_match "$FINAL.h5ad.meta.json" "${VCC_DE_MAG:--}" "${VCC_DE_MDE_SIGN:--}" \\
    || die "the final meta does not carry the requested DE riders (VCC_DE_MAG=${VCC_DE_MAG:-} VCC_DE_MDE_SIGN=${VCC_DE_MDE_SIGN:-})"
fi
'''
C_NEW = C_OLD + '''if [ -n "${VCC_DE_MDE_SIGN_KEYGATE:-}" ]; then   # 612dd: key gate pass-through -- the final meta must carry the gate it was built with
  python3 - "$FINAL.h5ad.meta.json" "$VCC_DE_MDE_SIGN_KEYGATE" <<'PY' || die "the final meta does not carry VCC_DE_MDE_SIGN_KEYGATE=$VCC_DE_MDE_SIGN_KEYGATE"
import json, sys
g = json.load(open(sys.argv[1]))["config"].get("de_mde_sign_keygate") or {}
sys.exit(0 if abs(float(g.get("q", -1.0)) - float(sys.argv[2])) < 1e-12 else 1)
PY
fi
'''
F_OLD = '''[ -n "${VCC_DE_MDE_SIGN:-}" ] && printf 'VCC_DE_MDE_SIGN=%s\\n' "$VCC_DE_MDE_SIGN" >> "$W/def_day.final.env"
'''
F_NEW = F_OLD + '''[ -n "${VCC_DE_MDE_SIGN_KEYGATE:-}" ] && printf 'VCC_DE_MDE_SIGN_KEYGATE=%s\\n' "$VCC_DE_MDE_SIGN_KEYGATE" >> "$W/def_day.final.env"   # 612dd: key gate pass-through
'''
A_OLD = '''[ -n "${VCC_DE_MDE_SIGN:-}" ] && ALLOW="$ALLOW,de_mde_sign"
'''
A_NEW = A_OLD + '''[ -n "${VCC_DE_MDE_SIGN_KEYGATE:-}" ] && ALLOW="$ALLOW,de_mde_sign_keygate"   # 612dd: key gate pass-through
'''
EDITS = [("DOC", DOC_OLD, DOC_NEW), ("V validation", V_OLD, V_NEW), ("E export", E_OLD, E_NEW),
         ("C meta check", C_OLD, C_NEW), ("F final env", F_OLD, F_NEW), ("A allow", A_OLD, A_NEW)]


def patched_text(day: str, builder: str) -> str:
    if MARK in day:
        raise SystemExit(f"{DAY}: already carries {MARK!r} (a second apply is refused)")
    if MARK598 not in day:
        raise SystemExit(f"{DAY}: lacks 598's mark: apply 598 first (538 -> 554 -> 590 -> 598 -> 612 -> 612dd)")
    if MARK612 not in builder:
        raise SystemExit(f"{BUILDER}: lacks 612's mark {MARK612!r}: apply 612 first")
    for name, old, _ in EDITS:
        n = day.count(old)
        if n != 1:
            raise SystemExit(f"{DAY}: anchor {name} found {n} times (expected exactly 1)")
    for _, old, new in EDITS:
        day = day.replace(old, new, 1)
    return day


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    p = a.root / DAY
    new = patched_text(p.read_text(), (a.root / BUILDER).read_text())
    if a.check:
        print(f"612dd check OK: {DAY} carries 598's mark, the builder carries 612's, all {len(EDITS)} anchors found once")
        return 0
    fd, tmp = tempfile.mkstemp(dir=p.parent, prefix=".patch612dd-")
    with os.fdopen(fd, "w") as f:
        f.write(new)
    os.chmod(tmp, os.stat(p).st_mode)
    os.replace(tmp, p)
    print(f"612dd applied: {DAY} ({len(EDITS)} insertions)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
