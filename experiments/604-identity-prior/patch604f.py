#!/usr/bin/env python3
"""604f -- the DAY-OF pass-through for 625's octile flip (604's flip option; [7772dc]'s module readout / D/E/F decode).
UNAPPLIED. It applies AFTER the whole 10-06 chain (641 -> 598 -> 612 -> 612dd -> 630 -> 630dd -> 642 -> 604 -> 644), so its
earliest apply is 10-08 (the 10-06 chain froze 10-05 18:00Z). Two files:

scripts/def_sameline.py
  --flip-octiles <8 x 0/1> (octile order 0..7, 625's `octile` array): passed to def_identity_prior.identity_prior(...,
  flip_octiles=...) in 604's identity branch. Refused without an identity programme source (`identity`, `identity:<lam>`,
  `identity@<lam>`): a flip of a BANK programme is not defined. Unset (or 00000000) = today's artifact, byte for byte.
scripts/def_day.sh
  DEF_DAY_DE_PRIOR_FLIP=<8 x 0/1> (default off):
    V  validated: 8 characters of 0/1; it must reach the prior the recipe actually uses, so it needs an identity
       programme source (DEF_DAY_PROGRAMME_SOURCE=identity...) or an identity hook bank (DEF_DAY_MDE_SIGN=
       hookbank:identity...); and a CLEAN4 hook (DEF_DAY_MDE_SIGN=clean4:...) is REFUSED with a flip, because that hook
       reads CLEAN4 directly and would not see the flip (use programme:<tau> or hookbank:identity...:<tau>);
    S1 the main def_sameline call gets --flip-octiles when its programme source is identity...;
    S2 the hookbank call gets it when the hook bank is identity...;
    F  def_day.final.env records it.

    uv run python experiments/604-identity-prior/patch604f.py --root <repo> --check
    uv run python experiments/604-identity-prior/patch604f.py --root <repo>
"""
from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path

MARK = "604f: octile flip pass-through"
SAME = Path("scripts/def_sameline.py")
DAY = Path("scripts/def_day.sh")
NEED_SAME = ('startswith(("identity:", "identity@"))', "--skip-contexts")   # 604 and 644 applied
NEED_DAY = ("630dd: rank-neutral UP pass-through", "644: DEF_DAY_SAMELINE_SKIP", "SLSKIP")   # 644 v3

# ------------------------------------------------------------------------------------------------ def_sameline.py
SL_ARG_OLD = '''    ap.add_argument("--programme-source", default="",
                    help="460: a bank whose programme becomes EVERY context's block prior (h1_2025 = Arc's H1, 458); "
                         "matched contexts keep their same-line rows; default off")
'''
SL_ARG_NEW = SL_ARG_OLD + '''    ap.add_argument("--flip-octiles", default=None,   # 604f: octile flip pass-through
                    help="604f: 625's octile flip (8 x 0/1, octile order 0..7) for an identity programme source; default none")
'''
SL_VAL_OLD = '''    a = ap.parse_args()
    t0 = time.time()
'''
SL_VAL_NEW = SL_VAL_OLD + '''    if a.flip_octiles is not None and not (a.programme_source == "identity"   # 604f: octile flip pass-through
                                           or a.programme_source.startswith(("identity:", "identity@"))):
        raise SystemExit(f"--flip-octiles needs an identity programme source (got --programme-source {a.programme_source!r})")
'''
SL_CALL_OLD = '''        ires = DIP.identity_prior(official, icpm, {c: str(v) for c, v in call.items()}, rcpm, rcalls, ilam, P_c4, ps_c4)
'''
SL_CALL_NEW = '''        ires = DIP.identity_prior(official, icpm, {c: str(v) for c, v in call.items()}, rcpm, rcalls, ilam, P_c4, ps_c4,
                                  flip_octiles=a.flip_octiles)   # 604f: octile flip pass-through (None = no flip)
'''
SL_EDITS = [("SL arg", SL_ARG_OLD, SL_ARG_NEW), ("SL validation", SL_VAL_OLD, SL_VAL_NEW), ("SL call", SL_CALL_OLD, SL_CALL_NEW)]

# ------------------------------------------------------------------------------------------------ def_day.sh
DAY_V_OLD = '''# DEF_DAY_FORCE="D=<CCLE sample>,..." forces the same-line call (rehearsals only; never on the day)
'''
DAY_V_NEW = '''# 604f: octile flip pass-through -- DEF_DAY_DE_PRIOR_FLIP=<8 x 0/1> (625's octiles of CLEAN4 strength, order 0..7)
FLIP_MAIN=""; FLIP_HOOK=""
if [ -n "${DEF_DAY_DE_PRIOR_FLIP:-}" ]; then
  case "$DEF_DAY_DE_PRIOR_FLIP" in [01][01][01][01][01][01][01][01]) ;; *) die "DEF_DAY_DE_PRIOR_FLIP must be 8 characters of 0/1 (got $DEF_DAY_DE_PRIOR_FLIP)" ;; esac
  case "${DEF_DAY_PROGRAMME_SOURCE:-}" in identity|identity:*|identity@*) FLIP_MAIN="$DEF_DAY_DE_PRIOR_FLIP" ;; esac
  case "${DEF_DAY_MDE_SIGN:-}" in
    hookbank:identity:*|hookbank:identity@*) FLIP_HOOK="$DEF_DAY_DE_PRIOR_FLIP" ;;
    clean4:*) die "DEF_DAY_DE_PRIOR_FLIP with DEF_DAY_MDE_SIGN=$DEF_DAY_MDE_SIGN: that hook reads CLEAN4 directly and would not see the flip (use programme:<tau> or hookbank:identity...:<tau>)" ;;
  esac
  [ -n "$FLIP_MAIN$FLIP_HOOK" ] || die "DEF_DAY_DE_PRIOR_FLIP needs DEF_DAY_PROGRAMME_SOURCE=identity... or DEF_DAY_MDE_SIGN=hookbank:identity...: nothing else can carry the flip"
  grep -q -- "--flip-octiles" scripts/def_sameline.py || die "def_sameline.py has no --flip-octiles (604f)"
  log "604f: prior flip $DEF_DAY_DE_PRIOR_FLIP on${FLIP_MAIN:+ the block prior}${FLIP_HOOK:+ the hook}"
fi
''' + DAY_V_OLD
DAY_S1_OLD = '''    ${SLSKIP:+--skip-contexts "$SLSKIP"} \\
    >> "$LOG" 2>&1 || die "def_sameline.py failed"
'''   # 644 v3 (d4a3d275) rewrote the call's skip flag to its whitespace-stripped $SLSKIP
DAY_S1_NEW = '''    ${SLSKIP:+--skip-contexts "$SLSKIP"} \\
    ${FLIP_MAIN:+--flip-octiles "$FLIP_MAIN"} \\
    >> "$LOG" 2>&1 || die "def_sameline.py failed"
'''
DAY_S2_OLD = '''      ${DEF_DAY_HOOK_ESTIMATOR:+--programme-estimator "$DEF_DAY_HOOK_ESTIMATOR"} >> "$LOG" 2>&1 \\
'''
DAY_S2_NEW = '''      ${DEF_DAY_HOOK_ESTIMATOR:+--programme-estimator "$DEF_DAY_HOOK_ESTIMATOR"} ${FLIP_HOOK:+--flip-octiles "$FLIP_HOOK"} >> "$LOG" 2>&1 \\
'''
DAY_F_OLD = '''[ -n "${DEF_DAY_HOOK_ESTIMATOR:-}" ] && printf 'DEF_DAY_HOOK_ESTIMATOR=%s\\n' "$DEF_DAY_HOOK_ESTIMATOR" >> "$W/def_day.final.env"
'''
DAY_F_NEW = DAY_F_OLD + '''[ -n "${DEF_DAY_DE_PRIOR_FLIP:-}" ] && printf 'DEF_DAY_DE_PRIOR_FLIP=%s\\n' "$DEF_DAY_DE_PRIOR_FLIP" >> "$W/def_day.final.env"   # 604f: octile flip pass-through
'''
DAY_EDITS = [("DAY validation", DAY_V_OLD, DAY_V_NEW), ("DAY main call", DAY_S1_OLD, DAY_S1_NEW),
             ("DAY hookbank call", DAY_S2_OLD, DAY_S2_NEW), ("DAY final env", DAY_F_OLD, DAY_F_NEW)]


def plan(root: Path) -> dict[Path, str]:
    texts = {SAME: (root / SAME).read_text(), DAY: (root / DAY).read_text()}
    if any(MARK in t for t in texts.values()):
        raise SystemExit("604f already applied (a second apply is refused)")
    for f, needs in ((SAME, NEED_SAME), (DAY, NEED_DAY)):
        missing = [n for n in needs if n not in texts[f]]
        if missing:
            raise SystemExit(f"{f}: lacks {missing}: 604f applies after the whole 10-06 chain")
    out = {}
    for f, edits in ((SAME, SL_EDITS), (DAY, DAY_EDITS)):
        t = texts[f]
        for name, old, new in edits:
            n = t.count(old)
            if n != 1:
                raise SystemExit(f"{f}: anchor {name} found {n} times (expected exactly 1)")
            t = t.replace(old, new, 1)
        if f.suffix == ".py":
            compile(t, str(f), "exec")
        out[f] = t
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    new = plan(a.root)
    if a.check:
        print(f"604f check OK: {SAME} ({len(SL_EDITS)} anchors) and {DAY} ({len(DAY_EDITS)} anchors) found once each")
        return 0
    for f, t in new.items():
        p = a.root / f
        fd, tmp = tempfile.mkstemp(dir=p.parent, prefix=".patch604f-")
        with os.fdopen(fd, "w") as fh:
            fh.write(t)
        os.chmod(tmp, os.stat(p).st_mode)
        os.replace(tmp, p)
    print(f"604f applied: {SAME}, {DAY}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
