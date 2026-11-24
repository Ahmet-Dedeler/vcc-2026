#!/usr/bin/env python3
"""690a — make the late DE riders reachable from the October recipe. UNAPPLIED until the coordinator applies it.

Written 2026-10-15 by [17adab] (D/E/F readiness). Why: def_day.sh already passes DEF_DAY_DE_UPRN (630dd) and
DEF_DAY_MDE_SIGN_KEYGATE (612dd) through, but `scripts/october_recipe.py` has no switch for either, so a TRANSFERS read
could not reach the day's command. And 631's DOWN-magnitude knob (VCC_DE_DNMAG, builder 8f8f6d5c + dn_mag.py 0a8cc560)
has no def_day pass-through at all. Every insertion is inert while its knob is unset: today's resolved command is
unchanged byte for byte.

Four files:
  scripts/def_day.sh — DEF_DAY_DE_DNMAG=d_m<m*100>[_q<n>] -> VCC_DE_DNMAG, six anchored insertions mirroring 630dd:
    DOC  documented next to DEF_DAY_DE_UPRN;
    V    validation after DE_UPRN's: needs DEF_DAY_DE_MAG (the builder refuses VCC_DE_DNMAG without VCC_DE_MAG), the
         builder's own grammar `d_m(\\d{3})(?:_q([1-9]\\d*))?` with m*100 in [125, 400], and the builder must carry
         VCC_DE_DNMAG;
    E    export for the level builds next to VCC_DE_UPRN's (never for the stack dump: the riders are exported after it);
    C    the final meta must record config.de_dnmag.variant == the request AND module_sha256 == 0a8cc560… (the
         reviewed module); inserted in front of 630dd's meta check;
    F    def_day.final.env records it;
    A    postbuild's allowed config diff gains de_dnmag when set.
  scripts/october_recipe.py — switches de_uprn / mde_sign_keygate / de_dnmag in ENV_OF and BUILDER_NEEDS, plus
    their dependencies (de_uprn and de_dnmag need de_mag; mde_sign_keygate needs mde_sign) and grammar checks.
  docs/october-decisions.json — the three switches, PENDING, with decided_by / if_on / needs / held_out like the
    others. PENDING keeps `command` fail-closed until the coordinator decides each one (a board read plus a held-out
    read, AGENTS 11.1c), "off" included.
  scripts/test_october_recipe.py — the fixture's list of PENDING switches decided off (OFF538) gains the three, so the
    tests that resolve without --allow-pending keep resolving; and the live-state test becomes state-aware (rc 2 iff a
    switch is PENDING, else rc 0), so the all-decided Oct 21 state is green (coordinator's request, 10-15).

The patch refuses a def_day.sh without 630dd / 612dd, a builder without 631's pinned module sha, and a second apply.

    uv run python experiments/690-def-readiness-1015/patch690a.py --root <repo> --check
    uv run python experiments/690-def-readiness-1015/patch690a.py --root <repo>
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

MARK = "690a: DOWN magnitude pass-through"
DNMAG_MODULE_SHA256 = "0a8cc560231a16f0fff7fbbf316e3a5e190cbdf20490239ef354b49b3ff60dac"  # src/vcc2026/dn_mag.py (631)
DAY = Path("scripts/def_day.sh")
RECIPE = Path("scripts/october_recipe.py")
DECISIONS = Path("docs/october-decisions.json")
BUILDER = Path("scripts/make_submission_pm_cd4.py")
TESTS = Path("scripts/test_october_recipe.py")

# ---------------------------------------------------------------- def_day.sh
DOC_OLD = '''#   DEF_DAY_DE_UPRN=<s|b>_k<K*100>[_q<n>r<n>p<n>] (630dd: rank-neutral UP pass-through; default off; e.g. b_k100_q20r3p10):
#    630's rank-neutral UP pass after zs (VCC_DE_UPRN). Needs DEF_DAY_DE_MAG. Builder >= 630 (rank_up module 7e6b3e52).
'''
DOC_NEW = DOC_OLD + '''#   DEF_DAY_DE_DNMAG=d_m<m*100>[_q<n>] (690a: DOWN magnitude pass-through; default off; e.g. d_m200): 631's DOWN-magnitude
#    pass after zs and 630's UP pass (VCC_DE_DNMAG). Needs DEF_DAY_DE_MAG. Builder >= 631 (dn_mag module 0a8cc560).
'''
V_OLD = r'''  grep -q "630: rank-neutral UP" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_DE_UPRN (630)"
fi
'''
V_NEW = V_OLD + r'''# 690a: DOWN magnitude pass-through -- DEF_DAY_DE_DNMAG=<variant> (validated here, exported with the other DE riders)
if [ -n "${DEF_DAY_DE_DNMAG:-}" ]; then
  [ -n "${DEF_DAY_DE_MAG:-}" ] || die "DEF_DAY_DE_DNMAG runs after the zs pass: it needs DEF_DAY_DE_MAG"
  python3 -c 'import re, sys; m = re.fullmatch(r"d_m(\d{3})(?:_q([1-9]\d*))?", sys.argv[1]); sys.exit(0 if m and 125 <= int(m.group(1)) <= 400 else 1)' "$DEF_DAY_DE_DNMAG" 2>/dev/null \
    || die "DEF_DAY_DE_DNMAG must be d_m<m*100>[_q<n>] with m in [1.25, 4], e.g. d_m200 (got $DEF_DAY_DE_DNMAG)"
  grep -q "VCC_DE_DNMAG" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_DE_DNMAG (631)"
fi
'''
E_OLD = '''[ -n "${DEF_DAY_DE_UPRN:-}" ] && { export VCC_DE_UPRN="$DEF_DAY_DE_UPRN"; log "DE rider for the level builds: VCC_DE_UPRN=$VCC_DE_UPRN"; }   # 630dd: rank-neutral UP pass-through
'''
E_NEW = E_OLD + '''[ -n "${DEF_DAY_DE_DNMAG:-}" ] && { export VCC_DE_DNMAG="$DEF_DAY_DE_DNMAG"; log "DE rider for the level builds: VCC_DE_DNMAG=$VCC_DE_DNMAG"; }   # 690a: DOWN magnitude pass-through
'''
C_OLD = '''if [ -n "${VCC_DE_UPRN:-}" ]; then   # 630dd: rank-neutral UP pass-through -- the final meta must carry the variant and the reviewed module
'''
C_NEW = f'''if [ -n "${{VCC_DE_DNMAG:-}}" ]; then   # 690a: DOWN magnitude pass-through -- the final meta must carry the variant and the reviewed module
  python3 - "$FINAL.h5ad.meta.json" "$VCC_DE_DNMAG" "{DNMAG_MODULE_SHA256}" <<'PY' || die "the final meta does not carry VCC_DE_DNMAG=$VCC_DE_DNMAG with dn_mag module {DNMAG_MODULE_SHA256[:8]}"
import json, sys
u = json.load(open(sys.argv[1]))["config"].get("de_dnmag") or {{}}
sys.exit(0 if u.get("variant") == sys.argv[2] and u.get("module_sha256") == sys.argv[3] else 1)
PY
fi
''' + C_OLD
F_OLD = '''[ -n "${VCC_DE_UPRN:-}" ] && printf 'VCC_DE_UPRN=%s\\n' "$VCC_DE_UPRN" >> "$W/def_day.final.env"   # 630dd: rank-neutral UP pass-through
'''
F_NEW = F_OLD + '''[ -n "${VCC_DE_DNMAG:-}" ] && printf 'VCC_DE_DNMAG=%s\\n' "$VCC_DE_DNMAG" >> "$W/def_day.final.env"   # 690a: DOWN magnitude pass-through
'''
A_OLD = '''[ -n "${VCC_DE_UPRN:-}" ] && ALLOW="$ALLOW,de_uprn"   # 630dd: rank-neutral UP pass-through
'''
A_NEW = A_OLD + '''[ -n "${VCC_DE_DNMAG:-}" ] && ALLOW="$ALLOW,de_dnmag"   # 690a: DOWN magnitude pass-through
'''
DAY_EDITS = [("DOC", DOC_OLD, DOC_NEW), ("V", V_OLD, V_NEW), ("E", E_OLD, E_NEW), ("C", C_OLD, C_NEW),
             ("F", F_OLD, F_NEW), ("A", A_OLD, A_NEW)]

# ---------------------------------------------------------------- october_recipe.py
R_ENV_OLD = '''    "callset_block": "DEF_DAY_CALLSET_BLOCK",
}
'''
R_ENV_NEW = '''    "callset_block": "DEF_DAY_CALLSET_BLOCK",
    # 690a: the late DE riders (def_day.sh carries all three; PENDING until a board read and a held-out read decide them)
    "de_uprn": "DEF_DAY_DE_UPRN",
    "mde_sign_keygate": "DEF_DAY_MDE_SIGN_KEYGATE",
    "de_dnmag": "DEF_DAY_DE_DNMAG",
}
'''
R_NEEDS_OLD = '''                 "DEF_DAY_CALLSET_BLOCK": "VCC_RC_CALLSET_BLOCK"}
'''
R_NEEDS_NEW = '''                 "DEF_DAY_CALLSET_BLOCK": "VCC_RC_CALLSET_BLOCK",
                 # 690a
                 "DEF_DAY_DE_UPRN": "630: rank-neutral UP", "DEF_DAY_MDE_SIGN_KEYGATE": "VCC_DE_MDE_SIGN_KEYGATE",
                 "DEF_DAY_DE_DNMAG": "VCC_DE_DNMAG"}
'''
R_CALL_OLD = '''    _problems_598(env, problems)
'''
R_CALL_NEW = R_CALL_OLD + '''    _problems_690a(env, problems)
'''
R_FN_OLD = '''def installed_problems(env: dict, flags: str) -> list:
'''
R_FN_NEW = '''def _problems_690a(env: dict, problems: list) -> None:
    """690a: dependencies and grammar of the late DE riders, matching def_day.sh's own validation. PENDING stays out of env."""
    import re  # noqa: PLC0415
    if "DEF_DAY_DE_UPRN" in env:
        if "DEF_DAY_DE_MAG" not in env:
            problems.append("de_uprn runs after the zs pass: it needs de_mag")
        m = re.fullmatch(r"(s|b)_k(\\d{2,3})(?:_q([1-9]\\d*)r([1-9]\\d*)p([1-9]\\d*))?", env["DEF_DAY_DE_UPRN"])
        if not (m and 50 <= int(m.group(2)) <= 300):
            problems.append(f"de_uprn {env['DEF_DAY_DE_UPRN']!r} must be <s|b>_k<K*100>[_q<n>r<n>p<n>] with K in [0.5, 3]")
    if "DEF_DAY_MDE_SIGN_KEYGATE" in env:
        if "DEF_DAY_MDE_SIGN" not in env:
            problems.append("mde_sign_keygate gates the hook's sign: it needs mde_sign")
        try:
            q = float(env["DEF_DAY_MDE_SIGN_KEYGATE"])
        except ValueError:
            q = -1.0
        if not 0.0 < q < 1.0:
            problems.append(f"mde_sign_keygate {env['DEF_DAY_MDE_SIGN_KEYGATE']!r} must be a quantile 0 < q < 1")
    if "DEF_DAY_DE_DNMAG" in env:
        if "DEF_DAY_DE_MAG" not in env:
            problems.append("de_dnmag runs after the zs pass: it needs de_mag")
        m = re.fullmatch(r"d_m(\\d{3})(?:_q([1-9]\\d*))?", env["DEF_DAY_DE_DNMAG"])
        if not (m and 125 <= int(m.group(1)) <= 400):
            problems.append(f"de_dnmag {env['DEF_DAY_DE_DNMAG']!r} must be d_m<m*100>[_q<n>] with m in [1.25, 4]")


''' + R_FN_OLD
RECIPE_EDITS = [("ENV_OF", R_ENV_OLD, R_ENV_NEW), ("BUILDER_NEEDS", R_NEEDS_OLD, R_NEEDS_NEW),
                ("resolve", R_CALL_OLD, R_CALL_NEW), ("fn", R_FN_OLD, R_FN_NEW)]

# ---------------------------------------------------------------- test_october_recipe.py
T_OLD = '''          ("lodo_hi", "off"), ("callset_block", "off")]
'''
T_NEW = '''          ("lodo_hi", "off"), ("callset_block", "off"),
          # 690a: three more PENDING switches, decided off for the same tests
          ("de_uprn", "off"), ("mde_sign_keygate", "off"), ("de_dnmag", "off")]
'''
T2_OLD = '''    def test_live_state_refuses_pending_then_resolves_with_allow(self):
        rc, out, err = self.cmd()
        self.assertEqual(rc, 2)
        self.assertIn("PENDING", err)
'''
T2_NEW = '''    def test_live_state_refuses_pending_then_resolves_with_allow(self):
        # 690a: state-aware. Plain `command` refuses (rc 2) iff a switch is PENDING; in the all-decided state (Oct 21) it
        # resolves (rc 0), and the --allow-pending checks below hold either way.
        live = json.loads((Path(__file__).resolve().parents[1] / "docs/october-decisions.json").read_text())
        pending = [k for k, v in live.get("switches", live).items() if v.get("status") == "PENDING"]
        rc, out, err = self.cmd()
        if pending:
            self.assertEqual(rc, 2)
            self.assertIn("PENDING", err)
        else:
            self.assertEqual(rc, 0, err)
'''
TEST_EDITS = [("OFF538", T_OLD, T_NEW), ("live_state", T2_OLD, T2_NEW)]

# ---------------------------------------------------------------- october-decisions.json
NEW_SWITCHES = {
    "de_uprn": {
        "status": "PENDING", "value": None,
        "decided_by": "630arm's board read (446bj2z + VCC_DE_UPRN=b_k100_q20r3p10, entry QJRJK6OC6Zw1kjwnZqid, 10-08) "
                      "against its preregistered bins, plus a held-out read",
        "if_on": "<s|b>_k<K*100>[_q<n>r<n>p<n>] with K in [0.5, 3] (e.g. b_k100_q20r3p10)",
        "otherwise": "off",
        "needs": "de_mag on (the UP pass runs after zs); def_day 630dd; builder 630 (rank_up module 7e6b3e52). Meta "
                 "config.de_uprn carries variant and module_sha256; def_day checks both and reports FAILSAFE_TOTAL",
        "held_out": "a prior-only DE arm: price it on H1's real emitter (618 / f 0.35) before it is DECIDED on",
        "added_by": "690a (2026-10-15)"},
    "mde_sign_keygate": {
        "status": "PENDING", "value": None,
        "decided_by": "a board read of 612's key-gated hook against the plain CLEAN4 hook (496a), plus a held-out read",
        "if_on": "a quantile 0 < q < 1 (e.g. 0.5)",
        "otherwise": "off",
        "needs": "mde_sign on; def_day 612dd; builder VCC_DE_MDE_SIGN_KEYGATE (612). Meta config.de_mde_sign_keygate.q",
        "held_out": "612's H1 / fold reads; the board read decides (prior-only arms read near-exactly, rule 13)",
        "added_by": "690a (2026-10-15)"},
    "de_dnmag": {
        "status": "PENDING", "value": None,
        "decided_by": "the DN671 arm's board read (631's DOWN magnitude on the October base) against its preregistered "
                      "bins, plus a held-out read",
        "if_on": "d_m<m*100>[_q<n>] with m in [1.25, 4] (e.g. d_m200)",
        "otherwise": "off",
        "needs": "de_mag on (the DOWN pass runs after zs and 630's UP pass); def_day 690a; builder 631 (dn_mag module "
                 "0a8cc560). Meta config.de_dnmag carries variant and module_sha256; def_day checks both",
        "held_out": "631's H1 price on the real emitter with production's zero repair (AGENTS 6.5, 654: count lfc < -30 "
                    "cells on any DOWN-magnitude price)",
        "added_by": "690a (2026-10-15)"},
}


def apply_edits(text: str, edits: list, label: str) -> tuple[str, list]:
    errs = []
    for name, old, new in edits:
        n = text.count(old)
        if n != 1:
            errs.append(f"{label} {name}: anchor found {n} times (want 1)")
            continue
        text = text.replace(old, new)
    return text, errs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--check", action="store_true", help="verify every anchor and precondition; write nothing")
    a = ap.parse_args()
    root = a.root.resolve()
    day, rec, dec, bu, tst = (root / p for p in (DAY, RECIPE, DECISIONS, BUILDER, TESTS))
    errs = []
    dd, rr = day.read_text(), rec.read_text()
    dj = json.loads(dec.read_text())
    if MARK in dd or "_problems_690a" in rr or any(k in dj["switches"] for k in NEW_SWITCHES):
        print("REFUSED: 690a is already applied (marker, recipe function or decisions entry present)", file=sys.stderr)
        return 1
    if "630dd: rank-neutral UP pass-through" not in dd or "612dd: key gate pass-through" not in dd:
        errs.append("def_day.sh lacks 630dd / 612dd (apply the 10-06 bundle first)")
    if DNMAG_MODULE_SHA256 not in bu.read_text() or "VCC_DE_DNMAG" not in bu.read_text():
        errs.append(f"the builder does not pin 631's dn_mag module {DNMAG_MODULE_SHA256[:8]} (apply patch_builder631 first)")
    dd2, e1 = apply_edits(dd, DAY_EDITS, "def_day.sh")
    rr2, e2 = apply_edits(rr, RECIPE_EDITS, "october_recipe.py")
    tt2, e3 = apply_edits(tst.read_text(), TEST_EDITS, "test_october_recipe.py")
    errs += e1 + e2 + e3
    if errs:
        print("REFUSED:\n  " + "\n  ".join(errs), file=sys.stderr)
        return 1
    sw = dj["switches"]
    order = list(sw)
    tail = [k for k in order if sw[k].get("status") == "DAY-OF"]
    new_sw = {k: sw[k] for k in order if k not in tail}
    new_sw.update(NEW_SWITCHES)
    new_sw.update({k: sw[k] for k in tail})
    dj["switches"] = new_sw
    dj2 = json.dumps(dj, indent=1) + "\n"
    if a.check:
        print("CHECK OK: 6 def_day.sh anchors, 4 october_recipe.py anchors, 2 test anchors, 3 new PENDING switches; "
              "nothing written")
        return 0
    day.write_text(dd2)
    rec.write_text(rr2)
    dec.write_text(dj2)
    tst.write_text(tt2)
    print(f"APPLIED 690a to {root}: def_day.sh +DEF_DAY_DE_DNMAG; october_recipe.py +de_uprn/mde_sign_keygate/de_dnmag; "
          f"october-decisions.json +3 PENDING; test_october_recipe.py OFF538 +3")
    return 0


if __name__ == "__main__":
    sys.exit(main())
