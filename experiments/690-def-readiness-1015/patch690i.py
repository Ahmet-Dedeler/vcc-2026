#!/usr/bin/env python3
"""690i — the stack dump's H1 and CD4 doses as day-path knobs (for 708). UNAPPLIED until R-ABC's def_day run has ended.

Why: 708 (prereg 69bd6105, science) re-prices the stack weights on D/E/F's REAL panel after R0, with 438 §5's rule
(dayof438.py verbatim). The rule may move H1 off x2 or CD4 off x1 / x1.5, and a move ships only as a later final entry
built by def_day. def_day hard-codes VCC_H1_UNION_G=0.56 (H1 x2) and knows only CD4 x1 / x1.5, so the other grid points
need a pass-through.

`DEF_DAY_H1_STACK=x1|x1.5|x2|x2.5|x3` (new; default x2) and `DEF_DAY_CD4_STACK=x0.75|x1|x1.25|x1.5|x1.75|x2` (530's
knob, grid widened): both act inside the dump subshell only, with recommend438's values (H1 0.28 x dose printed with
%.10g; CD4 0.3443655278649641 x dose with %.16g). The dump is named dump_h1x<dose>_cd4x<dose>.npz with the dot dropped
(x2.5 -> h1x25, x0.75 -> cd4x075), so the level meta's bulk_load names both doses and riders_match can check it.
Unset (H1 x2, CD4 x1) = today: the same file name dump_h1x2_cd4x1.npz, VCC_H1_UNION_G=0.56, no VCC_CD4_TARGET_G export;
x1.5 keeps 530's name and value. test690i checks this on the dump subshell's exported environment, against the
unpatched script, for every combination it shares.

Files: scripts/def_day.sh (8 anchors), scripts/def_day_helpers.py (4: the dump name, the spec's optional h1, the grid
checks), scripts/october_recipe.py (3: cd4_stack grid, h1_stack switch, installed checks), docs/october-decisions.json
(+ h1_stack DECIDED x2, inert: plain `command` still resolves to today's env).

    uv run python experiments/690-def-readiness-1015/patch690i.py --root <repo> [--check]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

MARK = "690i: the doses come from 438's grid"

DAY_EDITS = [
    ("doc", '''#   DEF_DAY_CD4_STACK=x1|x1.5 (530, default x1): x1.5 exports VCC_CD4_TARGET_G inside the dump subshell only (457).
''', '''#   DEF_DAY_CD4_STACK=x0.75|x1|x1.25|x1.5|x1.75|x2 (530; 690i widens it to 438's grid; default x1): a dose other than x1
#    exports VCC_CD4_TARGET_G = <dose> x 0.3443655278649641 inside the dump subshell only (457); x1 exports nothing (DEF_ENV).
#   DEF_DAY_H1_STACK=x1|x1.5|x2|x2.5|x3 (690i: H1 dose; default x2): the dump subshell's VCC_H1_UNION_G = <dose> x 0.28 (438's
#    grid, recommend438's values; 708 re-prices it on D/E/F's panel after R0). The dump is dump_h1x<dose>_cd4x<dose>.npz with
#    the dot dropped (x2.5 -> h1x25, x0.75 -> cd4x075), so the level meta's bulk_load names both doses.
'''),
    ("validate", '''case "${DEF_DAY_CD4_STACK:-x1}" in
  x1|x1.5) ;;
  *) die "DEF_DAY_CD4_STACK must be x1 or x1.5 (got ${DEF_DAY_CD4_STACK})" ;;
esac
''', '''case "${DEF_DAY_CD4_STACK:-x1}" in
  x0.75|x1|x1.25|x1.5|x1.75|x2) ;;   # 690i: 438's CD4 grid
  *) die "DEF_DAY_CD4_STACK must be one of x0.75 x1 x1.25 x1.5 x1.75 x2 (got ${DEF_DAY_CD4_STACK})" ;;
esac
case "${DEF_DAY_H1_STACK:-x2}" in   # 690i: 438's H1 grid
  x1|x1.5|x2|x2.5|x3) ;;
  *) die "DEF_DAY_H1_STACK must be one of x1 x1.5 x2 x2.5 x3 (got ${DEF_DAY_H1_STACK})" ;;
esac
'''),
    ("dump", '''case "${DEF_DAY_CD4_STACK:-x1}" in
  x1.5) DUMP="$W/dump_h1x2_cd4x15.npz" ;;
  *) DUMP="$W/dump_h1x2_cd4x1.npz" ;;
esac
CD4G=""
[ "${DEF_DAY_CD4_STACK:-x1}" = "x1.5" ] && CD4G="0.5165482917974462"
''', '''# 690i: the doses come from 438's grid with recommend438's values. Unset (H1 x2, CD4 x1) is today's dump: the same file
# name, VCC_H1_UNION_G=0.56 and no VCC_CD4_TARGET_G export; x1.5 keeps 530's name and value.
case "${DEF_DAY_H1_STACK:-x2}" in
  x1) H1T=h1x1; H1G=0.28 ;;
  x1.5) H1T=h1x15; H1G=0.42 ;;
  x2.5) H1T=h1x25; H1G=0.7 ;;
  x3) H1T=h1x3; H1G=0.84 ;;
  *) H1T=h1x2; H1G=0.56 ;;
esac
CD4G=""
case "${DEF_DAY_CD4_STACK:-x1}" in
  x0.75) C4T=cd4x075; CD4G="0.2582741458987231" ;;
  x1.25) C4T=cd4x125; CD4G="0.4304569098312051" ;;
  x1.5) C4T=cd4x15; CD4G="0.5165482917974462" ;;
  x1.75) C4T=cd4x175; CD4G="0.6026396737636872" ;;
  x2) C4T=cd4x2; CD4G="0.6887310557299282" ;;
  *) C4T=cd4x1 ;;
esac
DUMP="$W/dump_${H1T}_${C4T}.npz"
'''),
    ("log", '''[ -n "$CD4G" ] && log "CD4 stack x1.5: dump ${DUMP##*/} VCC_CD4_TARGET_G=$CD4G (level builds keep DEF_ENV)"
''', '''[ -n "$CD4G" ] && log "CD4 stack ${DEF_DAY_CD4_STACK}: dump ${DUMP##*/} VCC_CD4_TARGET_G=$CD4G (level builds keep DEF_ENV)"
[ "$H1G" != 0.56 ] && log "H1 stack ${DEF_DAY_H1_STACK}: dump ${DUMP##*/} VCC_H1_UNION_G=$H1G (level builds keep DEF_ENV)"   # 690i
'''),
    ("export", '''  export VCC_H1_UNION_G=0.56 VCC_BULK_DUMP="$DUMP"
''', '''  export VCC_H1_UNION_G="$H1G" VCC_BULK_DUMP="$DUMP"
'''),
    ("r530", '''[ -n "${DEF_DAY_CALLSET_BLOCK:-}" ] && _r530=1
if [ -n "$_r530" ]; then
''', '''[ -n "${DEF_DAY_CALLSET_BLOCK:-}" ] && _r530=1
[ "${DEF_DAY_H1_STACK:-x2}" != "x2" ] && _r530=1   # 690i
if [ -n "$_r530" ]; then
'''),
    ("spec", '''  [ -n "${DEF_DAY_CALLSET_BLOCK:-}" ] && _spec="${_spec},callset=${DEF_DAY_CALLSET_BLOCK}"
''', '''  [ -n "${DEF_DAY_CALLSET_BLOCK:-}" ] && _spec="${_spec},callset=${DEF_DAY_CALLSET_BLOCK}"
  [ "${DEF_DAY_H1_STACK:-x2}" != "x2" ] && _spec="${_spec},h1=${DEF_DAY_H1_STACK}"   # 690i
'''),
    ("final_env", '''[ "${DEF_DAY_CD4_STACK:-x1}" != "x1" ] && printf 'DEF_DAY_CD4_STACK=%s\\n' "$DEF_DAY_CD4_STACK" >> "$W/def_day.final.env"
''', '''[ "${DEF_DAY_CD4_STACK:-x1}" != "x1" ] && printf 'DEF_DAY_CD4_STACK=%s\\n' "$DEF_DAY_CD4_STACK" >> "$W/def_day.final.env"
[ "${DEF_DAY_H1_STACK:-x2}" != "x2" ] && printf 'DEF_DAY_H1_STACK=%s\\n' "$DEF_DAY_H1_STACK" >> "$W/def_day.final.env"   # 690i
'''),
]

HELPER_EDITS = [
    ("doc", '''stack=<x1|x1.5>,rho=<-|s>,strength=<-|theta>[,lodo_hi=<h>][,callset=<prop[:swap[:steer]]>]]''',
     '''stack=<x0.75..x2>,rho=<-|s>,strength=<-|theta>[,lodo_hi=<h>][,callset=<prop[:swap[:steer]]>][,h1=<x1..x3>]]'''),
    ("name", '''def _dump_name530(stack: str, strength: str) -> str:
    """The stack-dump filename def_day.sh writes. Unset stack and strength -> today's name."""
    base = "dump_h1x2_cd4x15" if stack == "x1.5" else "dump_h1x2_cd4x1"
''', '''# 690i: 438's dose grid; the dump-name token drops the dot (def_day.sh's case tables)
H1_TOK690I = {"x1": "h1x1", "x1.5": "h1x15", "x2": "h1x2", "x2.5": "h1x25", "x3": "h1x3"}
CD4_TOK690I = {"x0.75": "cd4x075", "x1": "cd4x1", "x1.25": "cd4x125", "x1.5": "cd4x15", "x1.75": "cd4x175", "x2": "cd4x2"}


def _dump_name530(stack: str, strength: str, h1: str = "x2") -> str:
    """The stack-dump filename def_day.sh writes. Unset stack, strength and h1 -> today's name (690i: any 438 grid dose)."""
    base = f"dump_{H1_TOK690I[h1]}_{CD4_TOK690I[stack]}"
'''),
    ("optional", '''    optional = {"lodo_hi", "callset"}
''', '''    optional = {"lodo_hi", "callset", "h1"}   # 690i: h1
'''),
    ("check", '''    if stack not in ("x1", "x1.5"):
        return f"stack {stack!r} is not x1 or x1.5"
    want_name = _dump_name530(stack, strength)
''', '''    if stack not in CD4_TOK690I:
        return f"stack {stack!r} is not one of {sorted(CD4_TOK690I)}"
    h1 = spec.get("h1", "x2")   # 690i
    if h1 not in H1_TOK690I:
        return f"h1 {h1!r} is not one of {sorted(H1_TOK690I)}"
    want_name = _dump_name530(stack, strength, h1)
'''),
]

RECIPE_EDITS = [
    ("cd4", '''        if cd4_val == "x1.5":
            env["DEF_DAY_CD4_STACK"] = "x1.5"
        elif cd4_val not in {"x1", "off", "0", "none"}:
            problems.append(f"cd4_stack {cd4_val!r} is not x1 or x1.5")
''', '''        if cd4_val in {"x0.75", "x1.25", "x1.5", "x1.75", "x2"}:   # 690i: 438's CD4 grid (530 had x1.5 only)
            env["DEF_DAY_CD4_STACK"] = cd4_val
        elif cd4_val not in {"x1", "off", "0", "none"}:
            problems.append(f"cd4_stack {cd4_val!r} is not one of x0.75 x1 x1.25 x1.5 x1.75 x2")
    # 690i: the stack dump's H1 dose; x2 (438's October default) is today's dump and is not exported
    h1s = s.get("h1_stack", {})
    if h1s.get("status") == "DECIDED":
        h1_val = str(h1s.get("value"))
        if h1_val in {"x1", "x1.5", "x2.5", "x3"}:
            env["DEF_DAY_H1_STACK"] = h1_val
        elif h1_val not in {"x2", "off", "none"}:
            problems.append(f"h1_stack {h1_val!r} is not one of x1 x1.5 x2 x2.5 x3")
'''),
    ("export_of", '''        "DEF_DAY_CD4_STACK": 'VCC_CD4_TARGET_G="$CD4G"',
''', '''        "DEF_DAY_CD4_STACK": 'VCC_CD4_TARGET_G="$CD4G"',
        "DEF_DAY_H1_STACK": 'VCC_H1_UNION_G="$H1G"',   # 690i
'''),
    ("installed", '''        out.append("cd4_stack x1.5: scripts/def_day.sh has no DEF_DAY_CD4_STACK dump override (530)")
''', '''        out.append(f"cd4_stack {env['DEF_DAY_CD4_STACK']}: scripts/def_day.sh has no DEF_DAY_CD4_STACK dump override (530)")
    if env.get("DEF_DAY_CD4_STACK", "x1.5") != "x1.5" and "690i: the doses come from 438's grid" not in dd:
        out.append(f"cd4_stack {env['DEF_DAY_CD4_STACK']}: scripts/def_day.sh predates 690i (530 has x1.5 only)")
'''),
]

SWITCH = {"h1_stack": {
    "status": "DECIDED", "value": "x2", "at": "2026-10-15T22:00:00Z",
    "evidence": "438 §5: H1 x2 / CD4 x1 is the October default (q 1.75); def_day has always dumped at VCC_H1_UNION_G=0.56. "
                "690i makes the dose a knob (DEF_DAY_H1_STACK) so 708 (prereg 69bd6105) can re-price it on D/E/F's real "
                "panel after R0 and ship a move only as a LATER final entry with its own def_day build and gates.",
    "held_out": "438's October-like panels (x2 kept); on A/B/C's own panel 438's rule picks H1 x2.5 / CD4 x1.5 (+0.0022, "
                "708/work/expected_abcdef.txt): panel-dependent, hence 708",
    "added_by": "690i (2026-10-15)"}}


def apply(text: str, edits: list, label: str) -> tuple[str, list]:
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
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    r = a.root.resolve()
    day, hel, rec, dec = (r / "scripts/def_day.sh", r / "scripts/def_day_helpers.py", r / "scripts/october_recipe.py",
                          r / "docs/october-decisions.json")
    dd, hh, rr = day.read_text(), hel.read_text(), rec.read_text()
    dj = json.loads(dec.read_text())
    if MARK in dd or "H1_TOK690I" in hh or "h1_stack" in rr or "h1_stack" in dj["switches"]:
        print("REFUSED: 690i is already applied (or partly)", file=sys.stderr)
        return 1
    dd2, e1 = apply(dd, DAY_EDITS, "def_day.sh")
    hh2, e2 = apply(hh, HELPER_EDITS, "def_day_helpers.py")
    rr2, e3 = apply(rr, RECIPE_EDITS, "october_recipe.py")
    errs = e1 + e2 + e3
    if errs:
        print("REFUSED:\n  " + "\n  ".join(errs), file=sys.stderr)
        return 1
    sw = dj["switches"]
    tail = [k for k in sw if sw[k].get("status") == "DAY-OF"]
    new = {k: v for k, v in sw.items() if k not in tail}
    new.update(SWITCH)
    new.update({k: sw[k] for k in tail})
    dj["switches"] = new
    if a.check:
        print(f"CHECK OK: {len(DAY_EDITS)} def_day.sh anchors, {len(HELPER_EDITS)} helper anchors, {len(RECIPE_EDITS)} recipe "
              "anchors, +h1_stack (DECIDED x2)")
        return 0
    day.write_text(dd2)
    hel.write_text(hh2)
    rec.write_text(rr2)
    dec.write_text(json.dumps(dj, indent=1) + "\n")   # october_recipe.py's own format
    print("APPLIED 690i")
    return 0


if __name__ == "__main__":
    sys.exit(main())
