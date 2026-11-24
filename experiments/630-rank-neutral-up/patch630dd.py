#!/usr/bin/env python3
"""630dd -- the def_day.sh pass-through for 630's rank-neutral UP pass: DEF_DAY_DE_UPRN=<variant> -> VCC_DE_UPRN.
UNAPPLIED. Written by [4d0207] at [d4cd5f]'s request (council 600, 2026-10-02 ~13:30Z). Order: 598 -> 612 -> 630 -> 630dd
(after the 10-05 bundle). It composes with 612dd in EITHER order to the same bytes: every anchor is an ORIGINAL line that
612dd neither edits nor appends to (612dd appends after the MDE_SIGN lines, 630dd after the DE_MAG lines, and its meta
check goes in front of the RC_CALLSET_BLOCK check, while 612dd's goes after the riders_match block).

ONE file, scripts/def_day.sh (the POST-598 file), seven anchored insertions. Unset knob = no export, no check, no meta key,
no final.env line, no ALLOW change, so today's command is unchanged.
  DOC  the knob documented next to DEF_DAY_DE_MAG;
  V    validation after DEF_DAY_DE_MAG's: needs DEF_DAY_DE_MAG (the builder runs UPRN after the zs pass and refuses
       VCC_DE_UPRN without VCC_DE_MAG); the variant must match the builder's own grammar
       `(s|b)_k(\\d{2,3})(?:_q([1-9]\\d*)r([1-9]\\d*)p([1-9]\\d*))?` (607: q / r / p >= 1, no leading zero) with K = k/100 in
       [0.5, 3]; the builder must carry "630: rank-neutral UP";
  E    export VCC_DE_UPRN for the level builds, next to VCC_DE_MAG's export;
  C    the final meta must record config.de_uprn.variant == the requested variant AND config.de_uprn.module_sha256 ==
       UPRN_MODULE_SHA256 (the reviewed module; patch time refuses a builder pinned to any other module); then a REPORT
       (never fatal) per context from the final rc_audit dump: mean frac_got / genes moved, guard_failed and
       invariant_failed target counts (a failed target keeps its zs block), the prereg's capacity floor; logged and
       written to <work>/uprn_audit.json, ending ALWAYS with the greppable `630dd: uprn FAILSAFE_TOTAL <n>`: n = guard +
       invariant fail-safe triggers, -1 = the dump was missing, unreadable or had no ur_ rows (607b). A KA630 build
       passes only with the line present and n == 0 (coordinator 10-02). The dump is 554's per-build
       `$W/$N-$tag.rc_audit.jsonl` (`build()` sets VCC_RC_AUDIT_DUMP and removes a stale one; the post-bundle file);
  F    def_day.final.env records it;
  A    postbuild's allowed config diff gains de_uprn when set;
  S    a context with the pass effectively OFF (fail-safe on > 10 % of its targets; C writes <work>/uprn_effectively_off.txt)
       is repeated as a WARNING right before def_day's DONE line (d4cd5f: loud, not fatal).

    uv run python experiments/630-rank-neutral-up/patch630dd.py --root <repo> --check
    uv run python experiments/630-rank-neutral-up/patch630dd.py --root <repo>
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import tempfile
from pathlib import Path

MARK = "630dd: rank-neutral UP pass-through"
MARK598 = "598: tercilec"
MARK630 = "630: rank-neutral UP"
UPRN_MODULE_SHA256 = "7e6b3e52e2ad374b7b6e21e92327fcf65771433693eb3846c09c0189f3f2959f"   # c30d64ce (605 fix + fail-safe invariants)
DAY = Path("scripts/def_day.sh")
BUILDER = Path("scripts/make_submission_pm_cd4.py")

DOC_OLD = '''#   DEF_DAY_DE_MAG=<variant> (502, default off; e.g. dn_k250_zs): the zs exact-U pass (VCC_DE_MAG) on the level builds.
'''
DOC_NEW = DOC_OLD + '''#   DEF_DAY_DE_UPRN=<s|b>_k<K*100>[_q<n>r<n>p<n>] (630dd: rank-neutral UP pass-through; default off; e.g. b_k100_q20r3p10):
#    630's rank-neutral UP pass after zs (VCC_DE_UPRN). Needs DEF_DAY_DE_MAG. Builder >= 630 (rank_up module 7e6b3e52).
'''
V_OLD = '''if [ -n "${DEF_DAY_DE_MAG:-}" ]; then
  case "$DEF_DAY_DE_MAG" in *[!a-z0-9_]*) die "DEF_DAY_DE_MAG must be a variant name like dn_k250_zs (got $DEF_DAY_DE_MAG)" ;; esac
  grep -q "VCC_DE_MAG" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_DE_MAG"
fi
'''
V_NEW = V_OLD + r'''# 630dd: rank-neutral UP pass-through -- DEF_DAY_DE_UPRN=<variant> (validated here, exported with the other DE riders)
if [ -n "${DEF_DAY_DE_UPRN:-}" ]; then
  [ -n "${DEF_DAY_DE_MAG:-}" ] || die "DEF_DAY_DE_UPRN runs after the zs pass: it needs DEF_DAY_DE_MAG"
  python3 -c 'import re, sys; m = re.fullmatch(r"(s|b)_k(\d{2,3})(?:_q([1-9]\d*)r([1-9]\d*)p([1-9]\d*))?", sys.argv[1]); sys.exit(0 if m and 50 <= int(m.group(2)) <= 300 else 1)' "$DEF_DAY_DE_UPRN" 2>/dev/null \
    || die "DEF_DAY_DE_UPRN must be <s|b>_k<K*100>[_q<n>r<n>p<n>] with K in [0.5, 3], e.g. b_k100_q20r3p10 (got $DEF_DAY_DE_UPRN)"
  grep -q "630: rank-neutral UP" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_DE_UPRN (630)"
fi
'''
E_OLD = '''[ -n "${DEF_DAY_DE_MAG:-}" ] && { export VCC_DE_MAG="$DEF_DAY_DE_MAG"; log "DE rider for the level builds: VCC_DE_MAG=$VCC_DE_MAG"; }
'''
E_NEW = E_OLD + '''[ -n "${DEF_DAY_DE_UPRN:-}" ] && { export VCC_DE_UPRN="$DEF_DAY_DE_UPRN"; log "DE rider for the level builds: VCC_DE_UPRN=$VCC_DE_UPRN"; }   # 630dd: rank-neutral UP pass-through
'''
C_OLD = '''if [ -n "${VCC_RC_CALLSET_BLOCK:-}" ]; then
  python3 - "$FINAL.h5ad.meta.json" "$VCC_RC_CALLSET_BLOCK" <<'PY' || die "the final meta does not carry VCC_RC_CALLSET_BLOCK=$VCC_RC_CALLSET_BLOCK"
'''
C_NEW = f'''if [ -n "${{VCC_DE_UPRN:-}}" ]; then   # 630dd: rank-neutral UP pass-through -- the final meta must carry the variant and the reviewed module
  python3 - "$FINAL.h5ad.meta.json" "$VCC_DE_UPRN" "{UPRN_MODULE_SHA256}" <<'PY' || die "the final meta does not carry VCC_DE_UPRN=$VCC_DE_UPRN with rank_up module {UPRN_MODULE_SHA256[:8]}"
import json, sys
u = json.load(open(sys.argv[1]))["config"].get("de_uprn") or {{}}
sys.exit(0 if u.get("variant") == sys.argv[2] and u.get("module_sha256") == sys.argv[3] else 1)
PY
  # REPORT (never fatal; d4cd5f): per context, the pass's delivery and its fail-safe outcomes. A target whose guard or
  # invariants failed keeps its zs block, which is the safe outcome. The prereg's capacity floor (mean frac_got >= 0.127,
  # moved >= 500 per target) is printed for reading, not enforced here.
  python3 - "$FINAL.rc_audit.jsonl" "$W/uprn_audit.json" "$W/uprn_effectively_off.txt" <<'PY' 2>&1 | while IFS= read -r _l; do log "630dd: $_l"; done || true
import collections, json, sys
# ONE greppable line always ends the report (coordinator 10-02; 607b): `uprn FAILSAFE_TOTAL <n>`, where n = the fail-safe
# triggers (guard + invariant) summed over contexts, and -1 = the audit could not be read or held no ur_ rows. A KA630
# build passes only with the line PRESENT and n == 0.
try:
    rows = collections.defaultdict(list)
    for line in open(sys.argv[1]):
        r = json.loads(line)
        if "ur_n_sel" in r:
            rows[str(r.get("context"))].append(r)
    rep = {{}}
    for c, L in sorted(rows.items()):
        n = len(L)
        fg = sum(r.get("ur_frac_got", 0.0) for r in L) / n
        nm = sum(r.get("ur_n_moved", 0.0) for r in L) / n
        rep[c] = dict(targets=n, mean_frac_got=fg, mean_n_moved=nm,
                      guard_failed=int(sum(r.get("ur_guard_failed", 0.0) > 0 for r in L)),
                      invariant_failed=int(sum(r.get("ur_invariant_failed", 0.0) > 0 for r in L)),
                      no_pos_donor_genes=int(sum(r.get("ur_n_no_pos_donor", 0.0) for r in L)),
                      capacity_floor_met=bool(fg >= 0.127 and nm >= 500))
        print(f"uprn {{c}}: " + json.dumps(rep[c]))
    json.dump(rep, open(sys.argv[2], "w"), indent=1)
    # LOUD, not fatal (d4cd5f): a context where > 10 % of the targets fell back to their zs block has the pass effectively off
    off = [f"{{c}} ({{v['guard_failed'] + v['invariant_failed']}}/{{v['targets']}})" for c, v in rep.items()
           if v["guard_failed"] + v["invariant_failed"] > 0.10 * v["targets"]]
    for o in off:
        print(f"630 EFFECTIVELY OFF in context {{o}}")
    if off:
        open(sys.argv[3], "w").write("\\n".join(off) + "\\n")
    else:
        import os
        if os.path.exists(sys.argv[3]):
            os.remove(sys.argv[3])
    if not rep:
        print("uprn: WARNING no ur_ diagnostics in the final rc_audit dump")
        print("uprn FAILSAFE_TOTAL -1")
    else:
        print(f"uprn FAILSAFE_TOTAL {{sum(v['guard_failed'] + v['invariant_failed'] for v in rep.values())}}")
except Exception as e:   # noqa: BLE001 -- never fatal, never silent
    print(f"uprn: ERROR reading {{sys.argv[1]}}: {{type(e).__name__}}: {{e}}")
    print("uprn FAILSAFE_TOTAL -1")
PY
fi
''' + C_OLD
F_OLD = '''[ -n "${VCC_DE_MAG:-}" ] && printf 'VCC_DE_MAG=%s\\n' "$VCC_DE_MAG" >> "$W/def_day.final.env"
'''
F_NEW = F_OLD + '''[ -n "${VCC_DE_UPRN:-}" ] && printf 'VCC_DE_UPRN=%s\\n' "$VCC_DE_UPRN" >> "$W/def_day.final.env"   # 630dd: rank-neutral UP pass-through
'''
A_OLD = '''[ -n "${VCC_DE_MAG:-}" ] && ALLOW="$ALLOW,de_mag"
'''
S_OLD = '''if [ "${DEF_DAY_KEEP_PINS:-0}" = "1" ]; then
  log "DONE. Next, by hand: vcc prep (real) on $W/${N}_final.h5ad, commit the pins (builder left pinned), vcc submit (-m joke, -d empty)."
'''
S_NEW = '''[ -s "$W/uprn_effectively_off.txt" ] && log "630dd: WARNING 630 EFFECTIVELY OFF (fail-safe on > 10 % of targets) in: $(tr '\\n' ' ' < "$W/uprn_effectively_off.txt")"   # 630dd: rank-neutral UP pass-through
''' + S_OLD
A_NEW = A_OLD + '''[ -n "${VCC_DE_UPRN:-}" ] && ALLOW="$ALLOW,de_uprn"   # 630dd: rank-neutral UP pass-through
'''
EDITS = [("DOC", DOC_OLD, DOC_NEW), ("V validation", V_OLD, V_NEW), ("E export", E_OLD, E_NEW),
         ("C meta check", C_OLD, C_NEW), ("F final env", F_OLD, F_NEW), ("A allow", A_OLD, A_NEW),
         ("S summary", S_OLD, S_NEW)]


def patched_text(day: str, builder: str) -> str:
    if MARK in day:
        raise SystemExit(f"{DAY}: already carries {MARK!r} (a second apply is refused)")
    if MARK598 not in day:
        raise SystemExit(f"{DAY}: lacks 598's mark: apply 598 first (538 -> 554 -> 590 -> 598 -> 612 -> 630 -> 630dd)")
    if MARK630 not in builder:
        raise SystemExit(f"{BUILDER}: lacks 630's mark {MARK630!r}: apply 630 first")
    m = re.search(r'^DE_UPRN_MODULE_SHA256 = "([0-9a-f]{64})"', builder, re.M)
    if m is None or m.group(1) != UPRN_MODULE_SHA256:
        raise SystemExit(f"{BUILDER}: pins rank_up module {m.group(1)[:12] if m else None}, not the reviewed "
                         f"{UPRN_MODULE_SHA256[:12]}: regenerate 630dd for that module")
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
        print(f"630dd check OK: {DAY} carries 598's mark, the builder carries 630's (module {UPRN_MODULE_SHA256[:12]}), "
              f"all {len(EDITS)} anchors found once")
        return 0
    fd, tmp = tempfile.mkstemp(dir=p.parent, prefix=".patch630dd-")
    with os.fdopen(fd, "w") as f:
        f.write(new)
    os.chmod(tmp, os.stat(p).st_mode)
    os.replace(tmp, p)
    print(f"630dd applied: {DAY} ({len(EDITS)} insertions)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
