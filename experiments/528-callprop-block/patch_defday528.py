#!/usr/bin/env python3
"""528: October pass-through for the call-propensity CPM block, on scripts/def_day.sh. NOT APPLIED. Needs the builder with
patch_builder528.py (VCC_RC_CALLSET_BLOCK) and the emitter with patch_emitter528.py.

  DEF_DAY_CALLSET_BLOCK=prop[:<swap>[:<steer>]]
                                        (default off) after preflight, compute the day's artifact from D/E/F's OWN control
                                        cells (experiments/528-callprop-block/callblock528.py --bundle $B --contexts $CTXS;
                                        beta 22 and H1's pinned median dropout: nothing fitted on the day) and export
                                        VCC_RC_CALLSET_BLOCK=prop:$W/callblock_def.npz:<sha256>[:<swap>[:<steer>]] (steer
                                        kept by default: see patch_builder528.py) for the lvl0 and
                                        re-level builds (a DE-only knob: pds/mse and so the re-level are untouched). The final
                                        meta must carry it (config.rc_callset_block: source, sha256, swap, steer); the rule-10 diff
                                        allows `rc_callset_block`; def_day.final.env records it.
Unset: byte-identical behaviour (every added branch is guarded by the variable). The edits anchor on 502's lines and on
lines 517's patch_defday517.py leaves intact, so the two patches apply in either order.
usage: uv run python experiments/528-callprop-block/patch_defday528.py <def_day.sh>
"""
import sys
from pathlib import Path

dd = Path(sys.argv[1])
s = dd.read_text()
if "DEF_DAY_CALLSET_BLOCK" in s:
    raise SystemExit("already patched")
VALIDATE = (
    '# 528: the call-propensity CPM block (DE-only; validated here, computed after preflight, exported for the level builds)\n'
    'if [ -n "${DEF_DAY_CALLSET_BLOCK:-}" ]; then\n'
    '  case "$DEF_DAY_CALLSET_BLOCK" in\n'
    '    prop) ;;\n'
    '    prop:*) python3 -c "import sys; p = sys.argv[1].split(\':\'); v = float(p[0]); sys.exit(0 if (0.0 < v <= 1.0 and len(p) <= 2 and (len(p) == 1 or p[1] in (\'kept\', \'key\', \'cpm\'))) else 1)" "${DEF_DAY_CALLSET_BLOCK#prop:}" 2>/dev/null \\\n'
    '              || die "DEF_DAY_CALLSET_BLOCK=prop:<swap>[:<steer>] needs 0 < swap <= 1 and steer kept|key|cpm (got $DEF_DAY_CALLSET_BLOCK)" ;;\n'
    '    *) die "DEF_DAY_CALLSET_BLOCK must be prop or prop:<swap>[:<steer>] (got $DEF_DAY_CALLSET_BLOCK)" ;;\n'
    '  esac\n'
    '  grep -q "VCC_RC_CALLSET_BLOCK" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_RC_CALLSET_BLOCK (528)"\n'
    '  grep -q "block_steer=" src/vcc2026/rc_emit_k.py || die "the emitter has no block_key / block_steer (528)"\n'
    '  [ -f experiments/528-callprop-block/callblock528.py ] || die "experiments/528-callprop-block/callblock528.py is missing"\n'
    'fi\n'
)
COMPUTE = (
    '# 2c. 528: DEF_DAY_CALLSET_BLOCK=prop[:<swap>] -- the day\'s block keys from D/E/F\'s own controls (nothing fitted today)\n'
    'CSB=""\n'
    'if [ -n "${DEF_DAY_CALLSET_BLOCK:-}" ]; then\n'
    '  uv run python experiments/528-callprop-block/callblock528.py --bundle "$B" --contexts "$CTXS" --out "$W/callblock_def.npz" \\\n'
    '      >> "$LOG" 2>&1 || die "callblock528.py (the call-propensity block keys) failed"\n'
    '  CSB="prop:$W/callblock_def.npz:$(shasum -a 256 "$W/callblock_def.npz" | cut -d\' \' -f1)"\n'
    '  case "$DEF_DAY_CALLSET_BLOCK" in prop:*) CSB="$CSB:${DEF_DAY_CALLSET_BLOCK#prop:}" ;; esac\n'
    '  log "call-propensity block keys for $CTXS: VCC_RC_CALLSET_BLOCK=prop:${CSB#prop:$W/}"\n'
    'fi\n\n'
)
EXPORT = '[ -n "$CSB" ] && { export VCC_RC_CALLSET_BLOCK="$CSB"; log "DE rider for the level builds: VCC_RC_CALLSET_BLOCK=${CSB#prop:$W/}"; }\n'
FINALCHECK = (
    'if [ -n "${VCC_RC_CALLSET_BLOCK:-}" ]; then\n'
    '  python3 - "$FINAL.h5ad.meta.json" "$VCC_RC_CALLSET_BLOCK" <<\'PY\' || die "the final meta does not carry VCC_RC_CALLSET_BLOCK=$VCC_RC_CALLSET_BLOCK"\n'
    'import json, sys\n'
    'c = json.load(open(sys.argv[1]))["config"].get("rc_callset_block") or {}\n'
    'p = sys.argv[2].split(":")\n'
    'want = float(p[3]) if len(p) >= 4 else 1.0\n'
    'steer = p[4].lower() if len(p) == 5 else "kept"\n'
    'sys.exit(0 if (c.get("source") == "prop" and c.get("sha256") == p[2] and abs(float(c.get("swap", -1.0)) - want) < 1e-12\n'
    '               and c.get("steer") == steer) else 1)\n'
    'PY\n'
    'fi\n'
)
E = [
    ("#   DEF_DAY_KEEP_PINS=1 (default off): skip the exit-trap restore so `git diff` on the builder is the pin to commit.\n",
     "#   DEF_DAY_CALLSET_BLOCK=prop[:<swap>[:<steer>]] (528, default off): the DE call set's CPM block chosen per context by a\n"
     "#    call-propensity key from the day's OWN controls (log10 control CPM + 22 x excess dropout, median-rescaled to H1's\n"
     "#    dropout scale; callblock528.py), exported as VCC_RC_CALLSET_BLOCK for the level builds. Builder/emitter >= 528's.\n"
     "#   DEF_DAY_KEEP_PINS=1 (default off): skip the exit-trap restore so `git diff` on the builder is the pin to commit.\n"),
    ('  grep -q "VCC_DE_MDE_SIGN" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_DE_MDE_SIGN (496p)"\nfi\n',
     '  grep -q "VCC_DE_MDE_SIGN" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_DE_MDE_SIGN (496p)"\nfi\n'
     + VALIDATE),
    ('# 2b. DEF_DAY_LODO=1 (460; default off)', COMPUTE + '# 2b. DEF_DAY_LODO=1 (460; default off)'),
    ('[ -n "$SL" ] && export VCC_SAMELINE="$SL"\n', EXPORT + '[ -n "$SL" ] && export VCC_SAMELINE="$SL"\n'),
    ('for P in "${POINTS[@]}"; do [ "${P%%.ntclevel.json:*}" != "$FINAL" ] && rm -f',
     FINALCHECK + 'for P in "${POINTS[@]}"; do [ "${P%%.ntclevel.json:*}" != "$FINAL" ] && rm -f'),
    ('log "final: ${FINAL##*/} factors $FF -> $W/${N}_final.h5ad"\n',
     '[ -n "${VCC_RC_CALLSET_BLOCK:-}" ] && printf \'VCC_RC_CALLSET_BLOCK=%s\\n\' "$VCC_RC_CALLSET_BLOCK" >> "$W/def_day.final.env"\n'
     'log "final: ${FINAL##*/} factors $FF -> $W/${N}_final.h5ad"\n'),
    ('[ -n "${VCC_DE_MDE_SIGN:-}" ] && ALLOW="$ALLOW,de_mde_sign"\n',
     '[ -n "${VCC_DE_MDE_SIGN:-}" ] && ALLOW="$ALLOW,de_mde_sign"\n'
     '[ -n "${VCC_RC_CALLSET_BLOCK:-}" ] && ALLOW="$ALLOW,rc_callset_block"\n'),
]
for old, new in E:
    n = s.count(old)
    if n != 1:
        raise SystemExit(f"anchor found {n} times (want 1): {old[:90]!r}")
    s = s.replace(old, new)
dd.write_text(s)
print(f"patched {dd}: DEF_DAY_CALLSET_BLOCK")
