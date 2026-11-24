"""502: October pass-throughs for the two DE-channel riders, on top of 497's def_day.sh / def_day_helpers.py.
  DEF_DAY_DE_MAG=<variant>          (e.g. dn_k250_zs, 423/429's exact-U pass; board read 09-28 as 446bj2z)
  DEF_DAY_MDE_SIGN=clean4:<tau>     (496's uncalled-gene sign hook; builder >= 24b95dec; board read 10-05 as 496a)
Both are exported for the level builds only (never the stack dump), added to 5b's allowed config diff (de_mag,
de_mde_sign), written to def_day.final.env, and checked on the final meta (`riders_match`). Both unset = no change.
Apply only after the rider's board read says transfers, with an inert known answer (both unset == 497's retrace).
usage: uv run python experiments/502-october-de-riders/patch_defday502.py <def_day.sh> <def_day_helpers.py> [--check]
"""
import os
import sys
from pathlib import Path

dd, hp = Path(sys.argv[1]), Path(sys.argv[2]); s, h = dd.read_text(), hp.read_text()
E_DD = [
    ("#   DEF_DAY_KEEP_PINS=1 (default off): skip the exit-trap restore so `git diff` on the builder is the pin to commit.\n",
     "#   DEF_DAY_DE_MAG=<variant> (502, default off; e.g. dn_k250_zs): the zs exact-U pass (VCC_DE_MAG) on the level builds.\n"
     "#   DEF_DAY_MDE_SIGN=clean4:<tau> (502, default off; e.g. clean4:0.05): 496's uncalled-gene sign hook (VCC_DE_MDE_SIGN).\n"
     "#   DEF_DAY_KEEP_PINS=1 (default off): skip the exit-trap restore so `git diff` on the builder is the pin to commit.\n"),
    ('SQ="$Q"\n',
     '# 502: the DE-channel riders (validated here, exported for the level builds after the dump)\n'
     'if [ -n "${DEF_DAY_DE_MAG:-}" ]; then\n'
     '  case "$DEF_DAY_DE_MAG" in *[!a-z0-9_]*) die "DEF_DAY_DE_MAG must be a variant name like dn_k250_zs (got $DEF_DAY_DE_MAG)" ;; esac\n'
     '  grep -q "VCC_DE_MAG" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_DE_MAG"\n'
     'fi\n'
     'if [ -n "${DEF_DAY_MDE_SIGN:-}" ]; then\n'
     '  case "$DEF_DAY_MDE_SIGN" in clean4:*) ;; *) die "DEF_DAY_MDE_SIGN must be clean4:<tau> (got $DEF_DAY_MDE_SIGN)" ;; esac\n'
     '  grep -q "VCC_DE_MDE_SIGN" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_DE_MDE_SIGN (496p)"\n'
     'fi\n'
     'SQ="$Q"\n'),
    ('export VCC_BULK_LOAD="$DUMP:$DSHA" VCC_BULK_SHARPEN_Q="$SQ"\n',
     'export VCC_BULK_LOAD="$DUMP:$DSHA" VCC_BULK_SHARPEN_Q="$SQ"\n'
     '[ -n "${DEF_DAY_DE_MAG:-}" ] && { export VCC_DE_MAG="$DEF_DAY_DE_MAG"; log "DE rider for the level builds: VCC_DE_MAG=$VCC_DE_MAG"; }\n'
     '[ -n "${DEF_DAY_MDE_SIGN:-}" ] && { export VCC_DE_MDE_SIGN="$DEF_DAY_MDE_SIGN"; log "DE rider for the level builds: VCC_DE_MDE_SIGN=$VCC_DE_MDE_SIGN"; }\n'),
    ('[ -n "${VCC_DE_PRIOR_GATE:-}" ] && printf \'VCC_DE_PRIOR_GATE=%s\\n\' "$VCC_DE_PRIOR_GATE" >> "$W/def_day.final.env"\n',
     '[ -n "${VCC_DE_PRIOR_GATE:-}" ] && printf \'VCC_DE_PRIOR_GATE=%s\\n\' "$VCC_DE_PRIOR_GATE" >> "$W/def_day.final.env"\n'
     '[ -n "${VCC_DE_MAG:-}" ] && printf \'VCC_DE_MAG=%s\\n\' "$VCC_DE_MAG" >> "$W/def_day.final.env"\n'
     '[ -n "${VCC_DE_MDE_SIGN:-}" ] && printf \'VCC_DE_MDE_SIGN=%s\\n\' "$VCC_DE_MDE_SIGN" >> "$W/def_day.final.env"\n'),
    ('  || die "final meta does not match factors $FF q $SQ (the h5ad is not the build this run measured)"\n',
     '  || die "final meta does not match factors $FF q $SQ (the h5ad is not the build this run measured)"\n'
     'uv run python scripts/def_day_helpers.py riders_match "$FINAL.h5ad.meta.json" "${VCC_DE_MAG:--}" "${VCC_DE_MDE_SIGN:--}" \\\n'
     '  || die "the final meta does not carry the requested DE riders (VCC_DE_MAG=${VCC_DE_MAG:-} VCC_DE_MDE_SIGN=${VCC_DE_MDE_SIGN:-})"\n'),
    ('[ "$SQ" != "$Q" ] && ALLOW="$ALLOW,bulk_sharpen_q,bulk_sharpen_pt"\n',
     '[ "$SQ" != "$Q" ] && ALLOW="$ALLOW,bulk_sharpen_q,bulk_sharpen_pt"\n'
     '[ -n "${VCC_DE_MAG:-}" ] && ALLOW="$ALLOW,de_mag"\n'
     '[ -n "${VCC_DE_MDE_SIGN:-}" ] && ALLOW="$ALLOW,de_mde_sign"\n'),
]
E_HP = [
    ('  def_day_helpers.py meta_matches <meta.json> <factors> <q | tercile:q,q,q | stercile:q,q,q>\n',
     '  def_day_helpers.py meta_matches <meta.json> <factors> <q | tercile:q,q,q | stercile:q,q,q>\n'
     '  def_day_helpers.py riders_match <meta.json> <de_mag variant | -> <clean4:tau | ->   (502; "-" = must be absent)\n'),
    ('    if cmd == "meta_matches":\n',
     '    if cmd == "riders_match":\n'
     '        meta_path, mag, sign = args[0], args[1], args[2]\n'
     '        try:\n'
     '            cfg = json.load(open(meta_path))["config"]\n'
     '        except Exception as exc:\n'
     '            print(f"riders_match failed: {exc}", file=sys.stderr)\n'
     '            return 2\n'
     '        got_mag = (cfg.get("de_mag") or {}).get("variant") if isinstance(cfg.get("de_mag"), dict) else None\n'
     '        if (mag == "-") != (got_mag is None) or (mag != "-" and got_mag != mag):\n'
     '            print(f"de_mag mismatch: meta {got_mag!r} vs {mag!r}", file=sys.stderr)\n'
     '            return 1\n'
     '        ms = cfg.get("de_mde_sign")\n'
     '        if sign == "-":\n'
     '            if ms is not None:\n'
     '                print(f"de_mde_sign present but not requested: {ms}", file=sys.stderr)\n'
     '                return 1\n'
     '        else:\n'
     '            src, _, tau = sign.partition(":")\n'
     '            if not isinstance(ms, dict) or ms.get("source") != src or abs(float(ms.get("tau", -1)) - float(tau)) >= 1e-12:\n'
     '                print(f"de_mde_sign mismatch: meta {ms} vs {sign}", file=sys.stderr)\n'
     '                return 1\n'
     '        return 0\n'
     '    if cmd == "meta_matches":\n'),
]
for src, edits, name in ((s, E_DD, "def_day.sh"), (h, E_HP, "def_day_helpers.py")):
    for old, new in edits:
        if src.count(old) != 1:
            raise SystemExit(f"{name}: anchor count {src.count(old)} != 1: {old[:80]!r}")
if "DEF_DAY_DE_MAG" in s or "riders_match" in h:
    raise SystemExit("already patched")
for old, new in E_DD:
    s = s.replace(old, new)
for old, new in E_HP:
    h = h.replace(old, new)
compile(h, str(hp), "exec")
if "--check" in sys.argv:
    print("patches apply cleanly (not written)")
else:
    # 519 finding 2: write both temps, then swap them in, so a kill leaves either the old or the new pair
    tmp_dd, tmp_hp = dd.with_name(dd.name + ".502tmp"), hp.with_name(hp.name + ".502tmp")
    tmp_dd.write_text(s); tmp_hp.write_text(h)
    os.replace(tmp_dd, dd); os.replace(tmp_hp, hp)
    print(f"patched {dd} and {hp}")
