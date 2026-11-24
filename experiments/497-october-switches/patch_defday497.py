"""497: the October pass-throughs for the two pending bulk switches, on top of 476's def_day.sh / def_day_helpers.py.
  DEF_DAY_LODO_HI_PCT=<P>  (needs DEF_DAY_LODO=1): def_lodo_key.py --hi-pct P; pin id lodo_lo0.5_pct_def.
  DEF_DAY_SHARPEN_PT=<tercile|stercile>:<q>,<q>,<q>: the level builds take VCC_BULK_SHARPEN_Q's per-target form instead
      of the scalar q; meta_matches checks bulk_sharpen_pt; 5b's allowed diff gains bulk_sharpen_q, bulk_sharpen_pt.
Both unset = the script's behaviour is unchanged (same exports, same argv, same checks).
Also (510): an empty flags argument ("") now means no flags (`${6-pt:noprog}`); unset keeps pt:noprog.
Also (499): F2, de_sign_emitter_sha256 is allowed in 5b only when the builder's emitter pin differs from the incumbent's;
F3, the closing meta line prints bulk_sharpen_pt.
usage: uv run python experiments/497-october-switches/patch_defday497.py <def_day.sh> <def_day_helpers.py> [--check]
"""
import sys
from pathlib import Path

dd, hp = Path(sys.argv[1]), Path(sys.argv[2]); s, h = dd.read_text(), hp.read_text()
E_DD = [
    # 487/510: an EMPTY flags argument must mean "no flags" (w 0.5 with the programme source and no pt); `:-` turned
    # it into the default pt:noprog. Unset keeps the default.
    ('SLW="${5:-auto}"; SLF="${6:-pt:noprog}"; MAXL="${7:-4}"',
     'SLW="${5:-auto}"; SLF="${6-pt:noprog}"; MAXL="${7:-4}"'),
    ("#   DEF_DAY_KEEP_PINS=1 (default off): skip the exit-trap restore so `git diff` on the builder is the pin to commit.\n",
     "#   DEF_DAY_LODO_HI_PCT=<P> (497, default off; needs DEF_DAY_LODO=1): the LODO key's cap is the panel's own P-th\n"
     "#    percentile of carried key / median (def_lodo_key.py --hi-pct; 543/546), pinned as lodo_lo0.5_pct_def.\n"
     "#    P 96.64697522769474 is exactly cap 6 on A/B/C (495's artifact); the rounded 96.647 is not.\n"
     "#   DEF_DAY_SHARPEN_PT=<tercile|stercile>:<q>,<q>,<q> (497, default off; e.g. stercile:2.5,2.0,1.75 = 479's S_P2): the\n"
     "#    level builds sharpen each target's bulk with its tercile's q (VCC_BULK_SHARPEN_Q's per-target form, builder >=\n"
     "#    d09a986c) instead of the scalar q. The stack dump is unchanged (q is a render knob).\n"
     "#   DEF_DAY_KEEP_PINS=1 (default off): skip the exit-trap restore so `git diff` on the builder is the pin to commit.\n"),
    ('  *) die "DEF_DAY_CD4_THRESHOLDS must be relaxed or unset (got $CD4T)" ;;\nesac\n',
     '  *) die "DEF_DAY_CD4_THRESHOLDS must be relaxed or unset (got $CD4T)" ;;\nesac\n'
     '# 497: the October forms of 495 (percentile cap) and 479 (per-target sharpening); both default off\n'
     'if [ -n "${DEF_DAY_LODO_HI_PCT:-}" ]; then\n'
     '  [ "${DEF_DAY_LODO:-0}" = "1" ] || die "DEF_DAY_LODO_HI_PCT=$DEF_DAY_LODO_HI_PCT needs DEF_DAY_LODO=1"\n'
     '  python3 -c "import sys; p = float(sys.argv[1]); sys.exit(0 if 50.0 < p < 100.0 else 1)" "$DEF_DAY_LODO_HI_PCT" 2>/dev/null \\\n'
     '    || die "DEF_DAY_LODO_HI_PCT must be a percentile in (50, 100) (got $DEF_DAY_LODO_HI_PCT)"\n'
     '  grep -q -- "--hi-pct" scripts/def_lodo_key.py || die "def_lodo_key.py has no --hi-pct (497)"\n'
     'fi\n'
     'SQ="$Q"\n'
     'if [ -n "${DEF_DAY_SHARPEN_PT:-}" ]; then\n'
     '  case "$DEF_DAY_SHARPEN_PT" in tercile:*,*,*|stercile:*,*,*) ;; *) die "DEF_DAY_SHARPEN_PT must be tercile:<q>,<q>,<q> or stercile:<q>,<q>,<q> (got $DEF_DAY_SHARPEN_PT)" ;; esac\n'
     '  grep -q "def sharpen_q_by_tercile" scripts/make_submission_pm_cd4.py || die "the builder has no per-target sharpening (479)"\n'
     '  SQ="$DEF_DAY_SHARPEN_PT"\n'
     '  log "per-target sharpening for the level builds: VCC_BULK_SHARPEN_Q=$SQ (scalar q $Q not used)"\n'
     'fi\n'),
    ('  LK="$W/pert_scale_lodo_def_lo0.5_hi3.npz"\n'
     '  uv run python scripts/def_lodo_key.py --bundle "$B" --op "$W/whitener_union_pruned_def.npz" --cd4-root "$W/cd4_compact" \\\n'
     '      --out "$LK" >> "$LOG" 2>&1 || die "def_lodo_key.py failed"\n',
     '  LK="$W/pert_scale_lodo_def_lo0.5_hi3.npz"; LID=lodo_lo0.5_hi3_def\n'
     '  [ -n "${DEF_DAY_LODO_HI_PCT:-}" ] && { LK="$W/pert_scale_lodo_def_lo0.5_pct.npz"; LID=lodo_lo0.5_pct_def; }\n'
     '  uv run python scripts/def_lodo_key.py --bundle "$B" --op "$W/whitener_union_pruned_def.npz" --cd4-root "$W/cd4_compact" \\\n'
     '      ${DEF_DAY_LODO_HI_PCT:+--hi-pct "$DEF_DAY_LODO_HI_PCT"} --out "$LK" >> "$LOG" 2>&1 || die "def_lodo_key.py failed"\n'),
    ('      --pert-scale "$LK" --pert-scale-id lodo_lo0.5_hi3_def --apply --out "$W/pin-lodo.json" >> "$LOG" 2>&1 \\\n'
     '      || die "LODO key pin failed"\n'
     '  LODO=lodo_lo0.5_hi3_def\n'
     '  log "LODO key ON for the stack dump: $LODO $(shasum -a 256 "$LK" | cut -c1-16)"\n',
     '      --pert-scale "$LK" --pert-scale-id "$LID" --apply --out "$W/pin-lodo.json" >> "$LOG" 2>&1 \\\n'
     '      || die "LODO key pin failed"\n'
     '  LODO="$LID"\n'
     '  log "LODO key ON for the stack dump: $LODO $(shasum -a 256 "$LK" | cut -c1-16); cap $(python3 -c \'import json, sys; print(json.load(open(sys.argv[1]))["hi"])\' "${LK%.npz}.json")"\n'),
    ('export VCC_BULK_LOAD="$DUMP:$DSHA" VCC_BULK_SHARPEN_Q="$Q"\n',
     'export VCC_BULK_LOAD="$DUMP:$DSHA" VCC_BULK_SHARPEN_Q="$SQ"\n'),
    ('uv run python scripts/def_day_helpers.py meta_matches "$FINAL.h5ad.meta.json" "$FF" "$Q" \\\n'
     '  || die "final meta does not match factors $FF q $Q (the h5ad is not the build this run measured)"\n',
     'uv run python scripts/def_day_helpers.py meta_matches "$FINAL.h5ad.meta.json" "$FF" "$SQ" \\\n'
     '  || die "final meta does not match factors $FF q $SQ (the h5ad is not the build this run measured)"\n'),
    ("printf 'VCC_BULK_LOAD=%s\\nVCC_BULK_SHARPEN_Q=%s\\nVCC_BULK_CTX_SCALE=%s\\nVCC_SAMELINE=%s\\n' \"$VCC_BULK_LOAD\" \"$Q\" \"$FF\" \"${SL:-}\" > \"$W/def_day.final.env\"\n",
     "printf 'VCC_BULK_LOAD=%s\\nVCC_BULK_SHARPEN_Q=%s\\nVCC_BULK_CTX_SCALE=%s\\nVCC_SAMELINE=%s\\n' \"$VCC_BULK_LOAD\" \"$SQ\" \"$FF\" \"${SL:-}\" > \"$W/def_day.final.env\"\n"),
    ('[ -n "${VCC_DE_PRIOR_GATE:-}" ] && ALLOW="$ALLOW,de_prior_gate"\n',
     '[ -n "${VCC_DE_PRIOR_GATE:-}" ] && ALLOW="$ALLOW,de_prior_gate"\n'
     '[ "$SQ" != "$Q" ] && ALLOW="$ALLOW,bulk_sharpen_q,bulk_sharpen_pt"\n'
     '# 499 F2: 496\'s emitter pin moves config.de_sign_emitter_sha256 on every build, knob unset included. Allow that one key\n'
     '# only when the builder source pin differs from the incumbent\'s recorded sha (today both are edb5b8dc, so nothing is\n'
     '# added). A build whose recorded sha moved while the source pin still matches the incumbent still fails the gate.\n'
     '_inc_em="$(python3 -c \'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("de_sign_emitter_sha256",""))\' "$INC" 2>/dev/null || true)"\n'
     '_src_em="$(python3 -c \'import re,sys; t=open(sys.argv[1]).read(); m=re.search(r"^DE_SIGN_EMITTER_SHA256 = .([0-9a-f]{64}).", t, re.M); print(m.group(1) if m else "")\' scripts/make_submission_pm_cd4.py 2>/dev/null || true)"\n'
     '[ -n "$_src_em" ] && [ -n "$_inc_em" ] && [ "$_src_em" != "$_inc_em" ] && ALLOW="$ALLOW,de_sign_emitter_sha256"\n'),
    # 499 F3: the closing meta line shows the per-target form
    ('print("meta: bulk_ctx_scale", c.get("bulk_ctx_scale"), "| q", c.get("bulk_sharpen_q"), "| sameline", json.dumps(c.get("sameline"))[:300])\n',
     'print("meta: bulk_ctx_scale", c.get("bulk_ctx_scale"), "| q", c.get("bulk_sharpen_q"), "| pt", c.get("bulk_sharpen_pt"), "| sameline", json.dumps(c.get("sameline"))[:300])\n'),
]
E_HP = [
    ('  def_day_helpers.py meta_matches <meta.json> <factors> <q>\n',
     '  def_day_helpers.py meta_matches <meta.json> <factors> <q | tercile:q,q,q | stercile:q,q,q>\n'),
    ('        meta_path, fact, q = args[0], args[1], float(args[2])\n'
     '        try:\n'
     '            cfg = json.load(open(meta_path))["config"]\n'
     '            want = factors(fact)\n'
     '            got_q = cfg.get("bulk_sharpen_q")\n'
     '            if got_q is None or abs(float(got_q) - q) >= 5e-5:\n'
     '                print(f"q mismatch: meta {got_q} vs {q}", file=sys.stderr)\n'
     '                return 1\n',
     '        meta_path, fact, qarg = args[0], args[1], args[2]\n'
     '        try:\n'
     '            cfg = json.load(open(meta_path))["config"]\n'
     '            want = factors(fact)\n'
     '            got_q = cfg.get("bulk_sharpen_q")\n'
     '            if qarg.split(":", 1)[0] in ("tercile", "stercile"):\n'
     '                # 497: the per-target form (479): the scalar q is absent and bulk_sharpen_pt carries mode and qs\n'
     '                mode, qs = qarg.split(":", 1)\n'
     '                want_q = [float(x) for x in qs.split(",")]\n'
     '                pt = cfg.get("bulk_sharpen_pt")\n'
     '                if (got_q is not None or not isinstance(pt, dict) or pt.get("mode") != mode or len(pt.get("q", [])) != len(want_q)\n'
     '                        or any(abs(float(a) - b) >= 5e-5 for a, b in zip(pt["q"], want_q))):\n'
     '                    print(f"per-target sharpening mismatch: meta q {got_q}, pt {pt} vs {qarg}", file=sys.stderr)\n'
     '                    return 1\n'
     '            else:\n'
     '                q = float(qarg)\n'
     '                if got_q is None or abs(float(got_q) - q) >= 5e-5:\n'
     '                    print(f"q mismatch: meta {got_q} vs {q}", file=sys.stderr)\n'
     '                    return 1\n'),
]
for src, edits, name in ((s, E_DD, "def_day.sh"), (h, E_HP, "def_day_helpers.py")):
    for old, new in edits:
        if src.count(old) != 1:
            raise SystemExit(f"{name}: anchor count {src.count(old)} != 1: {old[:80]!r}")
if "DEF_DAY_SHARPEN_PT" in s or "tercile" in h:
    raise SystemExit("already patched")
for old, new in E_DD:
    s = s.replace(old, new)
for old, new in E_HP:
    h = h.replace(old, new)
compile(h, str(hp), "exec")
if "--check" in sys.argv:
    print("patches apply cleanly (not written)")
else:
    dd.write_text(s); hp.write_text(h)
    print(f"patched {dd} and {hp}")
