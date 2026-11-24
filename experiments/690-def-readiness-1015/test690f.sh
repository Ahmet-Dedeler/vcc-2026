#!/usr/bin/env bash
# 690f tests: arm1015f.sh's FINAL_REUPLOAD branch with a stub `vcc` on PATH, a scratch ledger copy and a scratch
# board-entries dir. The package is a real ledgered one (630arm's, entry QJRJK6OC6Zw1kjwnZqid; only its .vcc sha is
# hashed, no content digest). Never touches the network; restores every file it could touch.
#   bash experiments/690-def-readiness-1015/test690f.sh
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
A=experiments/670-coord-1015/arm1015f.sh
PKG=submissions/446bj2z-uprn-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-uprnk100.vcc
OF=QJRJK6OC6Zw1kjwnZqid; NEWID=TESTREUPLOADTESTREUP
TD="$(mktemp -d /tmp/claude-501/t690f.XXXX)"; mkdir -p "$TD/bin" "$TD/be"
SHA="$(shasum -a 256 "$PKG" | cut -d' ' -f1)"
RES="$(python3 -c 'import re,sys; m=re.search(r"\*\*Final Nov 5 entry:\*\* `([^`]+)`", open(sys.argv[1]).read()); print(m.group(1))' docs/MODEL-NAMES.md)"
cat > "$TD/bin/vcc" <<EOF
#!/usr/bin/env bash
case "\$*" in
  "--json whoami") echo '{"identity": {"can_submit": true}}' ;;
  "--version") echo "vcc 0.2.2 (690f stub)" ;;
  submit*) echo "\$*" >> "$TD/submits.txt"; echo "entry $NEWID"; echo "scored: 0.1 (stub)" ;;
  "--json status $NEWID") python3 -c 'import json; e={"entry_id": "$NEWID", "status": "published", "is_final": True, "partition": "test", "panel_id": "vcc2026-test-1", "anchor_version": "t", "rank": None, "score_avg": 0.1}
for k in ("pds_cosine","expr_mse_unbiased_capped_norm","de_wilcoxon_lfc_nmae","de_wilcoxon_direction_fidelity_yield_raw","de_wilcoxon_direction_reach_raw","de_wilcoxon_sig_jaccard","score_pds","score_mse","score_nmae","score_fid","score_reach","score_jac"): e[k]=0.1
print(json.dumps(e))' ;;
  *) echo "stub: unhandled \$*" >&2; exit 9 ;;
esac
EOF
chmod +x "$TD/bin/vcc"
cp submissions/SUBMITTED_LEDGER.json "$TD/ledger.json"
LED0="$(shasum -a 256 submissions/SUBMITTED_LEDGER.json | cut -d' ' -f1)"
python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); d["entry"]["is_final"]=True; d["entry"]["partition"]="test"; json.dump(d, open(sys.argv[2],"w"))' \
  "data/board-entries/$OF.json" "$TD/be/$OF.json"
printf '%s READY t690f %s\n' "$(date -u +%FT%TZ)" "$SHA" > "$TD/ready.log"
printf '[t690f] PASS\n[t690f] CLIFF OK\nnovelty rc=0\n' > "$TD/run.log"
printf '%s END\n' "$(date -u +%FT%TZ)" > "$TD/prior.log"
NOWM="$(date -u +%Y%m%d%H%M)"; CAPM="$(date -u -v+20M +%Y%m%d%H%M)"
pass=0; fail=0
case_(){  # $1 label, $2 expect regex in the log, then env assignments
  local L="$1" X="$2"; shift 2; local SL="t690f-$L"; local LG="experiments/670-coord-1015/work/submit-$SL.log"
  rm -f "$LG" "$TD/submits.txt"
  env PATH="$TD/bin:$PATH" LEDGER691="$TD/ledger.json" BOARD_ENTRIES691="$TD/be" FINALPHASE_FROM=202601010000 "$@" \
    bash "$A" "$SL" t690f "$TD/ready.log" "$TD/run.log" "$PKG" "${NAME690F:-$RES}" "690f test" "$TD/prior.log" "$NOWM" "$CAPM" > /dev/null 2>&1
  if /usr/bin/grep -qE "$X" "$LG"; then pass=$((pass+1)); echo "PASS $L: $(/usr/bin/grep -oE "$X" "$LG" | head -1 | cut -c1-120)"
  else fail=$((fail+1)); echo "FAIL $L (want /$X/): $(tail -3 "$LG" | tr '\n' ' ' | cut -c1-300)"; fi
  cp "$LG" "$TD/$SL.log"; rm -f "$LG" "experiments/670-coord-1015/work/submit-$SL.pid"
}
# positive: final capture, reserved name, final phase -> submits; records a #reupload key; the original record is kept
case_ reupload "ledger: re-upload of $OF recorded as $NEWID" FINAL_REUPLOAD=1
python3 - "$TD/ledger.json" "$SHA" "$NEWID" "$OF" <<'PY' && { pass=$((pass+1)); echo "PASS ledger: #reupload key added, original kept"; } || { fail=$((fail+1)); echo "FAIL ledger"; }
import json, sys
L, sha, new, of = sys.argv[1:5]; d = json.load(open(L))
r = d.get(f"{sha}#reupload-{new}")
sys.exit(0 if r and r["reupload_of"] == of and d[sha]["entry_id"] == of else 1)
PY
/usr/bin/grep -q -- "-m $RES" "$TD/submits.txt" 2>/dev/null; true
rm -f "data/board-entries/$NEWID.json"
# negatives
case_ abcpkg "NOT SUBMITTING: FINAL_REUPLOAD: ledgered entry $OF has no FINAL-phase capture" FINAL_REUPLOAD=1 BOARD_ENTRIES691=data/board-entries
NAME690F="some queue name" case_ wrongname "NOT SUBMITTING: FINAL_REUPLOAD needs the reserved final name" FINAL_REUPLOAD=1
case_ early "NOT SUBMITTING: FINAL_REUPLOAD before the final phase" FINAL_REUPLOAD=1 FINALPHASE_FROM=209912312359
NAME690F="some queue name" case_ testflag "FINAL_REUPLOAD_TEST=1: the early identical-bytes test" FINAL_REUPLOAD=1 FINAL_REUPLOAD_TEST=1
rm -f "data/board-entries/$NEWID.json"
case_ normalpath "NOT SUBMITTING: novelty check refused" FINAL_REUPLOAD=0
# the real ledger must be untouched
[ "$(shasum -a 256 submissions/SUBMITTED_LEDGER.json | cut -d' ' -f1)" = "$LED0" ] && { pass=$((pass+1)); echo "PASS the real ledger is untouched"; } \
  || { fail=$((fail+1)); echo "FAIL the real ledger CHANGED: restore it from git"; }
echo "$pass PASS, $fail FAIL"
rm -rf "$TD"
[ "$fail" = 0 ]
