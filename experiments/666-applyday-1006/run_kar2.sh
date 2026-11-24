#!/usr/bin/env bash
# 666: re-run of KAr only (coordinator, 2026-10-05 ~21:55Z). The first KAr (KA1006r) died on a FULL DISK: OSError Errno 28 while
# writing its h5ad at 21:46:50Z (another session's 9.7 GB model download + a git repack + six 1.3 GB whitener copies, since cloned).
# It is not a 630 defect. Same env and checks as run_kas.sh's KAr block, under a fresh tag and a fresh audit dump.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/666-applyday-1006; LOG=$X/work/kas.log; K=experiments/527-arc-hook-arm; KW=$K/work
NEWB=fe4678c31b3c1f282976ad4fbbc23b44ffe1efc43f151fd832e0b139d1a44a89
DL=$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz; DLS=00ccfcba1b67fa24ff743ee0c64067e01158147ffe49c5b8ca37cef67600b667
F0="A=1.0962,B=1.1461,C=1.1681"
REFZ=$REPO/submissions/446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250.h5ad.meta.json
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/kar2.pid
free_gib(){ df -k "$REPO" | awk 'NR==2{printf "%.0f", $4/1048576}'; }
[ "$(free_gib)" -ge 25 ] || { log "KAr2: disk $(free_gib) GiB < 25 -- not starting"; log "KAr2 FAIL disk"; log "END"; exit 1; }
[ "$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)" = "$NEWB" ] || { log "KAr2: builder != NEWB"; log "KAr2 FAIL builder"; log "END"; exit 1; }
DUMP=$REPO/$X/work/ka1006r2.rc_audit.jsonl; rm -f "$DUMP"
log "KAr2: re-run of KAr after the full-disk failure (disk $(free_gib) GiB)"
( unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_RC_AUDIT_DUMP VCC_DE_MDE_SIGN VCC_DE_MDE_SIGN_KEYGATE VCC_RC_CALLSET_CPM \
        VCC_RC_CALLSET_BLOCK VCC_DE_MAG VCC_BULK_CTX_SCALE VCC_DE_UPRN
  export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F0" VCC_DE_MAG=dn_k250_zs \
         VCC_DE_UPRN=b_k100_q20r3p10 VCC_RC_AUDIT_DUMP="$DUMP" KA_Q=1.75
  bash $K/run_ka_527.sh KA1006r2 "$REFZ" "$NEWB" ) >> $X/work/KA1006r2.nohup 2>&1
log "KAr2: $(/usr/bin/grep -E 'driver rc=|KA PASS|KA FAIL|builder after restore' $KW/build_KA1006r2/runKA1006r2.log | tr '\n' ' ' | cut -c1-300)"
D630="$(python3 - $KW/build_KA1006r2/ka527-KA1006r2.h5ad.meta.json $KW/build_KA1006u/ka527-KA1006u.h5ad.meta.json <<'PY'
import json, sys
try:
    a = json.load(open(sys.argv[1]))["config"]; b = json.load(open(sys.argv[2]))["config"]
except Exception as e:
    print(f"unreadable {e}"); sys.exit()
print(",".join(sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))))
PY
)"
log "KAr2 config diff vs KAu: [$D630]"
if /usr/bin/grep -q "driver rc=0" $KW/build_KA1006r2/runKA1006r2.log && [ "$D630" = "de_uprn" ] \
   && uv run python experiments/630-rank-neutral-up/ka630_audit.py "$DUMP" >> $LOG 2>&1; then log "KAr PASS (re-run KA1006r2: diff == de_uprn; ka630_audit rc 0)"; log "KAS ALL PASS (KAu, KAb, KAc, KAh, KAg; KAr via KA1006r2)"
else log "KAr2 FAIL"; fi
log "END"
