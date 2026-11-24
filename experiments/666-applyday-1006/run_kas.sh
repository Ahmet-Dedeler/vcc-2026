#!/usr/bin/env bash
# 666 (coordinator, written 2026-10-05 ~16:40Z): "10-06 APPLY DAY v2" step 6, the six known-answer builds, serial, on the tree
# with every APPLY DAY patch applied (steps 1-4b) and step 5's builder tests green. Run it only after those, with
#   NEWB=<the post-apply builder sha> bash experiments/666-applyday-1006/run_kas.sh
# The commands are 598 RESULTS "KA plan" and 612 RESULTS's KA block verbatim (fresh tags so no stale build dir is reused), plus
# KA630 per APPLY DAY v2 step 6 and 606 RESULTS (VCC_DE_UPRN on the all-unset recipe with a fresh VCC_RC_AUDIT_DUMP; config diff
# vs the all-unset KA's meta exactly de_uprn; ka630_audit.py exit 0).
#   KAu  all-unset == 446bj2z ec3a5ae6        KAb tercile == 479p2g dd82740d      KAc tercilec -> ka598.py PASS
#   KAh  clean4:0.05 == 496a e21ff752          KAg clean4:0.05 + KEYGATE 0.5 -> ka612.py PASS
#   KAr  VCC_DE_UPRN=b_k100_q20r3p10 on KAu's env -> diff == {de_uprn} and ka630_audit.py rc 0
# Memory: each KA waits for free+inactive >= 25.5 GiB twice 60 s apart; def_pipeline's own 23 GiB floor refuses cleanly, and a
# refused KA is retried after 5 min in a fresh dir. KAu failing stops everything (the bundle would be broken). Cap 10-06 12:00Z.
# Ends with "KAS ALL PASS" or "KAS FAIL <tags>". It never decides or commits anything.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/666-applyday-1006; mkdir -p $X/work; LOG=$X/work/kas.log; K=experiments/527-arc-hook-arm; KW=$K/work
NEWB="${NEWB:?NEWB = the post-apply builder sha}"; CAP=202610061200
DL=$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz; DLS=00ccfcba1b67fa24ff743ee0c64067e01158147ffe49c5b8ca37cef67600b667
F0="A=1.0962,B=1.1461,C=1.1681"; FP="A=1.0957,B=1.1505,C=1.1807"
REFZ=$REPO/submissions/446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250.h5ad.meta.json
REFP=$REPO/submissions/479p2g-cap15-partial50-clean4q80-pertscale3-ballast-sharpenp2g-stacklodo05load-ctxlevelr-mde010coarse-zr5.h5ad.meta.json
REF4=$REPO/submissions/496a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-as05.h5ad.meta.json
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
now(){ date -u +%Y%m%d%H%M; }
fi_gib(){ vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}'; }
busy(){ uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; raise SystemExit(0 if P.other_builders() else 1)' \
        || ps -axo args= | grep -vE 'claude|grep' | grep -qE '(^|/| )vcc prep( |$)'; }
clean(){ unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_RC_AUDIT_DUMP VCC_DE_MDE_SIGN VCC_DE_MDE_SIGN_KEYGATE VCC_RC_CALLSET_CPM \
               VCC_RC_CALLSET_BLOCK VCC_DE_MAG VCC_BULK_CTX_SCALE VCC_DE_UPRN; }
for f in "$DL" "$REFZ" "$REFP" "$REF4"; do [ -s "$f" ] || { log "missing input $f -- stopping"; log "KAS FAIL inputs"; log "END"; exit 1; }; done
[ "$(shasum -a 256 "$DL" | cut -d' ' -f1)" = "$DLS" ] || { log "dump sha mismatch -- stopping"; log "KAS FAIL inputs"; log "END"; exit 1; }
echo $$ > $X/work/kas.pid
log "armed (pid $$): NEWB ${NEWB:0:16}; six KAs serial; memory >= 25.5 GiB twice; cap $CAP"
wait_mem(){ while [ "$(now)" -lt $CAP ]; do a="$(fi_gib)"
    if awk -v v="$a" 'BEGIN{exit !(v>=25.5)}' && ! busy; then sleep 60; b="$(fi_gib)"
      awk -v v="$b" 'BEGIN{exit !(v>=25.5)}' && ! busy && { log "  memory $a then $b GiB"; return 0; }; fi
    sleep 60; done; return 1; }
# run_ka TAG REF : the caller has exported the env; retries on def_pipeline's memory refusal
run_ka(){ local tag="$1" ref="$2" n=0
  while :; do
    wait_mem || { log "$tag: cap reached waiting for memory"; return 1; }
    [ "$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)" = "$NEWB" ] || { log "$tag: builder != NEWB before the build -- stopping"; return 1; }
    [ -e "$KW/build_$tag" ] && mv "$KW/build_$tag" "$KW/build_$tag.aborted-$(date -u +%Y%m%dT%H%M%SZ)"
    n=$((n + 1)); log "$tag: attempt $n (KA_Q=${KA_Q:-1.75})"
    bash $K/run_ka_527.sh "$tag" "$ref" "$NEWB" >> $X/work/$tag.nohup 2>&1
    if /usr/bin/grep -q "refusing" "$KW/build_$tag/driver.log" 2>/dev/null; then log "$tag: memory refusal -- retry in 5 min"; sleep 300; continue; fi
    log "$tag: $(/usr/bin/grep -E 'driver rc=|KA PASS|KA FAIL|builder after restore' $KW/build_$tag/run$tag.log | tr '\n' ' ' | cut -c1-300)"
    return 0
  done; }
ok(){ /usr/bin/grep -q "driver rc=0" "$KW/build_$1/run$1.log" 2>/dev/null; }
same(){ ok "$1" && /usr/bin/grep -q "KA PASS (identical)" "$KW/build_$1/run$1.log"; }
FAILS=""
# 1. KAu: all-unset == 446bj2z (ec3a5ae6)
( clean; export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F0" VCC_DE_MAG=dn_k250_zs
  KA_Q=1.75; export KA_Q; run_ka KA1006u "$REFZ" )
if same KA1006u; then log "KAu PASS (== ec3a5ae6)"; else log "KAu FAIL: the all-unset build does not reproduce 446bj2z -- stopping (bundle suspect)"; log "KAS FAIL KAu"; log "END"; exit 1; fi
# 2. KAb: tercile == 479p2g (dd82740d)
( clean; export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$FP"
  KA_Q=tercile:2.5,2.0,1.75; export KA_Q; run_ka KA1006b "$REFP" )
same KA1006b && log "KAb PASS (== dd82740d)" || { log "KAb FAIL"; FAILS="$FAILS KAb"; }
# 3. KAc: tercilec on 479p2g's env -> ka598.py
( clean; export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$FP"
  KA_Q=tercilec:2.5,2.0,1.75; export KA_Q; run_ka KA1006c "$REFP" )
if ok KA1006c && uv run python experiments/598-tercilec-knob/ka598.py $KW/build_KA1006c/ka527-KA1006c.h5ad.meta.json >> $LOG 2>&1; then log "KAc PASS (ka598.py)"; else log "KAc FAIL"; FAILS="$FAILS KAc"; fi
# 4. KAh: clean4:0.05 == 496a (e21ff752)
( clean; export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F0" VCC_DE_MDE_SIGN=clean4:0.05
  KA_Q=1.75; export KA_Q; run_ka KA1006h "$REF4" )
same KA1006h && log "KAh PASS (== e21ff752)" || { log "KAh FAIL"; FAILS="$FAILS KAh"; }
# 5. KAg: clean4:0.05 + KEYGATE 0.5 -> ka612.py
( clean; export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F0" VCC_DE_MDE_SIGN=clean4:0.05 VCC_DE_MDE_SIGN_KEYGATE=0.5
  KA_Q=1.75; export KA_Q; run_ka KA1006g "$REF4" )
if ok KA1006g && uv run python experiments/612-key-gated-hook/ka612.py $KW/build_KA1006g/ka527-KA1006g.h5ad.meta.json >> $LOG 2>&1; then log "KAg PASS (ka612.py)"; else log "KAg FAIL"; FAILS="$FAILS KAg"; fi
# 6. KAr: 630's knob on KAu's env, fresh audit dump -> config diff == {de_uprn} and ka630_audit.py rc 0
DUMP=$REPO/$X/work/ka1006r.rc_audit.jsonl; rm -f "$DUMP"
( clean; export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F0" VCC_DE_MAG=dn_k250_zs \
         VCC_DE_UPRN=b_k100_q20r3p10 VCC_RC_AUDIT_DUMP="$DUMP"
  KA_Q=1.75; export KA_Q; run_ka KA1006r "$REFZ" )
D630="$(python3 - $KW/build_KA1006r/ka527-KA1006r.h5ad.meta.json $KW/build_KA1006u/ka527-KA1006u.h5ad.meta.json <<'PY'
import json, sys
try:
    a = json.load(open(sys.argv[1]))["config"]; b = json.load(open(sys.argv[2]))["config"]
except Exception as e:
    print(f"unreadable {e}"); sys.exit()
print(",".join(sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))))
PY
)"
log "KAr config diff vs KAu: [$D630]"
if ok KA1006r && [ "$D630" = "de_uprn" ] && uv run python experiments/630-rank-neutral-up/ka630_audit.py "$DUMP" >> $LOG 2>&1; then log "KAr PASS (diff == de_uprn; ka630_audit rc 0)"; else log "KAr FAIL"; FAILS="$FAILS KAr"; fi
[ -z "$FAILS" ] && log "KAS ALL PASS" || log "KAS FAIL$FAILS"
log "END"
