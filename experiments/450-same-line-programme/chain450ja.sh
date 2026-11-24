#!/usr/bin/env bash
# 446bja = 446b with VCC_DE_SIGN_PRIOR=clean4q80ja (context A's block prior = Nadig Jurkat's P_down, 396's artifact d5dee640;
# B/C keep CLEAN4). Config diff vs 446b: exactly de_sign_prior, de_sign_prior_sha256, de_sign_prior_context_A. Serial:
# after the 450 H1 queues (drive450c). Builder b6c3f131, emitter edb5b8dc.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/450-same-line-programme; LOG="$X/work/chain450ja.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
busy(){ ps -axo comm,args | awk '$1 ~ /^bash$/ && /(drive450b?c?\.sh|run448b\.sh|chain449\.sh)/' | grep -q .; }
sleep 5; until ! busy; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609251800 ] && { log "450 queues never ended by 18:00Z"; exit 1; }; sleep 30; done
log "H1 queues and earlier builds are finished: $(tail -1 $X/work/drive450c.log 2>/dev/null)"
B=b6c3f13141e707012588cc354f62758889db9c38c3757f9765ebea37a8a733f8; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
[ "$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)" = "$B" ] || { log "builder is not b6c3f131 -- stopping"; exit 1; }
REF="$REPO/submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
log "446bja start"
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export ARM_DE_SIGN_PRIOR=clean4q80ja VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.1004,B=1.1517,C=1.1720"
  bash $X/run_arm_450.sh 446bja 1.75 "$REF" the-vortex-mixer-wrote-a-sonnet \
    446bja-cap15-partial50-clean4q80ja-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5 "$B" "$EM" \
    de_sign_prior,de_sign_prior_sha256,de_sign_prior_context_A ) >> "$X/work/build446bja.nohup" 2>&1
log "446bja: $(grep -E '^\[446bja\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $X/work/build_446bja/run446bja.log | tr '\n' ' ' | cut -c1-700)"
A="$X/work/build_446bja/the-vortex-mixer-wrote-a-sonnet.h5ad"
[ -s "$A" ] && cp -c "$A" "$X/work/446bja-ka-copy.h5ad" 2>/dev/null && log "kept an APFS clone for 451's known-answer row"
rm -f "$A"
log "END"
