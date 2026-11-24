#!/usr/bin/env bash
# Resume after the 10:39Z 09-23 reboot: chain433hk.sh finished 433h (e7495644) and died at the start of 433k (builder
# was never pinned: sha 2c90d9e1 after reboot). The aborted work dir was moved to build_433k.aborted-reboot-1039Z.
# 433k = 428k with coarse (417k: q 1.5 + its own re-level, coarse MDE, zero repair): vs 428k exactly de_mde + emitter sha.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/428-zero-gene-audit; X=experiments/433-coarse-mde-arm; LOG="$X/work/chain433hk.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
B=2c90d9e1cc2c19fd427da35ffe6a6d1306d3097382ff7a343233423cb710573f; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
log "433k start (resume after reboot)"
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.0464,B=1.092,C=1.114"
  bash $E/run_arm_generic.sh 433k 1.5 "$REPO/submissions/428k-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stackload-ctxlevel-mde010-zr5.h5ad.meta.json" \
    the-freezer-alarm-wrote-poetry 433k-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stackload-ctxlevel-mde010coarse-zr5 "$B" "$EM" de_mde,de_sign_emitter_sha256 ) >> "$X/work/build433k.nohup" 2>&1
log "433k: $(grep -E '^\[433k\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_433k/run433k.log | tr '\n' ' ' | cut -c1-700)"
rm -f $E/work/build_433k/the-freezer-alarm-wrote-poetry.h5ad
log "END"
