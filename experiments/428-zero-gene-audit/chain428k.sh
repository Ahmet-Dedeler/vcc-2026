#!/usr/bin/env bash
# 428k = 417k (417b + stack LOAD, q 1.5, its own re-level) + MDE + zero repair: vs 417k exactly de_mde, de_zero_repair,
# emitter sha; vs 428h the bulk step q 1.25 -> 1.5 (with 417k's re-level). Waits for chain428hg to END.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/428-zero-gene-audit; LOG="$E/work/chain428k.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
until grep -q " END$" "$E/work/chain428hg.log" 2>/dev/null; do [ "$(date -u +%Y%m%d%H%M)" -ge 202609232200 ] && { log "chain428hg never ended by 22:00Z"; exit 1; }; sleep 30; done
B=986e8fad8e4a09bb9e698757070157db01479a44c8e437d20fbe828d558cc7b4; EM=72be7e166d8673b64ecebec81f28234c348b02c16879d5ce6937c18785c0f61d
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
log "428k start"
( export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3 VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.0464,B=1.092,C=1.114"
  bash $E/run_arm_generic.sh 428k 1.5 "$REPO/submissions/417k-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stackload-ctxlevel.h5ad.meta.json" \
    the-mitochondria-filed-for-divorce 428k-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stackload-ctxlevel-mde010-zr5 "$B" "$EM" de_mde,de_zero_repair,de_sign_emitter_sha256 ) >> "$E/work/build428k.nohup" 2>&1
log "428k: $(grep -E '^\[428k\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_428k/run428k.log | tr '\n' ' ' | cut -c1-700)"
rm -f $E/work/build_428k/the-mitochondria-filed-for-divorce.h5ad
log "END"
