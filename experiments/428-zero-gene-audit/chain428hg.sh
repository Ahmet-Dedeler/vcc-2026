#!/usr/bin/env bash
# 428h = 428a + 417h's re-level (one variable vs 428a: bulk_ctx_scale);
# 428g = 417g (417a + re-level) + MDE + zero repair (vs 417g: de_mde, de_zero_repair, emitter sha). Sequential builds.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/428-zero-gene-audit; LOG="$E/work/chain428hg.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
B=986e8fad8e4a09bb9e698757070157db01479a44c8e437d20fbe828d558cc7b4; EM=72be7e166d8673b64ecebec81f28234c348b02c16879d5ce6937c18785c0f61d
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
log "428h start"
( export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3 VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.0088,B=1.05,C=1.0736"
  bash $E/run_arm_generic.sh 428h 1.25 "$REPO/submissions/428a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-mde010-zr5.h5ad.meta.json" \
    the-nucleolus-wants-a-window-seat 428h-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-mde010-zr5-ctxlevel "$B" "$EM" bulk_ctx_scale ) >> "$E/work/build428h.nohup" 2>&1
log "428h: $(grep -E '^\[428h\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_428h/run428h.log | tr '\n' ' ' | cut -c1-600)"
rm -f $E/work/build_428h/the-nucleolus-wants-a-window-seat.h5ad
log "428g start"
( unset VCC_BULK_LOAD; export VCC_DE_MDE=0.1,0.1,0.3 VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.002,B=1.0408,C=1.0648"
  bash $E/run_arm_generic.sh 428g 1.25 "$REPO/submissions/417g-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-ctxlevel.h5ad.meta.json" \
    the-golgi-requested-remote-work 428g-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-ctxlevel-mde010-zr5 "$B" "$EM" de_mde,de_zero_repair,de_sign_emitter_sha256 ) >> "$E/work/build428g.nohup" 2>&1
log "428g: $(grep -E '^\[428g\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_428g/run428g.log | tr '\n' ' ' | cut -c1-600)"
rm -f $E/work/build_428g/the-golgi-requested-remote-work.h5ad
log "END"
