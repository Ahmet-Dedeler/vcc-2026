#!/usr/bin/env bash
# 09-25 candidates rebuilt with the count-balanced MDE (426h `coarse`), after 433a's build:
#   433h = 428h with coarse   (417e bulk + stack LOAD, coarse MDE, zero repair, 417h's re-level): vs 428h exactly de_mde + emitter sha
#   433k = 428k with coarse   (417k: q 1.5 + its own re-level, coarse MDE, zero repair):         vs 428k exactly de_mde + emitter sha
# Stops unless 433a's rc dump shows context A's mde_lambda median >= 0.5 (the reason coarse exists).
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/428-zero-gene-audit; X=experiments/433-coarse-mde-arm; LOG="$X/work/chain433hk.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
until grep -qE " DONE$|stopping|NO ARTIFACT|GATE never" "$E/work/build_433a/run433a.log" 2>/dev/null && ! pgrep -f "build433a.sh" >/dev/null; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609232200 ] && { log "433a never finished by 22:00Z 09-23"; exit 1; }; sleep 30; done
grep -q " DONE$" "$E/work/build_433a/run433a.log" || { log "433a build did not finish clean: $(tail -1 $E/work/build_433a/run433a.log)"; exit 1; }
AMED="$(uv run python -c '
import json, numpy as np
R=[json.loads(l) for l in open("'"$X"'/work/rc_audit_433a.jsonl")]
print(round(float(np.median([r["mde_lambda"] for r in R if r["context"]=="A" and "mde_lambda" in r])),4))')"
log "433a dump: context A mde_lambda median $AMED"
awk -v v="$AMED" 'BEGIN{exit !(v>=0.5)}' || { log "context A dose not restored ($AMED < 0.5) -- stopping the chain"; exit 1; }
B=2c90d9e1cc2c19fd427da35ffe6a6d1306d3097382ff7a343233423cb710573f; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
log "433h start"
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.0088,B=1.05,C=1.0736"
  bash $E/run_arm_generic.sh 433h 1.25 "$REPO/submissions/428h-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-mde010-zr5-ctxlevel.h5ad.meta.json" \
    the-pipette-took-a-gap-year 433h-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-mde010coarse-zr5-ctxlevel "$B" "$EM" de_mde,de_sign_emitter_sha256 ) >> "$X/work/build433h.nohup" 2>&1
log "433h: $(grep -E '^\[433h\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_433h/run433h.log | tr '\n' ' ' | cut -c1-700)"
rm -f $E/work/build_433h/the-pipette-took-a-gap-year.h5ad
log "433k start"
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.0464,B=1.092,C=1.114"
  bash $E/run_arm_generic.sh 433k 1.5 "$REPO/submissions/428k-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stackload-ctxlevel-mde010-zr5.h5ad.meta.json" \
    the-freezer-alarm-wrote-poetry 433k-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stackload-ctxlevel-mde010coarse-zr5 "$B" "$EM" de_mde,de_sign_emitter_sha256 ) >> "$X/work/build433k.nohup" 2>&1
log "433k: $(grep -E '^\[433k\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_433k/run433k.log | tr '\n' ' ' | cut -c1-700)"
rm -f $E/work/build_433k/the-freezer-alarm-wrote-poetry.h5ad
log "END"
