#!/usr/bin/env bash
# 433a = 417e + VCC_DE_MDE=0.1,0.1,0.3,coarse (426h's count-balanced MDE), one variable vs 417e: config diff exactly
# de_mde + de_sign_emitter_sha256. It is 426a with context A's MDE dose restored (429a's rc dump: 426a-style MDE gives
# mde_lambda 0 in all 300 A blocks). Dumps per-block rc diagnostics (431, no config key) to check A's dose on the real data.
# usage: build433a.sh BUILDER_SHA
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
B="${1:?builder sha}"; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
E=experiments/428-zero-gene-audit; W=experiments/433-coarse-mde-arm/work
mkdir -p "$W"; rm -f "$W/rc_audit_433a.jsonl"
unset VCC_DE_ZERO_REPAIR VCC_DE_MAG VCC_BULK_CTX_SCALE
export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_RC_AUDIT_DUMP="$REPO/$W/rc_audit_433a.jsonl"
bash $E/run_arm_generic.sh 433a 1.25 "$REPO/submissions/417e-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload.h5ad.meta.json" \
  the-centrifuge-wrote-a-memoir 433a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-mde010coarse "$B" "$EM" \
  de_mde,de_sign_emitter_sha256
rm -f $E/work/build_433a/the-centrifuge-wrote-a-memoir.h5ad
