#!/usr/bin/env bash
# 429a = 428h (428a + 417h's re-level: stack LOAD, MDE, zero repair) + VCC_DE_MAG=dn_k250_zs, one variable vs 428h
# (config diff exactly de_mag). Also dumps per-block rc diagnostics (431, no config key) for the A/B/C MDE-dose question.
# usage: build429a.sh BUILDER_SHA
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
B="${1:?builder sha}"; EM=72be7e166d8673b64ecebec81f28234c348b02c16879d5ce6937c18785c0f61d
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
E=experiments/428-zero-gene-audit; W=experiments/429-zs-production/work
rm -f "$W/rc_audit_429a.jsonl"
export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3 VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.0088,B=1.05,C=1.0736" \
       VCC_DE_MAG=dn_k250_zs VCC_RC_AUDIT_DUMP="$REPO/$W/rc_audit_429a.jsonl"
bash $E/run_arm_generic.sh 429a 1.25 "$REPO/submissions/428h-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-mde010-zr5-ctxlevel.h5ad.meta.json" \
  the-ribosome-joined-a-book-club 429a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-mde010-zr5-ctxlevel-zs250 "$B" "$EM" de_mag
rm -f $E/work/build_429a/the-ribosome-joined-a-book-club.h5ad
