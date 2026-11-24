#!/usr/bin/env bash
# 417 overnight chain (PI, 2026-09-22 21:25Z): after the 417d build restores the builder,
#  1. KA-LOAD: 417a's env + VCC_BULK_LOAD of the same-env dump (dump_ka417a_g065.npz, from window2_chain) + carrier gamma
#     0.65 must reproduce 417a's X exactly (run417ka_build.sh). If it fails, STOP: the LOAD path is not trusted.
#  2. bulk DUMPs for 418's winners (q-agnostic; sharpening is applied at render time in the arm build):
#     h1x2 (VCC_H1_UNION_G=0.56), cd4x15 (VCC_CD4_TARGET_F=0.51), and the stack h1x2_cd4x15.
# Arms built from these dumps are the PI's next call (after 09-23's readings fix the sharpening dose).
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
E="$REPO/experiments/417-bulk-sharpen"; W="$E/work"; LOG="$W/chain_overnight417.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
log "WAITING for 417d (or its launcher's refusal) (pid $$)"
n=0; until grep -qE "builder after restore|stopping|GATE never" "$W/build_d/run417d.log" 2>/dev/null || grep -q "NOT launched" "$W/chain417d.log" 2>/dev/null; do
  n=$((n+1)); [ $n -ge 1440 ] && { log "417d never finished in 12 h"; exit 1; }; sleep 30; done
sleep 60   # let 417d's prep start and its log settle; the builds below wait on the memory/prep gate themselves
B="$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -d' ' -f1)"
[ "$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)" = "$B" ] || { log "builder file != HEAD -- stopping"; exit 1; }
D="$W/dump_ka417a_g065.npz"; [ -s "$D" ] || { log "no KA dump at $D -- stopping"; exit 1; }
DS="$(shasum -a 256 "$D" | cut -d' ' -f1)"
log "1. KA-LOAD with $D ($DS), builder $B"
bash "$E/run417ka_build.sh" "$D:$DS" "$B" >> "$LOG" 2>&1
grep -q "^\[417ka\] PASS" "$W/build_ka/run417ka.log" || { log "KA-LOAD did not PASS -- STOP (see build_ka/run417ka.log)"; exit 1; }
log "KA-LOAD PASS: $(grep '^\[417ka\] PASS' "$W/build_ka/run417ka.log" | tail -1)"
rm -f "$W/build_ka/ka-load-417a.h5ad" && log "removed the KA artifact (12 GB; its X was identical to 417a's)"
base_env(){ export VCC_NORMREL_EXP=-1.5 VCC_CD4_PARTIAL_G=0.5 VCC_NORMREL_MODE=cap VCC_DE_SIGN_PRIOR=clean4q80 VCC_PERT_SCALE=wshare_lo0.3_hi3 VCC_PERT_SCALE_BALLAST=1 VCC_BULK_SHARPEN_Q=1.25
            unset VCC_NORMREL_PM VCC_RC_ORDER VCC_CTX_A_PRIOR VCC_CTX_A_PRIOR_BETA VCC_SIGN_REPAIR VCC_BULK_GENE_GAIN VCC_BULK_CTX_SCALE VCC_BULK_LOAD VCC_CARRIER_GAMMA VCC_H1_UNION_G VCC_CD4_TARGET_F; }
for spec in "h1x2:VCC_H1_UNION_G=0.56" "cd4x15:VCC_CD4_TARGET_F=0.51" "h1x2_cd4x15:VCC_H1_UNION_G=0.56,VCC_CD4_TARGET_F=0.51"; do
  name="${spec%%:*}"; kv="${spec#*:}"
  base_env; IFS=',' read -ra pairs <<< "$kv"; for p in "${pairs[@]}"; do export "$p"; done
  log "2. DUMP $name ($kv)"
  bash "$E/run_dump417.sh" "$name" "$W/dump_$name.npz" "$B" >> "$LOG" 2>&1
  if [ -s "$W/dump_$name.npz" ]; then log "   dump_$name.npz sha256 $(shasum -a 256 "$W/dump_$name.npz" | cut -d' ' -f1)"; else log "   dump $name FAILED"; fi
done
log "END"
