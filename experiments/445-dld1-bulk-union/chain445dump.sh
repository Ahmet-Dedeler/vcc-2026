#!/usr/bin/env bash
# 445 dumps, after 446's build ends: (1) KA -- the stack dump remade with TODAY's builder/driver must equal 417's
# dump_h1x2_cd4x15.npz array for array (bulk, carriers, perts, genes); (2) the DLD-1 rider dump = the same env +
# VCC_UNION2_SOURCES=dld1_lowmoi VCC_UNION2_G=0.2 (248's second-union slot; 445's PLOLO dose). No arm is built here.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/445-dld1-bulk-union; W="$REPO/$X/work"; LOG="$W/chain445dump.log"; E=experiments/417-bulk-sharpen
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
B=63dfd7d0e402738854072e57c81777a71fd773212877a6f26540a363b2aa34d3
until grep -q " END$" experiments/446-q175-level/work/build446.log 2>/dev/null; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609250600 ] && { log "446 never ended by 06:00Z 09-25"; exit 1; }; sleep 30; done
log "446 ended; starting the KA dump"
base_env(){ export VCC_NORMREL_EXP=-1.5 VCC_CD4_PARTIAL_G=0.5 VCC_NORMREL_MODE=cap VCC_DE_SIGN_PRIOR=clean4q80 VCC_PERT_SCALE=wshare_lo0.3_hi3 VCC_PERT_SCALE_BALLAST=1 VCC_BULK_SHARPEN_Q=1.25 VCC_H1_UNION_G=0.56 VCC_CD4_TARGET_F=0.51
            unset VCC_NORMREL_PM VCC_RC_ORDER VCC_CTX_A_PRIOR VCC_CTX_A_PRIOR_BETA VCC_SIGN_REPAIR VCC_BULK_GENE_GAIN VCC_BULK_CTX_SCALE VCC_BULK_LOAD VCC_CARRIER_GAMMA VCC_UNION2_SOURCES VCC_UNION2_G VCC_DE_MAG VCC_RC_AUDIT_DUMP; }
( base_env; bash $E/run_dump417.sh ka445 "$W/dump_ka445_h1x2_cd4x15.npz" "$B" ) >> "$W/dump_ka445.nohup" 2>&1
uv run python - "$W/dump_ka445_h1x2_cd4x15.npz" "$REPO/$E/work/dump_h1x2_cd4x15.npz" <<'PY' | tee -a "$LOG"
import sys, numpy as np
a, b = np.load(sys.argv[1]), np.load(sys.argv[2])
ok = all(np.array_equal(a[k], b[k]) for k in ("bulk", "carriers", "perts", "genes"))
print(("KA445 PASS: the remade stack dump equals 417's array for array" if ok else "KA445 FAIL: " +
       str({k: float(np.max(np.abs(a[k].astype(float) - b[k].astype(float)))) for k in ("bulk", "carriers")})))
sys.exit(0 if ok else 1)
PY
[ ${PIPESTATUS[0]} -eq 0 ] || { log "KA445 failed -- no DLD-1 dump"; exit 1; }
rm -f "$W/dump_ka445_h1x2_cd4x15.npz"
log "KA passed; starting the DLD-1 union-2 dump (g 0.2)"
( base_env; export VCC_UNION2_SOURCES=dld1_lowmoi VCC_UNION2_G=0.2
  bash $E/run_dump417.sh dld1u2g020 "$W/dump_h1x2_cd4x15_dld1u2g020.npz" "$B" ) >> "$W/dump_dld1.nohup" 2>&1
[ -s "$W/dump_h1x2_cd4x15_dld1u2g020.npz" ] && log "DLD-1 dump sha256 $(shasum -a 256 "$W/dump_h1x2_cd4x15_dld1u2g020.npz" | cut -d' ' -f1)" || log "DLD-1 dump FAILED"
log "END"
