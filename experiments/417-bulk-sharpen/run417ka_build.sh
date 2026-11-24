#!/usr/bin/env bash
# 417-KA -- known answer for BOTH new paths: 417a's env + VCC_BULK_LOAD of a same-env DUMP + VCC_CARRIER_GAMMA=0.65 (the
# carriers' second pass at the shipped gamma) must reproduce 417a's artifact X exactly.
# Config diff vs 417a's meta must be exactly bulk_load + bulk_load_sha256 + carrier_gamma; then X is compared array by array.
# Carriers = 395's DE design on 400-K cells; K ballast cells (~999k counts) put each block's per-gene totals on the scaled profile.
# Nothing is submitted here.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
W="$HERE/work/build_ka"; mkdir -p "$W"; LOG="$W/run417ka.log"
BUNDLE="$REPO/data/challenge"
NAME="ka-load-417a"
KA_H5AD="$HERE/work/build_a/the-pigeon-passed-peer-review.h5ad"
REF_META="$HERE/work/build_a/the-pigeon-passed-peer-review.h5ad.meta.json"
LOADSPEC="${1:?<dump.npz>:<sha256>}"; WANT_B="${2:?builder sha}"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
free_gib(){ df -k "$REPO" | awk 'NR==2{printf "%.1f", $4/1048576}'; }
gate_ok(){
  local pct fi others preps disk
  pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  others="$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null || echo 99)"
  preps="$(ps -axo args= | grep -vE 'claude|grep' | grep -cE '(^|/| )vcc prep( |$)' | tr -d ' ')"
  disk="$(free_gib)"
  log "gate: free ${pct}%  free+inactive ${fi} GiB  other builders ${others}  vcc preps ${preps}  disk ${disk} GiB"
  [ "${pct:-0}" -ge 40 ] && awk -v v="$fi" 'BEGIN{exit !(v>=23)}' && [ "$others" = "0" ] \
    && [ "$preps" = "0" ] && awk -v v="$disk" 'BEGIN{exit !(v>=20)}'
}
wait_gate(){ local n=0; until gate_ok; do n=$((n+1)); [ $n -ge 240 ] && { log "GATE never opened"; exit 3; }; sleep 30; done; }
restore(){ [ -f "$W/pin-artifacts.json" ] && uv run python scripts/def_pin_panel_artifacts.py --restore --record "$W/pin-artifacts.json" | tee -a "$LOG"; }

export VCC_NORMREL_EXP=-1.5          # 319's dose
export VCC_CD4_PARTIAL_G=0.5         # 319p2's partial-rows dose
export VCC_NORMREL_MODE=cap          # 319's cap form
export VCC_DE_SIGN_PRIOR=clean4q80    # 395's DE-channel prior (unchanged)
export VCC_PERT_SCALE=wshare_lo0.3_hi3  # 410b key (unchanged)
export VCC_PERT_SCALE_BALLAST=1          # 410b's unwelding (unchanged)
export VCC_BULK_SHARPEN_Q=1.25           # 417a (unchanged)
export VCC_BULK_LOAD="$LOADSPEC"         # KNOWN ANSWER 1: a same-env dump
export VCC_CARRIER_GAMMA=0.65            # KNOWN ANSWER 2: the carriers' second pass at the shipped gamma
unset VCC_NORMREL_PM VCC_RC_ORDER VCC_CTX_A_PRIOR VCC_CTX_A_PRIOR_BETA VCC_SIGN_REPAIR VCC_BULK_GENE_GAIN
log "ARMED pid=$$  VCC_NORMREL_EXP=$VCC_NORMREL_EXP  VCC_NORMREL_MODE=$VCC_NORMREL_MODE  VCC_CD4_PARTIAL_G=$VCC_CD4_PARTIAL_G  VCC_DE_SIGN_PRIOR=$VCC_DE_SIGN_PRIOR  VCC_PERT_SCALE=$VCC_PERT_SCALE  VCC_PERT_SCALE_BALLAST=$VCC_PERT_SCALE_BALLAST  VCC_BULK_SHARPEN_Q=$VCC_BULK_SHARPEN_Q  VCC_BULK_GENE_GAIN/VCC_SIGN_REPAIR/VCC_NORMREL_PM/VCC_RC_ORDER/VCC_CTX_A_PRIOR=<unset>"
log "builder: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)  driver: $(shasum -a 256 scripts/def_pipeline.py | cut -d' ' -f1)  emitter: $(shasum -a 256 src/vcc2026/rc_emit_k.py | cut -d' ' -f1)  git: $(git rev-parse --short HEAD)"
B="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"; D="$(shasum -a 256 scripts/def_pipeline.py | cut -d' ' -f1)"; EM="$(shasum -a 256 src/vcc2026/rc_emit_k.py | cut -d' ' -f1)"
[ "$B" = "$WANT_B" ] || { log "BUILDER hash $B != $WANT_B (417-KA) -- stopping"; exit 4; }
[ "$D" = "27edc486232afff2e19c46bdd87c8c84322cdc4b7781f55c6e5d40adcd8c180f" ] || { log "DRIVER hash $D != 27edc486... -- stopping"; exit 4; }
[ "$EM" = "dbbb91210a2d808fd03ae291ae3dc1479a30e58c675c5dc9536da6c266e5c09c" ] || { log "EMITTER hash $EM != dbbb9121... -- stopping"; exit 4; }
wait_gate
log "DRIVER invocation 1: preflight ... postbuild  (panel A/B/C, real challenge bundle)"
/usr/bin/time -l uv run python scripts/def_pipeline.py \
    --bundle "$BUNDLE" --work "$W" --name "$NAME" \
    --whitener-id union_pruned_def --emitter rc_callset \
    --frozen-inputs "$W/frozen_inputs.json" --min-free-gib 23 \
    --stages preflight,identity,cd4-extract,census,whitener,repin-frozen,pin-artifacts,repin-panel-support,build,postbuild \
    > "$W/driver.log" 2>&1
RC=$?; log "invocation 1 rc=$RC  builder now: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)  disk $(free_gib) GiB"
log 'restoring builder (the pin is finished; the prep does not need it)'
restore
log "builder after restore: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
PB="$W/logs/postbuild.log"
NFAIL="$(grep -c '^  FAIL' "$PB" 2>/dev/null || echo 99)"
log "postbuild: $(grep -hE 'checks pass' "$PB" 2>/dev/null | tail -1)  FAIL lines: $NFAIL"
grep -h '^  FAIL' "$PB" 2>/dev/null | tee -a "$LOG"
if [ ! -s "$W/$NAME.h5ad" ]; then log "NO ARTIFACT -- stopping"; exit 1; fi
if [ "$NFAIL" != "1" ]; then log "UNEXPECTED postbuild failure count ($NFAIL != 1) -- stopping before prep"; exit 1; fi
if ! grep -q '^  FAIL  config differences' "$PB"; then log "the single FAIL is NOT the config-drift line -- stopping"; exit 1; fi
META="$W/$NAME.h5ad.meta.json"
log "meta pert_scale=$(uv run python -c 'import json,sys; m=json.load(open(sys.argv[1])); c=m["config"]; o=m["output"].get("rc_callset",{}); print(c.get("pert_scale_ballast","<ABSENT -- ARM DID NOT LAND>"), {k: o[k] for k in o if k.startswith("ballast")})' "$META")"
log "meta rc_callset prior diagnostics: $(uv run python -c 'import json,sys; m=json.load(open(sys.argv[1])); s=m["output"].get("rc_callset",{}); print({k: s[k] for k in s if k.startswith(("prior_","dose_","steer_","n_classes","sign_flip"))})' "$META")"
log "artifact: $(ls -la "$W/$NAME.h5ad" | awk '{print $5, $9}')  sha256 $(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' "$META")"

# the config diff: exactly bulk_load + bulk_load_sha256 vs 417a's meta; then the arrays.
uv run python - "$META" "$REF_META" "$W/$NAME.h5ad" "$KA_H5AD" <<'PY' | tee -a "$LOG"
import json, sys, h5py, numpy as np
built=json.load(open(sys.argv[1]))["config"]; ref=json.load(open(sys.argv[2]))["config"]
d={k:(ref.get(k,"<absent>"),built.get(k,"<absent>")) for k in sorted(set(ref)|set(built)) if ref.get(k,"<absent>")!=built.get(k,"<absent>")}
for k,(a,b) in d.items(): print(f"[417ka]     {k}: {a!r} -> {b!r}")
okc = set(d)=={"bulk_load","bulk_load_sha256","carrier_gamma"} and d["carrier_gamma"][1]==0.65
m = json.load(open(sys.argv[1])); bl = m.get("effects", {}).get("bulk_load") or {}
print("[417ka] max|loaded-own| =", bl.get("max_abs_diff_vs_own_bulk"))
cg = m.get("effects", {}).get("carrier_gamma") or {}
print("[417ka] carrier row cos vs bulk q =", cg.get("row_cos_vs_bulk_q"), "norm scale q =", cg.get("norm_scale_q"))
okg = cg.get("norm_scale_q") == [1.0, 1.0, 1.0] and cg.get("row_cos_vs_bulk_q", [0])[0] >= 1.0 - 1e-12
same = True
with h5py.File(sys.argv[3], "r") as a, h5py.File(sys.argv[4], "r") as b:
    for key in ("X/indptr", "X/indices", "X/data"):
        n = a[key].shape[0]
        if n != b[key].shape[0]: same = False; print("[417ka] length differs", key); break
        for s in range(0, n, 1 << 26):
            if not np.array_equal(a[key][s:s + (1 << 26)], b[key][s:s + (1 << 26)]): same = False; print("[417ka] differs in", key, "at chunk", s); break
        if not same: break
ok = okc and okg and same and bl.get("max_abs_diff_vs_own_bulk") == 0.0
print("[417ka] PASS: config diff exactly bulk_load(+sha)+carrier_gamma, max|loaded-own| = 0, carriers == bulk rows, X identical to 417a" if ok else f"[417ka] FAIL: config {okc} carriers {okg} same {same}")
sys.exit(0 if ok else 1)
PY
log "DONE"
