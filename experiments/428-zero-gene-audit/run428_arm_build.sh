#!/usr/bin/env bash
# 428 zero-repair arm -- REBUILD 426a (VCC_BULK_LOAD=<dump> + VCC_DE_MDE) + VCC_DE_ZERO_REPAIR=<floor_log2>.
# usage: run428_arm_build.sh ARM Q DUMP.npz REF_META NAME PKGSTEM BUILDER_SHA EMITTER_SHA MDE FLOOR [CTX]
# One variable against REF_META: exactly de_zero_repair. The audit (audit428.py) must find 0 zeroed gated genes.
# Carriers = 395's DE design on 400-K cells; K ballast cells (~999k counts) put each block's per-gene totals on the scaled profile.
# Nothing is submitted here.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
ARM="${1:?arm}"; QV="${2:?q}"; DUMP="${3:?dump}"; REF_META="${4:?ref meta}"; NAME="${5:?name}"; STEM="${6:?pkg stem}"; WANT_B="${7:?builder sha}"; WANT_EM="${8:?emitter sha}"; MDE="${9:?mde c,cap,min}"; ZR="${10:?zero-repair floor}"; CTX="${11:-}"
DSHA="$(shasum -a 256 "$DUMP" | cut -d' ' -f1)"
W="$HERE/work/build_$ARM"; mkdir -p "$W"; LOG="$W/run$ARM.log"
BUNDLE="$REPO/data/challenge"
PKG="$REPO/submissions/$STEM.vcc"
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
export VCC_BULK_SHARPEN_Q="$QV"          # the reference arm's q (unchanged)
export VCC_BULK_LOAD="$DUMP:$DSHA"       # the reference arm's loaded bulk (unchanged)
export VCC_DE_MDE="$MDE"                 # the reference arm's MDE (unchanged)
export VCC_DE_ZERO_REPAIR="$ZR"          # THE ARM: 428's zero-gene repair
[ -n "$CTX" ] && export VCC_BULK_CTX_SCALE="$CTX"   # optional: the reference arm's re-level (unchanged)
unset VCC_NORMREL_PM VCC_RC_ORDER VCC_CTX_A_PRIOR VCC_CTX_A_PRIOR_BETA VCC_SIGN_REPAIR VCC_BULK_GENE_GAIN VCC_CARRIER_GAMMA VCC_H1_UNION_G VCC_CD4_TARGET_F
[ -z "$CTX" ] && unset VCC_BULK_CTX_SCALE
log "ARMED pid=$$  VCC_DE_ZERO_REPAIR=$VCC_DE_ZERO_REPAIR  VCC_DE_MDE=$VCC_DE_MDE  VCC_BULK_LOAD=$VCC_BULK_LOAD  VCC_BULK_CTX_SCALE=${VCC_BULK_CTX_SCALE:-<unset>}  VCC_NORMREL_EXP=$VCC_NORMREL_EXP  VCC_NORMREL_MODE=$VCC_NORMREL_MODE  VCC_CD4_PARTIAL_G=$VCC_CD4_PARTIAL_G  VCC_DE_SIGN_PRIOR=$VCC_DE_SIGN_PRIOR  VCC_PERT_SCALE=$VCC_PERT_SCALE  VCC_PERT_SCALE_BALLAST=$VCC_PERT_SCALE_BALLAST  VCC_BULK_SHARPEN_Q=$VCC_BULK_SHARPEN_Q  VCC_BULK_GENE_GAIN/VCC_SIGN_REPAIR/VCC_NORMREL_PM/VCC_RC_ORDER/VCC_CTX_A_PRIOR=<unset>"
log "builder: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)  driver: $(shasum -a 256 scripts/def_pipeline.py | cut -d' ' -f1)  emitter: $(shasum -a 256 src/vcc2026/rc_emit_k.py | cut -d' ' -f1)  git: $(git rev-parse --short HEAD)"
B="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"; D="$(shasum -a 256 scripts/def_pipeline.py | cut -d' ' -f1)"; EM="$(shasum -a 256 src/vcc2026/rc_emit_k.py | cut -d' ' -f1)"
[ "$B" = "$WANT_B" ] || { log "BUILDER hash $B != $WANT_B ($ARM) -- stopping"; exit 4; }
[ "$D" = "27edc486232afff2e19c46bdd87c8c84322cdc4b7781f55c6e5d40adcd8c180f" ] || { log "DRIVER hash $D != 27edc486... -- stopping"; exit 4; }
[ "$EM" = "$WANT_EM" ] || { log "EMITTER hash $EM != $WANT_EM -- stopping"; exit 4; }
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

# the config diff, BEFORE the prep: exactly de_mde (+ the emitter pin) vs the reference arm's meta.
uv run python - "$META" "$REF_META" "$ARM" <<'PY' | tee -a "$LOG"
import json,sys
built=json.load(open(sys.argv[1]))["config"]; ref=json.load(open(sys.argv[2]))["config"]
d={k:(ref.get(k,"<absent>"),built.get(k,"<absent>")) for k in sorted(set(ref)|set(built)) if ref.get(k,"<absent>")!=built.get(k,"<absent>")}
for k,(a,b) in d.items(): print(f"[{sys.argv[3]}]     {k}: {a!r} -> {b!r}")
ok = set(d)=={"de_zero_repair"}
print(f"[{sys.argv[3]}] PASS: exactly de_zero_repair" if ok else f"[{sys.argv[3]}] FAIL: got {d}")
sys.exit(0 if ok else 1)
PY
RC=${PIPESTATUS[0]}; [ $RC -ne 0 ] && { log "CONFIG DIFF is not exactly de_zero_repair -- stopping before prep"; exit 1; }
log "zero-gene audit (428) on the built artifact: $(uv run python experiments/428-zero-gene-audit/audit428.py "$W/$NAME.h5ad" 2>&1 | head -1)"
log "meta zero-repair diagnostics: $(uv run python -c 'import json,sys; m=json.load(open(sys.argv[1])); s=m["output"].get("rc_callset",{}); print({k: s[k] for k in s if k.startswith(("zero_","mde_lambda","call_set_unsupp"))})' "$META")"

log "cliffcheck (414): the bulk vs 410b -- must equal the reference arm's (the MDE does not touch the bulk)"
uv run python experiments/414-bulk-reprice/cliffcheck414.py "$W/$NAME.h5ad" 2>&1 | tail -3 | tee -a "$LOG"
uv run python - "$NAME.h5ad" "$ARM" <<'PY' | tee -a "$LOG"
import json,sys,os
r=json.load(open("experiments/414-bulk-reprice/work/cliffcheck414.json"))[sys.argv[1]]
a=r["A_vs_410b"]; lo,hi=(0.85,1.15) if os.environ.get("VCC_BULK_CTX_SCALE") else (0.80,1.03); ok=all(lo<=a[c]<=hi for c in "ABC")
print((f"[{sys.argv[2]}] CLIFF OK " if ok else f"[{sys.argv[2]}] CLIFF FAIL ")+json.dumps({c: round(a[c],4) for c in "ABC"})+f" amp-only d raw vs 410b {r['amp_only_d_raw_mean']:+.5f}")
PY
wait_gate
log "INVOCATION 2: vcc prep -> $PKG"
/usr/bin/time -l vcc prep "$W/$NAME.h5ad" -g "$BUNDLE/gene_names.csv" \
    --perts "$BUNDLE/pert_counts.csv" --contexts A,B,C \
    -o "$PKG" > "$W/prep.log" 2>&1
RC=$?; log "prep rc=$RC  disk $(free_gib) GiB"
tail -20 "$W/prep.log" | tee -a "$LOG"
log "vcc: $(ls -la "$PKG" 2>/dev/null | awk '{print $5, $9}')"
cp "$META" "$REPO/submissions/$STEM.h5ad.meta.json"
log "meta copied beside the package: $(shasum -a 256 "$REPO/submissions/$STEM.h5ad.meta.json" | cut -d' ' -f1)"
uv run python scripts/check_submission_novelty.py check "$PKG" 2>&1 | tee -a "$LOG"
log "novelty rc=${PIPESTATUS[0]}"
log "package sha256: $(shasum -a 256 "$PKG" | cut -d' ' -f1)"
log "git status scripts/ src/: $(git status --short scripts/ src/ | tr '\n' ' ')"
log "DONE"
