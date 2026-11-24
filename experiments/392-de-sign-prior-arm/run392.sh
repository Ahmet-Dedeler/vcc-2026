#!/usr/bin/env bash
# 392 -- BUILD the DE-channel block sign prior arm (VCC_DE_SIGN_PRIOR=clean4, 390/391) on the 319p2 composition, REAL A/B/C panel.
# Modelled on run387.sh. One knob against entry lqWNXX5YaYArdglu6Stn (319p2): de_sign_* keys absent -> clean4 (five keys).
# sign_repair absent -> cons. Nothing is submitted here and `vcc submit` is never run.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
W="$HERE/work/build"; mkdir -p "$W"; LOG="$W/run392.log"
BUNDLE="$REPO/data/challenge"
NAME="the-stapler-went-on-sabbatical"
PKG="$REPO/submissions/392-cap15-partial50-designclean4.vcc"
REF_META="$REPO/submissions/319p2-cap15-partial50.h5ad.meta.json"
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
export VCC_DE_SIGN_PRIOR=clean4       # THE ARM (390 C5b / 391): CPM-block signs from the 4-bank cross-line prior, DE channel only
unset VCC_NORMREL_PM VCC_RC_ORDER VCC_CTX_A_PRIOR VCC_CTX_A_PRIOR_BETA VCC_SIGN_REPAIR
log "ARMED pid=$$  VCC_NORMREL_EXP=$VCC_NORMREL_EXP  VCC_NORMREL_MODE=$VCC_NORMREL_MODE  VCC_CD4_PARTIAL_G=$VCC_CD4_PARTIAL_G  VCC_DE_SIGN_PRIOR=$VCC_DE_SIGN_PRIOR  VCC_SIGN_REPAIR/VCC_NORMREL_PM/VCC_RC_ORDER/VCC_CTX_A_PRIOR=<unset>"
log "builder: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)  driver: $(shasum -a 256 scripts/def_pipeline.py | cut -d' ' -f1)  git: $(git rev-parse --short HEAD)"
B="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"; D="$(shasum -a 256 scripts/def_pipeline.py | cut -d' ' -f1)"
# builder: 319p2 ran 1318cd3c...; since then c47d2343 (VCC_CTX_A_PRIOR), 123df616 (VCC_RC_ORDER) and the 392 knob
# VCC_DE_SIGN_PRIOR (inert unset; pins src/vcc2026/rc_emit_k.py a12cf8f9... and the clean4 prior 7e532f84...).
[ "$B" = "260e03aaef3984d81f4b074ca528590d60b4d0bd8c02e06eb799cc8b68321bae" ] || { log "BUILDER hash $B != 260e03aa... (VCC_DE_SIGN_PRIOR knob, inert unset) -- stopping"; exit 4; }
[ "$D" = "27edc486232afff2e19c46bdd87c8c84322cdc4b7781f55c6e5d40adcd8c180f" ] || { log "DRIVER hash $D != 27edc486... -- stopping"; exit 4; }
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
PB="$W/logs/postbuild.log"
NFAIL="$(grep -c '^  FAIL' "$PB" 2>/dev/null || echo 99)"
log "postbuild: $(grep -hE 'checks pass' "$PB" 2>/dev/null | tail -1)  FAIL lines: $NFAIL"
grep -h '^  FAIL' "$PB" 2>/dev/null | tee -a "$LOG"
if [ ! -s "$W/$NAME.h5ad" ]; then log "NO ARTIFACT -- stopping"; exit 1; fi
if [ "$NFAIL" != "1" ]; then log "UNEXPECTED postbuild failure count ($NFAIL != 1) -- stopping before prep"; exit 1; fi
if ! grep -q '^  FAIL  config differences' "$PB"; then log "the single FAIL is NOT the config-drift line -- stopping"; exit 1; fi
META="$W/$NAME.h5ad.meta.json"
log "meta de_sign_prior=$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("de_sign_prior","<ABSENT -- ARM DID NOT LAND>"))' "$META")"
log "meta rc_callset prior diagnostics: $(uv run python -c 'import json,sys; m=json.load(open(sys.argv[1])); s=m["output"].get("rc_callset",{}); print({k: s[k] for k in s if k.startswith(("prior_","dose_","steer_","n_classes"))})' "$META")"
log "artifact: $(ls -la "$W/$NAME.h5ad" | awk '{print $5, $9}')  sha256 $(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' "$META")"

# the config diff, BEFORE the prep: exactly the five de_sign_* keys appear, everything else equals the 319p2 package's meta.
uv run python - "$META" "$REF_META" <<'PY' | tee -a "$LOG"
import json,sys
built=json.load(open(sys.argv[1]))["config"]; ref=json.load(open(sys.argv[2]))["config"]
d={k:(ref.get(k,"<absent>"),built.get(k,"<absent>")) for k in sorted(set(ref)|set(built)) if ref.get(k,"<absent>")!=built.get(k,"<absent>")}
for k,(a,b) in d.items(): print(f"[392]     {k}: {a!r} -> {b!r}")
want={"de_sign_prior","de_sign_prior_sha256","de_sign_emitter_sha256","de_sign_classes","de_sign_shrink_floor"}
ok = set(d)==want and all(a=="<absent>" for a,_ in d.values()) and d["de_sign_prior"][1]=="clean4"
print("[392] PASS: exactly the five de_sign_* keys, de_sign_prior -> 'clean4'" if ok else f"[392] FAIL: got {sorted(d)}")
sys.exit(0 if ok else 1)
PY
RC=${PIPESTATUS[0]}; [ $RC -ne 0 ] && { log "CONFIG DIFF is not exactly the five de_sign_* keys -- stopping before prep"; exit 1; }

wait_gate
log "INVOCATION 2: vcc prep -> $PKG"
/usr/bin/time -l vcc prep "$W/$NAME.h5ad" -g "$BUNDLE/gene_names.csv" \
    --perts "$BUNDLE/pert_counts.csv" --contexts A,B,C \
    -o "$PKG" > "$W/prep.log" 2>&1
RC=$?; log "prep rc=$RC  disk $(free_gib) GiB"
tail -20 "$W/prep.log" | tee -a "$LOG"
log "vcc: $(ls -la "$PKG" 2>/dev/null | awk '{print $5, $9}')"
cp "$META" "$REPO/submissions/392-cap15-partial50-designclean4.h5ad.meta.json"
log "meta copied beside the package: $(shasum -a 256 "$REPO/submissions/392-cap15-partial50-designclean4.h5ad.meta.json" | cut -d' ' -f1)"
uv run python scripts/check_submission_novelty.py check "$PKG" 2>&1 | tee -a "$LOG"
log "novelty rc=${PIPESTATUS[0]}"
log "package sha256: $(shasum -a 256 "$PKG" | cut -d' ' -f1)"
log "git status scripts/: $(git status --short scripts/ | tr '\n' ' ')"
log "DONE"
