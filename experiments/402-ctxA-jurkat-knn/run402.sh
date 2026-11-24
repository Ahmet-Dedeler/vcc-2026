#!/usr/bin/env bash
# 399 -- BUILD the per-target neighbour-conditioned tilt arm (VCC_DE_SIGN_PRIOR=knn25o400q80, 397/398) on the 395 composition,
# REAL A/B/C panel. Modelled on run392.sh. One knob value against entry TJTOj9qrbiobJBpqLSa3 (392): the config diff vs
# 392's meta must be exactly de_sign_protect_q (absent -> 0.8) and de_sign_emitter_sha256 (a12cf8f9 -> 6ff0b866; the
# emitter's clean4 path is byte-identical, the diff is the protect_q lines). Nothing is submitted here.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
W="$HERE/work/build"; mkdir -p "$W"; LOG="$W/run402.log"
BUNDLE="$REPO/data/challenge"
NAME="the-paperclip-demanded-a-raise"
PKG="$REPO/submissions/402-cap15-partial50-knn25o400jaq80.vcc"
REF_META="$REPO/submissions/399-cap15-partial50-knn25o400q80.h5ad.meta.json"
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
export VCC_DE_SIGN_PRIOR=knn25o400jaq80    # THE ARM (399): 395 with the per-target neighbour-conditioned tilt, six-bank family, shrink 400 (397/398/398c/398d) instead of the marginal CLEAN4 tilt
unset VCC_NORMREL_PM VCC_RC_ORDER VCC_CTX_A_PRIOR VCC_CTX_A_PRIOR_BETA VCC_SIGN_REPAIR
log "ARMED pid=$$  VCC_NORMREL_EXP=$VCC_NORMREL_EXP  VCC_NORMREL_MODE=$VCC_NORMREL_MODE  VCC_CD4_PARTIAL_G=$VCC_CD4_PARTIAL_G  VCC_DE_SIGN_PRIOR=$VCC_DE_SIGN_PRIOR  VCC_SIGN_REPAIR/VCC_NORMREL_PM/VCC_RC_ORDER/VCC_CTX_A_PRIOR=<unset>"
log "builder: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)  driver: $(shasum -a 256 scripts/def_pipeline.py | cut -d' ' -f1)  emitter: $(shasum -a 256 src/vcc2026/rc_emit_k.py | cut -d' ' -f1)  git: $(git rev-parse --short HEAD)"
B="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"; D="$(shasum -a 256 scripts/def_pipeline.py | cut -d' ' -f1)"; EM="$(shasum -a 256 src/vcc2026/rc_emit_k.py | cut -d' ' -f1)"
[ "$B" = "7d3f2b4536015c44b7f1de5aeb72697894ce6ba88a25411abd0d3f8a0fe8676a" ] || { log "BUILDER hash $B != 7d3f2b4536015c44b7f1de5aeb72697894ce6ba88a25411abd0d3f8a0fe8676a (per-target prior knob) -- stopping"; exit 4; }
[ "$D" = "27edc486232afff2e19c46bdd87c8c84322cdc4b7781f55c6e5d40adcd8c180f" ] || { log "DRIVER hash $D != 27edc486... -- stopping"; exit 4; }
[ "$EM" = "dbbb91210a2d808fd03ae291ae3dc1479a30e58c675c5dc9536da6c266e5c09c" ] || { log "EMITTER hash $EM != dbbb9121... (399 guard) -- stopping"; exit 4; }
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
log "meta de_sign_prior=$(uv run python -c 'import json,sys; c=json.load(open(sys.argv[1]))["config"]; print(c.get("de_sign_prior","<ABSENT -- ARM DID NOT LAND>"), c.get("de_sign_protect_q","<protect_q ABSENT>"), c.get("de_sign_prior_targets_with_row","<per-target rows ABSENT>"))' "$META")"
log "meta rc_callset prior diagnostics: $(uv run python -c 'import json,sys; m=json.load(open(sys.argv[1])); s=m["output"].get("rc_callset",{}); print({k: s[k] for k in s if k.startswith(("prior_","dose_","steer_","n_classes","sign_flip"))})' "$META")"
log "artifact: $(ls -la "$W/$NAME.h5ad" | awk '{print $5, $9}')  sha256 $(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' "$META")"

# the config diff, BEFORE the prep: exactly de_sign_prior (knn25o400q80 -> knn25o400jaq80), de_sign_prior_sha256 (399 -> 402
# artifact) and de_sign_prior_context_A (absent -> nadig_jurkat_knn25_shrink50); everything else equals 399's meta.
uv run python - "$META" "$REF_META" <<'PY' | tee -a "$LOG"
import json,sys
built=json.load(open(sys.argv[1]))["config"]; ref=json.load(open(sys.argv[2]))["config"]
d={k:(ref.get(k,"<absent>"),built.get(k,"<absent>")) for k in sorted(set(ref)|set(built)) if ref.get(k,"<absent>")!=built.get(k,"<absent>")}
for k,(a,b) in d.items(): print(f"[402]     {k}: {a!r} -> {b!r}")
want={"de_sign_prior","de_sign_prior_sha256","de_sign_prior_context_A"}
sha402=json.load(open("experiments/402-ctxA-jurkat-knn/work/knn25o400_A_jurkat50_prior_challenge.json"))["artifact_sha256"]
sha399=json.load(open("experiments/399-knn-tilt-arm/work/knn25o400_prior_challenge.json"))["artifact_sha256"]
ok = (set(d)==want and d["de_sign_prior"]==("knn25o400q80","knn25o400jaq80") and d["de_sign_prior_sha256"]==(sha399,sha402)
      and d["de_sign_prior_context_A"]==("<absent>","nadig_jurkat_knn25_shrink50"))
print("[402] PASS: exactly the knob value, the artifact sha and the context-A tag" if ok else f"[402] FAIL: got {d}")
sys.exit(0 if ok else 1)
PY
RC=${PIPESTATUS[0]}; [ $RC -ne 0 ] && { log "CONFIG DIFF is not exactly the prior identity keys -- stopping before prep"; exit 1; }

wait_gate
log "INVOCATION 2: vcc prep -> $PKG"
/usr/bin/time -l vcc prep "$W/$NAME.h5ad" -g "$BUNDLE/gene_names.csv" \
    --perts "$BUNDLE/pert_counts.csv" --contexts A,B,C \
    -o "$PKG" > "$W/prep.log" 2>&1
RC=$?; log "prep rc=$RC  disk $(free_gib) GiB"
tail -20 "$W/prep.log" | tee -a "$LOG"
log "vcc: $(ls -la "$PKG" 2>/dev/null | awk '{print $5, $9}')"
cp "$META" "$REPO/submissions/402-cap15-partial50-knn25o400jaq80.h5ad.meta.json"
log "meta copied beside the package: $(shasum -a 256 "$REPO/submissions/402-cap15-partial50-knn25o400jaq80.h5ad.meta.json" | cut -d' ' -f1)"
uv run python scripts/check_submission_novelty.py check "$PKG" 2>&1 | tee -a "$LOG"
log "novelty rc=${PIPESTATUS[0]}"
log "package sha256: $(shasum -a 256 "$PKG" | cut -d' ' -f1)"
log "git status scripts/ src/: $(git status --short scripts/ src/ | tr '\n' ' ')"
log "DONE"
