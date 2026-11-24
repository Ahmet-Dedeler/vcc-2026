#!/usr/bin/env bash
# 319 -- BUILD the SUPPRESSION-ONLY normrel arm (`VCC_NORMREL_MODE=cap`) on the REAL A/B/C panel.
# Modelled on run316.sh, which is run312.sh's gated two-invocation pattern.  ONE variable against
# entry Q519I2pLWRCW7PAEx5ul (`312-normrel15`): the same exponent -1.5, `VCC_NORMREL_PM` UNSET, only
# the ratio clamp added.  Nothing is submitted here and `vcc submit` is never run.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
W="$HERE/work/build_p"; mkdir -p "$W"; LOG="$W/run319p.log"
BUNDLE="$REPO/data/challenge"
NAME="the-stapler-has-tenure"
PKG="$REPO/submissions/319p-cap15-partial.vcc"
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
  # the build writes 13.2 GiB (h5ad 11.8 + whitener 1.2 + cd4_compact 0.2) and the prep 2.6 GiB
  # beside it, measured on 316; 20 GiB free is that with ~4 GiB of slack.
  [ "${pct:-0}" -ge 40 ] && awk -v v="$fi" 'BEGIN{exit !(v>=23)}' && [ "$others" = "0" ] \
    && [ "$preps" = "0" ] && awk -v v="$disk" 'BEGIN{exit !(v>=20)}'
}
wait_gate(){ local n=0; until gate_ok; do n=$((n+1)); [ $n -ge 240 ] && { log "GATE never opened"; exit 3; }; sleep 30; done; }
restore(){ [ -f "$W/pin-artifacts.json" ] && uv run python scripts/def_pin_panel_artifacts.py --restore --record "$W/pin-artifacts.json" | tee -a "$LOG"; }

export VCC_NORMREL_EXP=-1.5          # 319's dose
export VCC_CD4_PARTIAL_G=0.34        # THE ARM: 254's partial-condition CD4 rows as their own residual (PREREGISTRATION-partial.md)
export VCC_NORMREL_MODE=cap          # THE ARM: ratio clamped at 1, so only OVERSIZED rows are suppressed
unset VCC_NORMREL_PM                 # the shipped weight-first rendering, i.e. 312's, not 316's
log "ARMED pid=$$  VCC_NORMREL_EXP=$VCC_NORMREL_EXP  VCC_NORMREL_MODE=$VCC_NORMREL_MODE  VCC_CD4_PARTIAL_G=$VCC_CD4_PARTIAL_G  VCC_NORMREL_PM=<unset>"
log "builder: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)  driver: $(shasum -a 256 scripts/def_pipeline.py | cut -d' ' -f1)  git: $(git rev-parse --short HEAD)"
B="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"; D="$(shasum -a 256 scripts/def_pipeline.py | cut -d' ' -f1)"
[ "$B" = "1318cd3c7dbe8ed129a4f45dbf7e1e4ce7bed75c5f73af7174d3dd3530f118f4" ] || { log "BUILDER hash $B != 319d3's 1318cd3c... -- not one variable away, stopping"; exit 4; }
# driver: 319d3 ran 3fb0c4bc... (commit 04950e4b); current 27edc486... (commit 230a7fe5) differs ONLY in one
# comment line where commit-hash citations were remapped after the 136MB blob strip (git diff 04950e4b 230a7fe5
# -- scripts/def_pipeline.py: zero non-comment lines). Logic identical; gate on the current hash.
[ "$D" = "27edc486232afff2e19c46bdd87c8c84322cdc4b7781f55c6e5d40adcd8c180f" ] || { log "DRIVER hash $D != 27edc486... (230a7fe5, comment-only delta from 319d3's 3fb0c4bc) -- not one variable away, stopping"; exit 4; }
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
log "meta normrel_exp=$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("normrel_exp","<ABSENT -- ARM DID NOT LAND>"))' "$META")"
log "meta normrel_mode=$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("normrel_mode","<ABSENT -- ARM DID NOT LAND>"))' "$META")"
log "meta normrel_pm=$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("normrel_pm","<absent, as intended>"))' "$META")"
log "artifact: $(ls -la "$W/$NAME.h5ad" | awk '{print $5, $9}')  sha256 $(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' "$META")"

# the one-key config diff, BEFORE the prep: exactly `cd4_partial_g` <absent> -> 0.34 against the 319 package's own meta.
uv run python - "$META" "$REPO/submissions/319-normrel15-cap.h5ad.meta.json" <<'PY' | tee -a "$LOG"
import json,sys
built=json.load(open(sys.argv[1]))["config"]; ref=json.load(open(sys.argv[2]))["config"]
d={k:(ref.get(k,"<absent>"),built.get(k,"<absent>")) for k in sorted(set(ref)|set(built)) if ref.get(k,"<absent>")!=built.get(k,"<absent>")}
for k,(a,b) in d.items(): print(f"[319p]     {k}: {a!r} -> {b!r}")
ok = set(d)=={"cd4_partial_g"} and d["cd4_partial_g"]==("<absent>",0.34)
print("[319p] PASS: exactly one key, cd4_partial_g <absent> -> 0.34" if ok else f"[319p] FAIL: got {sorted(d)}")
sys.exit(0 if ok else 1)
PY
RC=${PIPESTATUS[0]}; [ $RC -ne 0 ] && { log "CONFIG DIFF is not exactly cd4_partial_g <absent> -> 0.34 -- stopping before prep"; exit 1; }

wait_gate
log "INVOCATION 2: vcc prep -> $PKG  (the driver's prep stage is --dry-run, i.e. validation only,"
log "  so the package is written by this real prep; ~22-24 GiB RSS, gated above)"
/usr/bin/time -l vcc prep "$W/$NAME.h5ad" -g "$BUNDLE/gene_names.csv" \
    --perts "$BUNDLE/pert_counts.csv" --contexts A,B,C \
    -o "$PKG" > "$W/prep.log" 2>&1
RC=$?; log "prep rc=$RC  disk $(free_gib) GiB"
tail -20 "$W/prep.log" | tee -a "$LOG"
log "vcc: $(ls -la "$PKG" 2>/dev/null | awk '{print $5, $9}')"
cp "$META" "$REPO/submissions/319p-cap15-partial.h5ad.meta.json"
log "meta copied beside the package: $(shasum -a 256 "$REPO/submissions/319p-cap15-partial.h5ad.meta.json" | cut -d' ' -f1)"
uv run python scripts/check_submission_novelty.py check "$PKG" 2>&1 | tee -a "$LOG"
log "novelty rc=${PIPESTATUS[0]}"
log "git status scripts/: $(git status --short scripts/ | tr '\n' ' ')"
log "DONE"
