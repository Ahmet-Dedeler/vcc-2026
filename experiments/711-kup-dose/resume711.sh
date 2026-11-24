#!/usr/bin/env bash
# 711 ku38 resume (science [d80a10], 2026-10-16; coordinator [2d6f46] asked): finish ku38 from PREP on the existing h5ad
# (no rebuild), under experiments/711-kup-dose/OVERRIDE-554.md (554's call-set-missed floor, this information arm only).
# Everything after postbuild is run_arm_generic.sh's own code, extracted from it verbatim (gate_ok, the config diff, the
# zero-gene audit, the cliff check, vcc prep, the meta copy, novelty, the package sha), appending to the same run log
# (build_ku38/runku38.log, which armready 1018s2 reads). Two changes only:
#   - run_arm_generic's "postbuild FAIL count must be 1" becomes "the FAIL lines are exactly the known two" (the config-drift
#     line and 554's "pooled: median call_set_missed 141.5 past its floor 100"); any other FAIL still stops;
#   - its 240-iteration prep wait becomes a timestamp cap (AGENTS 6.2): 10-17 23:00Z.
# Then build711's own READY checks (PASS exactly, CLIFF OK, zero-gene audit 0, novelty rc=0, package sha) and its READY line
# in experiments/701-kup-callset/work/READY-science.log. The builder / emitter are checked as the BUILT artifact's recorded
# provenance (run log + meta), not as the live files: the build is over, and a later build's pin may hold the builder file.
# Environment as build711 exported for the build where the post-build code reads it (VCC_BULK_CTX_SCALE sets the cliff band).
#   nohup bash scripts/detach.sh bash experiments/711-kup-dose/resume711.sh > experiments/711-kup-dose/work/resume711.nohup 2>&1 &
# Terminal line in its own log (work/resume711.log), always: "<ts> READY ku38 <sha>" | "<ts> NOT READY <why>", then "<ts> END".
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
HERE="$REPO/experiments/428-zero-gene-audit"
ARM=ku38; NAME=kup701; QV=1.75
STEM="446bj2zku38-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cpm3800"
REF_META="$REPO/submissions/446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250.h5ad.meta.json"
EXPECT="rc_callset_cpm,de_sign_emitter_sha256"
W="$HERE/work/build_$ARM"; LOG="$W/run$ARM.log"; BUNDLE="$REPO/data/challenge"; PKG="$REPO/submissions/$STEM.vcc"
RLOG="$REPO/experiments/711-kup-dose/work/resume711.log"; S=submissions
BUS=e99bbd4a63ad834c5ac7a042d890513f7c1e098f09df4173ab6843021c296cc5; EMS=0f56209d57e8385abd08a7a15692ae3241b9373b436b4921546e546d571e5876
export VCC_BULK_CTX_SCALE="A=1.0962,B=1.1461,C=1.1681"
CAPTS=202610172300
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
rlog(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$RLOG"; }
sha(){ shasum -a 256 "$1" | cut -d' ' -f1; }
TERM_=0; trap 'rc=$?; [ "$TERM_" = 1 ] || rlog "NOT READY unexpected exit rc $rc"; rlog "END"' EXIT
nr(){ TERM_=1; rlog "NOT READY $1"; log "RESUME stopped: $1"; exit 1; }
# ---- run_arm_generic.sh's gate (verbatim)
free_gib(){ df -k "$REPO" | awk 'NR==2{printf "%.1f", $4/1048576}'; }
gate_ok(){
  # $1 = build | prep. Both need no other builder, no vcc prep, disk >= 20 GiB, and one memory branch:
  #   fi23        memory_pressure free >= 40 % AND vm_stat free+inactive >= 23 GiB (the original gate; both modes)
  #   pressure55  memory_pressure free >= 55 % AND swap free >= 400 MiB (swap total 0 = no swap file = passes, 670 10-16) (591; PREP ONLY: other apps' pages held in the
  #               compressor keep free+inactive under 23 GiB while the machine has room; 549 was prepped by hand 09-30)
  # The build keeps fi23 alone: def_pipeline.py re-checks free+inactive(+speculative+purgeable) >= --min-free-gib before
  # its build stages and REFUSES instead of waiting, so a looser gate there would turn a wait into a failed build.
  local mode="${1:-build}" pct fi swapfree others preps disk branch
  pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  swapfree="$(sysctl -n vm.swapusage 2>/dev/null | awk '{for(i=1;i<NF;i++) if($i=="total" && $(i+2)+0==0){print 999999; exit}; for(i=1;i<NF;i++) if($i=="free"){v=$(i+2); u=substr(v,length(v),1); n=substr(v,1,length(v)-1)+0; if(u=="G") n=n*1024; else if(u=="K") n=n/1024; printf "%.0f", n; exit}}')"
  others="$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null || echo 99)"
  preps="$(ps -axo args= | grep -vE 'claude|grep' | grep -cE '(^|/| )vcc prep( |$)' | tr -d ' ')"
  disk="$(free_gib)"
  branch=none
  if [ "${pct:-0}" -ge 40 ] && awk -v v="$fi" 'BEGIN{exit !(v>=23)}'; then branch=fi23
  elif [ "$mode" = prep ] && [ "${pct:-0}" -ge 55 ] && awk -v v="${swapfree:-0}" 'BEGIN{exit !(v>=400)}'; then branch=pressure55
  fi
  log "gate: free ${pct}%  free+inactive ${fi} GiB  other builders ${others}  vcc preps ${preps}  disk ${disk} GiB  swap free ${swapfree:-?} MiB  mode ${mode}  memory branch ${branch}"
  [ "$branch" != none ] && [ "$others" = "0" ] && [ "$preps" = "0" ] && awk -v v="$disk" 'BEGIN{exit !(v>=20)}'
}
wait_gate_ts(){ until gate_ok "${1:-build}"; do [ "$(date -u +%Y%m%d%H%M)" -ge "$CAPTS" ] && { log "GATE never opened by $CAPTS ($1)"; nr "prep gate never opened by $CAPTS"; }; sleep 30; done; }
rlog "armed (pid $$): resume ku38 from prep under OVERRIDE-554.md; cap $CAPTS"
[ -e experiments/691-release-watcher/work/RELEASE-DETECTED ] && nr "release detected"
[ -s "$W/$NAME.h5ad" ] || nr "no artifact $W/$NAME.h5ad"
[ -e "$PKG" ] && nr "a package already exists at $PKG"
# provenance of the BUILT artifact (the build is over; a later build's pin may hold the live builder file, as R-NEW's does
# while it builds, and vcc prep / novelty read only the h5ad): the run log's ARMED-time builder and its post-restore builder
# are the build's (e99bbd4a), and the meta's emitter pin is the build's (0f56209d)
grep -q "^.* builder: $BUS  driver: " "$LOG" || nr "the run log's ARMED builder is not $BUS"
grep -q "builder after restore: $BUS$" "$LOG" || nr "the run log's post-restore builder is not $BUS"
[ "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"]["de_sign_emitter_sha256"])' "$W/$NAME.h5ad.meta.json")" = "$EMS" ] \
  || nr "the meta's emitter pin is not $EMS"
log "RESUME (experiments/711-kup-dose/resume711.sh, OVERRIDE-554.md): from prep on the existing artifact; no rebuild"
PB="$W/logs/postbuild.log"
NFAIL="$(grep -c '^  FAIL' "$PB" 2>/dev/null || echo 99)"
grep -q '^  FAIL  config differences' "$PB" || nr "the config-drift FAIL line is missing"
grep -qE '^  FAIL  emitter steering dose per context: .* \(554\): pooled: median call_set_missed 141\.5 past its floor 100$' "$PB" || nr "554's FAIL line is not the overridden one"
[ "$NFAIL" = "2" ] || nr "postbuild FAIL count $NFAIL (the override covers exactly the two known lines)"
log "postbuild FAIL lines: exactly the known two (config drift; 554 call_set_missed 141.5 > 100, overridden for this arm only)"
# ---- run_arm_generic.sh after postbuild (verbatim; the prep wait is timestamp-capped)
META="$W/$NAME.h5ad.meta.json"
log "meta pert_scale=$(uv run python -c 'import json,sys; m=json.load(open(sys.argv[1])); c=m["config"]; o=m["output"].get("rc_callset",{}); print(c.get("pert_scale_ballast","<ABSENT -- ARM DID NOT LAND>"), {k: o[k] for k in o if k.startswith("ballast")})' "$META")"
log "meta rc_callset prior diagnostics: $(uv run python -c 'import json,sys; m=json.load(open(sys.argv[1])); s=m["output"].get("rc_callset",{}); print({k: s[k] for k in s if k.startswith(("prior_","dose_","steer_","n_classes","sign_flip"))})' "$META")"
log "artifact: $(ls -la "$W/$NAME.h5ad" | awk '{print $5, $9}')  sha256 $(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' "$META")"

# the config diff, BEFORE the prep: exactly de_mde (+ the emitter pin) vs the reference arm's meta.
uv run python - "$META" "$REF_META" "$ARM" "$EXPECT" <<'PY' | tee -a "$LOG"
import json,sys
built=json.load(open(sys.argv[1]))["config"]; ref=json.load(open(sys.argv[2]))["config"]
d={k:(ref.get(k,"<absent>"),built.get(k,"<absent>")) for k in sorted(set(ref)|set(built)) if ref.get(k,"<absent>")!=built.get(k,"<absent>")}
for k,(a,b) in d.items(): print(f"[{sys.argv[3]}]     {k}: {a!r} -> {b!r}")
want=set(x for x in sys.argv[4].split(",") if x)
ok = set(d)==want
print(f"[{sys.argv[3]}] PASS: exactly {sorted(want)}" if ok else f"[{sys.argv[3]}] FAIL: got {sorted(d)} want {sorted(want)}")
sys.exit(0 if ok else 1)
PY
RC=${PIPESTATUS[0]}; [ $RC -ne 0 ] && { log "CONFIG DIFF is not exactly $EXPECT -- stopping before prep"; exit 1; }
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
wait_gate_ts prep
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
# ---- build711's READY checks (verbatim logic)
RL="$LOG"
# build711's "builder or emitter moved during the build" guarded a build in progress; this resume does not build (see above)
/usr/bin/grep -q "^\[$ARM\] PASS: exactly" "$RL" 2>/dev/null || nr "config diff is not exactly $EXPECT"
/usr/bin/grep -q "^\[$ARM\] CLIFF OK" "$RL" 2>/dev/null || nr "cliff check"
ZL="$(/usr/bin/grep -E 'zero-gene audit' "$RL" | tail -1)"
echo "$ZL" | /usr/bin/grep -qE 'ZERO block total \(target excluded\): 0 in 0 blocks' || nr "zero-gene audit: ${ZL:-missing}"
NOV="$(/usr/bin/grep -E 'novelty rc=' "$RL" | tail -1 | /usr/bin/grep -oE 'rc=[0-9]+')"
[ "$NOV" = "rc=0" ] || nr "novelty ${NOV:-missing}"
PK="$(/usr/bin/grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"
[ -s "$S/$STEM.vcc" ] && [ -n "$PK" ] && [ "$(sha "$S/$STEM.vcc")" = "$PK" ] || nr "no package or sha differs"
printf '%s READY ku38 %s %s\n' "$(date -u +%FT%TZ)" "$PK" "$S/$STEM.vcc" >> experiments/701-kup-callset/work/READY-science.log
printf '%s READY ku38 %s (resume711 under OVERRIDE-554)\n' "$(date -u +%FT%TZ)" "$PK" >> experiments/711-kup-dose/work/build701.log
TERM_=1; rlog "READY ku38 $PK"
