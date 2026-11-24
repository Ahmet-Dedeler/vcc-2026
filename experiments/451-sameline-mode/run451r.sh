#!/usr/bin/env bash
# 451r — mechanics rehearsal of VCC_SAMELINE on ef-bundle-r with a FORCED identity (D = K562 -> k562_gwps 300/300 covered,
# E = HCT116 -> orion_hct116 300/300, F = A549 -> jiang_A549/myllia rows only, programme_ok False). No truth: this checks the path end to end
# (pins, blend, per-target prior, audits) and measures how far the blend moves each context's spread (rho) at 448's lvl2
# factors, i.e. whether the day-of re-level must run after the same-line blend. Runs after chain451 (patched builder).
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/451-sameline-mode; W=experiments/448-def-rehearsal-q175/work/efr; LOG=$X/work/run451r.log
B=experiments/238-essential-free-panel/work/ef-bundle-r
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
until grep -q " END builder\| RESTORED\|stopping" $X/work/chain451.log 2>/dev/null && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain451\.sh/' | grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609252100 ] && { log "chain451 never ended by 21:00Z"; exit 1; }; sleep 30; done
grep -q " END builder" $X/work/chain451.log || { log "chain451 did not end with the patched builder -- stopping"; exit 1; }
for v in $(env | grep -o '^VCC_[A-Z0-9_]*'); do unset "$v"; done
BSHA0="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"; log "start: builder ${BSHA0:0:16}"
DRV=(uv run python scripts/def_pipeline.py --bundle "$B" --work "$W" --whitener-id union_pruned_def --emitter rc_callset
     --frozen-inputs "$W/frozen_inputs.json" --min-free-gib 23)
restore(){ [ -f "$W/pin-artifacts.json" ] && uv run python scripts/def_pin_panel_artifacts.py --restore --record "$W/pin-artifacts.json" >> "$LOG" 2>&1
           log "builder after restore: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -c1-16) (want ${BSHA0:0:16})"; }
trap restore EXIT
"${DRV[@]}" --name 451r-setup --stages preflight,identity,cd4-extract,census,whitener,key,repin-frozen,pin-artifacts,repin-panel-support > "$X/work/driver-setup-451r.log" 2>&1
grep -q "^--- \[repin-panel-support\].*rc=0" "$X/work/driver-setup-451r.log" || { log "setup did not complete -- stopping"; exit 1; }
SL=$X/work/sameline_efr_forced.npz; SS=$(shasum -a 256 $SL | cut -d' ' -f1)
DUMP="$W/dump_h1x2_cd4x1.npz"; DSHA="$(shasum -a 256 "$DUMP" | cut -d' ' -f1)"
export VCC_BULK_LOAD="$DUMP:$DSHA" VCC_BULK_SHARPEN_Q=1.75 VCC_BULK_CTX_SCALE="D=1.2218,E=1.2721,F=1.2356" VCC_SAMELINE="$SL:$SS:1:pt"
log "build: VCC_SAMELINE=$(basename $SL):${SS:0:16}:1:pt at 448 lvl2 factors"
( "${DRV[@]}" --name 451r-sl --stages build > "$X/work/driver-451r.log" 2>&1 )
A="$W/451r-sl.h5ad"
[ -s "$A" ] || { log "build failed: $(tail -3 $X/work/driver-451r.log)"; exit 1; }
grep -h "\[sameline\]" "$X/work/driver-451r.log" "$W/logs/build.log" 2>/dev/null | sort -u | tee -a "$LOG"
LV="$(uv run python experiments/417-bulk-sharpen/ntclevel417.py "$A" --bundle "$B" --contexts D,E,F --sd 0.01 2>&1)"; echo "$LV" >> "$LOG"
log "rho with the same-line blend at lvl2 factors: $(echo "$LV" | grep -E '^  [DEF]:' | awk '{print $1, $13}' | tr '\n' ' ') (without: D 0.9798 E 0.9798 F 0.9463)"
uv run python experiments/428-zero-gene-audit/audit428.py "$A" --bundle "$B" 2>&1 | head -1 | tee -a "$LOG"
uv run python - "$A.meta.json" <<'PY' | tee -a "$LOG"
import json, sys
m = json.load(open(sys.argv[1])); print("meta sameline:", json.dumps(m["config"].get("sameline"))[:600])
PY
rm -f "$A"
log "END"
