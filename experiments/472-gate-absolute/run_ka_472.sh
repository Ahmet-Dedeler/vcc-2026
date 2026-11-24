#!/usr/bin/env bash
# 472 known-answer build (copied from 468 <- 451): the shipped base + the caller's VCC_* exports, build only (no prep), then the built artifact's
# output sha256 must equal REF_META's. usage: run_ka_451.sh TAG REF_META BUILDER_SHA
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
TAG="${1:?tag}"; REF_META="${2:?ref meta}"; WANT_B="${3:?builder sha}"
W="$HERE/work/build_$TAG"; mkdir -p "$W"; LOG="$W/run$TAG.log"; NAME="ka472-$TAG"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
restore(){ [ -f "$W/pin-artifacts.json" ] && uv run python scripts/def_pin_panel_artifacts.py --restore --record "$W/pin-artifacts.json" | tee -a "$LOG"; }
export VCC_NORMREL_EXP=-1.5 VCC_CD4_PARTIAL_G=0.5 VCC_NORMREL_MODE=cap VCC_DE_SIGN_PRIOR="${ARM_DE_SIGN_PRIOR:-clean4q80}" \
       VCC_PERT_SCALE=wshare_lo0.3_hi3 VCC_PERT_SCALE_BALLAST=1 VCC_BULK_SHARPEN_Q="${KA_Q:-1.75}"
unset VCC_NORMREL_PM VCC_RC_ORDER VCC_CTX_A_PRIOR VCC_CTX_A_PRIOR_BETA VCC_SIGN_REPAIR VCC_BULK_GENE_GAIN VCC_CARRIER_GAMMA VCC_H1_UNION_G VCC_CD4_TARGET_F
# 472 (2026-09-26 21:05Z): wait for memory before the driver (whose preflight refuses below 23 GiB free+inactive); other
# sessions' analysis workers may be loading data. free% >= 40 and free+inactive >= 25 GiB, checked every 60 s for up to 3 h.
memok(){ local pct fi; pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  log "memory gate: free ${pct}%  free+inactive ${fi} GiB"; [ "${pct:-0}" -ge 40 ] && awk -v v="$fi" 'BEGIN{exit !(v>=25)}'; }
n=0; until memok; do n=$((n+1)); [ $n -ge 180 ] && { log "memory gate never opened in 3 h -- stopping"; exit 3; }; sleep 60; done
B="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
log "KA $TAG: builder ${B:0:16} (want ${WANT_B:0:16})  VCC_SAMELINE=${VCC_SAMELINE:-<unset>}  prior=$VCC_DE_SIGN_PRIOR"
[ "$B" = "$WANT_B" ] || { log "builder hash mismatch -- stopping"; exit 4; }
/usr/bin/time -l uv run python scripts/def_pipeline.py --bundle "$REPO/data/challenge" --work "$W" --name "$NAME" \
    --whitener-id union_pruned_def --emitter rc_callset --frozen-inputs "$W/frozen_inputs.json" --min-free-gib 23 \
    --stages preflight,identity,cd4-extract,census,whitener,repin-frozen,pin-artifacts,repin-panel-support,build \
    > "$W/driver.log" 2>&1
log "driver rc=$?"; restore
log "builder after restore: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -c1-16)"
META="$W/$NAME.h5ad.meta.json"
[ -s "$META" ] || { log "no meta -- build failed: $(tail -3 $W/driver.log)"; exit 1; }
uv run python - "$META" "$REF_META" "$TAG" <<'PY' | tee -a "$LOG"
import json, sys
a = json.load(open(sys.argv[1])); b = json.load(open(sys.argv[2]))
sa, sb = a["output"]["sha256"], b["output"]["sha256"]
print(f"[{sys.argv[3]}] output sha256 {sa[:16]} vs reference {sb[:16]}: " + ("KA PASS (identical)" if sa == sb else "KA FAIL"))
print(f"[{sys.argv[3]}] sameline meta: {a['config'].get('sameline', '<absent>')}")
PY
rm -f "$W/$NAME.h5ad"
log "END"
