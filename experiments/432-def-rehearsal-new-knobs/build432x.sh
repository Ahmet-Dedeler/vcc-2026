#!/usr/bin/env bash
# 432 EXTRAS (build432x.sh, a copy of build432.sh made 2026-09-23 ~12:20Z while the main chain was executing build432.sh;
# never edit a script a live bash is reading). Adds: dump_b0env (KA for LOAD: b0's own-env dump), and on SET=ef
# phase b1 also exports VCC_DE_MAG=dn_k250_zs (one low-coverage build exercising every knob but the re-level).
# 432: D/E/F dry rehearsal of the September knobs on 199's realistic synthetic bundle (def-bundle-r).
# usage: build432.sh <phase>      phases, in order:
#   whitener   the driver's operator re-fit for the synthetic panel (scratch removed afterwards)
#   key        410's key rebuilt for the panel on that operator (key432.py) -> work/key/pert_scale_def432_lo0.3_hi3.npz
#   dump_stack | dump_lit | dump_h1       VCC_BULK_DUMP: the D/E/F stack (H1 x2 + CD4 x1.5 as VCC_CD4_TARGET_G x1.5),
#              417's literal recipe (VCC_CD4_TARGET_F=0.51), and H1 x2 alone
#   neg_abcdump | neg_abckey | neg_abcctx  negative controls: A/B/C dump LOADed, A/B/C key, A/B/C re-level string
#   b0 | b1 | b2                           base (433h - LOAD - CTX), + LOAD, + LOAD + CTX + VCC_DE_MAG
#   level_b0 | level_b1 | level_b2         ntclevel417.py on the artifact (b1's prints the CTX string b2 uses)
#   audit_b0 | audit_b1 | audit_b2         audit428 + analyse431 + postbuild (strict allow-list) on the artifact
#   post_b0 | post_b1 | post_b2            the driver's postbuild stage + def_postbuild.py with a strict allow-list vs 433h
#   prep_b0 | prep_b2                      the driver's `prep` stage (vcc prep --dry-run --contexts D,E,F)
# Every builder invocation: memory/builder/prep/disk gate -> driver repin-frozen + pin-artifacts (arm env UNSET) ->
# run_builder432.py (the builder + ONE in-memory PERT_SCALE_ARTIFACTS entry; the file is never written) for
# repin-panel-support and build -> byte-restore from the pin record -> assert the builder is HEAD's 2c90d9e1 again.
# Nothing is submitted; nothing is written outside experiments/432-def-rehearsal-new-knobs/.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=$REPO/experiments/432-def-rehearsal-new-knobs; W=$X/work
BUNDLE=$REPO/experiments/199-def-runbook/work/def-bundle-r
# SET=ef: the low-coverage extension on 238's essential-free pessimistic bundle (K562 146/300, CD4 union 186), reusing
# 244's operator fitted for it (373df206..., axis f060e926...; read-only) and 244's CD4 banks; work dirs get an ef_ prefix.
SET="${SET:-r}"; P=""
HEADB=2c90d9e1cc2c19fd427da35ffe6a6d1306d3097382ff7a343233423cb710573f
EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
DRV=27edc486232afff2e19c46bdd87c8c84322cdc4b7781f55c6e5d40adcd8c180f
WOUT=$W/base/whitener_union_pruned_def.npz; WMETA=$W/base/whitener_union_pruned_def.json; CD4=$W/base/cd4_compact
KEY=$W/key/pert_scale_def432_lo0.3_hi3.npz; KEYID=wshare_lo0.3_hi3_def
ABCDUMP=$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz
STACKDUMP=$W/dump_stack/dump_h1x2_cd4x15.npz
if [ "$SET" = ef ]; then
  P=ef_; BUNDLE=$REPO/experiments/238-essential-free-panel/work/ef-bundle
  WOUT=$REPO/experiments/244-october-rehearsal/work/run1/whitener_union_pruned_def.npz
  WMETA=$REPO/experiments/244-october-rehearsal/work/run1/whitener_union_pruned_def.json
  CD4=$REPO/experiments/244-october-rehearsal/work/run1/cd4_compact
  KEY=$W/key/pert_scale_ef432_lo0.3_hi3.npz; STACKDUMP=$W/ef_dump_stack/dump_h1x2_cd4x15.npz
fi
G_DEF=0.3443655278649641; G_DEF_X15=0.5165482917974462      # DEF_ENV's carrying-row CD4 dose, and x1.5
PHASE="${1:?phase}"
LOG=$W/state432.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "[${P:-}$PHASE] $*" | tee -a "$LOG"; }
bsha(){ shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1; }
disk_gib(){ df -k "$REPO" | awk 'NR==2{printf "%.1f", $4/1048576}'; }
gate_ok(){
  local pct fi others preps disk need="${1:-20}"
  pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  others="$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null || echo 99)"
  preps="$(ps -axo args= | grep -vE 'claude|grep' | grep -cE '(^|/| )vcc prep( |$)' | tr -d ' ')"
  disk="$(disk_gib)"
  log "gate: free ${pct}%  free+inactive ${fi} GiB  other builders ${others}  vcc preps ${preps}  disk ${disk} GiB (need ${need})"
  [ "${pct:-0}" -ge 40 ] && awk -v v="$fi" 'BEGIN{exit !(v>=23)}' && [ "$others" = "0" ] \
    && [ "$preps" = "0" ] && awk -v v="$disk" -v n="$need" 'BEGIN{exit !(v>=n)}'
}
wait_gate(){ local n=0; until gate_ok "${1:-20}"; do n=$((n+1)); [ $n -ge 240 ] && { log "GATE never opened"; exit 3; }; sleep 30; done; }
clean_env(){ local v; for v in $(env | grep -o '^VCC_[A-Z0-9_]*'); do unset "$v"; done; }
# 433h's fixed base (run_arm_generic.sh) with the D/E/F key id; the caller adds the arm's own knobs on top.
base_env(){ clean_env; export VCC_NORMREL_EXP=-1.5 VCC_CD4_PARTIAL_G=0.5 VCC_NORMREL_MODE=cap VCC_DE_SIGN_PRIOR=clean4q80 \
  VCC_PERT_SCALE=$KEYID VCC_PERT_SCALE_BALLAST=1 VCC_BULK_SHARPEN_Q=1.25; }
# the driver's day defaults (SHIPPED_ENV + DEF_ENV), resolved exactly as def_pipeline.resolve_env does: exported wins
driver_env(){ eval "$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P
r,_=P.resolve_env(P.day_env(False)); print("\n".join(f"export {k}={v}" for k,v in r.items()))')"; }
check_hashes(){
  local e d; e="$(shasum -a 256 src/vcc2026/rc_emit_k.py | cut -d' ' -f1)"; d="$(shasum -a 256 scripts/def_pipeline.py | cut -d' ' -f1)"
  [ "$e" = "$EM" ] || { log "EMITTER $e != $EM -- stopping"; exit 4; }
  [ "$d" = "$DRV" ] || { log "DRIVER $d != $DRV -- stopping"; exit 4; }
  [ "$(bsha)" = "$HEADB" ] || { log "BUILDER $(bsha) != HEAD $HEADB (someone else's pin live?) -- stopping"; exit 4; }
}
WX=""
restore_pin(){
  [ -n "$WX" ] && [ -f "$WX/pin-artifacts.json" ] || return 0
  if [ "$(bsha)" != "$HEADB" ]; then
    uv run python scripts/def_pin_panel_artifacts.py --restore --record "$WX/pin-artifacts.json" 2>&1 | tee -a "$LOG"
  fi
  log "builder after restore: $(bsha) ($([ "$(bsha)" = "$HEADB" ] && echo 'HEAD, byte-identical' || echo 'NOT HEAD -- INVESTIGATE'))"
}
trap restore_pin EXIT
trap 'exit 143' TERM INT
# pin -> repin-panel-support -> build (both through the wrapper) -> restore. $1 = work name, $2 = expect build to fail (0/1)
builder_run(){
  local name="$1" expect_fail="${2:-0}" pinned rc common out
  WX=$W/$name; mkdir -p "$WX/logs"
  check_hashes; wait_gate 20; check_hashes
  log "builder before: $(bsha)  env: $(env | grep '^VCC_' | sort | tr '\n' ' ')"
  # the pin stages run with the ARM env unset: the unpatched builder refuses VCC_PERT_SCALE=$KEYID at import
  ( clean_env; /usr/bin/time -l uv run python scripts/def_pipeline.py --bundle "$BUNDLE" --work "$WX" --name "$name" \
      --frozen-inputs "$WX/frozen_inputs.json" --whitener-id union_pruned_def --whitener-out "$WOUT" --cd4-out "$CD4" \
      --min-free-gib 23 --stages repin-frozen,pin-artifacts > "$WX/driver_pin.log" 2>&1 )
  rc=$?; log "driver repin-frozen,pin-artifacts rc=$rc"
  [ $rc -eq 0 ] || { log "PIN FAILED -- stopping"; exit 1; }
  pinned="$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["builder_sha256_after"])' "$WX/pin-artifacts.json")"
  [ "$(bsha)" = "$pinned" ] || { log "pinned sha mismatch $(bsha) != $pinned -- stopping"; exit 1; }
  log "pinned: $pinned"
  common=(--bundle "$BUNDLE" --frozen-inputs "$WX/frozen_inputs.json" --expect-contexts D,E,F --whitener "$WOUT"
          --whitener-metadata "$WMETA" --whitener-id union_pruned_def --cd4-root "$CD4" --max-rss-gib 24.0)
  wrap(){ /usr/bin/time -l uv run python "$X/run_builder432.py" --builder scripts/make_submission_pm_cd4.py --id "$KEYID" \
            --artifact "$KEY" --expect-builder-sha "$pinned" --record "$WX/builder432_$1.json" -- "${@:2}"; }
  driver_env
  log "resolved env: $(env | grep '^VCC_' | sort | tr '\n' ' ')"
  wrap repin "${common[@]}" --repin-panel-support --min-cpm 10 > "$WX/logs/repin-panel-support.log" 2>&1
  rc=$?; log "repin-panel-support rc=$rc  $(grep -hE 'repin\] panel_support pinned|panel-support\]' "$WX/logs/repin-panel-support.log" | tail -2 | tr '\n' ' ')"
  if [ $rc -ne 0 ]; then
    restore_pin
    if [ "$expect_fail" = "1" ]; then
      log "EXPECTED REFUSAL observed at repin-panel-support: $(grep -hE 'SystemExit|ValueError|Error|REFUS|differ|must' "$WX/logs/repin-panel-support.log" | grep -v '^\s*#' | tail -2 | tr '\n' ' ')"
      return 0
    fi
    log "REPIN FAILED -- stopping"; exit 1
  fi
  out="$WX/$name.h5ad"
  local t0=$SECONDS
  wrap build "${common[@]}" --min-cpm 10 --tail-spread 1.0 --keep-top 8000 --keep-spread 0.25 --overdispersion 4.0 \
      --emitter rc_callset -o "$out" > "$WX/logs/build.log" 2>&1
  rc=$?
  log "build rc=$rc  wall $((SECONDS - t0)) s  peak RSS $(awk '/maximum resident set size/{printf "%.2f GiB", $1/2^30}' "$WX/logs/build.log")"
  if [ $rc -eq 0 ] && [ "$expect_fail" = "0" ] && [ -s "$out.meta.json" ]; then postbuild_now "$name"; fi
  restore_pin
  [ "$(bsha)" = "$HEADB" ] || { log "BUILDER NOT RESTORED -- stopping"; exit 9; }
  if [ "$expect_fail" = "1" ]; then
    [ $rc -ne 0 ] && log "EXPECTED REFUSAL observed: $(grep -hE 'SystemExit|ValueError|Error:|REFUS|differ|must' "$WX/logs/build.log" | grep -v '^\s*#' | tail -2 | tr '\n' ' ')" \
                  || log "UNEXPECTED: the negative control BUILT (rc 0)"
    rm -f "$out" "$out.meta.json.tmp"
    return 0
  fi
  [ $rc -eq 0 ] || { log "BUILD FAILED: $(tail -3 "$WX/logs/build.log" | tr '\n' ' ')"; exit 1; }
}
postbuild_now(){   # $1 = work name. Run WHILE PINNED: def_postbuild.py reads WHITENER_REGISTRY from the builder source.
  local a="$1" REF ALLOW
  REF="$(ls "$REPO"/submissions/433h-*.h5ad.meta.json 2>/dev/null | head -1)"
  ALLOW="pert_scale,pert_scale_sha256,bulk_load,bulk_load_sha256,bulk_ctx_scale"
  { [ "$a" = "${P}b2" ] || [ "$a" = "ef_b1" ]; } && ALLOW="$ALLOW,de_mag"
  if [ -z "$REF" ]; then REF="$(ls "$REPO"/submissions/428h-*.h5ad.meta.json | head -1)"; ALLOW="$ALLOW,de_mde"; fi
  log "postbuild (pinned $(bsha | cut -c1-16)) reference meta $(basename "$REF"); strict allow-list $ALLOW (+ the script's own new-panel keys)"
  ( clean_env; uv run python scripts/def_pipeline.py --bundle "$BUNDLE" --work "$W/$a" --name "$a" \
      --frozen-inputs "$W/$a/frozen_inputs.json" --whitener-id union_pruned_def --whitener-out "$WOUT" --cd4-out "$CD4" \
      --incumbent-meta "$REF" --stages postbuild > "$W/$a/driver_postbuild.log" 2>&1 )
  log "driver postbuild rc=$?: $(grep -hE 'checks pass|^  FAIL' "$W/$a/logs/postbuild.log" | tr '\n' ' ' | cut -c1-900)"
  ( clean_env; uv run python scripts/def_postbuild.py --meta "$W/$a/$a.h5ad.meta.json" --bundle "$BUNDLE" --incumbent "$REF" \
      --cd4-root "$CD4" --allow-config-diff "whitener_id,$ALLOW" --out "$W/$a/postbuild_strict.json" > "$W/$a/postbuild_strict.log" 2>&1 )
  log "strict postbuild rc=$?: $(grep -hE 'checks pass|^  FAIL' "$W/$a/postbuild_strict.log" | tr '\n' ' ' | cut -c1-900)"
}
arms_and_meta(){   # 240's arm assertion off the written meta + the headline meta fields
  local meta="$W/$1/$1.h5ad.meta.json"
  ( clean_env; uv run python - "$meta" <<'PY'
import json, sys
sys.path.insert(0, "scripts")
from pathlib import Path
import def_pipeline as P
r, _ = P.resolve_env(P.day_env(False))
P.assert_build_carries_arms(Path(sys.argv[1]), r)
m = json.loads(Path(sys.argv[1]).read_text())
c = m["config"]; o = m["output"]
print("config keys:", {k: c.get(k) for k in ("contexts", "pert_scale", "pert_scale_sha256", "pert_scale_ballast", "bulk_sharpen_q",
      "bulk_load", "bulk_load_sha256", "bulk_ctx_scale", "de_mde", "de_zero_repair", "de_mag", "cd4_sizing", "cd4_target_g",
      "cd4_deflate_basis", "h1_union_g", "normrel_mode", "normrel_exp", "cd4_partial_g", "whitener_id")})
print("operator:", {k: m["operator"].get(k) for k in ("sha256", "axis_sha256", "panel_sha256") if k in m["operator"]})
print("output sha256:", o.get("sha256"), "rows", o.get("n_rows"), "peak_rss_gib", m.get("peak_rss_gib"))
rc = o.get("rc_callset", {})
print("rc diagnostics:", {k: rc[k] for k in rc if k.startswith(("ballast", "mde_", "zero_", "zs_", "prior_lambda", "call_set_unsupp"))})
e = m.get("effects", {})
print("effects:", {k: (e[k] if not isinstance(e[k], dict) else {kk: e[k][kk] for kk in list(e[k])[:12]}) for k in e if k in ("pert_scale", "bulk_load", "h1_union", "cd4")})
PY
  ) 2>&1 | tee -a "$W/$1/meta_summary.txt"
}

case "$PHASE" in
  whitener)
    check_hashes; wait_gate 20
    ( clean_env; /usr/bin/time -l uv run python scripts/def_pipeline.py --bundle "$BUNDLE" --work "$W/base" --name base432 \
        --frozen-inputs "$W/base/frozen_inputs.json" --whitener-id union_pruned_def --stages whitener > "$W/base/driver_whitener.log" 2>&1 )
    rc=$?; log "whitener rc=$rc  $(grep -hE 'axis sha256|panel sha256|union axis|dropped_in_panel|sha256 [0-9a-f]{64}' "$W/base/logs/whitener.log" | head -6 | tr '\n' ' ')"
    rm -rf "$W/base/.whitener_scratch"; log "scratch removed; operator $(shasum -a 256 "$WOUT" | cut -d' ' -f1)"
    [ $rc -eq 0 ] || exit 1 ;;
  key)
    /usr/bin/time -l uv run python "$X/key432.py" --bundle "$BUNDLE" --op "$WOUT" --lo 0.3 --hi 3 --out "$KEY" > "$W/key/key_$SET.log" 2>&1
    rc=$?; log "key rc=$rc sha $(shasum -a 256 "$KEY" 2>/dev/null | cut -d' ' -f1)  $(grep -E 'maximum resident|axis check' "$W/key/key_$SET.log" | tr '\n' ' ')"
    [ $rc -eq 0 ] || exit 1 ;;
  dump_stack|dump_lit|dump_h1|dump_b0env)
    base_env
    case "$PHASE" in
      dump_b0env) : ;;   # b0's own bulk env (H1 x1, CD4 x1): the LOAD known answer
      dump_stack) export VCC_H1_UNION_G=0.56 VCC_CD4_TARGET_G=$G_DEF_X15 ;;
      dump_lit)   export VCC_H1_UNION_G=0.56 VCC_CD4_TARGET_F=0.51 ;;
      dump_h1)    export VCC_H1_UNION_G=0.56 ;;
    esac
    out="$W/$P$PHASE/dump_h1x2_cd4x15.npz"; [ "$PHASE" = dump_h1 ] && out="$W/$P$PHASE/dump_h1x2.npz"
    [ "$PHASE" = dump_lit ] && out="$W/$P$PHASE/dump_h1x2_cd4f051.npz"
    [ "$PHASE" = dump_b0env ] && out="$W/$P$PHASE/dump_b0env.npz"
    mkdir -p "$W/$P$PHASE"; rm -f "$out"; export VCC_BULK_DUMP="$out"
    builder_run "$P$PHASE" 0 || true
    L="$(grep -h '\[bulk-dump\] wrote' "$W/$P$PHASE/logs/build.log" | tail -1)"
    log "cd4/h1 lines: $(grep -hE '^\[cd4\] sizing|^\[h1-union\]|^\[cd4-partial|^\[pert-scale\]' "$W/$P$PHASE/logs/build.log" | tr '\n' ' ' | cut -c1-900)"
    [ -s "$out" ] && [ -n "$L" ] && log "DUMP OK sha256 $(shasum -a 256 "$out" | cut -d' ' -f1)" || { log "DUMP FAILED"; exit 1; } ;;
  neg_abcdump)
    base_env; export VCC_BULK_LOAD="$ABCDUMP:$(shasum -a 256 "$ABCDUMP" | cut -d' ' -f1)"
    builder_run "${P}neg_abcdump" 1 ;;
  neg_abckey)
    base_env; export VCC_PERT_SCALE=wshare_lo0.3_hi3
    builder_run "${P}neg_abckey" 1 ;;
  neg_abcctx)
    base_env; export VCC_BULK_CTX_SCALE="A=1.0088,B=1.05,C=1.0736"
    builder_run "${P}neg_abcctx" 1 ;;
  b0|b1|b2)
    base_env; export VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_RC_AUDIT_DUMP="$W/$P$PHASE/rc_audit_$P$PHASE.jsonl"
    if [ "$PHASE" != b0 ]; then export VCC_BULK_LOAD="$STACKDUMP:$(shasum -a 256 "$STACKDUMP" | cut -d' ' -f1)"; fi
    if [ "$SET" = ef ] && [ "$PHASE" = b1 ]; then export VCC_DE_MAG=dn_k250_zs; fi
    if [ "$PHASE" = b2 ]; then
      CTX="$(grep -h '^VCC_BULK_CTX_SCALE' "$W/${P}b1/${P}b1.ntclevel.out" | sed 's/.*= //')"
      [ -n "$CTX" ] || { log "no CTX string from level_b1 -- stopping"; exit 1; }
      export VCC_BULK_CTX_SCALE="$CTX" VCC_DE_MAG=dn_k250_zs
    fi
    mkdir -p "$W/$P$PHASE"; rm -f "$W/$P$PHASE/rc_audit_$P$PHASE.jsonl"
    builder_run "$P$PHASE" 0
    arms_and_meta "$P$PHASE"
    log "artifact $(ls -la "$W/$P$PHASE/$P$PHASE.h5ad" | awk '{print $5}') bytes; disk $(disk_gib) GiB" ;;
  level_b0|level_b1|level_b2)
    a="$P${PHASE#level_}"; wait_gate 20
    /usr/bin/time -l uv run python experiments/417-bulk-sharpen/ntclevel417.py "$W/$a/$a.h5ad" --bundle "$BUNDLE" \
        --contexts D,E,F --sd 0.01 > "$W/$a/$a.ntclevel.out" 2>&1
    rc=$?; log "ntclevel rc=$rc  $(grep -E '^  [DEF]:|VCC_BULK_CTX_SCALE|maximum resident' "$W/$a/$a.ntclevel.out" | tr '\n' ' ')" ;;
  audit_b0|audit_b1|audit_b2)
    a="$P${PHASE#audit_}"
    uv run python experiments/428-zero-gene-audit/audit428.py "$W/$a/$a.h5ad" --bundle "$BUNDLE" > "$W/$a/audit428.out" 2>&1
    log "audit428 rc=$?: $(head -3 "$W/$a/audit428.out" | tr '\n' ' ' | cut -c1-600)"
    uv run python experiments/431-rc-audit-dump/analyse431.py "$W/$a/rc_audit_$a.jsonl" > "$W/$a/analyse431.out" 2>&1
    log "analyse431 rc=$?" ;;
  post_b0|post_b1|post_b2)
    a="$P${PHASE#post_}"
    if [ -s "$W/$a/postbuild_strict.json" ]; then
      log "postbuild ran while pinned: $(grep -hE 'checks pass' "$W/$a/logs/postbuild.log" "$W/$a/postbuild_strict.log" | tr '\n' ' ')"
    else
      log "no pinned postbuild on record; running unpinned (the registry row reads the RESTORED builder and will FAIL)"
      postbuild_now "$a"
    fi ;;
  prep_b0|prep_b2)
    a="$P${PHASE#prep_}"; wait_gate 20
    ( clean_env; /usr/bin/time -l uv run python scripts/def_pipeline.py --bundle "$BUNDLE" --work "$W/$a" --name "$a" \
        --frozen-inputs "$W/$a/frozen_inputs.json" --whitener-id union_pruned_def --whitener-out "$WOUT" --cd4-out "$CD4" \
        --min-free-gib 23 --stages prep > "$W/$a/driver_prep.log" 2>&1 )
    log "prep rc=$?  $(grep -hE 'dry run|cells|D=|targets|counts|maximum resident' "$W/$a/logs/prep.log" | tr '\n' ' ' | cut -c1-600)" ;;
  *) echo "unknown phase $PHASE"; exit 2 ;;
esac
log "phase done"
