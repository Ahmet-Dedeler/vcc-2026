#!/usr/bin/env bash
# 658a: build ALLSW (experiments/658-prereg-1007-builds/PREREGISTRATION-ALLSW.md). It is 498p2g's round-2 production build
# (498-combined-bulk/chain498.sh build()) with VCC_DE_MAG=dn_k250_zs, plus VCC_DE_MDE_SIGN=clean4:0.05 iff
# docs/october-decisions.json has mde_sign DECIDED clean4:0.05 when this starts. tercilec is never added (prereg 1.1).
# Then the pre-slot bulk identity check (bulkid658 vs 498p2g's restored package). Never submits; unarmed (a chain runs it).
#
# The LAST line of the log is always ONE terminal line:
#   <ts> READY allsw658 <package sha256>
#   <ts> NOT READY <reason>
# The exact env and arguments of the build are recorded in work/allsw/env.json.
#
# usage: ALLSW_BUILDER_SHA=<post-10-06 builder sha> ALLSW_EMITTER_SHA=<post-10-06 emitter sha> \
#          [ALLSW_FLOOR=YYYYmmddHHMM] [ALLSW_CAP=YYYYmmddHHMM] [ALLSW_EXTRA_EXPECT=k1,k2 ALLSW_EXTRA_REASON="..."] \
#          bash experiments/658-prereg-1007-builds/build_allsw658.sh
#   DRY=1: check every input and print the plan (env, EXPECT, gates); build nothing. The shas default to the current files.
# Waits are capped by ALLSW_CAP (a timestamp), never by an iteration count (AGENTS 6.2). bash 3.2 / BSD grep safe.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/658-prereg-1007-builds; W=$X/work/allsw; LOG=$W/build_allsw658.log; E=experiments/428-zero-gene-audit
BU=scripts/make_submission_pm_cd4.py; EMP=src/vcc2026/rc_emit_k.py; DEC=docs/october-decisions.json
TAG=allsw658; NAME=allsw658
R498=submissions/498p2g-cap15-partial50-clean4q80-lodo6-ballast-sharpenp2g-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5.h5ad.meta.json
V498=submissions/498p2g-cap15-partial50-clean4q80-lodo6-ballast-sharpenp2g-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5.vcc
N2=experiments/495-lodo-cap6/work/dump_h1x2_cd4x15_lodo05hi6.npz
DS2=d58750b83915d40d488125926614c83793bd6d7260ceef1192b4b10c19eeb45c
F1="A=1.0687,B=1.1223,C=1.1259"; QV="tercile:2.5,2.0,1.75"
DRY="${DRY:-0}"; FLOOR="${ALLSW_FLOOR:-202610070000}"; CAP="${ALLSW_CAP:-202610081200}"
PFX=""; [ "$DRY" = 1 ] && PFX="dry-"          # a DRY run never writes the real run's log or env.json
LOG=$W/${PFX}build_allsw658.log
mkdir -p "$W"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
TERMINAL=0
finish(){ TERMINAL=1; log "$1"; exit "${2:-0}"; }
trap 'rc=$?; [ "$TERMINAL" = 1 ] || log "NOT READY unexpected exit (rc $rc) -- see the lines above"' EXIT
eq_head(){ [ "$(shasum -a 256 "$1" | cut -d' ' -f1)" = "$(git show HEAD:"$1" | shasum -a 256 | cut -d' ' -f1)" ]; }
free_gib(){ df -k "$REPO" | awk 'NR==2{printf "%d", $4/1048576}'; }
builders(){
  # recipe_dr524.sh's detector: bash launchers by argv word, plus def_pipeline.other_builders()
  ps -axo comm,args | awk '$1 ~ /^bash$/ { for (i = 2; i <= NF; i++) if ($i ~ /(^|\/)(chain4[0-9][0-9][a-z0-9]*\.sh|chain5[0-9][0-9][a-z0-9]*\.sh|relevel_sp2\.sh|run_ka_[0-9][^[:space:]]*|run_arm_generic\.sh|def_day\.sh|run_dump4[^[:space:]]*|recipe_dr524[a-z0-9]*\.sh)$/) found = 1 } found { print; exit }' | /usr/bin/grep -q . && return 0
  uv run python -c 'import sys; sys.path.insert(0, "scripts"); import def_pipeline; raise SystemExit(0 if def_pipeline.other_builders() else 1)'
}
echo $$ > "$W/${PFX}build_allsw658.pid"
log "start (pid $$, DRY=$DRY): ALLSW = 498p2g's round-2 build + zs (+ the hook iff DECIDED); floor $FLOOR cap $CAP"

# ---- fixed inputs: 498p2g's meta, its package, the cap-6 dump
[ -s "$R498" ] && [ -s "$V498" ] || finish "NOT READY 498p2g's meta or package missing"
[ "$(shasum -a 256 "$N2" | cut -d' ' -f1)" = "$DS2" ] || finish "NOT READY the cap-6 dump's sha is not $DS2"
CHK="$(uv run python - "$R498" "$DS2" "$F1" <<'PY'
import json, sys
c = json.load(open(sys.argv[1]))["config"]
f = {k: float(v) for k, v in (p.split("=") for p in sys.argv[3].split(","))}
bad = []
if c.get("bulk_load_sha256") != sys.argv[2]: bad.append("bulk_load_sha256")
if any(abs(c["bulk_ctx_scale"][k] - f[k]) > 1e-12 for k in "ABC"): bad.append("bulk_ctx_scale")
pt = c.get("bulk_sharpen_pt") or {}
if pt.get("mode") != "tercile" or pt.get("q") != [2.5, 2.0, 1.75]: bad.append("bulk_sharpen_pt")
if c.get("de_mag") is not None or c.get("de_mde_sign") is not None: bad.append("de_mag/de_mde_sign already set")
print("ok" if not bad else "bad " + ",".join(bad))
PY
)"
[ "$CHK" = ok ] || finish "NOT READY 498p2g's meta does not carry the expected recipe ($CHK)"

# ---- builder / emitter: the committed post-10-06 tree
B="$(shasum -a 256 $BU | cut -d' ' -f1)"; EM="$(shasum -a 256 $EMP | cut -d' ' -f1)"
if [ "$DRY" = 1 ]; then WANT_B="${ALLSW_BUILDER_SHA:-$B}"; WANT_EM="${ALLSW_EMITTER_SHA:-$EM}"
else
  [ -n "${ALLSW_BUILDER_SHA:-}" ] && [ -n "${ALLSW_EMITTER_SHA:-}" ] || finish "NOT READY ALLSW_BUILDER_SHA / ALLSW_EMITTER_SHA not given (the post-10-06 shas from APPLY DAY step 8)"
  WANT_B="$ALLSW_BUILDER_SHA"; WANT_EM="$ALLSW_EMITTER_SHA"
fi

# ---- the hook: only if DECIDED clean4:0.05
HOOK="$(uv run python - "$DEC" <<'PY'
import json, sys
s = json.load(open(sys.argv[1]))["switches"]["mde_sign"]
print("on" if s.get("status") == "DECIDED" and s.get("value") == "clean4:0.05" else "off")
print(json.dumps({k: s.get(k) for k in ("status", "value", "decided_by")})[:300])
PY
)"
HOOK_REC="$(printf '%s\n' "$HOOK" | sed -n 2p)"; HOOK="$(printf '%s\n' "$HOOK" | sed -n 1p)"
case "$HOOK" in on|off) ;; *) finish "NOT READY could not read mde_sign from $DEC";; esac
log "mde_sign record: $HOOK_REC -> hook $HOOK"

# ---- EXPECT: de_mag (+ de_mde_sign) (+ de_sign_emitter_sha256 iff the builder's pinned emitter differs from 498p2g's)
PIN_NOW="$(/usr/bin/grep -E '^DE_SIGN_EMITTER_SHA256 = "' $BU | head -1 | sed -E 's/^DE_SIGN_EMITTER_SHA256 = "([0-9a-f]{64})".*/\1/')"
PIN_498="$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("de_sign_emitter_sha256",""))' "$R498")"
EXPECT="de_mag"
[ "$HOOK" = on ] && EXPECT="$EXPECT,de_mde_sign"
[ -n "$PIN_NOW" ] && [ "$PIN_NOW" != "$PIN_498" ] && EXPECT="$EXPECT,de_sign_emitter_sha256"
if [ -n "${ALLSW_EXTRA_EXPECT:-}" ]; then
  [ -n "${ALLSW_EXTRA_REASON:-}" ] || finish "NOT READY ALLSW_EXTRA_EXPECT given without ALLSW_EXTRA_REASON"
  EXPECT="$EXPECT,$ALLSW_EXTRA_EXPECT"; log "EXTRA EXPECT $ALLSW_EXTRA_EXPECT: $ALLSW_EXTRA_REASON"
fi
SUF="-zs250"; [ "$HOOK" = on ] && SUF="$SUF-hk005"
STEM="658allsw-cap15-partial50-clean4q80-lodo6-ballast-sharpenp2g-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5$SUF"
log "builder now ${B:0:16} (want ${WANT_B:0:16}), emitter now ${EM:0:16} (want ${WANT_EM:0:16}), builder emitter pin ${PIN_NOW:0:16} vs 498p2g's ${PIN_498:0:16}"
log "EXPECT=$EXPECT  STEM=$STEM"

# ---- the env, recorded before anything runs
uv run python - "$W/${PFX}env.json" "$REPO/$N2:$DS2" "$F1" "$HOOK" "$TAG" "$QV" "$REPO/$R498" "$NAME" "$STEM" "$WANT_B" "$WANT_EM" "$EXPECT" <<'PY'
import json, sys
o, load, f1, hook, *a = sys.argv[1:]
env = {"VCC_BULK_LOAD": load, "VCC_DE_MDE": "0.1,0.1,0.3,coarse", "VCC_DE_ZERO_REPAIR": "-5", "VCC_BULK_CTX_SCALE": f1,
       "VCC_DE_MAG": "dn_k250_zs"}
if hook == "on":
    env["VCC_DE_MDE_SIGN"] = "clean4:0.05"
unset = ["VCC_SAMELINE", "VCC_DE_PRIOR_GATE", "VCC_UNION_SOURCES", "VCC_RC_AUDIT_DUMP"] + ([] if hook == "on" else ["VCC_DE_MDE_SIGN"])
json.dump({"exports": env, "unset": unset, "run_arm_generic_args": a,
           "fixed_by_run_arm_generic": "VCC_NORMREL_EXP -1.5, VCC_CD4_PARTIAL_G 0.5, VCC_NORMREL_MODE cap, VCC_DE_SIGN_PRIOR clean4q80, "
                                       "VCC_PERT_SCALE wshare_lo0.3_hi3, VCC_PERT_SCALE_BALLAST 1, VCC_BULK_SHARPEN_Q = arg 2"},
          open(o, "w"), indent=1)
PY
log "env recorded: $W/${PFX}env.json"

if [ "$DRY" = 1 ]; then
  log "DRY plan: builders now: $(builders && echo busy || echo none); disk $(free_gib) GiB (needs >= 35 to start)"
  eq_head $BU && eq_head $EMP && log "DRY builder and emitter == HEAD" || log "DRY builder or emitter differs from HEAD"
  finish "NOT READY DRY=1 (plan only; nothing built)"
fi

# ---- wait: floor, no other build, builder/emitter == HEAD == the post-10-06 shas, disk >= 35 GiB (capped by CAP)
log "waiting for >= $FLOOR, no other build, builder/emitter == HEAD and == the given shas, disk >= 35 GiB"
until [ "$(date -u +%Y%m%d%H%M)" -ge "$FLOOR" ] && ! builders && eq_head $BU && eq_head $EMP \
      && [ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$WANT_B" ] && [ "$(shasum -a 256 $EMP | cut -d' ' -f1)" = "$WANT_EM" ] \
      && [ "$(free_gib)" -ge 35 ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && finish "NOT READY prerequisites not met by $CAP UTC"
  sleep 60; done
[ -e "$E/work/build_$TAG" ] && mv "$E/work/build_$TAG" "$E/work/build_$TAG.aborted-$(date -u +%Y%m%dT%H%MZ)"
log "released: building"
( unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_RC_AUDIT_DUMP VCC_DE_MDE_SIGN VCC_DE_MAG
  export VCC_BULK_LOAD="$REPO/$N2:$DS2" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F1" VCC_DE_MAG=dn_k250_zs
  [ "$HOOK" = on ] && export VCC_DE_MDE_SIGN=clean4:0.05
  bash $E/run_arm_generic.sh "$TAG" "$QV" "$REPO/$R498" "$NAME" "$STEM" "$WANT_B" "$WANT_EM" "$EXPECT" ) >> "$W/build.nohup" 2>&1
RL="$E/work/build_$TAG/run$TAG.log"; H5="$E/work/build_$TAG/$NAME.h5ad"; PKG="submissions/$STEM.vcc"
[ -s "$RL" ] || finish "NOT READY run_arm_generic left no log ($RL)"
log "build: $(/usr/bin/grep -E "^\[$TAG\] (PASS|FAIL|CLIFF)|package sha256|novelty rc|stopping|GATE never" "$RL" | tr '\n' ' ' | cut -c1-700)"
/usr/bin/grep -qE "^\[$TAG\] PASS: exactly" "$RL" || finish "NOT READY config diff is not exactly $EXPECT (see $RL)"
/usr/bin/grep -qE "^\[$TAG\] CLIFF OK" "$RL" || finish "NOT READY cliff check did not pass (see $RL)"
/usr/bin/grep -qE "novelty rc=0" "$RL" || finish "NOT READY novelty rc != 0 (see $RL)"
PKGSHA="$(/usr/bin/grep -E "package sha256: [0-9a-f]{64}" "$RL" | tail -1 | sed -E 's/.*package sha256: ([0-9a-f]{64}).*/\1/')"
[ -s "$PKG" ] && [ -n "$PKGSHA" ] && [ "$(shasum -a 256 "$PKG" | cut -d' ' -f1)" = "$PKGSHA" ] || finish "NOT READY package missing or its sha does not match the build log"
log "package $PKG sha256 $PKGSHA"

# ---- pre-slot bulk identity vs 498p2g (F1 before the slot)
[ -s "$H5" ] || finish "NOT READY the built h5ad is gone before the bulk check ($H5)"
[ "$(free_gib)" -ge 20 ] || finish "NOT READY disk $(free_gib) GiB < 20 for the restored 498p2g copy (12 GB); free space and rerun the check by hand"
RST="$W/r498p2g.h5ad"; rm -f "$RST"
uv run python scripts/restore_prediction_h5ad.py "$V498" "$RST" >> "$LOG" 2>&1 || { rm -f "$RST"; finish "NOT READY could not restore 498p2g's package"; }
uv run python $X/bulkid658.py "$H5" "$RST" --out "$W/bulkid_allsw.json" >> "$LOG" 2>&1; BID=$?
rm -f "$RST"
[ "$BID" = 0 ] || finish "NOT READY bulk is NOT 498p2g's (bulkid658 rc $BID; $W/bulkid_allsw.json): F1 would fire -- diagnose, do not submit"
log "bulk identity vs 498p2g: IDENTICAL"
rm -f "$H5"     # the package exists, its sha is logged and the bulk check passed
finish "READY $TAG $PKGSHA"
