#!/usr/bin/env bash
# 690 regate_r0 (readiness [401d9b], 2026-10-22; the coordinator's call 05:3xZ: apply 554's preregistered override
# DEF_DAY_ALLOW_LOW_DOSE=1 by RE-GATING the finished R0 build, the 697/698/673 precedent; never rebuilds).
# R0's def_day (experiments/692-def-r0/work/r0) CONVERGED at lvl3 (rho D 0.9762 / E 0.9764 / F 0.9749; output fa832988…) and
# STOPPED on one postbuild row of 42: 554's steering dose, D median prior_lambda 0.8975 < 0.9 (MDE 0.954, missed 1; E/F 1.0).
# 554 RESULTS: a context at 0.85-0.9 is a ~10-15 % dose loss; "stop, look, and override with DEF_DAY_ALLOW_LOW_DOSE=1 if
# nothing else is wrong". def_day stops at postbuild, so its later audits never ran. This script runs them on the finished build:
#   0. preconditions: r0's def_day.log has the postbuild STOP, the restored builder (want == got), the final line for
#      r0def-lvl3 at D=1.0737,E=0.9063,F=1.0294, FAILSAFE_TOTAL 0; no def_day runs; the final meta's output is fa832988…;
#      r0/builder.pinned.py (690c: the pin def_day's postbuild read) is the meta's builder sha; postbuild.json fails ONLY 554's
#      row with every lambda >= 0.85; the fresh evidence dir does not exist; submissions/r0def.vcc does not exist;
#   1. KA: def_postbuild re-run exactly as def_day ran it (same meta, bundle, cd4-root, incumbent, allow list, rc-audit; the
#      pinned builder via --builder instead of the then-pinned production file), no override: rc 1 and the SAME 42 rows
#      (check, pass, value) as r0/postbuild.json. A different row = the reconstruction is wrong: NOT READY;
#   2. the override: the same call with DEF_DAY_ALLOW_LOW_DOSE=1: rc 0, and versus the KA exactly ONE row differs, 554's, now
#      passing with "ALLOWED by DEF_DAY_ALLOW_LOW_DOSE=1";
#   3. def_day's step 6 audit: the zero-gene audit (audit428) must read "ZERO block total (target excluded): 0 in ";
#      (def_day's driver prep stage is `vcc prep --dry-run`; step 5's real prep covers it);
#   4-7. r0ready692's steps 3-7 unchanged: FAILSAFE 0, DNMAG report, past_cliff on lvl3 == 1 (safe), the real `vcc prep`
#      under the same gate (no builder/operator fit, no other vcc prep, free+inactive >= 26 GiB: it waits out day708's fits),
#      novelty, provenance; the raw "[r0def] PASS" / "[r0def] CLIFF OK" lines arm1015 greps;
#   8. "READY r0def <sha>" into the FRESH ready log; commit the small records (explicit pathspec).
# Evidence dir (default experiments/692-def-r0/work/r0g-ready): r0def.run.log, r0def.ready.log, postbuild-ka.json,
# postbuild.json, audit428.txt, r0def.prep.log. The coordinator re-arms arm1015 R0 on its ready/run logs
# (rearm_r0b.sh with R0B_READY_DIR=<dir> R0B_NO_692=1). Terminal line, always: "<ts> END" in the run log.
#   nohup bash scripts/detach.sh bash experiments/690-def-readiness-1015/regate_r0.sh > experiments/690-def-readiness-1015/work/regate_r0.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
if [ -z "${REGATE_R0_FROZEN:-}" ]; then
  FZ="experiments/690-def-readiness-1015/work/.frozen-$$-regate_r0.sh"; cp experiments/690-def-readiness-1015/regate_r0.sh "$FZ" || exit 1
  REGATE_R0_FROZEN=1 exec bash "$FZ" "$@"
fi
# 10-22 06:2xZ: parametrised for X1 (6a (d), the coordinator's ruling d98c1c48: a 554-only postbuild FAIL with every lambda
# >= 0.85 is re-gated exactly as R0 was). Defaults are R0's run; X1: REGATE_W=experiments/692-def-r0/work/x1 REGATE_N=x1def
# REGATE_OUT=experiments/692-def-r0/work/x1g-ready REGATE_KA_SHA=meta (the final meta's own output; FT/FF from its final line).
B=data/challenge-def; W="${REGATE_W:-experiments/692-def-r0/work/r0}"; N="${REGATE_N:-r0def}"; DL=$W/def_day.log
OUT="${REGATE_OUT:-experiments/692-def-r0/work/r0g-ready}"; PKGDIR=submissions; PKG=$PKGDIR/$N.vcc; CAP="${REGATE_CAP:-202610270000}"
KA_SHA="${REGATE_KA_SHA:-fa83298813978cd8210d08f320f62d63d9501a630ca102e6310f1dfb95c07799}"
FTL="$(/usr/bin/grep -oE ' final: [^ ]+ factors [^ ]+ -> ' "$DL" 2>/dev/null | tail -1)"; FT="$(echo "$FTL" | awk '{print $2}')"; FF="$(echo "$FTL" | awk '{print $4}')"
[ "$KA_SHA" = meta ] && KA_SHA="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' "$W/${N}_final.h5ad.meta.json" 2>/dev/null)"
INC=submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json
H5="$W/${N}_final.h5ad"; META="$H5.meta.json"; PIN="$W/builder.pinned.py"
[ -e "$OUT" ] && { echo "STOP: $OUT exists (a re-gate writes a FRESH evidence dir)" >&2; rm -f "experiments/690-def-readiness-1015/work/.frozen-$$-regate_r0.sh"; exit 1; }
mkdir -p "$OUT"; RL="$OUT/$N.run.log"; RDY="$OUT/$N.ready.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$RL"; }
raw(){ printf '%s\n' "$*" >> "$RL"; }
notready(){ log "NOT READY: $*"; printf '%s NOT READY %s %s\n' "$(date -u +%FT%TZ)" "$N" "$*" >> "$RDY"; exit 1; }
now(){ date -u +%Y%m%d%H%M; }
trap 'log "END"; rm -f "experiments/690-def-readiness-1015/work/.frozen-$$-regate_r0.sh"' EXIT
log "regate_r0 armed (pid $$): re-gate of $H5 ($FT, $FF) with DEF_DAY_ALLOW_LOW_DOSE=1 (554); evidence $OUT; cap $CAP"

# 0. preconditions
/usr/bin/grep -q "STOP: postbuild FAILED" "$DL" || notready "r0's def_day.log has no postbuild STOP"
RB="$(/usr/bin/grep 'builder after restore:' "$DL" | tail -1)"
echo "$RB" | /usr/bin/grep -qE 'restore: ([0-9a-f]{16}) \(want \1\)' || notready "builder not restored to its pre-pin sha: $RB"
[ -n "$FT" ] && [ -n "$FF" ] && /usr/bin/grep -q " final: $FT factors $FF -> " "$DL" || notready "no 'final: <tag> factors <f>' line in $DL"
[ -n "$KA_SHA" ] || notready "no expected output sha"
[ "$(ps -axo comm=,args= | awk '$1 ~ /bash$/ && /scripts\/def_day\.sh/' | wc -l | tr -d ' ')" = 0 ] || notready "a def_day.sh is running"
[ -s "$H5" ] && [ -s "$META" ] || notready "missing $H5 or its meta"
[ "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' "$META")" = "$KA_SHA" ] || notready "final meta output is not $KA_SHA"
PB="$(shasum -a 256 "$PIN" | cut -d' ' -f1)"
MB="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["hashes"]["implementation"].get("scripts/make_submission_pm_cd4.py",""))' "$META")"
[ "$PB" = "$MB" ] || notready "builder.pinned.py ${PB:0:16} != the meta's builder ${MB:0:16}"
python3 - "$W/postbuild.json" >> "$RL" 2>&1 <<'PY' || notready "r0's postbuild fails something besides 554's dose row at >= 0.85"
import json, sys
d = json.load(open(sys.argv[1]))
bad = [c for c in d["checks"] if not c["pass"]]
assert len(bad) == 1 and "(554)" in bad[0]["check"], [c["check"] for c in bad]
g = d["dose_guard"]
for ctx, v in g["per_context"].items():
    for k in ("prior_lambda", "mde_lambda"):
        assert v[k] >= 0.85, (ctx, k, v[k])
    assert v["call_set_missed"] <= g["floors"]["call_set_missed"], (ctx, v["call_set_missed"])
print("precondition: r0's only failing postbuild row is 554's; every lambda >= 0.85; missed calls under the floor")
PY
[ -e "$PKG" ] && notready "$PKG exists (never overwrite a package)"
log "preconditions PASS: $RB; meta output ${KA_SHA:0:16}; pinned builder ${PB:0:16} == meta"

# 1. KA: def_day's postbuild call, reconstructed (def_day.sh 5b with today's env: same-line off, no prior gate, tercilec sharpen,
#    DE_MAG / UPRN / DNMAG / MDE_SIGN set, no keygate / callset block / cpm; the emitter key only if the pin moved it)
ALLOW="whitener_id,pert_scale,pert_scale_sha256,bulk_load,bulk_load_sha256,bulk_ctx_scale,bulk_sharpen_q,bulk_sharpen_pt,de_mag,de_uprn,de_dnmag,de_mde_sign"
[ -s "$W/dumphook/bulk_ctx.npz" ] && ALLOW="$ALLOW,bulk_load_ctx,bulk_load_ctx_sha256"   # 690g: a hook that wrote a per-context load
_inc_em="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("de_sign_emitter_sha256",""))' "$INC" 2>/dev/null || true)"
_src_em="$(python3 -c 'import re,sys; t=open(sys.argv[1]).read(); m=re.search(r"^DE_SIGN_EMITTER_SHA256 = .([0-9a-f]{64}).", t, re.M); print(m.group(1) if m else "")' "$PIN" 2>/dev/null || true)"
[ -n "$_src_em" ] && [ -n "$_inc_em" ] && [ "$_src_em" != "$_inc_em" ] && ALLOW="$ALLOW,de_sign_emitter_sha256"
PBARGS=(--meta "$META" --bundle "$B" --cd4-root "$W/cd4_compact" --builder "$PIN" --incumbent "$INC" --allow-config-diff "$ALLOW" \
        --rc-audit "$W/${N}_final.rc_audit.jsonl")
env -u DEF_DAY_ALLOW_LOW_DOSE uv run python scripts/def_postbuild.py "${PBARGS[@]}" --out "$OUT/postbuild-ka.json" > "$OUT/postbuild-ka.txt" 2>&1; K=$?
log "KA postbuild (no override) rc $K; allow $ALLOW"
[ "$K" = 1 ] || notready "KA postbuild rc $K (want 1, as def_day's)"
python3 - "$W/postbuild.json" "$OUT/postbuild-ka.json" >> "$RL" 2>&1 <<'PY' || notready "KA postbuild rows differ from def_day's: the reconstruction is wrong"
import json, sys
a = [(c["check"], c["pass"], c["value"]) for c in json.load(open(sys.argv[1]))["checks"]]
b = [(c["check"], c["pass"], c["value"]) for c in json.load(open(sys.argv[2]))["checks"]]
assert a == b, [x for x in zip(a, b) if x[0] != x[1]][:3]
print(f"KA PASS: the re-run reproduces def_day's {len(a)} postbuild rows exactly ({sum(p for _, p, _ in a)} pass)")
PY
# 2. the override
DEF_DAY_ALLOW_LOW_DOSE=1 uv run python scripts/def_postbuild.py "${PBARGS[@]}" --out "$OUT/postbuild.json" > "$OUT/postbuild.txt" 2>&1; P=$?
log "postbuild with DEF_DAY_ALLOW_LOW_DOSE=1 rc $P"
[ "$P" = 0 ] || notready "postbuild with the override rc $P"
python3 - "$OUT/postbuild-ka.json" "$OUT/postbuild.json" >> "$RL" 2>&1 <<'PY' || notready "the override changed more than 554's row"
import json, sys
a = json.load(open(sys.argv[1]))["checks"]; b = json.load(open(sys.argv[2]))["checks"]
assert [c["check"] for c in a] == [c["check"] for c in b]
diff = [(x, y) for x, y in zip(a, b) if (x["pass"], x["value"]) != (y["pass"], y["value"])]
assert len(diff) == 1 and "(554)" in diff[0][0]["check"] and diff[0][1]["pass"] and "ALLOWED by DEF_DAY_ALLOW_LOW_DOSE=1" in diff[0][1]["value"], diff
assert all(c["pass"] for c in b)
print("override: exactly one row changed, 554's, now ALLOWED: " + diff[0][1]["value"][:200])
PY
log "[$N] postbuild PASS (re-gate, 554 override DEF_DAY_ALLOW_LOW_DOSE=1; KA reproduced def_day's rows)"

# 3. def_day's zero-gene audit on the final build
until [ "$(uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(int(d.free_gib()))" 2>/dev/null || echo 0)" -ge 20 ]; do
  [ "$(now)" -ge "$CAP" ] && notready "no memory for audit428 by CAP"; sleep 60; done
uv run python experiments/428-zero-gene-audit/audit428.py "$H5" --bundle "$B" > "$OUT/audit428.txt" 2>&1; A=$?
summary="$(head -1 "$OUT/audit428.txt" 2>/dev/null || true)"; log "audit: rc $A: $summary"
[ "$A" = 0 ] || notready "zero-gene audit failed rc=$A"
case "$summary" in *"ZERO block total (target excluded): 0 in "*) ;; *) notready "zero-gene audit is not clean: $summary" ;; esac
raw "[$N] PASS postbuild PASS (re-gate with the 554 override) + zero-gene audit clean + builder restored"

# 4-7. r0ready692's steps 3-7
FS="$(/usr/bin/grep -oE 'uprn FAILSAFE_TOTAL -?[0-9]+' "$DL" | tail -1 | awk '{print $3}')"
[ "$FS" = 0 ] || notready "630dd FAILSAFE_TOTAL ${FS:-missing} (KA630 needs 0)"; log "630dd FAILSAFE_TOTAL 0"
DNL="$(/usr/bin/grep -oE 'dnmag LOW_CONTEXTS .*' "$DL" | tail -1)"; [ -n "$DNL" ] && log "690h: $DNL (report only)"
NJ="$W/$FT.ntclevel.json"; [ -s "$NJ" ] || notready "missing $NJ"
uv run python scripts/def_day_helpers.py past_cliff "$NJ" >> "$RL" 2>&1; PC=$?
[ "$PC" = 1 ] || notready "past_cliff rc $PC on $FT (1 = safe)"
RHO="$(uv run python scripts/def_day_helpers.py rho "$NJ" 2>/dev/null)"; log "[$N] CLIFF OK final $FT: $RHO"; raw "[$N] CLIFF OK final $FT: $RHO"
gate(){
  [ "$(uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(len(d.other_builders()))" 2>/dev/null || echo 99)" = 0 ] || return 1
  [ "$(ps -axo args= | grep -vE 'claude|grep' | grep -cE '(^|/| )vcc prep( |$)' | tr -d ' ')" = 0 ] || return 1
  [ "$(uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(int(d.free_gib()))" 2>/dev/null || echo 0)" -ge 26 ] || return 1
}
last=""; until gate; do [ "$last" = w ] || log "prep gate closed (a builder / operator fit, a vcc prep, or < 26 GiB): waiting"; last=w
  [ "$(now)" -ge "$CAP" ] && notready "prep gate not open by CAP"; sleep 60; done
[ -e "$PKG" ] && notready "$PKG exists (never overwrite a package)"
CTXS="$(python3 -c 'import json,sys; print(",".join(json.load(open(sys.argv[1]))["contexts"]))' "$B/manifest.json")"
t0=$(date +%s)
/usr/bin/time -l vcc prep "$H5" -g "$B/gene_names.csv" --perts "$B/pert_counts.csv" --contexts "$CTXS" -o "$PKG" > "$OUT/$N.prep.log" 2>&1
PRC=$?; log "prep rc=$PRC in $(( $(date +%s) - t0 )) s; $(/usr/bin/grep 'maximum resident' "$OUT/$N.prep.log" | tr -s ' '); package $(stat -f %z "$PKG" 2>/dev/null) bytes"
[ "$PRC" = 0 ] && [ -s "$PKG" ] || notready "vcc prep failed rc $PRC (see $OUT/$N.prep.log)"
cp "$META" "$PKGDIR/$N.h5ad.meta.json"
uv run python scripts/check_submission_novelty.py check "$PKG" >> "$RL" 2>&1; NRC=$?
log "novelty rc=$NRC"
[ "$NRC" = 0 ] || notready "novelty check refused"
SHA="$(shasum -a 256 "$PKG" | cut -d' ' -f1)"; log "package sha256: $SHA"
log "provenance: builder.pinned.py == the meta's builder sha ${PB:0:16}"
printf '%s READY %s %s\n' "$(date -u +%FT%TZ)" "$N" "$SHA" >> "$RDY"
log "READY $N $SHA"
P=()
for f in def_day.log postbuild.json identity.json census.json preflight.json pin-artifacts.json pin-lodo.json def_day.final.env \
         builder.pinned.py "${N}_final.h5ad.meta.json" timings.json whitener_union_pruned_def.json "$FT.ntclevel.json"; do
  [ -s "$W/$f" ] && P+=("$W/$f")
done
P+=("$PKGDIR/$N.h5ad.meta.json" "$RL" "$RDY" "$OUT/$N.prep.log" "$OUT/postbuild-ka.json" "$OUT/postbuild.json" "$OUT/audit428.txt")
git add -f "${P[@]}" && git commit -q -m "692: $N re-gated with 554's override (package $SHA): KA postbuild == def_day's, override changes only 554's row, zero-gene clean (regate_r0)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" -- "${P[@]}" \
  && log "committed the records: $(git log -1 --format=%h)" || log "WARNING commit of the records failed (commit by hand)"
