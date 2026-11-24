#!/usr/bin/env bash
# 663 (coordinator, 2026-10-04): build idall3, the UNGATED v3 identity hook, as 10-05 slot 2's fill. rule1005's s2 list is
# exhausted (498p2g and 460dr are both scored).
# Arm per 612/PREREGISTRATION-draft-IDALL3.md (FINAL, cf3f8004): base = the best scored entry at slot time = 498p2g
# (JWSghc0o7MN5D6o8Iadc, 0.2749190) + VCC_DE_MDE_SIGN=file:<604 v3 identity artifact, λ 0.5>:<sha>:0, no KEYGATE. It reads
# against Δ_496a (496a = 446bj2 + the plain CLEAN4 hook), which is 10-05's slot 1.
# Env = 498p2g's round-2 build (498-combined-bulk/chain498.sh): cap-6 dump, coarse MDE, zr -5, round-2 factors, tercile
# sharpening. Live builder 4daae237 (== HEAD), emitter 0f56209d. EXPECT exactly de_mde_sign + de_sign_emitter_sha256
# (498p2g used edb5b8dc). Never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/663-idall3-arm; LOG=$X/work/build_idall3.log; E=experiments/428-zero-gene-audit; BU=scripts/make_submission_pm_cd4.py
mkdir -p $X/work
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
term(){ log "$*"; exit 0; }
echo $$ > $X/work/build_idall3.pid
R="$REPO/submissions/498p2g-cap15-partial50-clean4q80-lodo6-ballast-sharpenp2g-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5.h5ad.meta.json"
N2="$REPO/experiments/495-lodo-cap6/work/dump_h1x2_cd4x15_lodo05hi6.npz"
ID="$REPO/experiments/604-identity-prior/work/ident604v3_flex_lam0.5_ABC.npz"
IDSHA=6cd84dc258cb609d08d93ec4c6c6ec7c7134ff8e3093eccde5ac15bbd2ac618f
[ -s "$R" ] && [ -s "$N2" ] && [ -s "$ID" ] || term "NOT READY inputs missing (498p2g meta, cap-6 dump or identity artifact)"
[ "$(shasum -a 256 "$ID" | cut -d' ' -f1)" = "$IDSHA" ] || term "NOT READY identity artifact sha mismatch"
DS2=$(shasum -a 256 "$N2" | cut -d' ' -f1)
[ "$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"]["bulk_load_sha256"])' "$R")" = "$DS2" ] || term "NOT READY 498p2g's meta does not carry the cap-6 dump"
F="$(uv run python -c 'import json,sys; s=json.load(open(sys.argv[1]))["config"]["bulk_ctx_scale"]; print(",".join(f"{k}={s[k]:.4f}" for k in "ABC"))' "$R")"
B="$(shasum -a 256 $BU | cut -d' ' -f1)"
[ "$B" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ] || term "NOT READY builder != HEAD"
EM=0f56209d57e8385abd08a7a15692ae3241b9373b436b4921546e546d571e5876
# run_arm_generic arg 4 is the artifact name (file-safe); the leaderboard -m name (claimed 17:49Z) goes to the arm script.
NAME=idall3-ungated-identity-hook
STEM=idall3-cap15-partial50-clean4q80-lodo6-ballast-sharpenp2g-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5-idhook
log "released: builder ${B:0:16} emitter ${EM:0:16}; base 498p2g factors $F; hook file ${IDSHA:0:16}"
printf '{"base":"498p2g","VCC_BULK_LOAD":"%s:%s","VCC_DE_MDE":"0.1,0.1,0.3,coarse","VCC_DE_ZERO_REPAIR":"-5","VCC_BULK_CTX_SCALE":"%s","VCC_DE_MDE_SIGN":"file:%s:%s:0","sharpen":"tercile:2.5,2.0,1.75"}\n' \
  "$N2" "$DS2" "$F" "$ID" "$IDSHA" > $X/work/env.json
( unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_DE_MAG VCC_RC_AUDIT_DUMP VCC_DE_MDE_SIGN VCC_DE_UPRN
  export VCC_BULK_LOAD="$N2:$DS2" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F" \
         VCC_DE_MDE_SIGN="file:$ID:$IDSHA:0"
  bash $E/run_arm_generic.sh idall3 tercile:2.5,2.0,1.75 "$R" "$NAME" "$STEM" "$B" "$EM" de_mde_sign,de_sign_emitter_sha256 ) \
  >> $X/work/build_idall3.nohup 2>&1
RL=$E/work/build_idall3/runidall3.log
log "idall3: $(/usr/bin/grep -E '^\[idall3\] (PASS|FAIL|CLIFF)|package sha|novelty rc|stopping' $RL 2>/dev/null | tr '\n' ' ' | cut -c1-600)"
if /usr/bin/grep -q '^\[idall3\] PASS' $RL && /usr/bin/grep -q '^\[idall3\] CLIFF OK' $RL && /usr/bin/grep -q 'novelty rc=0' $RL; then
  SHA="$(/usr/bin/grep 'package sha256:' $RL | tail -1 | awk '{print $NF}')"
  [ -s "submissions/$STEM.vcc" ] && [ "$(shasum -a 256 submissions/$STEM.vcc | cut -d' ' -f1)" = "$SHA" ] && term "READY idall3 $SHA"
  term "NOT READY package missing or sha mismatch"
fi
term "NOT READY gates did not pass (see $RL)"
