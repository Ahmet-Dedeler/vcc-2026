#!/usr/bin/env bash
# def_day.sh — the October D/E/F day in one gated script (2026-09-25, coordinator; generalises 448's run448.sh).
# 476 (2026-09-27): 485's review fixes merged onto 474 (481's patches). Fail closed on driver rc, dump rc, audit and prep;
# a bad secant step ships the nearest safe build instead of aborting; lvlC must sit inside the 1% band; past_cliff crashes
# are not "safe". Needs the driver with 485's dump-assert skip (def_pipeline.py): without it the driver exits 1 after
# every successful dump and this script's dump rc check stops the day.
#
#   bash scripts/def_day.sh <bundle> <work_dir> <name> [q=1.75] [sameline_w=auto] [sameline_flags=pt:noprog] [max_levels=4]
#   (sameline_w `auto`, builder >= 4feae27d: each matched context gets w from its bank's quality class, 2 strong / 0.5 weak,
#    written by def_sameline.py; a number forces one global w)
#   DEF_DAY_LODO=1 (460, default off): the stack dump is scaled by 415's LODO-V1 key re-fitted on the day's panel
#    (446bj2's key on A/B/C) instead of 410's; see stage 2b.
#   DEF_DAY_PROGRAMME_SOURCE=<bank> (461, default off; e.g. h1_2025_full = Arc's H1 from all public cells, 458d): every
#    context's block prior comes from that bank (def_sameline.py --programme-source). Needs sameline_flags WITHOUT noprog
#    (pass `pt`); the script refuses `noprog` with it, because the builder would then ignore the programme.
#   DEF_DAY_PRIOR_GATE=q<NN> (468, default off; e.g. q33): with DEF_DAY_PROGRAMME_SOURCE, each context's design-weakest
#    NN % of targets keep CLEAN4 instead of the programme (the builder's VCC_DE_PRIOR_GATE, builder >= 693fde15; 467 on H1).
#    Refused without DEF_DAY_PROGRAMME_SOURCE. Exported for the lvl0 and re-level builds only (never for the stack dump).
#   DEF_DAY_LODO_HI_PCT=<P> (497, default off; needs DEF_DAY_LODO=1): the LODO key's cap is the panel's own P-th
#    percentile of carried key / median (def_lodo_key.py --hi-pct; 543/546), pinned as lodo_lo0.5_pct_def.
#    P 96.64697522769474 is exactly cap 6 on A/B/C (495's artifact); the rounded 96.647 is not.
#   DEF_DAY_LODO_HI=<h> (538, default off; needs DEF_DAY_LODO=1; refused with DEF_DAY_LODO_HI_PCT):
#    a fixed cap, passed to def_lodo_key.py as --hi. 1 < h <= 20, a plain decimal. Unset keeps cap 3,
#    or the percentile form. The pinned key is lodo_lo0.5_hi<h>_def and its json records hi.
#   DEF_DAY_CENSUS_MAX_ZERO=<n> (538, default off; integer 6..30; a fresh work dir only): raises the
#    census stop for zero-carrier targets from 5 to n. Logged, and recorded in def_day.final.env.
#   DEF_DAY_SHARPEN_PT=<tercile|stercile|tercilec>:<q>,<q>,<q> (497, default off; e.g. stercile:2.5,2.0,1.75 = 479's S_P2): the
#    level builds sharpen each target's bulk with its tercile's q (VCC_BULK_SHARPEN_Q's per-target form, builder >=
#    d09a986c) instead of the scalar q. The stack dump is unchanged (q is a render knob).
#    598: tercilec = tercile with 493h's carried-only cut points (the carried flags of the key that scaled the dump;
#    non-carried targets in MID). Builder >= 598's patch (dump_key_carried).
#   DEF_DAY_DUMP_HOOK=<committed script> (690g: dump hook; default off): after the stack dump, `bash <script> <bundle> <work>
#    <dump.npz> <identity.json> <work>/dumphook` writes a replacement dump (same name and format, env unchanged) and maybe
#    bulk_ctx.npz (672's VCC_BULK_LOAD_CTX); scripts/def_dumphook_check.py validates both before any level build.
#   DEF_DAY_DE_MAG=<variant> (502, default off; e.g. dn_k250_zs): the zs exact-U pass (VCC_DE_MAG) on the level builds.
#   DEF_DAY_DE_UPRN=<s|b>_k<K*100>[_q<n>r<n>p<n>] (630dd: rank-neutral UP pass-through; default off; e.g. b_k100_q20r3p10):
#    630's rank-neutral UP pass after zs (VCC_DE_UPRN). Needs DEF_DAY_DE_MAG. Builder >= 630 (rank_up module 7e6b3e52).
#   DEF_DAY_DE_DNMAG=d_m<m*100>[_q<n>] (690a: DOWN magnitude pass-through; default off; e.g. d_m200): 631's DOWN-magnitude
#    pass after zs and 630's UP pass (VCC_DE_DNMAG). Needs DEF_DAY_DE_MAG. Builder >= 631 + 693p (dn_mag2 module b9d6cfcf).
#   DEF_DAY_MDE_SIGN=clean4:<tau> (502, default off; e.g. clean4:0.05): 496's uncalled-gene sign hook (VCC_DE_MDE_SIGN).
#    517: or programme:<tau>[:<tau_c4>] (the hook reads the programme of DEF_DAY_PROGRAMME_SOURCE's artifact and follows
#    DEF_DAY_PRIOR_GATE), or hookbank:<bank>:<tau>[:<tau_c4>[:le<N>]] (the hook ALONE from <bank>'s programme, CLEAN4
#    block: a second def_sameline.py artifact, exported as VCC_DE_MDE_SIGN=file:<npz>:<sha>:...). Builder >= 517's patch.
#    DEF_DAY_HOOK_ESTIMATOR=meanz_ne (517, default pdown): the hook artifact's estimator (def_sameline.py >= 517's patch).
#    DEF_DAY_MDE_SIGN_KEYGATE=<q> (612dd: key gate pass-through; default off; e.g. 0.5): 612's key gate on the hook's sign
#    (VCC_DE_MDE_SIGN_KEYGATE): the hook's sign only on each target's weak-key genes. Needs DEF_DAY_MDE_SIGN. Builder >= 612.
#   DEF_DAY_RHO_SCALE=<s> (530, default off; 0.9 <= s <= 1.1): the re-level loop targets rho* x s.
#   DEF_DAY_CD4_STACK=x1|x1.5 (530, default x1): x1.5 exports VCC_CD4_TARGET_G inside the dump subshell only (457).
#   DEF_DAY_RC_CALLSET_CPM=<n> (530, default off; 100..5000): VCC_RC_CALLSET_CPM on the level builds (523; unset = 2200).
#   DEF_DAY_LODO_STRENGTH=<theta> (530, default off): 522's pct x exp(theta z_src) rider. Needs LODO and the percentile cap.
#   DEF_DAY_CALLSET_BLOCK=prop[:<swap>[:<steer>]] (528, default off): the DE call set's CPM block chosen per context by a
#    call-propensity key from the day's OWN controls (log10 control CPM + 22 x excess dropout, median-rescaled to H1's
#    dropout scale; callblock528.py), exported as VCC_RC_CALLSET_BLOCK for the level builds. Builder/emitter >= 528's.
#   DEF_DAY_KEEP_PINS=1 (default off): skip the exit-trap restore so `git diff` on the builder is the pin to commit.
#    The default restores, which is what the 452 rehearsal checks. Pins are re-applied by re-running the pin stages.
#
# Stages (each logged to <work>/def_day.log; a failed driver, audit, or prep stops the script):
#   1. the driver's light and pin stages with a bare env: preflight, identity, cd4-extract, census, the day's whitener,
#      the day's 410 key, repin-frozen, pin-artifacts, repin-panel-support (437/439);
#   2. the same-line artifact from the identity call (451: `scripts/def_sameline.py`). If any context is matched,
#      VCC_SAMELINE=<npz>:<sha>:<w>:<flags> is exported for every build below. Unmatched: nothing changes;
#   3. the stack dump, H1 x2 / CD4 x1 (438's October default: VCC_H1_UNION_G=0.56, CD4 at DEF_ENV's default);
#   4. lvl0, the un-levelled LOAD build at q, with SHIPPED_ENV's coarse MDE and zero repair (driver defaults);
#   5. the re-level loop (448): ntclevel417 (sd 0.01) on each build, then `relevel_iter.py` over every build so far,
#      until CONVERGED (every context within 1 % of rho* and none above 1) or max_levels levelled builds. A build with any
#      rho > 1 is never the final one. A proposal outside [0.7, 1.4], or a degenerate secant (relevel exit 2), ships the
#      nearest already-built safe artifact instead of stopping;
#   6. audits on the final build: ntclevel, zero-gene audit (stop unless the zero count is 0), and a `vcc prep --dry-run`
#      through the driver (stop unless that stage exits 0). lfc < -5 is logged and not gated (452 shipped with 7).
# It stops BEFORE the real `vcc prep` and `vcc submit`: those stay manual (DEF-RUNBOOK §2's tail, AGENTS §6.2, §14).
# The final artifact is <work>/<name>_final.h5ad (a hard link to the converged build). Its VCC_BULK_CTX_SCALE string is
# in <work>/def_day.final.env. The link is made only after the meta's q and scale match that string.
#
#   DEF_DAY_CD4_THRESHOLDS=relaxed (481 F1, default `default`): the CD4 extraction runs without 047's gates (275 per
#    condition, 290 union). Arc's other 300-target panel (VCC 2025 H1) covers 262/265/262, union 267, so the gates fire
#    on an Arc-style panel; even A/B/C cleared the union gate by one target. The census's own floors (union >= 120,
#    complete >= 100) still gate. Check the mapping first (244's mapping244.py), as DEF-RUNBOOK section 0.0 says.
#   DEF_DAY_INCUMBENT_META=<meta.json> (481 F5): the shipped A/B/C arm the final artifact's config is diffed against by
#    def_postbuild.py before the builder is restored (default: 446b's meta).
#
# Rehearse on a synthetic bundle first (e.g. experiments/238-essential-free-panel/work/ef-bundle-r).
# A rerun needs a fresh work dir: this script deletes a tag's h5ad before rebuilding it, and refuses to
# treat a pre-existing file as the build it just asked for.
set -uo pipefail
# 538: the repo is the directory that contains this script. On the production Mac this is
# ~/Code/virtual-cell-challenge-2026; in a worktree it is the worktree.
REPO="$(CDPATH= cd "$(dirname "$0")/.." && pwd)"
[ -f "$REPO/scripts/def_day.sh" ] || { echo "STOP: derived REPO=$REPO has no scripts/def_day.sh" >&2; exit 1; }
cd "$REPO" || exit 1
B="${1:?bundle}"; W="${2:?work dir}"; N="${3:?name}"; Q="${4:-1.75}"; SLW="${5:-auto}"; SLF="${6-pt:noprog}"; MAXL="${7:-4}"
mkdir -p "$W"; LOG="$W/def_day.log"
CTXS="$(python3 -c "import json, sys; print(','.join(json.load(open(sys.argv[1]))['contexts']))" "$B/manifest.json")" || exit 1
UNIT="$(python3 -c "import sys; print(','.join(c + '=1' for c in sys.argv[1].split(',')))" "$CTXS")"
CD4T="${DEF_DAY_CD4_THRESHOLDS:-default}"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
die(){ log "STOP: $*"; exit 1; }
# 538: the census override is a one-shot. Checked before log "start" writes def_day.log, so a fresh
# work dir is still empty here (mkdir -p above does not count). Unset: this whole block is skipped.
if [ -n "${DEF_DAY_CENSUS_MAX_ZERO:-}" ]; then
  python3 -c 'import sys; n=sys.argv[1]; sys.exit(0 if n.isdigit() and not (len(n)>1 and n[0]=="0") and 6<=int(n)<=30 else 1)' "$DEF_DAY_CENSUS_MAX_ZERO" 2>/dev/null \
    || die "DEF_DAY_CENSUS_MAX_ZERO must be an integer in [6, 30] with no leading zero (got $DEF_DAY_CENSUS_MAX_ZERO)"
  if [ -d "$W" ] && [ -n "$(ls -A "$W" 2>/dev/null)" ]; then
    die "DEF_DAY_CENSUS_MAX_ZERO=$DEF_DAY_CENSUS_MAX_ZERO is a one-shot override and needs a fresh work dir ($W is not empty)"
  fi
fi
for v in $(env | grep -o '^VCC_[A-Z0-9_]*'); do unset "$v"; done
BSHA0="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
log "start: bundle $B work $W name $N q $Q sameline w $SLW flags $SLF; builder ${BSHA0:0:16} driver $(shasum -a 256 scripts/def_pipeline.py | cut -c1-16)"
[ "${DEF_DAY_ALLOW_WIDE_LIBS:-}" = 1 ] && log "554: DEF_DAY_ALLOW_WIDE_LIBS=1: preflight's control-library spread floor is overridden"
[ "${DEF_DAY_ALLOW_LOW_DOSE:-}" = 1 ] && log "554: DEF_DAY_ALLOW_LOW_DOSE=1: postbuild's steering-dose floors are overridden"
[ "${DEF_DAY_ALLOW_DEPTH:-}" = 1 ] && log "641: the depth row: DEF_DAY_ALLOW_DEPTH=1: preflight's control-library median check is overridden"
[ -n "${DEF_DAY_ALLOW_NONCCLE_REF:-}" ] && log "642: the HEK293 contingency: DEF_DAY_ALLOW_NONCCLE_REF=$DEF_DAY_ALLOW_NONCCLE_REF (def_identity's non-CCLE second stage); decision: ${DEF_DAY_NONCCLE_DECISION:-<none recorded>}"
[ -n "${DEF_DAY_SAMELINE_SKIP:-}" ] && log "644: DEF_DAY_SAMELINE_SKIP=$DEF_DAY_SAMELINE_SKIP: these contexts take the unmatched same-line path"
# A reboot skips this script's EXIT trap, so the builder can still be the LODO pin. A fresh work
# dir would record that pin as BSHA0 and the trap would then "restore" it (526). Refuse.
_HEADB="$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -d' ' -f1)"
if [ "$BSHA0" != "$_HEADB" ]; then
  die "builder on disk ${BSHA0:0:16} is not HEAD ${_HEADB:0:16}. Restore before a fresh work dir: uv run python scripts/def_pin_panel_artifacts.py --restore --record <previous-work>/pin-artifacts.json"
fi
grep -q "VCC_SAMELINE" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_SAMELINE (apply 451's patch first)"
DRV=(uv run python scripts/def_pipeline.py --bundle "$B" --work "$W" --whitener-id union_pruned_def --emitter rc_callset
     --frozen-inputs "$W/frozen_inputs.json" --min-free-gib 23)
case "$CD4T" in
  default) ;;
  relaxed) log "CD4 extract thresholds relaxed (DEF_DAY_CD4_THRESHOLDS=relaxed): check the mapping first (DEF-RUNBOOK 0.0)" ;;
  *) die "DEF_DAY_CD4_THRESHOLDS must be relaxed or unset (got $CD4T)" ;;
esac
# 497: the October forms of 495 (percentile cap) and 479 (per-target sharpening); both default off
if [ -n "${DEF_DAY_LODO_HI_PCT:-}" ]; then
  [ "${DEF_DAY_LODO:-0}" = "1" ] || die "DEF_DAY_LODO_HI_PCT=$DEF_DAY_LODO_HI_PCT needs DEF_DAY_LODO=1"
  python3 -c "import sys; p = float(sys.argv[1]); sys.exit(0 if 50.0 < p < 100.0 else 1)" "$DEF_DAY_LODO_HI_PCT" 2>/dev/null \
    || die "DEF_DAY_LODO_HI_PCT must be a percentile in (50, 100) (got $DEF_DAY_LODO_HI_PCT)"
  grep -q -- "--hi-pct" scripts/def_lodo_key.py || die "def_lodo_key.py has no --hi-pct (497)"
fi
# 538: one cap form. Unset keeps today's --hi default (3) or DEF_DAY_LODO_HI_PCT.
if [ -n "${DEF_DAY_LODO_HI:-}" ]; then
  [ "${DEF_DAY_LODO:-0}" = "1" ] || die "DEF_DAY_LODO_HI=$DEF_DAY_LODO_HI needs DEF_DAY_LODO=1"
  [ -z "${DEF_DAY_LODO_HI_PCT:-}" ] || die "DEF_DAY_LODO_HI and DEF_DAY_LODO_HI_PCT are two cap forms; set only one (got HI=$DEF_DAY_LODO_HI PCT=$DEF_DAY_LODO_HI_PCT)"
  python3 -c 'import re,sys
t=sys.argv[1]
if re.fullmatch(r"[1-9][0-9]*(?:\.[0-9]+)?", t) is None:
    sys.exit(1)
h=float(t)
sys.exit(0 if 1.0 < h <= 20.0 else 1)' "$DEF_DAY_LODO_HI" 2>/dev/null \
    || die "DEF_DAY_LODO_HI must satisfy 1 < h <= 20, written as a plain decimal (got $DEF_DAY_LODO_HI)"
  grep -q 'add_argument("--hi"' scripts/def_lodo_key.py || die "def_lodo_key.py has no --hi"
fi
# 690g: dump hook -- DEF_DAY_DUMP_HOOK=<script> (validated here; it runs right after the stack dump)
if [ -n "${DEF_DAY_DUMP_HOOK:-}" ]; then
  [ -f "$DEF_DAY_DUMP_HOOK" ] || die "DEF_DAY_DUMP_HOOK=$DEF_DAY_DUMP_HOOK is not a file"
  git ls-files --error-unmatch "$DEF_DAY_DUMP_HOOK" >/dev/null 2>&1 || die "DEF_DAY_DUMP_HOOK must be a committed script (got $DEF_DAY_DUMP_HOOK)"
  [ -z "$(git status --porcelain -- "$DEF_DAY_DUMP_HOOK")" ] || die "DEF_DAY_DUMP_HOOK $DEF_DAY_DUMP_HOOK has uncommitted changes"
  [ -f scripts/def_dumphook_check.py ] || die "scripts/def_dumphook_check.py is missing (690g)"
fi
# 502: the DE-channel riders (validated here, exported for the level builds after the dump)
if [ -n "${DEF_DAY_DE_MAG:-}" ]; then
  case "$DEF_DAY_DE_MAG" in *[!a-z0-9_]*) die "DEF_DAY_DE_MAG must be a variant name like dn_k250_zs (got $DEF_DAY_DE_MAG)" ;; esac
  grep -q "VCC_DE_MAG" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_DE_MAG"
fi
# 630dd: rank-neutral UP pass-through -- DEF_DAY_DE_UPRN=<variant> (validated here, exported with the other DE riders)
if [ -n "${DEF_DAY_DE_UPRN:-}" ]; then
  [ -n "${DEF_DAY_DE_MAG:-}" ] || die "DEF_DAY_DE_UPRN runs after the zs pass: it needs DEF_DAY_DE_MAG"
  python3 -c 'import re, sys; m = re.fullmatch(r"(s|b)_k(\d{2,3})(?:_q([1-9]\d*)r([1-9]\d*)p([1-9]\d*))?", sys.argv[1]); sys.exit(0 if m and 50 <= int(m.group(2)) <= 300 else 1)' "$DEF_DAY_DE_UPRN" 2>/dev/null \
    || die "DEF_DAY_DE_UPRN must be <s|b>_k<K*100>[_q<n>r<n>p<n>] with K in [0.5, 3], e.g. b_k100_q20r3p10 (got $DEF_DAY_DE_UPRN)"
  grep -q "630: rank-neutral UP" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_DE_UPRN (630)"
fi
# 690a: DOWN magnitude pass-through -- DEF_DAY_DE_DNMAG=<variant> (validated here, exported with the other DE riders)
if [ -n "${DEF_DAY_DE_DNMAG:-}" ]; then
  [ -n "${DEF_DAY_DE_MAG:-}" ] || die "DEF_DAY_DE_DNMAG runs after the zs pass: it needs DEF_DAY_DE_MAG"
  python3 -c 'import re, sys; m = re.fullmatch(r"d_m(\d{3})(?:_q([1-9]\d*))?", sys.argv[1]); sys.exit(0 if m and 125 <= int(m.group(1)) <= 400 else 1)' "$DEF_DAY_DE_DNMAG" 2>/dev/null \
    || die "DEF_DAY_DE_DNMAG must be d_m<m*100>[_q<n>] with m in [1.25, 4], e.g. d_m200 (got $DEF_DAY_DE_DNMAG)"
  grep -q "VCC_DE_DNMAG" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_DE_DNMAG (631)"
fi
if [ -n "${DEF_DAY_MDE_SIGN:-}" ]; then
  grep -q "VCC_DE_MDE_SIGN" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_DE_MDE_SIGN (496p)"
  case "$DEF_DAY_MDE_SIGN" in
    clean4:*) ;;
    programme:*) [ -n "${DEF_DAY_PROGRAMME_SOURCE:-}" ] || die "DEF_DAY_MDE_SIGN=$DEF_DAY_MDE_SIGN reads the programme artifact: it needs DEF_DAY_PROGRAMME_SOURCE"
                 grep -q "programme:<tau>" scripts/make_submission_pm_cd4.py || die "the builder has no programme hook form (517)" ;;
    hookbank:*) python3 -c 'import re,sys
p=sys.argv[1].split(":")
ok=len(p) in (3,4,5) and bool(p[1])
try:
    t=float(p[2]); c=float(p[3]) if len(p)>=4 else 0.05
    ok = ok and 0.0 <= t < 0.5 and 0.0 <= c < 0.5
    if len(p)==5: ok = ok and re.fullmatch(r"le[0-9]{1,2}", p[4]) is not None
except Exception:
    ok=False
sys.exit(0 if ok else 1)' "$DEF_DAY_MDE_SIGN" 2>/dev/null \
                 || die "DEF_DAY_MDE_SIGN hookbank form must be hookbank:<bank>:<tau>[:<tau_c4>[:le<N>]] (got $DEF_DAY_MDE_SIGN)"
                 grep -q "file:<npz>" scripts/make_submission_pm_cd4.py || die "the builder has no file hook form (517)" ;;
    *) die "DEF_DAY_MDE_SIGN must be clean4:<tau>, programme:<tau>[:<tau_c4>] or hookbank:<bank>:<tau>[:<tau_c4>[:le<N>]] (got $DEF_DAY_MDE_SIGN)" ;;
  esac
fi
# 612dd: key gate pass-through -- DEF_DAY_MDE_SIGN_KEYGATE=<q> (validated here, exported with the other DE riders)
if [ -n "${DEF_DAY_MDE_SIGN_KEYGATE:-}" ]; then
  [ -n "${DEF_DAY_MDE_SIGN:-}" ] || die "DEF_DAY_MDE_SIGN_KEYGATE gates the hook's sign: it needs DEF_DAY_MDE_SIGN"
  python3 -c 'import sys; q=float(sys.argv[1]); sys.exit(0 if 0.0 < q < 1.0 else 1)' "$DEF_DAY_MDE_SIGN_KEYGATE" 2>/dev/null \
    || die "DEF_DAY_MDE_SIGN_KEYGATE must be a quantile 0 < q < 1 (got $DEF_DAY_MDE_SIGN_KEYGATE)"
  grep -q "VCC_DE_MDE_SIGN_KEYGATE" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_DE_MDE_SIGN_KEYGATE (612)"
fi
# 530: October pass-throughs. Unset = no export and no extra argument (today's command).
if [ -n "${DEF_DAY_RHO_SCALE:-}" ]; then
  python3 -c 'import sys; s=float(sys.argv[1]); sys.exit(0 if 0.9<=s<=1.1 else 1)' "$DEF_DAY_RHO_SCALE" 2>/dev/null \
    || die "DEF_DAY_RHO_SCALE must be a scale in [0.9, 1.1] (got $DEF_DAY_RHO_SCALE)"
fi
if [ -n "${DEF_DAY_RC_CALLSET_CPM:-}" ]; then
  python3 -c 'import sys; n=sys.argv[1]; sys.exit(0 if n.isdigit() and 100<=int(n)<=5000 else 1)' "$DEF_DAY_RC_CALLSET_CPM" 2>/dev/null \
    || die "DEF_DAY_RC_CALLSET_CPM must be an integer in [100, 5000] (got $DEF_DAY_RC_CALLSET_CPM)"
  grep -q "VCC_RC_CALLSET_CPM" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_RC_CALLSET_CPM (523)"
fi
case "${DEF_DAY_CD4_STACK:-x1}" in
  x1|x1.5) ;;
  *) die "DEF_DAY_CD4_STACK must be x1 or x1.5 (got ${DEF_DAY_CD4_STACK})" ;;
esac
if [ -n "${DEF_DAY_LODO_STRENGTH:-}" ]; then
  [ "${DEF_DAY_LODO:-0}" = "1" ] || die "DEF_DAY_LODO_STRENGTH=$DEF_DAY_LODO_STRENGTH needs DEF_DAY_LODO=1"
  [ -n "${DEF_DAY_LODO_HI_PCT:-}" ] || die "DEF_DAY_LODO_STRENGTH=$DEF_DAY_LODO_STRENGTH needs the percentile cap (DEF_DAY_LODO_HI_PCT)"
  python3 -c 'import re,sys
t=sys.argv[1]
if re.fullmatch(r"-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?", t) is None:
    sys.exit(1)
v=float(t)
sys.exit(0 if v != 0.0 and abs(v) <= 1.0 else 1)' "$DEF_DAY_LODO_STRENGTH" 2>/dev/null \
    || die "DEF_DAY_LODO_STRENGTH must be a non-zero theta in [-1, 1], written as a plain decimal (got $DEF_DAY_LODO_STRENGTH)"
fi
if [ -n "${DEF_DAY_HOOK_ESTIMATOR:-}" ]; then
  [ "$DEF_DAY_HOOK_ESTIMATOR" = "meanz_ne" ] || die "DEF_DAY_HOOK_ESTIMATOR must be meanz_ne (got $DEF_DAY_HOOK_ESTIMATOR)"
  case "${DEF_DAY_MDE_SIGN:-}" in
    hookbank:*) ;;
    *) die "DEF_DAY_HOOK_ESTIMATOR=$DEF_DAY_HOOK_ESTIMATOR needs DEF_DAY_MDE_SIGN=hookbank:<bank>:<tau>... (got ${DEF_DAY_MDE_SIGN:-})" ;;
  esac
  grep -q -- "--programme-estimator" scripts/def_sameline.py || die "def_sameline.py has no --programme-estimator (517)"
fi
# 528: the call-propensity CPM block (DE-only; validated here, computed after preflight, exported for the level builds)
if [ -n "${DEF_DAY_CALLSET_BLOCK:-}" ]; then
  case "$DEF_DAY_CALLSET_BLOCK" in
    prop) ;;
    prop:*) python3 -c "import sys; p = sys.argv[1].split(':'); v = float(p[0]); sys.exit(0 if (0.0 < v <= 1.0 and len(p) <= 2 and (len(p) == 1 or p[1] in ('kept', 'key', 'cpm'))) else 1)" "${DEF_DAY_CALLSET_BLOCK#prop:}" 2>/dev/null \
              || die "DEF_DAY_CALLSET_BLOCK=prop:<swap>[:<steer>] needs 0 < swap <= 1 and steer kept|key|cpm (got $DEF_DAY_CALLSET_BLOCK)" ;;
    *) die "DEF_DAY_CALLSET_BLOCK must be prop or prop:<swap>[:<steer>] (got $DEF_DAY_CALLSET_BLOCK)" ;;
  esac
  grep -q "VCC_RC_CALLSET_BLOCK" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_RC_CALLSET_BLOCK (528)"
  grep -q "block_steer=" src/vcc2026/rc_emit_k.py || die "the emitter has no block_key / block_steer (528)"
  [ -f experiments/528-callprop-block/callblock528.py ] || die "experiments/528-callprop-block/callblock528.py is missing"
fi
SQ="$Q"
if [ -n "${DEF_DAY_SHARPEN_PT:-}" ]; then
  case "$DEF_DAY_SHARPEN_PT" in tercile:*,*,*|stercile:*,*,*|tercilec:*,*,*) ;; *) die "DEF_DAY_SHARPEN_PT must be tercile:<q>,<q>,<q>, stercile:<q>,<q>,<q> or tercilec:<q>,<q>,<q> (got $DEF_DAY_SHARPEN_PT)" ;; esac
  grep -q "def sharpen_q_by_tercile" scripts/make_submission_pm_cd4.py || die "the builder has no per-target sharpening (479)"
  # 598: tercilec reads the carried flags of the dump's own key artifact (the builder's dump_key_carried)
  case "$DEF_DAY_SHARPEN_PT" in tercilec:*) grep -q "def dump_key_carried" scripts/make_submission_pm_cd4.py || die "the builder has no tercilec (598)" ;; esac
  SQ="$DEF_DAY_SHARPEN_PT"
  log "per-target sharpening for the level builds: VCC_BULK_SHARPEN_Q=$SQ (scalar q $Q not used)"
fi
restore(){
  local rc=$?
  # 690c: pinned builder kept -- the builder the day's builds ran on (after every pin; the final meta's
  # hashes.implementation sha), saved before the restore so the round can commit it next to its records
  if [ -f "$W/pin-artifacts.json" ] && cp scripts/make_submission_pm_cd4.py "$W/builder.pinned.py" 2>/dev/null; then
    log "690c: pinned builder saved as $W/builder.pinned.py ($(shasum -a 256 "$W/builder.pinned.py" | cut -c1-16))"
  fi
  if [ "${DEF_DAY_KEEP_PINS:-0}" = "1" ]; then
    log "DEF_DAY_KEEP_PINS=1: builder left pinned $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -c1-16) (pre-pin ${BSHA0:0:16}); commit the pins, then restore from $W/pin-artifacts.json"
    exit "$rc"
  fi
  if [ -f "$W/pin-artifacts.json" ]; then
    if ! uv run python scripts/def_pin_panel_artifacts.py --restore --record "$W/pin-artifacts.json" >> "$LOG" 2>&1; then
      log "STOP: builder restore command failed"
      rc=1
    fi
  fi
  local now
  now="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
  log "builder after restore: ${now:0:16} (want ${BSHA0:0:16})"
  if [ "$now" != "$BSHA0" ]; then
    log "STOP: restored builder sha does not match the pre-pin sha"
    rc=1
  fi
  exit "$rc"
}
trap restore EXIT

# 1. setup
"${DRV[@]}" --name "$N-setup" --stages preflight,identity,cd4-extract,census,whitener,key,repin-frozen,pin-artifacts,repin-panel-support \
    --cd4-thresholds "$CD4T" ${DEF_DAY_CENSUS_MAX_ZERO:+--census-max-zero "$DEF_DAY_CENSUS_MAX_ZERO"} > "$W/driver-setup.log" 2>&1
if ! grep -q "^--- \[repin-panel-support\].*rc=0" "$W/driver-setup.log"; then
  grep -qE "condition coverage below|panel targets map across any condition" "$W/logs/cd4-extract.log" 2>/dev/null \
    && die "CD4 extraction gate fired ($(grep -E 'panel coverage|only [0-9]+ panel targets|coverage below' "$W/logs/cd4-extract.log" | head -2 | tr '\n' ' ')). Check the mapping (experiments/244-october-rehearsal/mapping244.py), then rerun with DEF_DAY_CD4_THRESHOLDS=relaxed in a fresh work dir"
  die "setup did not complete (see driver-setup.log)"
fi
log "setup: $(grep -E '^--- \[' "$W/driver-setup.log" | awk '{print $2, $3}' | tr '\n' ' ')"
if [ -n "${DEF_DAY_CENSUS_MAX_ZERO:-}" ]; then
  ZC="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["consensus"]["targets_with_0_donors"])' "$W/census.json")" \
    || die "census override set but $W/census.json has no zero-carrier count"
  log "census override: max zero-carrier $DEF_DAY_CENSUS_MAX_ZERO (count $ZC)"
fi
[ -s "$W/identity.json" ] || die "no identity.json"
log "identity: $(python3 -c "import json; d=json.load(open('$W/identity.json')); print(d['call'], 'unanimous', d['call_unanimous_across_methods'])")"
# Context labels come from the preflight the setup stage just wrote. Hardcoding D,E,F ships a
# mismatch the moment the manifest uses other names: lvl0 has no scale and succeeds, then ntclevel
# looks for context_D.h5ad. isidentifier() rejects anything that is not a bare name.
CTXP="$(python3 -c 'import json,sys
cs=json.load(open(sys.argv[1]))["contexts"]
assert isinstance(cs, list) and len(cs)==3 and all(isinstance(c, str) and c.isidentifier() for c in cs), cs
print(",".join(cs))' "$W/preflight.json")" || die "could not read three context names from $W/preflight.json"
[ "$CTXP" = "$CTXS" ] || die "preflight contexts $CTXP differ from the manifest's $CTXS"
log "contexts $CTXS; unit factors $UNIT"

# 2. same-line
SL=""
if [ -n "${DEF_DAY_PROGRAMME_SOURCE:-}" ]; then
  case ":$SLF:" in *:noprog:*) die "DEF_DAY_PROGRAMME_SOURCE=$DEF_DAY_PROGRAMME_SOURCE needs sameline_flags without noprog (got $SLF)";; esac
  log "programme-source $DEF_DAY_PROGRAMME_SOURCE: every context's block prior from it (flags $SLF)"
fi
if [ -n "${DEF_DAY_PRIOR_GATE:-}" ]; then
  [ -n "${DEF_DAY_PROGRAMME_SOURCE:-}" ] || die "DEF_DAY_PRIOR_GATE=$DEF_DAY_PRIOR_GATE gates a programme: it needs DEF_DAY_PROGRAMME_SOURCE"
  grep -q "VCC_DE_PRIOR_GATE" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_DE_PRIOR_GATE (apply 468's patch first)"
  log "prior gate $DEF_DAY_PRIOR_GATE: each context's design-weakest targets keep CLEAN4"
fi
# DEF_DAY_FORCE="D=<CCLE sample>,..." forces the same-line call (rehearsals only; never on the day)
# 644: DEF_DAY_SAMELINE_SKIP=<ctx[,ctx]>: contexts that take the unmatched same-line path whatever their identity call
# (whitespace is ignored: "E, F" == "E,F"; 607)
SLSKIP="$(printf '%s' "${DEF_DAY_SAMELINE_SKIP:-}" | tr -d '[:space:]')"
if [ -n "$SLSKIP" ]; then
  for _c in $(printf '%s' "$SLSKIP" | tr ',' ' '); do
    case ",$CTXS," in *",$_c,"*) ;; *) die "644: DEF_DAY_SAMELINE_SKIP names $_c, which is not a manifest context ($CTXS)" ;; esac
  done
  log "644: same-line skipped for $SLSKIP (def_sameline --skip-contexts): treated as unmatched"
fi
uv run python scripts/def_sameline.py --bundle "$B" --identity "$W/identity.json" --out "$W/sameline.npz" \
    ${DEF_DAY_FORCE:+--force "$DEF_DAY_FORCE"} ${DEF_DAY_PROGRAMME_SOURCE:+--programme-source "$DEF_DAY_PROGRAMME_SOURCE"} \
    ${SLSKIP:+--skip-contexts "$SLSKIP"} \
    >> "$LOG" 2>&1 || die "def_sameline.py failed"
[ -n "${DEF_DAY_FORCE:-}" ] && log "REHEARSAL: same-line call forced to $DEF_DAY_FORCE"
NCTX="$(uv run python scripts/def_day_helpers.py nctx "$W/sameline.npz")" || die "nctx could not read $W/sameline.npz"
case "$NCTX" in
  ''|*[!0-9]*) die "nctx printed '$NCTX'" ;;
esac
if [ "$NCTX" -gt 0 ]; then
  SL="$W/sameline.npz:$(shasum -a 256 "$W/sameline.npz" | cut -d' ' -f1):$SLW${SLF:+:$SLF}"
  log "same-line mode ON for $NCTX context(s): VCC_SAMELINE=${SL##*/}"
else
  if [ -n "${DEF_DAY_PROGRAMME_SOURCE:-}" ]; then
    die "DEF_DAY_PROGRAMME_SOURCE=$DEF_DAY_PROGRAMME_SOURCE but the same-line artifact switches on 0 contexts"
  fi
  log "same-line mode off (no unanimous registered line)"
fi

# 2a. 517: DEF_DAY_MDE_SIGN=hookbank:<bank>:... -- the hook's own programme artifact (the block prior is untouched)
MDS="${DEF_DAY_MDE_SIGN:-}"
case "$MDS" in hookbank:*)
  HB="${MDS#hookbank:}"; HBANK="${HB%%:*}"; HREST="${HB#*:}"
  [ -n "$HBANK" ] && [ "$HREST" != "$HB" ] || die "DEF_DAY_MDE_SIGN=$MDS: expected hookbank:<bank>:<tau>[...]"
  uv run python scripts/def_sameline.py --bundle "$B" --identity "$W/identity.json" --out "$W/hookprog.npz" \
      ${DEF_DAY_FORCE:+--force "$DEF_DAY_FORCE"} --programme-source "$HBANK" \
      ${DEF_DAY_HOOK_ESTIMATOR:+--programme-estimator "$DEF_DAY_HOOK_ESTIMATOR"} >> "$LOG" 2>&1 \
      || die "def_sameline.py --programme-source $HBANK (the hook artifact) failed"
  MDS="file:$W/hookprog.npz:$(shasum -a 256 "$W/hookprog.npz" | cut -d' ' -f1):$HREST"
  log "hook programme artifact from $HBANK (estimator ${DEF_DAY_HOOK_ESTIMATOR:-pdown}): VCC_DE_MDE_SIGN=file:${MDS#file:$W/}" ;;
esac

# 2c. 528: DEF_DAY_CALLSET_BLOCK=prop[:<swap>] -- the day's block keys from D/E/F's own controls (nothing fitted today)
CSB=""
if [ -n "${DEF_DAY_CALLSET_BLOCK:-}" ]; then
  uv run python experiments/528-callprop-block/callblock528.py --bundle "$B" --contexts "$CTXS" --out "$W/callblock_def.npz" \
      >> "$LOG" 2>&1 || die "callblock528.py (the call-propensity block keys) failed"
  CSB="prop:$W/callblock_def.npz:$(shasum -a 256 "$W/callblock_def.npz" | cut -d' ' -f1)"
  case "$DEF_DAY_CALLSET_BLOCK" in prop:*) CSB="$CSB:${DEF_DAY_CALLSET_BLOCK#prop:}" ;; esac
  log "call-propensity block keys for $CTXS: VCC_RC_CALLSET_BLOCK=prop:${CSB#prop:$W/}"
fi

# 2b. DEF_DAY_LODO=1 (460; default off): the stack dump is scaled by 415's LODO-V1 key at clip (0.5, 3), re-fitted on
#     this panel, the day's operator and CD4 banks (`scripts/def_lodo_key.py`; its --ka reproduces A/B/C's 8c82bef3). A
#     second pin registers it as `lodo_lo0.5_hi3_def` next to the 410 key, and only the dump build selects it, exactly as
#     446bj2 did on A/B/C (config diff vs the default path: bulk_load, bulk_load_sha256). The restore trap undoes both pins.
LODO=""
if [ "${DEF_DAY_LODO:-0}" = "1" ]; then
  LK="$W/pert_scale_lodo_def_lo0.5_hi3.npz"; LID=lodo_lo0.5_hi3_def
  if [ -n "${DEF_DAY_LODO_HI:-}" ]; then
    LK="$W/pert_scale_lodo_def_lo0.5_hi${DEF_DAY_LODO_HI}.npz"; LID="lodo_lo0.5_hi${DEF_DAY_LODO_HI}_def"
  elif [ -n "${DEF_DAY_LODO_HI_PCT:-}" ]; then
    LK="$W/pert_scale_lodo_def_lo0.5_pct.npz"; LID=lodo_lo0.5_pct_def
  fi
  uv run python scripts/def_lodo_key.py --bundle "$B" --op "$W/whitener_union_pruned_def.npz" --cd4-root "$W/cd4_compact" \
      ${DEF_DAY_LODO_HI:+--hi "$DEF_DAY_LODO_HI"} ${DEF_DAY_LODO_HI_PCT:+--hi-pct "$DEF_DAY_LODO_HI_PCT"} --out "$LK" >> "$LOG" 2>&1 || die "def_lodo_key.py failed"
  LJSON="${LK%.npz}.json"
  # 530: 522's strength rider, between the LODO key and the pin. Unset: LK stays the key just written.
  if [ -n "${DEF_DAY_LODO_STRENGTH:-}" ]; then
    SK="$W/pert_scale_lodo_strength.npz"
    uv run python experiments/522-october-bulk-key/def_strength_key522.py --in "$LK" --bundle "$B" \
        --theta "$DEF_DAY_LODO_STRENGTH" --out "$SK" >> "$LOG" 2>&1 || die "def_strength_key522.py failed"
    LK="$SK"
    log "LODO strength rider: theta $DEF_DAY_LODO_STRENGTH"
  fi
  uv run python scripts/def_pin_panel_artifacts.py --whitener "$W/whitener_union_pruned_def.npz" \
      --whitener-metadata "$W/whitener_union_pruned_def.json" --whitener-id union_pruned_def --cd4-root "$W/cd4_compact" \
      --pert-scale "$LK" --pert-scale-id "$LID" --apply --out "$W/pin-lodo.json" >> "$LOG" 2>&1 \
      || die "LODO key pin failed"
  LODO="$LID"
  log "LODO key ON for the stack dump: $LODO $(shasum -a 256 "$LK" | cut -c1-16); cap $(python3 -c 'import json, sys; print(json.load(open(sys.argv[1]))["hi"])' "$LJSON")"
fi

# 3. the stack dump. The subshell keeps VCC_H1_UNION_G / VCC_BULK_DUMP / VCC_PERT_SCALE off the level builds.
#    Remove any previous npz first: a failed retry must not keep the last successful dump.
# 530: x1 (unset) keeps today's filename and does not export VCC_CD4_TARGET_G (DEF_ENV's g).
# x1.5 is 457: 1.5 x 0.3443655278649641 inside this subshell only. A strength rider renames the
# file so the level meta's bulk_load carries the theta; unset leaves dump_h1x2_cd4x1.npz.
case "${DEF_DAY_CD4_STACK:-x1}" in
  x1.5) DUMP="$W/dump_h1x2_cd4x15.npz" ;;
  *) DUMP="$W/dump_h1x2_cd4x1.npz" ;;
esac
CD4G=""
[ "${DEF_DAY_CD4_STACK:-x1}" = "x1.5" ] && CD4G="0.5165482917974462"
if [ -n "${DEF_DAY_LODO_STRENGTH:-}" ]; then
  STOKEN="$(printf '%s' "$DEF_DAY_LODO_STRENGTH" | tr '-' 'm')"
  DUMP="${DUMP%.npz}_s${STOKEN}.npz"
fi
[ -n "$CD4G" ] && log "CD4 stack x1.5: dump ${DUMP##*/} VCC_CD4_TARGET_G=$CD4G (level builds keep DEF_ENV)"
rm -f "$DUMP"
(
  export VCC_H1_UNION_G=0.56 VCC_BULK_DUMP="$DUMP"
  [ -n "$CD4G" ] && export VCC_CD4_TARGET_G="$CD4G"
  [ -n "$LODO" ] && export VCC_PERT_SCALE="$LODO"
  "${DRV[@]}" --name "$N-dump" --stages build > "$W/driver-dump.log" 2>&1
)
rc=$?
if [ "$rc" -ne 0 ] || [ ! -s "$DUMP" ]; then
  die "dump stage failed rc=$rc (see $W/driver-dump.log)"
fi
# The full builder log, not the driver's 25-line tail. time -l appends ~18 lines after the
# dump message, so the tail still contains it today; the stage log is the one that cannot drop it.
grep -q '\[bulk-dump\] wrote' "$W/logs/build.log" || die "dump stage log has no [bulk-dump] wrote line"
DSHA="$(shasum -a 256 "$DUMP" | cut -d' ' -f1)"; log "dump ${DSHA:0:16}"
# 690g: dump hook -- the hook's dump replaces the stack dump (kept on disk) under the same file name; validated first
BLC=""
if [ -n "${DEF_DAY_DUMP_HOOK:-}" ]; then
  HD="$W/dumphook"; rm -rf "$HD"; mkdir -p "$HD"; HB="${DUMP##*/}"
  log "690g: dump hook $DEF_DAY_DUMP_HOOK ($(shasum -a 256 "$DEF_DAY_DUMP_HOOK" | cut -c1-16)) on $HB"
  bash "$DEF_DAY_DUMP_HOOK" "$B" "$W" "$DUMP" "$W/identity.json" "$HD" > "$W/dumphook.log" 2>&1 || die "690g: the dump hook failed (see $W/dumphook.log)"
  [ -s "$HD/$HB" ] || die "690g: the dump hook wrote no $HD/$HB"
  CTXARG=""; [ -s "$HD/bulk_ctx.npz" ] && CTXARG="$HD/bulk_ctx.npz"
  uv run python scripts/def_dumphook_check.py "$DUMP" "$HD/$HB" "$B" $CTXARG >> "$LOG" 2>&1 || die "690g: the dump hook's output failed validation (see def_day.log)"
  DUMP="$HD/$HB"; DSHA="$(shasum -a 256 "$DUMP" | cut -d' ' -f1)"; log "690g: dump replaced by the hook's ${DSHA:0:16} (the stack dump stays at $W/$HB)"
  if [ -n "$CTXARG" ]; then
    grep -q "672: per-context bulk load" scripts/make_submission_pm_cd4.py || die "690g: the hook wrote bulk_ctx.npz but the builder has no VCC_BULK_LOAD_CTX (672)"
    BLC="$CTXARG:$(shasum -a 256 "$CTXARG" | cut -d' ' -f1)"
  fi
fi
export VCC_BULK_LOAD="$DUMP:$DSHA" VCC_BULK_SHARPEN_Q="$SQ"
[ -n "$BLC" ] && { export VCC_BULK_LOAD_CTX="$BLC"; log "690g: per-context bulk for the level builds: VCC_BULK_LOAD_CTX=${BLC##*/}"; }
[ -n "${DEF_DAY_DE_MAG:-}" ] && { export VCC_DE_MAG="$DEF_DAY_DE_MAG"; log "DE rider for the level builds: VCC_DE_MAG=$VCC_DE_MAG"; }
[ -n "${DEF_DAY_DE_UPRN:-}" ] && { export VCC_DE_UPRN="$DEF_DAY_DE_UPRN"; log "DE rider for the level builds: VCC_DE_UPRN=$VCC_DE_UPRN"; }   # 630dd: rank-neutral UP pass-through
[ -n "${DEF_DAY_DE_DNMAG:-}" ] && { export VCC_DE_DNMAG="$DEF_DAY_DE_DNMAG"; log "DE rider for the level builds: VCC_DE_DNMAG=$VCC_DE_DNMAG"; }   # 690a: DOWN magnitude pass-through
[ -n "${DEF_DAY_MDE_SIGN:-}" ] && { export VCC_DE_MDE_SIGN="$MDS"; log "DE rider for the level builds: VCC_DE_MDE_SIGN=$VCC_DE_MDE_SIGN"; }
[ -n "${DEF_DAY_MDE_SIGN_KEYGATE:-}" ] && { export VCC_DE_MDE_SIGN_KEYGATE="$DEF_DAY_MDE_SIGN_KEYGATE"; log "DE rider for the level builds: VCC_DE_MDE_SIGN_KEYGATE=$VCC_DE_MDE_SIGN_KEYGATE"; }   # 612dd: key gate pass-through
[ -n "${DEF_DAY_RC_CALLSET_CPM:-}" ] && { export VCC_RC_CALLSET_CPM="$DEF_DAY_RC_CALLSET_CPM"; log "DE rider for the level builds: VCC_RC_CALLSET_CPM=$VCC_RC_CALLSET_CPM"; }
[ -n "$CSB" ] && { export VCC_RC_CALLSET_BLOCK="$CSB"; log "DE rider for the level builds: VCC_RC_CALLSET_BLOCK=${CSB#prop:$W/}"; }
[ -n "$SL" ] && export VCC_SAMELINE="$SL"
if [ -n "${DEF_DAY_PRIOR_GATE:-}" ]; then
  [ -n "$SL" ] || die "DEF_DAY_PRIOR_GATE=$DEF_DAY_PRIOR_GATE, but same-line mode is off: the programme never reached the builder"
  export VCC_DE_PRIOR_GATE="$DEF_DAY_PRIOR_GATE"
fi

# 4-5. lvl0, then the re-level loop
POINTS=()
FACT="$UNIT"
build(){ local tag=$1 fact=$2 h5 rc pc
  if [ "$fact" = "$UNIT" ]; then unset VCC_BULK_CTX_SCALE; else export VCC_BULK_CTX_SCALE="$fact"; fi
  h5="$W/$N-$tag.h5ad"
  # The driver refuses to overwrite. Delete this tag's previous artifact so a retry rebuilds it,
  # then require the driver's own exit code. A pre-existing file is not evidence the new build ran.
  rm -f "$h5" "$h5.meta.json" "$W/$N-$tag.rc_audit.jsonl"
  # 554: the K-class emitter's per-block steering-dose record (diagnostics only; the emission is untouched), read by
  # def_postbuild.py --rc-audit for the final build
  VCC_RC_AUDIT_DUMP="$W/$N-$tag.rc_audit.jsonl" "${DRV[@]}" --name "$N-$tag" --stages build > "$W/driver-$tag.log" 2>&1
  rc=$?
  if [ "$rc" -ne 0 ] || [ ! -s "$h5" ]; then
    rm -f "$h5" "$h5.meta.json"
    die "$tag build failed rc=$rc: $(tail -3 "$W/driver-$tag.log" 2>/dev/null)"
  fi
  uv run python experiments/417-bulk-sharpen/ntclevel417.py "$h5" --bundle "$B" --contexts "$CTXS" --sd 0.01 >> "$LOG" 2>&1 \
    || die "ntclevel failed on $tag"
  POINTS+=("$W/$N-$tag.ntclevel.json:$fact")
  log "$tag ($fact) rho: $(uv run python scripts/def_day_helpers.py rho "$W/$N-$tag.ntclevel.json")"
  # past_cliff: 0 = a context is past its cliff (delete), 1 = safe, anything else = unreadable (stop).
  uv run python scripts/def_day_helpers.py past_cliff "$W/$N-$tag.ntclevel.json"
  pc=$?
  if [ "$pc" -eq 0 ]; then
    rm -f "$h5" "$h5.meta.json"; log "$tag has a context past its cliff: artifact deleted (never shippable)"
  elif [ "$pc" -ne 1 ]; then
    die "past_cliff failed on $tag rc=$pc; not treating the artifact as safe"
  fi; }
same_f(){ uv run python scripts/def_day_helpers.py same "$1" "$2"; }
build lvl0 "$FACT"
FINAL=""
for L in $(seq 1 "$MAXL"); do
  [ -n "${DEF_DAY_RHO_SCALE:-}" ] && log "re-level rho-scale=$DEF_DAY_RHO_SCALE"
  OUT="$(uv run python experiments/448-def-rehearsal-q175/relevel_iter.py "${POINTS[@]}" ${DEF_DAY_RHO_SCALE:+--rho-scale "$DEF_DAY_RHO_SCALE"} 2>&1)"
  rc=$?
  echo "$OUT" >> "$LOG"
  if [ "$rc" -ne 0 ]; then
    log "relevel_iter rc=$rc; shipping the nearest safe build"
    break
  fi
  NEXT="$(echo "$OUT" | grep -o 'VCC_BULK_CTX_SCALE=[A-Z0-9=.,]*' | cut -d= -f2-)"
  if [ -z "$NEXT" ]; then
    log "relevel_iter printed no factors; shipping the nearest safe build"
    break
  fi
  if echo "$OUT" | grep -q "^CONVERGED"; then
    # the converged factors are one of the builds already made: find it
    for P in "${POINTS[@]}"; do same_f "${P##*:}" "$NEXT" && [ -s "${P%%.ntclevel.json:*}.h5ad" ] && FINAL="${P%%.ntclevel.json:*}"; done
    if [ -z "$FINAL" ]; then
      # each context converged in a different build: build the combination once, and take it only
      # when the new measurement is inside the same 1% band relevel required of the pieces.
      log "CONVERGED per context at $NEXT, but no single build carries all of them: building it"
      build lvlC "$NEXT"
      if [ -s "$W/$N-lvlC.h5ad" ] && uv run python scripts/def_day_helpers.py inband "$W/$N-lvlC.ntclevel.json" ${DEF_DAY_RHO_SCALE:+"$DEF_DAY_RHO_SCALE"}; then
        FINAL="$W/$N-lvlC"
      else
        log "lvlC is missing or outside the 1% band; not labelling it CONVERGED"
        rm -f "$W/$N-lvlC.h5ad" "$W/$N-lvlC.h5ad.meta.json"
      fi
    fi
    if [ -n "$FINAL" ]; then
      log "CONVERGED at $NEXT after $((L - 1)) levelled build(s)"
    else
      log "CONVERGED factors were not a shippable build; shipping the nearest safe build"
    fi
    break
  fi
  if ! uv run python scripts/def_day_helpers.py inbound "$NEXT"; then
    log "next factors $NEXT leave the builder bound [0.7, 1.4]; shipping the nearest safe build"
    break
  fi
  build "lvl$L" "$NEXT"
done
if [ -z "$FINAL" ]; then
  # not converged: ship the last build whose contexts are all below their cliffs, the nearest to rho* among those.
  # rho* is computed in the helper from the same q and sd as relevel_iter.py.
  FINAL="$(uv run python scripts/def_day_helpers.py best_safe ${DEF_DAY_RHO_SCALE:+--scale "$DEF_DAY_RHO_SCALE"} "${POINTS[@]}")" || die "best_safe failed"
  [ -n "$FINAL" ] || die "no build with every rho <= 1; do not ship"
  log "NOT converged in $MAXL levels; shipping the nearest safe build: ${FINAL##*/}"
fi
FF="$(for P in "${POINTS[@]}"; do [ "${P%%.ntclevel.json:*}" = "$FINAL" ] && echo "${P##*:}"; done | head -1)"
[ -n "$FF" ] || die "final ${FINAL} has no factor string in the point list"
# Check the chosen artifact before deleting the others, so a meta mismatch still leaves them on disk.
uv run python scripts/def_day_helpers.py meta_matches "$FINAL.h5ad.meta.json" "$FF" "$SQ" \
  || die "final meta does not match factors $FF q $SQ (the h5ad is not the build this run measured)"
# 530's spec is appended only when one of its knobs is set, so the unset argv is today's two-argument call.
_r530=""
[ -n "${DEF_DAY_RC_CALLSET_CPM:-}" ] && _r530=1
[ "${DEF_DAY_CD4_STACK:-x1}" != "x1" ] && _r530=1
[ -n "${DEF_DAY_RHO_SCALE:-}" ] && _r530=1
[ -n "${DEF_DAY_LODO_STRENGTH:-}" ] && _r530=1
[ -n "${DEF_DAY_HOOK_ESTIMATOR:-}" ] && _r530=1
[ -n "${DEF_DAY_LODO_HI:-}" ] && _r530=1
[ -n "${DEF_DAY_CALLSET_BLOCK:-}" ] && _r530=1
if [ -n "$_r530" ]; then
  _rc="${DEF_DAY_RC_CALLSET_CPM:--}"; _st="${DEF_DAY_CD4_STACK:-x1}"; _rh="${DEF_DAY_RHO_SCALE:--}"; _ls="${DEF_DAY_LODO_STRENGTH:--}"
  _spec="rc=${_rc},stack=${_st},rho=${_rh},strength=${_ls}"
  [ -n "${DEF_DAY_LODO_HI:-}" ] && _spec="${_spec},lodo_hi=${DEF_DAY_LODO_HI}"
  [ -n "${DEF_DAY_CALLSET_BLOCK:-}" ] && _spec="${_spec},callset=${DEF_DAY_CALLSET_BLOCK}"
  uv run python scripts/def_day_helpers.py riders_match "$FINAL.h5ad.meta.json" "${VCC_DE_MAG:--}" "${VCC_DE_MDE_SIGN:--}" \
      "$_spec" \
    || die "the final meta does not carry the requested October riders (rc=${_rc} stack=${_st} rho=${_rh} strength=${_ls}${DEF_DAY_LODO_HI:+ lodo_hi=$DEF_DAY_LODO_HI}${DEF_DAY_CALLSET_BLOCK:+ callset=$DEF_DAY_CALLSET_BLOCK})"
else
  uv run python scripts/def_day_helpers.py riders_match "$FINAL.h5ad.meta.json" "${VCC_DE_MAG:--}" "${VCC_DE_MDE_SIGN:--}" \
    || die "the final meta does not carry the requested DE riders (VCC_DE_MAG=${VCC_DE_MAG:-} VCC_DE_MDE_SIGN=${VCC_DE_MDE_SIGN:-})"
fi
if [ -n "${VCC_DE_MDE_SIGN_KEYGATE:-}" ]; then   # 612dd: key gate pass-through -- the final meta must carry the gate it was built with
  python3 - "$FINAL.h5ad.meta.json" "$VCC_DE_MDE_SIGN_KEYGATE" <<'PY' || die "the final meta does not carry VCC_DE_MDE_SIGN_KEYGATE=$VCC_DE_MDE_SIGN_KEYGATE"
import json, sys
g = json.load(open(sys.argv[1]))["config"].get("de_mde_sign_keygate") or {}
sys.exit(0 if abs(float(g.get("q", -1.0)) - float(sys.argv[2])) < 1e-12 else 1)
PY
fi
if [ -n "${VCC_DE_DNMAG:-}" ]; then   # 690a: DOWN magnitude pass-through -- the final meta must carry the variant and the reviewed module
  python3 - "$FINAL.h5ad.meta.json" "$VCC_DE_DNMAG" "b9d6cfcf17ad783c1ea2607f95366173a0fc1f7920c4104604f9ee925ac55c87" <<'PY' || die "the final meta does not carry VCC_DE_DNMAG=$VCC_DE_DNMAG with dn_mag module b9d6cfcf (693 v2)"
import json, sys
u = json.load(open(sys.argv[1]))["config"].get("de_dnmag") or {}
sys.exit(0 if u.get("variant") == sys.argv[2] and u.get("module_sha256") == sys.argv[3] else 1)
PY
fi
if [ -n "${VCC_DE_DNMAG:-}" ]; then   # 690h: DNMAG delivery report -- per context from the final rc_audit dump; REPORT, never fatal
  python3 - "$FINAL.rc_audit.jsonl" "$W/dnmag_audit.json" "$W/dnmag_low.txt" <<'PY' 2>&1 | while IFS= read -r _l; do log "690h: $_l"; done || true
import collections, json, sys
src, out_json, out_low = sys.argv[1:4]
try:
    rows = [json.loads(l) for l in open(src) if l.strip()]
    rows = [r for r in rows if "dn_frac" in r and "dn_moved" in r]
except Exception as exc:
    print(f"dnmag: could not read {src}: {exc}")
    rows = []
if not rows:
    print("dnmag LOW_CONTEXTS -1")
    sys.exit(0)
by = collections.defaultdict(list)
for r in rows:
    by[str(r.get("context"))].append(r)
rep, low = {}, []
for c in sorted(by):
    rs = by[c]
    fr = sum(float(r["dn_frac"]) for r in rs) / len(rs)
    mv = sum(float(r["dn_moved"]) for r in rs) / len(rs)
    rec = {"targets": len(rs), "mean_dn_frac": round(fr, 4), "mean_dn_moved": round(mv, 1),
           "targets_moved_lt_200": sum(1 for r in rs if float(r["dn_moved"]) < 200),
           "invariant_failed": int(sum(float(r.get("dn_invariant_failed", 0)) for r in rs)),
           "guard_failed": int(sum(float(r.get("dn_guard_failed", 0)) for r in rs)),
           "low": bool(fr < 0.10 or mv < 200)}
    rep[c] = rec
    if rec["low"]:
        low.append(c)
    print(f"dnmag {c}: " + json.dumps(rec))
json.dump(rep, open(out_json, "w"), indent=1)
if low:
    open(out_low, "w").write("\n".join(low) + "\n")
print("dnmag LOW_CONTEXTS " + (",".join(low) if low else "none"))
PY
fi
if [ -n "${BLC:-}" ]; then   # 690g: dump hook -- the final meta must carry the per-context bulk it was built with
  python3 - "$FINAL.h5ad.meta.json" "${BLC##*:}" <<'PY' || die "690g: the final meta does not carry bulk_load_ctx_sha256 ${BLC##*:}"
import json, sys
c = json.load(open(sys.argv[1]))["config"]
sys.exit(0 if c.get("bulk_load_ctx_sha256") == sys.argv[2] else 1)
PY
  log "690g: per-context outcome (builder): $(grep -hE '672|bulk.load.ctx' "$W/logs/build.log" 2>/dev/null | tail -4 | tr '\n' ' ' | cut -c1-400)"
fi
if [ -n "${VCC_DE_UPRN:-}" ]; then   # 630dd: rank-neutral UP pass-through -- the final meta must carry the variant and the reviewed module
  python3 - "$FINAL.h5ad.meta.json" "$VCC_DE_UPRN" "7e6b3e52e2ad374b7b6e21e92327fcf65771433693eb3846c09c0189f3f2959f" <<'PY' || die "the final meta does not carry VCC_DE_UPRN=$VCC_DE_UPRN with rank_up module 7e6b3e52"
import json, sys
u = json.load(open(sys.argv[1]))["config"].get("de_uprn") or {}
sys.exit(0 if u.get("variant") == sys.argv[2] and u.get("module_sha256") == sys.argv[3] else 1)
PY
  # REPORT (never fatal; d4cd5f): per context, the pass's delivery and its fail-safe outcomes. A target whose guard or
  # invariants failed keeps its zs block, which is the safe outcome. The prereg's capacity floor (mean frac_got >= 0.127,
  # moved >= 500 per target) is printed for reading, not enforced here.
  python3 - "$FINAL.rc_audit.jsonl" "$W/uprn_audit.json" "$W/uprn_effectively_off.txt" <<'PY' 2>&1 | while IFS= read -r _l; do log "630dd: $_l"; done || true
import collections, json, sys
# ONE greppable line always ends the report (coordinator 10-02; 607b): `uprn FAILSAFE_TOTAL <n>`, where n = the fail-safe
# triggers (guard + invariant) summed over contexts, and -1 = the audit could not be read or held no ur_ rows. A KA630
# build passes only with the line PRESENT and n == 0.
try:
    rows = collections.defaultdict(list)
    for line in open(sys.argv[1]):
        r = json.loads(line)
        if "ur_n_sel" in r:
            rows[str(r.get("context"))].append(r)
    rep = {}
    for c, L in sorted(rows.items()):
        n = len(L)
        fg = sum(r.get("ur_frac_got", 0.0) for r in L) / n
        nm = sum(r.get("ur_n_moved", 0.0) for r in L) / n
        rep[c] = dict(targets=n, mean_frac_got=fg, mean_n_moved=nm,
                      guard_failed=int(sum(r.get("ur_guard_failed", 0.0) > 0 for r in L)),
                      invariant_failed=int(sum(r.get("ur_invariant_failed", 0.0) > 0 for r in L)),
                      no_pos_donor_genes=int(sum(r.get("ur_n_no_pos_donor", 0.0) for r in L)),
                      capacity_floor_met=bool(fg >= 0.127 and nm >= 500))
        print(f"uprn {c}: " + json.dumps(rep[c]))
    json.dump(rep, open(sys.argv[2], "w"), indent=1)
    # LOUD, not fatal (d4cd5f): a context where > 10 % of the targets fell back to their zs block has the pass effectively off
    off = [f"{c} ({v['guard_failed'] + v['invariant_failed']}/{v['targets']})" for c, v in rep.items()
           if v["guard_failed"] + v["invariant_failed"] > 0.10 * v["targets"]]
    for o in off:
        print(f"630 EFFECTIVELY OFF in context {o}")
    if off:
        open(sys.argv[3], "w").write("\n".join(off) + "\n")
    else:
        import os
        if os.path.exists(sys.argv[3]):
            os.remove(sys.argv[3])
    if not rep:
        print("uprn: WARNING no ur_ diagnostics in the final rc_audit dump")
        print("uprn FAILSAFE_TOTAL -1")
    else:
        print(f"uprn FAILSAFE_TOTAL {sum(v['guard_failed'] + v['invariant_failed'] for v in rep.values())}")
except Exception as e:   # noqa: BLE001 -- never fatal, never silent
    print(f"uprn: ERROR reading {sys.argv[1]}: {type(e).__name__}: {e}")
    print("uprn FAILSAFE_TOTAL -1")
PY
fi
if [ -n "${VCC_RC_CALLSET_BLOCK:-}" ]; then
  python3 - "$FINAL.h5ad.meta.json" "$VCC_RC_CALLSET_BLOCK" <<'PY' || die "the final meta does not carry VCC_RC_CALLSET_BLOCK=$VCC_RC_CALLSET_BLOCK"
import json, sys
c = json.load(open(sys.argv[1]))["config"].get("rc_callset_block") or {}
p = sys.argv[2].split(":")
want = float(p[3]) if len(p) >= 4 else 1.0
steer = p[4].lower() if len(p) == 5 else "kept"
sys.exit(0 if (c.get("source") == "prop" and c.get("sha256") == p[2] and abs(float(c.get("swap", -1.0)) - want) < 1e-12
               and c.get("steer") == steer) else 1)
PY
fi
for P in "${POINTS[@]}"; do [ "${P%%.ntclevel.json:*}" != "$FINAL" ] && rm -f "${P%%.ntclevel.json:*}.h5ad" "${P%%.ntclevel.json:*}.h5ad.meta.json" "${P%%.ntclevel.json:*}.rc_audit.jsonl"; done
ln -f "$FINAL.h5ad" "$W/${N}_final.h5ad" || die "ln of the final h5ad failed"
ln -f "$FINAL.h5ad.meta.json" "$W/${N}_final.h5ad.meta.json" || die "ln of the final meta failed"
# 554: the final build's per-block dose record (def_postbuild.py --rc-audit)
ln -f "$FINAL.rc_audit.jsonl" "$W/${N}_final.rc_audit.jsonl" || die "554: the final build wrote no rc_audit dump"
printf 'VCC_BULK_LOAD=%s\nVCC_BULK_SHARPEN_Q=%s\nVCC_BULK_CTX_SCALE=%s\nVCC_SAMELINE=%s\n' "$VCC_BULK_LOAD" "$SQ" "$FF" "${SL:-}" > "$W/def_day.final.env"
[ -n "${VCC_DE_PRIOR_GATE:-}" ] && printf 'VCC_DE_PRIOR_GATE=%s\n' "$VCC_DE_PRIOR_GATE" >> "$W/def_day.final.env"
[ -n "${DEF_DAY_DUMP_HOOK:-}" ] && printf 'DEF_DAY_DUMP_HOOK=%q\n' "$DEF_DAY_DUMP_HOOK" >> "$W/def_day.final.env"   # 690g: dump hook
[ -n "${BLC:-}" ] && printf 'VCC_BULK_LOAD_CTX=%s\n' "$BLC" >> "$W/def_day.final.env"   # 690g: dump hook
[ -n "${VCC_DE_MAG:-}" ] && printf 'VCC_DE_MAG=%s\n' "$VCC_DE_MAG" >> "$W/def_day.final.env"
[ -n "${VCC_DE_UPRN:-}" ] && printf 'VCC_DE_UPRN=%s\n' "$VCC_DE_UPRN" >> "$W/def_day.final.env"   # 630dd: rank-neutral UP pass-through
[ -n "${VCC_DE_DNMAG:-}" ] && printf 'VCC_DE_DNMAG=%s\n' "$VCC_DE_DNMAG" >> "$W/def_day.final.env"   # 690a: DOWN magnitude pass-through
[ -n "${VCC_DE_MDE_SIGN:-}" ] && printf 'VCC_DE_MDE_SIGN=%s\n' "$VCC_DE_MDE_SIGN" >> "$W/def_day.final.env"
[ -n "${VCC_DE_MDE_SIGN_KEYGATE:-}" ] && printf 'VCC_DE_MDE_SIGN_KEYGATE=%s\n' "$VCC_DE_MDE_SIGN_KEYGATE" >> "$W/def_day.final.env"   # 612dd: key gate pass-through
[ -n "${VCC_RC_CALLSET_CPM:-}" ] && printf 'VCC_RC_CALLSET_CPM=%s\n' "$VCC_RC_CALLSET_CPM" >> "$W/def_day.final.env"
[ -n "${DEF_DAY_RHO_SCALE:-}" ] && printf 'DEF_DAY_RHO_SCALE=%s\n' "$DEF_DAY_RHO_SCALE" >> "$W/def_day.final.env"
[ "${DEF_DAY_CD4_STACK:-x1}" != "x1" ] && printf 'DEF_DAY_CD4_STACK=%s\n' "$DEF_DAY_CD4_STACK" >> "$W/def_day.final.env"
[ -n "${DEF_DAY_LODO_STRENGTH:-}" ] && printf 'DEF_DAY_LODO_STRENGTH=%s\n' "$DEF_DAY_LODO_STRENGTH" >> "$W/def_day.final.env"
[ -n "${DEF_DAY_HOOK_ESTIMATOR:-}" ] && printf 'DEF_DAY_HOOK_ESTIMATOR=%s\n' "$DEF_DAY_HOOK_ESTIMATOR" >> "$W/def_day.final.env"
[ -n "${DEF_DAY_LODO_HI:-}" ] && printf 'DEF_DAY_LODO_HI=%s\n' "$DEF_DAY_LODO_HI" >> "$W/def_day.final.env"
[ -n "${DEF_DAY_CENSUS_MAX_ZERO:-}" ] && printf 'DEF_DAY_CENSUS_MAX_ZERO=%s\n' "$DEF_DAY_CENSUS_MAX_ZERO" >> "$W/def_day.final.env"
[ -n "${VCC_RC_CALLSET_BLOCK:-}" ] && printf 'VCC_RC_CALLSET_BLOCK=%s\n' "$VCC_RC_CALLSET_BLOCK" >> "$W/def_day.final.env"
# 554: the dose guard's overrides, recorded
[ "${DEF_DAY_ALLOW_WIDE_LIBS:-}" = 1 ] && printf 'DEF_DAY_ALLOW_WIDE_LIBS=1\n' >> "$W/def_day.final.env"
[ "${DEF_DAY_ALLOW_LOW_DOSE:-}" = 1 ] && printf 'DEF_DAY_ALLOW_LOW_DOSE=1\n' >> "$W/def_day.final.env"
# 641: the depth row's override, recorded
[ "${DEF_DAY_ALLOW_DEPTH:-}" = 1 ] && printf 'DEF_DAY_ALLOW_DEPTH=1\n' >> "$W/def_day.final.env"
# 642: the HEK293 contingency's knob and decision, recorded
[ -n "${DEF_DAY_ALLOW_NONCCLE_REF:-}" ] && printf 'DEF_DAY_ALLOW_NONCCLE_REF=%q\n' "$DEF_DAY_ALLOW_NONCCLE_REF" >> "$W/def_day.final.env"
[ -n "${DEF_DAY_NONCCLE_DECISION:-}" ] && printf 'DEF_DAY_NONCCLE_DECISION=%q\n' "$DEF_DAY_NONCCLE_DECISION" >> "$W/def_day.final.env"
# 644: the same-line skip, recorded
[ -n "${DEF_DAY_SAMELINE_SKIP:-}" ] && printf 'DEF_DAY_SAMELINE_SKIP=%q\n' "$DEF_DAY_SAMELINE_SKIP" >> "$W/def_day.final.env"
log "final: ${FINAL##*/} factors $FF -> $W/${N}_final.h5ad"

# 5b. the rule-10 diff (481 F5), BEFORE the EXIT trap restores the builder: def_postbuild.py reads the day's operator
#     pin out of the builder source, and the allowed diffs are exactly what this run changed on purpose.
INC="${DEF_DAY_INCUMBENT_META:-submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json}"
ALLOW="whitener_id,pert_scale,pert_scale_sha256,bulk_load,bulk_load_sha256,bulk_ctx_scale"
[ -n "$SL" ] && ALLOW="$ALLOW,sameline"
[ -n "${VCC_DE_PRIOR_GATE:-}" ] && ALLOW="$ALLOW,de_prior_gate"
[ "$SQ" != "$Q" ] && ALLOW="$ALLOW,bulk_sharpen_q,bulk_sharpen_pt"
[ -n "${VCC_DE_MAG:-}" ] && ALLOW="$ALLOW,de_mag"
[ -n "${BLC:-}" ] && ALLOW="$ALLOW,bulk_load_ctx,bulk_load_ctx_sha256"   # 690g: dump hook
[ -n "${VCC_DE_UPRN:-}" ] && ALLOW="$ALLOW,de_uprn"   # 630dd: rank-neutral UP pass-through
[ -n "${VCC_DE_DNMAG:-}" ] && ALLOW="$ALLOW,de_dnmag"   # 690a: DOWN magnitude pass-through
[ -n "${VCC_DE_MDE_SIGN:-}" ] && ALLOW="$ALLOW,de_mde_sign"
[ -n "${VCC_DE_MDE_SIGN_KEYGATE:-}" ] && ALLOW="$ALLOW,de_mde_sign_keygate"   # 612dd: key gate pass-through
[ -n "${VCC_RC_CALLSET_BLOCK:-}" ] && ALLOW="$ALLOW,rc_callset_block"
[ -n "${VCC_RC_CALLSET_CPM:-}" ] && ALLOW="$ALLOW,rc_callset_cpm"
# 499 F2: 496's emitter pin moves config.de_sign_emitter_sha256 on every build, knob unset included. Allow that one key
# only when the builder source pin differs from the incumbent's recorded sha (today both are edb5b8dc, so nothing is
# added). A build whose recorded sha moved while the source pin still matches the incumbent still fails the gate.
_inc_em="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("de_sign_emitter_sha256",""))' "$INC" 2>/dev/null || true)"
_src_em="$(python3 -c 'import re,sys; t=open(sys.argv[1]).read(); m=re.search(r"^DE_SIGN_EMITTER_SHA256 = .([0-9a-f]{64}).", t, re.M); print(m.group(1) if m else "")' scripts/make_submission_pm_cd4.py 2>/dev/null || true)"
[ -n "$_src_em" ] && [ -n "$_inc_em" ] && [ "$_src_em" != "$_inc_em" ] && ALLOW="$ALLOW,de_sign_emitter_sha256"
uv run python scripts/def_postbuild.py --meta "$W/${N}_final.h5ad.meta.json" --bundle "$B" --cd4-root "$W/cd4_compact" \
    --incumbent "$INC" --allow-config-diff "$ALLOW" --rc-audit "$W/${N}_final.rc_audit.jsonl" \
    --out "$W/postbuild.json" >> "$LOG" 2>&1 \
  && log "postbuild PASS vs ${INC##*/} (allowed: $ALLOW)" || die "postbuild FAILED vs ${INC##*/}: see $W/postbuild.json"

# 6. audits. Write the zero-gene audit to a file so its exit code is the audit's, not head's.
uv run python experiments/428-zero-gene-audit/audit428.py "$W/${N}_final.h5ad" --bundle "$B" > "$W/audit428.txt" 2>&1
rc=$?
summary="$(head -1 "$W/audit428.txt" 2>/dev/null || true)"
log "audit: $summary"
if [ "$rc" -ne 0 ]; then
  die "zero-gene audit failed rc=$rc (see $W/audit428.txt)"
fi
case "$summary" in
  *"ZERO block total (target excluded): 0 in "*) ;;
  *) die "zero-gene audit is not clean: $summary" ;;
esac
python3 - "$W/${N}_final.h5ad.meta.json" <<'PY' | tee -a "$LOG"
import json, sys
c = json.load(open(sys.argv[1]))["config"]
print("meta: bulk_ctx_scale", c.get("bulk_ctx_scale"), "| q", c.get("bulk_sharpen_q"), "| pt", c.get("bulk_sharpen_pt"), "| sameline", json.dumps(c.get("sameline"))[:300])
if "de_prior_gate" in c:
    print("meta: de_prior_gate", json.dumps(c["de_prior_gate"])[:400])
PY
if ! "${DRV[@]}" --name "${N}_final" --stages prep > "$W/driver-prep.log" 2>&1; then
  die "prep dry run failed (see $W/driver-prep.log): $(tail -3 "$W/driver-prep.log")"
fi
log "prep dry run: $(grep -iE 'dry run|cells|verified' "$W/driver-prep.log" | head -3 | tr '\n' ' ')"
[ -s "$W/uprn_effectively_off.txt" ] && log "630dd: WARNING 630 EFFECTIVELY OFF (fail-safe on > 10 % of targets) in: $(tr '\n' ' ' < "$W/uprn_effectively_off.txt")"   # 630dd: rank-neutral UP pass-through
[ -s "$W/dnmag_low.txt" ] && log "690h: WARNING DNMAG under-delivered (mean dn_frac < 0.10 or < 200 genes/target) in: $(tr '\n' ' ' < "$W/dnmag_low.txt")"   # 690h: DNMAG delivery report
if [ "${DEF_DAY_KEEP_PINS:-0}" = "1" ]; then
  log "DONE. Next, by hand: vcc prep (real) on $W/${N}_final.h5ad, commit the pins (builder left pinned), vcc submit (-m joke, -d empty)."
else
  log "DONE. Next, by hand: vcc prep (real) on $W/${N}_final.h5ad, vcc submit (-m joke, -d empty). The exit trap restores the builder; re-run the pin stages if you need that diff to commit."
fi
