#!/usr/bin/env python3
"""538: October extras on top of 530's bundle and 528's call-propensity block. NOT APPLIED to the live tree.

Apply with experiments/538-cloud-bundle-extras/apply_october538.py, which runs 530, then 528, then this
file. 528's own anchors do not match a 530-patched def_day.sh (the MDE grep is no longer immediately
followed by fi). This file re-inserts 528's snippets, byte-for-byte, on anchors that survive 530, and
skips that insertion when DEF_DAY_CALLSET_BLOCK is already present.

Knobs (each default off; unset keeps today's exports and today's argv):

  DEF_DAY_LODO_HI=<h>            1 < h <= 20, a plain decimal. Passed to def_lodo_key.py as --hi.
                                 Needs DEF_DAY_LODO=1. Refused together with DEF_DAY_LODO_HI_PCT.
                                 Recorded on the pinned key (filename, id, sidecar json "hi").
  DEF_DAY_CALLSET_BLOCK=prop[:<swap>[:<steer>]]
                                 528's pass-through, re-anchored here when 528's patch does not apply.
  DEF_DAY_CENSUS_MAX_ZERO=<n>    integer 6..30, no leading zero. Fresh work dir only. Forwarded to
                                 def_pipeline.py --census-max-zero, which forwards def_census.py
                                 --max-zero-carrier. Logged and written into def_day.final.env.

Also: def_day.sh derives REPO from its own path; def_postbuild.py --require-config-diff;
resolver switches lodo_hi and callset_block (PENDING); riders_match optional keys; gate524 rows.

Unset: the guarded branches add no export and no argv word. --require-config-diff empty adds no check.

usage:
  uv run python experiments/538-cloud-bundle-extras/patch_extras538.py --root <repo> [--check]
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

# 528's snippets, copied verbatim from experiments/528-callprop-block/patch_defday528.py.
# Only the insertion point of VALIDATE changes (before SQ="$Q", which is where 530's validators end).
VALIDATE = (
    '# 528: the call-propensity CPM block (DE-only; validated here, computed after preflight, exported for the level builds)\n'
    'if [ -n "${DEF_DAY_CALLSET_BLOCK:-}" ]; then\n'
    '  case "$DEF_DAY_CALLSET_BLOCK" in\n'
    '    prop) ;;\n'
    '    prop:*) python3 -c "import sys; p = sys.argv[1].split(\':\'); v = float(p[0]); sys.exit(0 if (0.0 < v <= 1.0 and len(p) <= 2 and (len(p) == 1 or p[1] in (\'kept\', \'key\', \'cpm\'))) else 1)" "${DEF_DAY_CALLSET_BLOCK#prop:}" 2>/dev/null \\\n'
    '              || die "DEF_DAY_CALLSET_BLOCK=prop:<swap>[:<steer>] needs 0 < swap <= 1 and steer kept|key|cpm (got $DEF_DAY_CALLSET_BLOCK)" ;;\n'
    '    *) die "DEF_DAY_CALLSET_BLOCK must be prop or prop:<swap>[:<steer>] (got $DEF_DAY_CALLSET_BLOCK)" ;;\n'
    '  esac\n'
    '  grep -q "VCC_RC_CALLSET_BLOCK" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_RC_CALLSET_BLOCK (528)"\n'
    '  grep -q "block_steer=" src/vcc2026/rc_emit_k.py || die "the emitter has no block_key / block_steer (528)"\n'
    '  [ -f experiments/528-callprop-block/callblock528.py ] || die "experiments/528-callprop-block/callblock528.py is missing"\n'
    'fi\n'
)
COMPUTE = (
    '# 2c. 528: DEF_DAY_CALLSET_BLOCK=prop[:<swap>] -- the day\'s block keys from D/E/F\'s own controls (nothing fitted today)\n'
    'CSB=""\n'
    'if [ -n "${DEF_DAY_CALLSET_BLOCK:-}" ]; then\n'
    '  uv run python experiments/528-callprop-block/callblock528.py --bundle "$B" --contexts "$CTXS" --out "$W/callblock_def.npz" \\\n'
    '      >> "$LOG" 2>&1 || die "callblock528.py (the call-propensity block keys) failed"\n'
    '  CSB="prop:$W/callblock_def.npz:$(shasum -a 256 "$W/callblock_def.npz" | cut -d\' \' -f1)"\n'
    '  case "$DEF_DAY_CALLSET_BLOCK" in prop:*) CSB="$CSB:${DEF_DAY_CALLSET_BLOCK#prop:}" ;; esac\n'
    '  log "call-propensity block keys for $CTXS: VCC_RC_CALLSET_BLOCK=prop:${CSB#prop:$W/}"\n'
    'fi\n\n'
)
EXPORT = '[ -n "$CSB" ] && { export VCC_RC_CALLSET_BLOCK="$CSB"; log "DE rider for the level builds: VCC_RC_CALLSET_BLOCK=${CSB#prop:$W/}"; }\n'
FINALCHECK = (
    'if [ -n "${VCC_RC_CALLSET_BLOCK:-}" ]; then\n'
    '  python3 - "$FINAL.h5ad.meta.json" "$VCC_RC_CALLSET_BLOCK" <<\'PY\' || die "the final meta does not carry VCC_RC_CALLSET_BLOCK=$VCC_RC_CALLSET_BLOCK"\n'
    'import json, sys\n'
    'c = json.load(open(sys.argv[1]))["config"].get("rc_callset_block") or {}\n'
    'p = sys.argv[2].split(":")\n'
    'want = float(p[3]) if len(p) >= 4 else 1.0\n'
    'steer = p[4].lower() if len(p) == 5 else "kept"\n'
    'sys.exit(0 if (c.get("source") == "prop" and c.get("sha256") == p[2] and abs(float(c.get("swap", -1.0)) - want) < 1e-12\n'
    '               and c.get("steer") == steer) else 1)\n'
    'PY\n'
    'fi\n'
)
CALLSET_HEADER = (
    "#   DEF_DAY_CALLSET_BLOCK=prop[:<swap>[:<steer>]] (528, default off): the DE call set's CPM block chosen per context by a\n"
    "#    call-propensity key from the day's OWN controls (log10 control CPM + 22 x excess dropout, median-rescaled to H1's\n"
    "#    dropout scale; callblock528.py), exported as VCC_RC_CALLSET_BLOCK for the level builds. Builder/emitter >= 528's.\n"
)
CALLSET_ENV = '[ -n "${VCC_RC_CALLSET_BLOCK:-}" ] && printf \'VCC_RC_CALLSET_BLOCK=%s\\n\' "$VCC_RC_CALLSET_BLOCK" >> "$W/def_day.final.env"\n'
CALLSET_ALLOW = '[ -n "${VCC_RC_CALLSET_BLOCK:-}" ] && ALLOW="$ALLOW,rc_callset_block"\n'


def _once(src: str, old: str, name: str) -> None:
    n = src.count(old)
    if n != 1:
        raise SystemExit(f"{name}: anchor count {n} != 1: {old[:140]!r}")


def _apply(src: str, edits: list[tuple[str, str]], name: str) -> str:
    for old, _new in edits:
        _once(src, old, name)
    for old, new in edits:
        src = src.replace(old, new, 1)
    return src


def _markers(texts: dict[str, str]) -> dict[str, bool]:
    day = texts["def_day"]
    return {
        # DEF_DAY_LODO_HI is a prefix of DEF_DAY_LODO_HI_PCT. Match the --hi expansion, not the name.
        "def_day": (
            '--hi "$DEF_DAY_LODO_HI"' in day
            and 'REPO="$(CDPATH=' in day
            and "DEF_DAY_CENSUS_MAX_ZERO" in day
            and "DEF_DAY_CALLSET_BLOCK" in day
        ),
        "helpers": '{"lodo_hi", "callset"}' in texts["helpers"],
        "recipe": '"lodo_hi": "DEF_DAY_LODO_HI"' in texts["recipe"] and '"callset_block": "DEF_DAY_CALLSET_BLOCK"' in texts["recipe"],
        "decisions": '"lodo_hi": {' in texts["decisions"] and '"callset_block": {' in texts["decisions"],
        "gate": '"DEF_DAY_LODO_HI":' in texts["gate"] and '"DEF_DAY_CALLSET_BLOCK":' in texts["gate"],
        "postbuild": "--require-config-diff" in texts["postbuild"],
        "pipeline": "--census-max-zero" in texts["pipeline"],
    }


def callset_edits() -> list[tuple[str, str]]:
    """528's seven insertions. VALIDATE is re-anchored onto the line 530 leaves in front of the sharpening block."""
    keep = "#   DEF_DAY_KEEP_PINS=1 (default off): skip the exit-trap restore so `git diff` on the builder is the pin to commit.\n"
    return [
        (keep, CALLSET_HEADER + keep),
        ('SQ="$Q"\n', VALIDATE + 'SQ="$Q"\n'),
        ("# 2b. DEF_DAY_LODO=1 (460; default off)", COMPUTE + "# 2b. DEF_DAY_LODO=1 (460; default off)"),
        ('[ -n "$SL" ] && export VCC_SAMELINE="$SL"\n', EXPORT + '[ -n "$SL" ] && export VCC_SAMELINE="$SL"\n'),
        (
            'for P in "${POINTS[@]}"; do [ "${P%%.ntclevel.json:*}" != "$FINAL" ] && rm -f',
            FINALCHECK + 'for P in "${POINTS[@]}"; do [ "${P%%.ntclevel.json:*}" != "$FINAL" ] && rm -f',
        ),
        (
            'log "final: ${FINAL##*/} factors $FF -> $W/${N}_final.h5ad"\n',
            CALLSET_ENV + 'log "final: ${FINAL##*/} factors $FF -> $W/${N}_final.h5ad"\n',
        ),
        (
            '[ -n "${VCC_DE_MDE_SIGN:-}" ] && ALLOW="$ALLOW,de_mde_sign"\n',
            '[ -n "${VCC_DE_MDE_SIGN:-}" ] && ALLOW="$ALLOW,de_mde_sign"\n' + CALLSET_ALLOW,
        ),
    ]


def patch_def_day(src: str) -> str:
    edits = [
        (
            "#    P 96.64697522769474 is exactly cap 6 on A/B/C (495's artifact); the rounded 96.647 is not.\n",
            "#    P 96.64697522769474 is exactly cap 6 on A/B/C (495's artifact); the rounded 96.647 is not.\n"
            "#   DEF_DAY_LODO_HI=<h> (538, default off; needs DEF_DAY_LODO=1; refused with DEF_DAY_LODO_HI_PCT):\n"
            "#    a fixed cap, passed to def_lodo_key.py as --hi. 1 < h <= 20, a plain decimal. Unset keeps cap 3,\n"
            "#    or the percentile form. The pinned key is lodo_lo0.5_hi<h>_def and its json records hi.\n"
            "#   DEF_DAY_CENSUS_MAX_ZERO=<n> (538, default off; integer 6..30; a fresh work dir only): raises the\n"
            "#    census stop for zero-carrier targets from 5 to n. Logged, and recorded in def_day.final.env.\n",
        ),
        (
            'REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1\n',
            '# 538: the repo is the directory that contains this script. On the production Mac this is\n'
            '# ~/Code/virtual-cell-challenge-2026; in a worktree it is the worktree.\n'
            'REPO="$(CDPATH= cd "$(dirname "$0")/.." && pwd)"\n'
            '[ -f "$REPO/scripts/def_day.sh" ] || { echo "STOP: derived REPO=$REPO has no scripts/def_day.sh" >&2; exit 1; }\n'
            'cd "$REPO" || exit 1\n',
        ),
        (
            'die(){ log "STOP: $*"; exit 1; }\n'
            "for v in $(env | grep -o '^VCC_[A-Z0-9_]*'); do unset \"$v\"; done\n",
            'die(){ log "STOP: $*"; exit 1; }\n'
            '# 538: the census override is a one-shot. Checked before log "start" writes def_day.log, so a fresh\n'
            '# work dir is still empty here (mkdir -p above does not count). Unset: this whole block is skipped.\n'
            'if [ -n "${DEF_DAY_CENSUS_MAX_ZERO:-}" ]; then\n'
            '  python3 -c \'import sys; n=sys.argv[1]; sys.exit(0 if n.isdigit() and not (len(n)>1 and n[0]=="0") and 6<=int(n)<=30 else 1)\' "$DEF_DAY_CENSUS_MAX_ZERO" 2>/dev/null \\\n'
            '    || die "DEF_DAY_CENSUS_MAX_ZERO must be an integer in [6, 30] with no leading zero (got $DEF_DAY_CENSUS_MAX_ZERO)"\n'
            '  if [ -d "$W" ] && [ -n "$(ls -A "$W" 2>/dev/null)" ]; then\n'
            '    die "DEF_DAY_CENSUS_MAX_ZERO=$DEF_DAY_CENSUS_MAX_ZERO is a one-shot override and needs a fresh work dir ($W is not empty)"\n'
            '  fi\n'
            'fi\n'
            "for v in $(env | grep -o '^VCC_[A-Z0-9_]*'); do unset \"$v\"; done\n",
        ),
        (
            '  grep -q -- "--hi-pct" scripts/def_lodo_key.py || die "def_lodo_key.py has no --hi-pct (497)"\n'
            "fi\n"
            "# 502: the DE-channel riders (validated here, exported for the level builds after the dump)\n",
            '  grep -q -- "--hi-pct" scripts/def_lodo_key.py || die "def_lodo_key.py has no --hi-pct (497)"\n'
            "fi\n"
            '# 538: one cap form. Unset keeps today\'s --hi default (3) or DEF_DAY_LODO_HI_PCT.\n'
            'if [ -n "${DEF_DAY_LODO_HI:-}" ]; then\n'
            '  [ "${DEF_DAY_LODO:-0}" = "1" ] || die "DEF_DAY_LODO_HI=$DEF_DAY_LODO_HI needs DEF_DAY_LODO=1"\n'
            '  [ -z "${DEF_DAY_LODO_HI_PCT:-}" ] || die "DEF_DAY_LODO_HI and DEF_DAY_LODO_HI_PCT are two cap forms; set only one (got HI=$DEF_DAY_LODO_HI PCT=$DEF_DAY_LODO_HI_PCT)"\n'
            "  python3 -c 'import re,sys\n"
            "t=sys.argv[1]\n"
            'if re.fullmatch(r"[1-9][0-9]*(?:\\.[0-9]+)?", t) is None:\n'
            "    sys.exit(1)\n"
            "h=float(t)\n"
            "sys.exit(0 if 1.0 < h <= 20.0 else 1)' \"$DEF_DAY_LODO_HI\" 2>/dev/null \\\n"
            '    || die "DEF_DAY_LODO_HI must satisfy 1 < h <= 20, written as a plain decimal (got $DEF_DAY_LODO_HI)"\n'
            '  grep -q \'add_argument("--hi"\' scripts/def_lodo_key.py || die "def_lodo_key.py has no --hi"\n'
            "fi\n"
            "# 502: the DE-channel riders (validated here, exported for the level builds after the dump)\n",
        ),
        (
            '"${DRV[@]}" --name "$N-setup" --stages preflight,identity,cd4-extract,census,whitener,key,repin-frozen,pin-artifacts,repin-panel-support \\\n'
            '    --cd4-thresholds "$CD4T" > "$W/driver-setup.log" 2>&1\n',
            '"${DRV[@]}" --name "$N-setup" --stages preflight,identity,cd4-extract,census,whitener,key,repin-frozen,pin-artifacts,repin-panel-support \\\n'
            '    --cd4-thresholds "$CD4T" ${DEF_DAY_CENSUS_MAX_ZERO:+--census-max-zero "$DEF_DAY_CENSUS_MAX_ZERO"} > "$W/driver-setup.log" 2>&1\n',
        ),
        (
            'log "setup: $(grep -E \'^--- \\[\' "$W/driver-setup.log" | awk \'{print $2, $3}\' | tr \'\\n\' \' \')"\n'
            '[ -s "$W/identity.json" ] || die "no identity.json"\n',
            'log "setup: $(grep -E \'^--- \\[\' "$W/driver-setup.log" | awk \'{print $2, $3}\' | tr \'\\n\' \' \')"\n'
            'if [ -n "${DEF_DAY_CENSUS_MAX_ZERO:-}" ]; then\n'
            '  ZC="$(python3 -c \'import json,sys; print(json.load(open(sys.argv[1]))["consensus"]["targets_with_0_donors"])\' "$W/census.json")" \\\n'
            '    || die "census override set but $W/census.json has no zero-carrier count"\n'
            '  log "census override: max zero-carrier $DEF_DAY_CENSUS_MAX_ZERO (count $ZC)"\n'
            'fi\n'
            '[ -s "$W/identity.json" ] || die "no identity.json"\n',
        ),
        (
            '  LK="$W/pert_scale_lodo_def_lo0.5_hi3.npz"; LID=lodo_lo0.5_hi3_def\n'
            '  [ -n "${DEF_DAY_LODO_HI_PCT:-}" ] && { LK="$W/pert_scale_lodo_def_lo0.5_pct.npz"; LID=lodo_lo0.5_pct_def; }\n'
            '  uv run python scripts/def_lodo_key.py --bundle "$B" --op "$W/whitener_union_pruned_def.npz" --cd4-root "$W/cd4_compact" \\\n'
            '      ${DEF_DAY_LODO_HI_PCT:+--hi-pct "$DEF_DAY_LODO_HI_PCT"} --out "$LK" >> "$LOG" 2>&1 || die "def_lodo_key.py failed"\n',
            '  LK="$W/pert_scale_lodo_def_lo0.5_hi3.npz"; LID=lodo_lo0.5_hi3_def\n'
            '  if [ -n "${DEF_DAY_LODO_HI:-}" ]; then\n'
            '    LK="$W/pert_scale_lodo_def_lo0.5_hi${DEF_DAY_LODO_HI}.npz"; LID="lodo_lo0.5_hi${DEF_DAY_LODO_HI}_def"\n'
            '  elif [ -n "${DEF_DAY_LODO_HI_PCT:-}" ]; then\n'
            '    LK="$W/pert_scale_lodo_def_lo0.5_pct.npz"; LID=lodo_lo0.5_pct_def\n'
            '  fi\n'
            '  uv run python scripts/def_lodo_key.py --bundle "$B" --op "$W/whitener_union_pruned_def.npz" --cd4-root "$W/cd4_compact" \\\n'
            '      ${DEF_DAY_LODO_HI:+--hi "$DEF_DAY_LODO_HI"} ${DEF_DAY_LODO_HI_PCT:+--hi-pct "$DEF_DAY_LODO_HI_PCT"} --out "$LK" >> "$LOG" 2>&1 || die "def_lodo_key.py failed"\n',
        ),
        (
            '_r530=""\n'
            '[ -n "${DEF_DAY_RC_CALLSET_CPM:-}" ] && _r530=1\n'
            '[ "${DEF_DAY_CD4_STACK:-x1}" != "x1" ] && _r530=1\n'
            '[ -n "${DEF_DAY_RHO_SCALE:-}" ] && _r530=1\n'
            '[ -n "${DEF_DAY_LODO_STRENGTH:-}" ] && _r530=1\n'
            '[ -n "${DEF_DAY_HOOK_ESTIMATOR:-}" ] && _r530=1\n'
            'if [ -n "$_r530" ]; then\n'
            '  _rc="${DEF_DAY_RC_CALLSET_CPM:--}"; _st="${DEF_DAY_CD4_STACK:-x1}"; _rh="${DEF_DAY_RHO_SCALE:--}"; _ls="${DEF_DAY_LODO_STRENGTH:--}"\n'
            '  uv run python scripts/def_day_helpers.py riders_match "$FINAL.h5ad.meta.json" "${VCC_DE_MAG:--}" "${VCC_DE_MDE_SIGN:--}" \\\n'
            '      "rc=${_rc},stack=${_st},rho=${_rh},strength=${_ls}" \\\n'
            '    || die "the final meta does not carry the requested October riders (rc=${_rc} stack=${_st} rho=${_rh} strength=${_ls})"\n',
            '_r530=""\n'
            '[ -n "${DEF_DAY_RC_CALLSET_CPM:-}" ] && _r530=1\n'
            '[ "${DEF_DAY_CD4_STACK:-x1}" != "x1" ] && _r530=1\n'
            '[ -n "${DEF_DAY_RHO_SCALE:-}" ] && _r530=1\n'
            '[ -n "${DEF_DAY_LODO_STRENGTH:-}" ] && _r530=1\n'
            '[ -n "${DEF_DAY_HOOK_ESTIMATOR:-}" ] && _r530=1\n'
            '[ -n "${DEF_DAY_LODO_HI:-}" ] && _r530=1\n'
            '[ -n "${DEF_DAY_CALLSET_BLOCK:-}" ] && _r530=1\n'
            'if [ -n "$_r530" ]; then\n'
            '  _rc="${DEF_DAY_RC_CALLSET_CPM:--}"; _st="${DEF_DAY_CD4_STACK:-x1}"; _rh="${DEF_DAY_RHO_SCALE:--}"; _ls="${DEF_DAY_LODO_STRENGTH:--}"\n'
            '  _spec="rc=${_rc},stack=${_st},rho=${_rh},strength=${_ls}"\n'
            '  [ -n "${DEF_DAY_LODO_HI:-}" ] && _spec="${_spec},lodo_hi=${DEF_DAY_LODO_HI}"\n'
            '  [ -n "${DEF_DAY_CALLSET_BLOCK:-}" ] && _spec="${_spec},callset=${DEF_DAY_CALLSET_BLOCK}"\n'
            '  uv run python scripts/def_day_helpers.py riders_match "$FINAL.h5ad.meta.json" "${VCC_DE_MAG:--}" "${VCC_DE_MDE_SIGN:--}" \\\n'
            '      "$_spec" \\\n'
            '    || die "the final meta does not carry the requested October riders (rc=${_rc} stack=${_st} rho=${_rh} strength=${_ls}${DEF_DAY_LODO_HI:+ lodo_hi=$DEF_DAY_LODO_HI}${DEF_DAY_CALLSET_BLOCK:+ callset=$DEF_DAY_CALLSET_BLOCK})"\n',
        ),
        (
            '[ -n "${DEF_DAY_HOOK_ESTIMATOR:-}" ] && printf \'DEF_DAY_HOOK_ESTIMATOR=%s\\n\' "$DEF_DAY_HOOK_ESTIMATOR" >> "$W/def_day.final.env"\n',
            '[ -n "${DEF_DAY_HOOK_ESTIMATOR:-}" ] && printf \'DEF_DAY_HOOK_ESTIMATOR=%s\\n\' "$DEF_DAY_HOOK_ESTIMATOR" >> "$W/def_day.final.env"\n'
            '[ -n "${DEF_DAY_LODO_HI:-}" ] && printf \'DEF_DAY_LODO_HI=%s\\n\' "$DEF_DAY_LODO_HI" >> "$W/def_day.final.env"\n'
            '[ -n "${DEF_DAY_CENSUS_MAX_ZERO:-}" ] && printf \'DEF_DAY_CENSUS_MAX_ZERO=%s\\n\' "$DEF_DAY_CENSUS_MAX_ZERO" >> "$W/def_day.final.env"\n',
        ),
    ]
    if "DEF_DAY_CALLSET_BLOCK" not in src:
        edits = callset_edits() + edits
    return _apply(src, edits, "def_day.sh")


def patch_helpers(src: str) -> str:
    edits = [
        (
            '      file:npz:sha:tau[:tau_c4[:leN]] | ->   [rc=<-|N>,stack=<x1|x1.5>,rho=<-|s>,strength=<-|theta>]\n'
            '      (502; 517 adds programme / file; 530\'s optional 4th arg checks the October meta rows; "-" = absent)\n',
            '      file:npz:sha:tau[:tau_c4[:leN]] | ->   [rc=<-|N>,stack=<x1|x1.5>,rho=<-|s>,strength=<-|theta>[,lodo_hi=<h>][,callset=<prop[:swap[:steer]]>]]\n'
            '      (502; 517 adds programme / file; 530\'s optional 4th arg checks the October meta rows; 538 adds\n'
            '      optional lodo_hi and callset. "-" = absent. A 4-key spec with no lodo_hi/callset is 530\'s spec.)\n',
        ),
        (
            "def _october530_spec(spec: str):\n"
            "    out = {}\n"
            "    for part in spec.split(\",\"):\n"
            "        k, sep, v = part.partition(\"=\")\n"
            "        if sep != \"=\" or k not in {\"rc\", \"stack\", \"rho\", \"strength\"} or v == \"\":\n"
            "            return None\n"
            "        out[k] = v\n"
            "    if set(out) != {\"rc\", \"stack\", \"rho\", \"strength\"}:\n"
            "        return None\n"
            "    return out\n",
            "def _callset538(raw: str):\n"
            '    """538: DEF_DAY_CALLSET_BLOCK grammar -> the level-meta fields riders_match requires."""\n'
            "    if raw == \"prop\":\n"
            "        return {\"source\": \"prop\", \"swap\": 1.0, \"steer\": \"kept\"}\n"
            "    p = raw.split(\":\")\n"
            "    if not p or p[0] != \"prop\" or len(p) not in (2, 3):\n"
            "        return None\n"
            "    try:\n"
            "        swap = float(p[1])\n"
            "    except ValueError:\n"
            "        return None\n"
            "    steer = p[2].lower() if len(p) == 3 else \"kept\"\n"
            "    if not (0.0 < swap <= 1.0) or steer not in {\"kept\", \"key\", \"cpm\"}:\n"
            "        return None\n"
            "    return {\"source\": \"prop\", \"swap\": swap, \"steer\": steer}\n"
            "\n"
            "\n"
            "def _october530_spec(spec: str):\n"
            "    required = {\"rc\", \"stack\", \"rho\", \"strength\"}\n"
            "    optional = {\"lodo_hi\", \"callset\"}\n"
            "    out = {}\n"
            "    for part in spec.split(\",\"):\n"
            "        k, sep, v = part.partition(\"=\")\n"
            "        if sep != \"=\" or k not in required | optional or v == \"\":\n"
            "            return None\n"
            "        out[k] = v\n"
            "    if not required <= set(out):\n"
            "        return None\n"
            "    return out\n",
        ),
        (
            "        if got_rho is not None and abs(float(got_rho) - s) >= 1e-12:\n"
            "            return f\"rho_scale meta {got_rho} vs {rho}\"\n"
            "    return None\n",
            "        if got_rho is not None and abs(float(got_rho) - s) >= 1e-12:\n"
            "            return f\"rho_scale meta {got_rho} vs {rho}\"\n"
            "    # 538: the level meta does not record the LODO cap (the pinned key's json does). A value that\n"
            "    # IS recorded must be the one requested. Absence passes. callset IS a level-meta row.\n"
            "    if \"lodo_hi\" in spec:\n"
            "        got = cfg.get(\"lodo_hi\")\n"
            "        if got is not None and str(got) != spec[\"lodo_hi\"]:\n"
            "            try:\n"
            "                same = abs(float(got) - float(spec[\"lodo_hi\"])) < 1e-12\n"
            "            except (TypeError, ValueError):\n"
            "                same = False\n"
            "            if not same:\n"
            "                return f\"lodo_hi meta {got!r} vs {spec['lodo_hi']}\"\n"
            "    if \"callset\" in spec:\n"
            "        want = _callset538(spec[\"callset\"])\n"
            "        got = cfg.get(\"rc_callset_block\")\n"
            "        if want is None:\n"
            "            return f\"callset {spec['callset']!r} is not prop[:swap[:steer]]\"\n"
            "        if not isinstance(got, dict) or got.get(\"source\") != \"prop\" or got.get(\"steer\") != want[\"steer\"]:\n"
            "            return f\"rc_callset_block {got!r} vs {spec['callset']}\"\n"
            "        sha = got.get(\"sha256\")\n"
            "        if not isinstance(sha, str) or len(sha) != 64:\n"
            "            return f\"rc_callset_block sha256 {sha!r}\"\n"
            "        try:\n"
            "            swap_ok = abs(float(got.get(\"swap\", -1.0)) - want[\"swap\"]) < 1e-12\n"
            "        except (TypeError, ValueError):\n"
            "            swap_ok = False\n"
            "        if not swap_ok:\n"
            "            return f\"rc_callset_block swap {got.get('swap')!r} vs {want['swap']}\"\n"
            "    return None\n",
        ),
    ]
    return _apply(src, edits, "def_day_helpers.py")


def patch_recipe(src: str) -> str:
    edits = [
        (
            '    "lodo_strength": "DEF_DAY_LODO_STRENGTH",\n'
            "}\n",
            '    "lodo_strength": "DEF_DAY_LODO_STRENGTH",\n'
            '    "lodo_hi": "DEF_DAY_LODO_HI",\n'
            '    "callset_block": "DEF_DAY_CALLSET_BLOCK",\n'
            "}\n",
        ),
        (
            '                 "DEF_DAY_RC_CALLSET_CPM": "VCC_RC_CALLSET_CPM"}\n'
            'SAMELINE = Path(os.environ.get("OCTOBER_SAMELINE", ROOT / "scripts/def_sameline.py"))\n',
            '                 "DEF_DAY_RC_CALLSET_CPM": "VCC_RC_CALLSET_CPM",\n'
            '                 "DEF_DAY_CALLSET_BLOCK": "VCC_RC_CALLSET_BLOCK"}\n'
            'SAMELINE = Path(os.environ.get("OCTOBER_SAMELINE", ROOT / "scripts/def_sameline.py"))\n'
            'EMITTER = Path(os.environ.get("OCTOBER_EMITTER", ROOT / "src/vcc2026/rc_emit_k.py"))\n'
            'CALLBLOCK = Path(os.environ.get("OCTOBER_CALLBLOCK", ROOT / "experiments/528-callprop-block/callblock528.py"))\n',
        ),
        (
            "    _problems_530(env, problems)\n",
            "    _problems_530(env, problems)\n"
            "    _problems_538(env, problems)\n",
        ),
        (
            "def installed_problems(env: dict, flags: str) -> list:\n"
            "    dd, bu = DEF_DAY.read_text(), BUILDER.read_text()\n"
            "    sameline = SAMELINE.read_text() if SAMELINE.is_file() else \"\"\n",
            "def _callset_ok538(raw: str) -> bool:\n"
            "    if raw == \"prop\":\n"
            "        return True\n"
            "    p = raw.split(\":\")\n"
            "    if not p or p[0] != \"prop\" or len(p) not in (2, 3):\n"
            "        return False\n"
            "    try:\n"
            "        v = float(p[1])\n"
            "    except ValueError:\n"
            "        return False\n"
            "    if not (0.0 < v <= 1.0):\n"
            "        return False\n"
            "    if len(p) == 3 and p[2].lower() not in {\"kept\", \"key\", \"cpm\"}:\n"
            "        return False\n"
            "    return True\n"
            "\n"
            "\n"
            "def _problems_538(env: dict, problems: list) -> None:\n"
            '    """Range and dependency checks for 538\'s knobs. PENDING stays out of env."""\n'
            '    if "DEF_DAY_LODO_HI" in env:\n'
            '        if env.get("DEF_DAY_LODO") != "1":\n'
            '            problems.append("lodo_hi needs lodo = 1")\n'
            '        if "DEF_DAY_LODO_HI_PCT" in env:\n'
            '            problems.append("lodo_hi and lodo_hi_pct are two cap forms; set only one")\n'
            "        import re as _re\n"
            '        t = env["DEF_DAY_LODO_HI"]\n'
            "        try:\n"
            "            h = float(t)\n"
            '            ok = _re.fullmatch(r"[1-9]\\d*(?:\\.\\d+)?", t) is not None and 1.0 < h <= 20.0\n'
            "        except ValueError:\n"
            "            ok = False\n"
            "        if not ok:\n"
            '            problems.append(f"lodo_hi {t!r} must satisfy 1 < h <= 20")\n'
            '    if "DEF_DAY_CALLSET_BLOCK" in env and not _callset_ok538(env["DEF_DAY_CALLSET_BLOCK"]):\n'
            '        problems.append(f"callset_block {env[\'DEF_DAY_CALLSET_BLOCK\']!r} must be prop or prop:<swap>[:<steer>]")\n'
            "\n"
            "\n"
            "def installed_problems(env: dict, flags: str) -> list:\n"
            "    dd, bu = DEF_DAY.read_text(), BUILDER.read_text()\n"
            "    sameline = SAMELINE.read_text() if SAMELINE.is_file() else \"\"\n",
        ),
        (
            '        "DEF_DAY_HOOK_ESTIMATOR": \'--programme-estimator "$DEF_DAY_HOOK_ESTIMATOR"\',\n'
            "    }\n",
            '        "DEF_DAY_HOOK_ESTIMATOR": \'--programme-estimator "$DEF_DAY_HOOK_ESTIMATOR"\',\n'
            '        "DEF_DAY_LODO_HI": \'--hi "$DEF_DAY_LODO_HI"\',\n'
            '        "DEF_DAY_CALLSET_BLOCK": \'export VCC_RC_CALLSET_BLOCK="$CSB"\',\n'
            "    }\n",
        ),
        (
            '    if "DEF_DAY_HOOK_ESTIMATOR" in env and "--programme-estimator" not in sameline:\n'
            '        out.append("scripts/def_sameline.py has no --programme-estimator (517)")\n',
            '    if "DEF_DAY_HOOK_ESTIMATOR" in env and "--programme-estimator" not in sameline:\n'
            '        out.append("scripts/def_sameline.py has no --programme-estimator (517)")\n'
            '    if "DEF_DAY_CALLSET_BLOCK" in env:\n'
            '        emitter = EMITTER.read_text() if EMITTER.is_file() else ""\n'
            '        if "block_steer=" not in emitter:\n'
            '            out.append("the emitter has no block_steer (528)")\n'
            "        if not CALLBLOCK.is_file():\n"
            '            out.append("experiments/528-callprop-block/callblock528.py is missing")\n',
        ),
    ]
    return _apply(src, edits, "october_recipe.py")


def patch_decisions(src: str) -> str:
    """JSON-level, same as 530: string anchors on this file broke once held_out lines landed."""
    d = json.loads(src)
    s = d["switches"]
    new = {
        "lodo_hi": {
            "status": "PENDING", "value": None,
            "decided_by": "FINAL-PHASE row 4a / 531 finding 12: a fixed cap (row 4a uses 6) if it transfers",
            "if_on": "a cap h with 1 < h <= 20",
            "otherwise": "off (cap 3, or the percentile form when lodo_hi_pct is on)",
            "needs": "DEF_DAY_LODO=1; refused together with lodo_hi_pct (one cap form only). Recorded on the pinned key json (hi) and its id lodo_lo0.5_hi<h>_def",
            "held_out": "a fixed cap still has to clear a held-out line before it is DECIDED (531 #12); the percentile form's held-out read does not stand in for it",
        },
        "callset_block": {
            "status": "PENDING", "value": None,
            "decided_by": "528's H1 price, then a board read; not DECIDED here",
            "if_on": "prop or prop:<swap>[:<steer>] with 0 < swap <= 1 and steer kept|key|cpm",
            "otherwise": "off",
            "needs": "528's pass-through (builder VCC_RC_CALLSET_BLOCK, emitter block_steer, callblock528.py). Meta config.rc_callset_block carries source, sha256, swap, steer",
            "held_out": "528 priced the block on H1; it enters the recipe only if that price is positive and a held-out read agrees",
        },
    }
    out_sw = {}
    for k, v in s.items():
        if k == "cd4_thresholds":
            out_sw.update(new)
        out_sw[k] = v
    for k, v in new.items():
        out_sw.setdefault(k, v)
    d["switches"] = out_sw
    out = json.dumps(d, indent=1, ensure_ascii=False) + "\n"
    json.loads(out)
    return out


def patch_gate(src: str) -> str:
    edits = [
        (
            '    "DEF_DAY_HOOK_ESTIMATOR": {"de_mde_sign"},\n'
            "}\n",
            '    "DEF_DAY_HOOK_ESTIMATOR": {"de_mde_sign"},\n'
            '    "DEF_DAY_LODO_HI": {"bulk_load", "bulk_load_sha256"},\n'
            '    "DEF_DAY_CALLSET_BLOCK": {"rc_callset_block"},\n'
            "}\n",
        ),
        (
            'if "DEF_DAY_HOOK_ESTIMATOR" in env:\n'
            '    est = env["DEF_DAY_HOOK_ESTIMATOR"]\n'
            '    chk("hook estimator", est == "meanz_ne" and f"estimator {est}" in dd\n'
            '        and str(env.get("DEF_DAY_MDE_SIGN", "")).startswith("hookbank:"),\n'
            '        f"log names estimator {est}; mde {env.get(\'DEF_DAY_MDE_SIGN\')}")\n',
            'if "DEF_DAY_HOOK_ESTIMATOR" in env:\n'
            '    est = env["DEF_DAY_HOOK_ESTIMATOR"]\n'
            '    chk("hook estimator", est == "meanz_ne" and f"estimator {est}" in dd\n'
            '        and str(env.get("DEF_DAY_MDE_SIGN", "")).startswith("hookbank:"),\n'
            '        f"log names estimator {est}; mde {env.get(\'DEF_DAY_MDE_SIGN\')}")\n'
            'if "DEF_DAY_LODO_HI" in env:\n'
            '    h = env["DEF_DAY_LODO_HI"]\n'
            '    side = RUN / f"pert_scale_lodo_def_lo0.5_hi{h}.json"\n'
            "    got_hi = None\n"
            "    if side.is_file():\n"
            "        try:\n"
            '            got_hi = json.loads(side.read_text()).get("hi")\n'
            "        except (OSError, json.JSONDecodeError):\n"
            "            got_hi = None\n"
            "    try:\n"
            "        hi_ok = got_hi is not None and abs(float(got_hi) - float(h)) < 1e-12\n"
            "    except (TypeError, ValueError):\n"
            "        hi_ok = False\n"
            '    chk("LODO fixed cap", hi_ok and f"lodo_lo0.5_hi{h}_def" in dd, f"sidecar hi {got_hi}")\n'
            'if "DEF_DAY_CALLSET_BLOCK" in env:\n'
            '    raw = env["DEF_DAY_CALLSET_BLOCK"]\n'
            '    parts = raw.split(":")\n'
            "    form_ok, want_swap, want_steer = False, None, None\n"
            '    if raw == "prop":\n'
            '        form_ok, want_swap, want_steer = True, 1.0, "kept"\n'
            '    elif parts and parts[0] == "prop" and len(parts) in (2, 3):\n'
            "        try:\n"
            "            want_swap = float(parts[1])\n"
            '            want_steer = parts[2].lower() if len(parts) == 3 else "kept"\n'
            '            form_ok = 0.0 < want_swap <= 1.0 and want_steer in {"kept", "key", "cpm"}\n'
            "        except ValueError:\n"
            "            form_ok = False\n"
            '    blk = c.get("rc_callset_block") or {}\n'
            "    sha = blk.get(\"sha256\")\n"
            "    try:\n"
            "        swap_ok = want_swap is not None and abs(float(blk.get(\"swap\", -1.0)) - want_swap) < 1e-12\n"
            "    except (TypeError, ValueError):\n"
            "        swap_ok = False\n"
            '    chk("call-propensity block",\n'
            '        form_ok and blk.get("source") == "prop" and isinstance(sha, str) and len(sha) == 64\n'
            '        and swap_ok and blk.get("steer") == want_steer,\n'
            '        f"meta {json.dumps(blk)[:160]} vs {raw}")\n',
        ),
    ]
    return _apply(src, edits, "gate524.py")


def patch_postbuild(src: str) -> str:
    edits = [
        (
            '    ap.add_argument("--allow-config-diff", default="",\n'
            '                    help="comma-separated extra config keys this build INTENDS to move")\n',
            '    ap.add_argument("--allow-config-diff", default="",\n'
            '                    help="comma-separated extra config keys this build INTENDS to move")\n'
            '    ap.add_argument("--require-config-diff", default="",\n'
            '                    help="538: comma-separated config keys that MUST move against the incumbent. "\n'
            '                         "Fails unless every listed key moved AND nothing outside the allow-list "\n'
            '                         "(built-in allows, --allow-config-diff, the new-panel keys) moved. "\n'
            '                         "Empty, the default, adds no check and does not change today\'s output")\n',
        ),
        (
            '    chk(f"config differences vs {args.incumbent.name} are only {sorted(allowed)}",\n'
            '        not unexpected, f"moved={moved}")\n'
            "\n"
            "    # --- the emitted artifact ------------------------------------------------------------\n",
            '    chk(f"config differences vs {args.incumbent.name} are only {sorted(allowed)}",\n'
            '        not unexpected, f"moved={moved}")\n'
            "    # 538: empty --require-config-diff adds nothing, so the checks above and the exit code\n"
            "    # stay today's. A listed key is not auto-allowed: --allow-config-diff keeps its meaning.\n"
            '    required = [k.strip() for k in args.require_config_diff.split(",") if k.strip()]\n'
            "    if required:\n"
            "        missing = [k for k in required if k not in moved]\n"
            "        outside = [k for k in moved if k not in allowed]\n"
            '        chk("every --require-config-diff key moved and nothing outside the allow-list did",\n'
            "            not missing and not outside,\n"
            '            f"required={required} missing={missing} moved={moved} outside={outside}")\n'
            "\n"
            "    # --- the emitted artifact ------------------------------------------------------------\n",
        ),
    ]
    return _apply(src, edits, "def_postbuild.py")


def patch_pipeline(src: str) -> str:
    edits = [
        (
            '    ap.add_argument("--cd4-thresholds", choices=["default", "relaxed"], default="default",\n'
            '                    help="047\'s gates (275 per condition, 290 union) or none. Relax ONLY for a synthetic panel")\n',
            '    ap.add_argument("--cd4-thresholds", choices=["default", "relaxed"], default="default",\n'
            '                    help="047\'s gates (275 per condition, 290 union) or none. Relax ONLY for a synthetic panel")\n'
            '    ap.add_argument("--census-max-zero", default=None,\n'
            '                    help="538: forward --max-zero-carrier N to def_census.py (integer 6..30). "\n'
            '                         "Omitted, the census keeps its floor of 5 and this argument adds no argv word")\n',
        ),
        (
            '                r.run(stage, [*py, "scripts/def_census.py", "--bundle", str(bundle),\n'
            '                              *(["--verify"] if is_abc else []),\n'
            '                              *(["--cd4-root", str(cd4_root)] if cd4_root else []),\n'
            '                              "--out", str(work / "census.json")])\n',
            "                zero_arg: list[str] = []\n"
            "                if args.census_max_zero is not None:\n"
            "                    tok = str(args.census_max_zero)\n"
            '                    if not (tok.isdigit() and not (len(tok) > 1 and tok[0] == "0") and 6 <= int(tok) <= 30):\n'
            '                        raise SystemExit(f"--census-max-zero must be an integer in [6, 30] (got {args.census_max_zero})")\n'
            '                    zero_arg = ["--max-zero-carrier", str(int(tok))]\n'
            '                r.run(stage, [*py, "scripts/def_census.py", "--bundle", str(bundle),\n'
            '                              *(["--verify"] if is_abc else []),\n'
            '                              *(["--cd4-root", str(cd4_root)] if cd4_root else []),\n'
            "                              *zero_arg,\n"
            '                              "--out", str(work / "census.json")])\n',
        ),
    ]
    return _apply(src, edits, "def_pipeline.py")


def build(paths: dict[str, Path]) -> dict[str, str]:
    texts = {k: p.read_text() for k, p in paths.items()}
    marks = _markers(texts)
    if all(marks.values()):
        raise SystemExit("already patched")
    if any(marks.values()):
        have = ", ".join(k for k, v in marks.items() if v)
        missing = ", ".join(k for k, v in marks.items() if not v)
        raise SystemExit(f"partial extras (have {have}; missing {missing}); restore and re-apply")
    out = {
        "def_day": patch_def_day(texts["def_day"]),
        "helpers": patch_helpers(texts["helpers"]),
        "recipe": patch_recipe(texts["recipe"]),
        "decisions": patch_decisions(texts["decisions"]),
        "gate": patch_gate(texts["gate"]),
        "postbuild": patch_postbuild(texts["postbuild"]),
        "pipeline": patch_pipeline(texts["pipeline"]),
    }
    if not all(_markers(out).values()):
        raise SystemExit(f"patch did not land every marker: {_markers(out)}")
    for key in ("helpers", "recipe", "gate", "postbuild", "pipeline"):
        compile(out[key], str(paths[key]), "exec")
    json.loads(out["decisions"])
    with tempfile.NamedTemporaryFile("w", suffix=".sh", delete=False) as fh:
        fh.write(out["def_day"])
        tmp = fh.name
    try:
        r = subprocess.run(["bash", "-n", tmp], capture_output=True, text=True)
        if r.returncode != 0:
            raise SystemExit(f"bash -n failed: {r.stderr}")
    finally:
        os.unlink(tmp)
    return out


def resolve_paths(a) -> dict[str, Path]:
    root = a.root
    rel = {
        "def_day": Path("scripts/def_day.sh"),
        "helpers": Path("scripts/def_day_helpers.py"),
        "recipe": Path("scripts/october_recipe.py"),
        "decisions": Path("docs/october-decisions.json"),
        "gate": Path("experiments/524-october-recipe-dr/gate524.py"),
        "postbuild": Path("scripts/def_postbuild.py"),
        "pipeline": Path("scripts/def_pipeline.py"),
    }
    if root is None:
        raise SystemExit("pass --root")
    return {k: Path(root) / p for k, p in rel.items()}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args(argv)
    paths = resolve_paths(a)
    out = build(paths)
    if a.check:
        print("patches apply cleanly (not written)")
        return 0
    temps = []
    for key, text in out.items():
        dest = paths[key]
        tmp = dest.with_name(dest.name + ".538tmp")
        tmp.write_text(text)
        temps.append((tmp, dest))
    for tmp, dest in temps:
        os.replace(tmp, dest)
    print("patched " + ", ".join(str(p) for p in paths.values()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
