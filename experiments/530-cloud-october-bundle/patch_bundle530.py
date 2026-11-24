"""530: October pass-through bundle, as one string-anchored patch. NOT APPLIED by this file's own tests.

Knobs (each default off; unset keeps today's exports and today's argv):

  DEF_DAY_RHO_SCALE=<s>             0.9 <= s <= 1.1. The re-level loop targets rho* x s
                                    (relevel_iter.py --rho-scale; inband / best_safe use the same scale).
  DEF_DAY_CD4_STACK=x1|x1.5         x1 is today's dump (H1 x2 / CD4 at DEF_ENV's g, file dump_h1x2_cd4x1.npz).
                                    x1.5 is 457's one-off: VCC_CD4_TARGET_G=0.5165482917974462 inside the
                                    dump subshell only (1.5 x 0.3443655278649641). Level builds keep DEF_ENV's g.
  DEF_DAY_RC_CALLSET_CPM=<n>        100 <= n <= 5000. Exports VCC_RC_CALLSET_CPM on the level builds only
                                    (builder knob from 523; unset leaves the builder's 2200).
  DEF_DAY_MDE_SIGN=hookbank:<bank>:<tau>[:<tau_c4>[:le<N>]]
  DEF_DAY_HOOK_ESTIMATOR=meanz_ne   517's Arc-sourced hook (patch_defday517.py + patch_sameline517.py),
                                    folded in here. The shell rewrites hookbank to
                                    VCC_DE_MDE_SIGN=file:<npz>:<sha>:<tau>... The meta records source "file".
  DEF_DAY_LODO_STRENGTH=<theta>     522's rider: pct key, then w' = max(floor, w x exp(theta z_src)) above the
                                    floor. Needs DEF_DAY_LODO=1 and DEF_DAY_LODO_HI_PCT.

Also: a PENDING switch per new knob in docs/october-decisions.json (cd4_stack stays DECIDED x1, which is the
default and is not exported), the mapping in scripts/october_recipe.py, a riders_match spec for the keys the
final meta actually carries, and value rows in experiments/524-october-recipe-dr/gate524.py.

526's 6-line "builder is not HEAD" insertion sits after the BSHA0= line. None of the anchors below touch that
gap, so the same patch applies whether or not that insertion is already present.

Idempotent: a second run exits before writing. --check parses and syntax-checks and does not replace the targets.
Writes are temp files then os.replace, after every temp is on disk.

usage:
  uv run python experiments/530-cloud-october-bundle/patch_bundle530.py --root <repo> [--check]
  uv run python experiments/530-cloud-october-bundle/patch_bundle530.py \\
      --def-day <sh> --helpers <py> --relevel <py> --sameline <py> \\
      --recipe <py> --decisions <json> --gate <py> [--check]
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

X1_G = "0.3443655278649641"
X15_G = "0.5165482917974462"  # 1.5 * X1_G, the literal 457's def_day_cd4x15.sh exports


def _once(src: str, old: str, name: str) -> None:
    n = src.count(old)
    if n != 1:
        raise SystemExit(f"{name}: anchor count {n} != 1: {old[:120]!r}")


def _apply(src: str, edits: list[tuple[str, str]], name: str) -> str:
    for old, _new in edits:
        _once(src, old, name)
    for old, new in edits:
        _once(src, old, name)
        src = src.replace(old, new, 1)
    return src


def _markers(texts: dict[str, str]) -> dict[str, bool]:
    return {
        "def_day": "DEF_DAY_RHO_SCALE" in texts["def_day"],
        "helpers": "def _october530_spec" in texts["helpers"],
        "relevel": "--rho-scale" in texts["relevel"],
        "sameline": "--programme-estimator" in texts["sameline"],
        "recipe": '"rho_scale": "DEF_DAY_RHO_SCALE"' in texts["recipe"],
        "decisions": '"rho_scale"' in texts["decisions"],
        "gate": "def _dump_name530" in texts["gate"],
    }


def patch_def_day(src: str) -> str:
    edits = [
        (
            "#   DEF_DAY_MDE_SIGN=clean4:<tau> (502, default off; e.g. clean4:0.05): 496's uncalled-gene sign hook (VCC_DE_MDE_SIGN).\n",
            "#   DEF_DAY_MDE_SIGN=clean4:<tau> (502, default off; e.g. clean4:0.05): 496's uncalled-gene sign hook (VCC_DE_MDE_SIGN).\n"
            "#    517: or programme:<tau>[:<tau_c4>] (the hook reads the programme of DEF_DAY_PROGRAMME_SOURCE's artifact and follows\n"
            "#    DEF_DAY_PRIOR_GATE), or hookbank:<bank>:<tau>[:<tau_c4>[:le<N>]] (the hook ALONE from <bank>'s programme, CLEAN4\n"
            "#    block: a second def_sameline.py artifact, exported as VCC_DE_MDE_SIGN=file:<npz>:<sha>:...). Builder >= 517's patch.\n"
            "#    DEF_DAY_HOOK_ESTIMATOR=meanz_ne (517, default pdown): the hook artifact's estimator (def_sameline.py >= 517's patch).\n"
            "#   DEF_DAY_RHO_SCALE=<s> (530, default off; 0.9 <= s <= 1.1): the re-level loop targets rho* x s.\n"
            "#   DEF_DAY_CD4_STACK=x1|x1.5 (530, default x1): x1.5 exports VCC_CD4_TARGET_G inside the dump subshell only (457).\n"
            "#   DEF_DAY_RC_CALLSET_CPM=<n> (530, default off; 100..5000): VCC_RC_CALLSET_CPM on the level builds (523; unset = 2200).\n"
            "#   DEF_DAY_LODO_STRENGTH=<theta> (530, default off): 522's pct x exp(theta z_src) rider. Needs LODO and the percentile cap.\n",
        ),
        (
            'if [ -n "${DEF_DAY_MDE_SIGN:-}" ]; then\n'
            '  case "$DEF_DAY_MDE_SIGN" in clean4:*) ;; *) die "DEF_DAY_MDE_SIGN must be clean4:<tau> (got $DEF_DAY_MDE_SIGN)" ;; esac\n'
            '  grep -q "VCC_DE_MDE_SIGN" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_DE_MDE_SIGN (496p)"\n'
            "fi\n",
            'if [ -n "${DEF_DAY_MDE_SIGN:-}" ]; then\n'
            '  grep -q "VCC_DE_MDE_SIGN" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_DE_MDE_SIGN (496p)"\n'
            '  case "$DEF_DAY_MDE_SIGN" in\n'
            "    clean4:*) ;;\n"
            '    programme:*) [ -n "${DEF_DAY_PROGRAMME_SOURCE:-}" ] || die "DEF_DAY_MDE_SIGN=$DEF_DAY_MDE_SIGN reads the programme artifact: it needs DEF_DAY_PROGRAMME_SOURCE"\n'
            '                 grep -q "programme:<tau>" scripts/make_submission_pm_cd4.py || die "the builder has no programme hook form (517)" ;;\n'
            "    hookbank:*) python3 -c 'import re,sys\n"
            "p=sys.argv[1].split(\":\")\n"
            "ok=len(p) in (3,4,5) and bool(p[1])\n"
            "try:\n"
            "    t=float(p[2]); c=float(p[3]) if len(p)>=4 else 0.05\n"
            "    ok = ok and 0.0 <= t < 0.5 and 0.0 <= c < 0.5\n"
            "    if len(p)==5: ok = ok and re.fullmatch(r\"le[0-9]{1,2}\", p[4]) is not None\n"
            "except Exception:\n"
            "    ok=False\n"
            "sys.exit(0 if ok else 1)' \"$DEF_DAY_MDE_SIGN\" 2>/dev/null \\\n"
            '                 || die "DEF_DAY_MDE_SIGN hookbank form must be hookbank:<bank>:<tau>[:<tau_c4>[:le<N>]] (got $DEF_DAY_MDE_SIGN)"\n'
            '                 grep -q "file:<npz>" scripts/make_submission_pm_cd4.py || die "the builder has no file hook form (517)" ;;\n'
            '    *) die "DEF_DAY_MDE_SIGN must be clean4:<tau>, programme:<tau>[:<tau_c4>] or hookbank:<bank>:<tau>[:<tau_c4>[:le<N>]] (got $DEF_DAY_MDE_SIGN)" ;;\n'
            "  esac\n"
            "fi\n",
        ),
        (
            'SQ="$Q"\n',
            "# 530: October pass-throughs. Unset = no export and no extra argument (today's command).\n"
            'if [ -n "${DEF_DAY_RHO_SCALE:-}" ]; then\n'
            "  python3 -c 'import sys; s=float(sys.argv[1]); sys.exit(0 if 0.9<=s<=1.1 else 1)' \"$DEF_DAY_RHO_SCALE\" 2>/dev/null \\\n"
            '    || die "DEF_DAY_RHO_SCALE must be a scale in [0.9, 1.1] (got $DEF_DAY_RHO_SCALE)"\n'
            "fi\n"
            'if [ -n "${DEF_DAY_RC_CALLSET_CPM:-}" ]; then\n'
            "  python3 -c 'import sys; n=sys.argv[1]; sys.exit(0 if n.isdigit() and 100<=int(n)<=5000 else 1)' \"$DEF_DAY_RC_CALLSET_CPM\" 2>/dev/null \\\n"
            '    || die "DEF_DAY_RC_CALLSET_CPM must be an integer in [100, 5000] (got $DEF_DAY_RC_CALLSET_CPM)"\n'
            '  grep -q "VCC_RC_CALLSET_CPM" scripts/make_submission_pm_cd4.py || die "the builder has no VCC_RC_CALLSET_CPM (523)"\n'
            "fi\n"
            'case "${DEF_DAY_CD4_STACK:-x1}" in\n'
            "  x1|x1.5) ;;\n"
            '  *) die "DEF_DAY_CD4_STACK must be x1 or x1.5 (got ${DEF_DAY_CD4_STACK})" ;;\n'
            "esac\n"
            'if [ -n "${DEF_DAY_LODO_STRENGTH:-}" ]; then\n'
            '  [ "${DEF_DAY_LODO:-0}" = "1" ] || die "DEF_DAY_LODO_STRENGTH=$DEF_DAY_LODO_STRENGTH needs DEF_DAY_LODO=1"\n'
            '  [ -n "${DEF_DAY_LODO_HI_PCT:-}" ] || die "DEF_DAY_LODO_STRENGTH=$DEF_DAY_LODO_STRENGTH needs the percentile cap (DEF_DAY_LODO_HI_PCT)"\n'
            "  python3 -c 'import re,sys\n"
            "t=sys.argv[1]\n"
            "if re.fullmatch(r\"-?(?:0|[1-9][0-9]*)(?:\\.[0-9]+)?\", t) is None:\n"
            "    sys.exit(1)\n"
            "v=float(t)\n"
            "sys.exit(0 if v != 0.0 and abs(v) <= 1.0 else 1)' \"$DEF_DAY_LODO_STRENGTH\" 2>/dev/null \\\n"
            '    || die "DEF_DAY_LODO_STRENGTH must be a non-zero theta in [-1, 1], written as a plain decimal (got $DEF_DAY_LODO_STRENGTH)"\n'
            "fi\n"
            'if [ -n "${DEF_DAY_HOOK_ESTIMATOR:-}" ]; then\n'
            '  [ "$DEF_DAY_HOOK_ESTIMATOR" = "meanz_ne" ] || die "DEF_DAY_HOOK_ESTIMATOR must be meanz_ne (got $DEF_DAY_HOOK_ESTIMATOR)"\n'
            '  case "${DEF_DAY_MDE_SIGN:-}" in\n'
            "    hookbank:*) ;;\n"
            '    *) die "DEF_DAY_HOOK_ESTIMATOR=$DEF_DAY_HOOK_ESTIMATOR needs DEF_DAY_MDE_SIGN=hookbank:<bank>:<tau>... (got ${DEF_DAY_MDE_SIGN:-})" ;;\n'
            "  esac\n"
            '  grep -q -- "--programme-estimator" scripts/def_sameline.py || die "def_sameline.py has no --programme-estimator (517)"\n'
            "fi\n"
            'SQ="$Q"\n',
        ),
        (
            "# 2b. DEF_DAY_LODO=1 (460; default off)",
            "# 2a. 517: DEF_DAY_MDE_SIGN=hookbank:<bank>:... -- the hook's own programme artifact (the block prior is untouched)\n"
            'MDS="${DEF_DAY_MDE_SIGN:-}"\n'
            'case "$MDS" in hookbank:*)\n'
            '  HB="${MDS#hookbank:}"; HBANK="${HB%%:*}"; HREST="${HB#*:}"\n'
            '  [ -n "$HBANK" ] && [ "$HREST" != "$HB" ] || die "DEF_DAY_MDE_SIGN=$MDS: expected hookbank:<bank>:<tau>[...]"\n'
            '  uv run python scripts/def_sameline.py --bundle "$B" --identity "$W/identity.json" --out "$W/hookprog.npz" \\\n'
            '      ${DEF_DAY_FORCE:+--force "$DEF_DAY_FORCE"} --programme-source "$HBANK" \\\n'
            '      ${DEF_DAY_HOOK_ESTIMATOR:+--programme-estimator "$DEF_DAY_HOOK_ESTIMATOR"} >> "$LOG" 2>&1 \\\n'
            '      || die "def_sameline.py --programme-source $HBANK (the hook artifact) failed"\n'
            '  MDS="file:$W/hookprog.npz:$(shasum -a 256 "$W/hookprog.npz" | cut -d\' \' -f1):$HREST"\n'
            '  log "hook programme artifact from $HBANK (estimator ${DEF_DAY_HOOK_ESTIMATOR:-pdown}): VCC_DE_MDE_SIGN=file:${MDS#file:$W/}" ;;\n'
            "esac\n\n"
            "# 2b. DEF_DAY_LODO=1 (460; default off)",
        ),
        (
            '      ${DEF_DAY_LODO_HI_PCT:+--hi-pct "$DEF_DAY_LODO_HI_PCT"} --out "$LK" >> "$LOG" 2>&1 || die "def_lodo_key.py failed"\n'
            '  uv run python scripts/def_pin_panel_artifacts.py --whitener "$W/whitener_union_pruned_def.npz" \\\n',
            '      ${DEF_DAY_LODO_HI_PCT:+--hi-pct "$DEF_DAY_LODO_HI_PCT"} --out "$LK" >> "$LOG" 2>&1 || die "def_lodo_key.py failed"\n'
            '  LJSON="${LK%.npz}.json"\n'
            "  # 530: 522's strength rider, between the LODO key and the pin. Unset: LK stays the key just written.\n"
            '  if [ -n "${DEF_DAY_LODO_STRENGTH:-}" ]; then\n'
            '    SK="$W/pert_scale_lodo_strength.npz"\n'
            '    uv run python experiments/522-october-bulk-key/def_strength_key522.py --in "$LK" --bundle "$B" \\\n'
            '        --theta "$DEF_DAY_LODO_STRENGTH" --out "$SK" >> "$LOG" 2>&1 || die "def_strength_key522.py failed"\n'
            '    LK="$SK"\n'
            '    log "LODO strength rider: theta $DEF_DAY_LODO_STRENGTH"\n'
            "  fi\n"
            '  uv run python scripts/def_pin_panel_artifacts.py --whitener "$W/whitener_union_pruned_def.npz" \\\n',
        ),
        (
            'print(json.load(open(sys.argv[1]))["hi"])\' "${LK%.npz}.json")"\n',
            'print(json.load(open(sys.argv[1]))["hi"])\' "$LJSON")"\n',
        ),
        (
            'DUMP="$W/dump_h1x2_cd4x1.npz"\n'
            'rm -f "$DUMP"\n'
            "(\n"
            '  export VCC_H1_UNION_G=0.56 VCC_BULK_DUMP="$DUMP"\n'
            '  [ -n "$LODO" ] && export VCC_PERT_SCALE="$LODO"\n'
            '  "${DRV[@]}" --name "$N-dump" --stages build > "$W/driver-dump.log" 2>&1\n'
            ")\n",
            "# 530: x1 (unset) keeps today's filename and does not export VCC_CD4_TARGET_G (DEF_ENV's g).\n"
            "# x1.5 is 457: 1.5 x 0.3443655278649641 inside this subshell only. A strength rider renames the\n"
            "# file so the level meta's bulk_load carries the theta; unset leaves dump_h1x2_cd4x1.npz.\n"
            'case "${DEF_DAY_CD4_STACK:-x1}" in\n'
            '  x1.5) DUMP="$W/dump_h1x2_cd4x15.npz" ;;\n'
            '  *) DUMP="$W/dump_h1x2_cd4x1.npz" ;;\n'
            "esac\n"
            'CD4G=""\n'
            '[ "${DEF_DAY_CD4_STACK:-x1}" = "x1.5" ] && CD4G="' + X15_G + '"\n'
            'if [ -n "${DEF_DAY_LODO_STRENGTH:-}" ]; then\n'
            '  STOKEN="$(printf \'%s\' "$DEF_DAY_LODO_STRENGTH" | tr \'-\' \'m\')"\n'
            '  DUMP="${DUMP%.npz}_s${STOKEN}.npz"\n'
            "fi\n"
            '[ -n "$CD4G" ] && log "CD4 stack x1.5: dump ${DUMP##*/} VCC_CD4_TARGET_G=$CD4G (level builds keep DEF_ENV)"\n'
            'rm -f "$DUMP"\n'
            "(\n"
            '  export VCC_H1_UNION_G=0.56 VCC_BULK_DUMP="$DUMP"\n'
            '  [ -n "$CD4G" ] && export VCC_CD4_TARGET_G="$CD4G"\n'
            '  [ -n "$LODO" ] && export VCC_PERT_SCALE="$LODO"\n'
            '  "${DRV[@]}" --name "$N-dump" --stages build > "$W/driver-dump.log" 2>&1\n'
            ")\n",
        ),
        (
            '[ -n "${DEF_DAY_MDE_SIGN:-}" ] && { export VCC_DE_MDE_SIGN="$DEF_DAY_MDE_SIGN"; log "DE rider for the level builds: VCC_DE_MDE_SIGN=$VCC_DE_MDE_SIGN"; }\n',
            '[ -n "${DEF_DAY_MDE_SIGN:-}" ] && { export VCC_DE_MDE_SIGN="$MDS"; log "DE rider for the level builds: VCC_DE_MDE_SIGN=$VCC_DE_MDE_SIGN"; }\n'
            '[ -n "${DEF_DAY_RC_CALLSET_CPM:-}" ] && { export VCC_RC_CALLSET_CPM="$DEF_DAY_RC_CALLSET_CPM"; log "DE rider for the level builds: VCC_RC_CALLSET_CPM=$VCC_RC_CALLSET_CPM"; }\n',
        ),
        (
            '  OUT="$(uv run python experiments/448-def-rehearsal-q175/relevel_iter.py "${POINTS[@]}" 2>&1)"\n',
            '  [ -n "${DEF_DAY_RHO_SCALE:-}" ] && log "re-level rho-scale=$DEF_DAY_RHO_SCALE"\n'
            '  OUT="$(uv run python experiments/448-def-rehearsal-q175/relevel_iter.py "${POINTS[@]}" ${DEF_DAY_RHO_SCALE:+--rho-scale "$DEF_DAY_RHO_SCALE"} 2>&1)"\n',
        ),
        (
            'uv run python scripts/def_day_helpers.py inband "$W/$N-lvlC.ntclevel.json"; then\n',
            'uv run python scripts/def_day_helpers.py inband "$W/$N-lvlC.ntclevel.json" ${DEF_DAY_RHO_SCALE:+"$DEF_DAY_RHO_SCALE"}; then\n',
        ),
        (
            'FINAL="$(uv run python scripts/def_day_helpers.py best_safe "${POINTS[@]}")" || die "best_safe failed"\n',
            'FINAL="$(uv run python scripts/def_day_helpers.py best_safe ${DEF_DAY_RHO_SCALE:+--scale "$DEF_DAY_RHO_SCALE"} "${POINTS[@]}")" || die "best_safe failed"\n',
        ),
        (
            'uv run python scripts/def_day_helpers.py riders_match "$FINAL.h5ad.meta.json" "${VCC_DE_MAG:--}" "${VCC_DE_MDE_SIGN:--}" \\\n'
            '  || die "the final meta does not carry the requested DE riders (VCC_DE_MAG=${VCC_DE_MAG:-} VCC_DE_MDE_SIGN=${VCC_DE_MDE_SIGN:-})"\n',
            "# 530's spec is appended only when one of its knobs is set, so the unset argv is today's two-argument call.\n"
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
            '    || die "the final meta does not carry the requested October riders (rc=${_rc} stack=${_st} rho=${_rh} strength=${_ls})"\n'
            "else\n"
            '  uv run python scripts/def_day_helpers.py riders_match "$FINAL.h5ad.meta.json" "${VCC_DE_MAG:--}" "${VCC_DE_MDE_SIGN:--}" \\\n'
            '    || die "the final meta does not carry the requested DE riders (VCC_DE_MAG=${VCC_DE_MAG:-} VCC_DE_MDE_SIGN=${VCC_DE_MDE_SIGN:-})"\n'
            "fi\n",
        ),
        (
            '[ -n "${VCC_DE_MDE_SIGN:-}" ] && printf \'VCC_DE_MDE_SIGN=%s\\n\' "$VCC_DE_MDE_SIGN" >> "$W/def_day.final.env"\n',
            '[ -n "${VCC_DE_MDE_SIGN:-}" ] && printf \'VCC_DE_MDE_SIGN=%s\\n\' "$VCC_DE_MDE_SIGN" >> "$W/def_day.final.env"\n'
            '[ -n "${VCC_RC_CALLSET_CPM:-}" ] && printf \'VCC_RC_CALLSET_CPM=%s\\n\' "$VCC_RC_CALLSET_CPM" >> "$W/def_day.final.env"\n'
            '[ -n "${DEF_DAY_RHO_SCALE:-}" ] && printf \'DEF_DAY_RHO_SCALE=%s\\n\' "$DEF_DAY_RHO_SCALE" >> "$W/def_day.final.env"\n'
            '[ "${DEF_DAY_CD4_STACK:-x1}" != "x1" ] && printf \'DEF_DAY_CD4_STACK=%s\\n\' "$DEF_DAY_CD4_STACK" >> "$W/def_day.final.env"\n'
            '[ -n "${DEF_DAY_LODO_STRENGTH:-}" ] && printf \'DEF_DAY_LODO_STRENGTH=%s\\n\' "$DEF_DAY_LODO_STRENGTH" >> "$W/def_day.final.env"\n'
            '[ -n "${DEF_DAY_HOOK_ESTIMATOR:-}" ] && printf \'DEF_DAY_HOOK_ESTIMATOR=%s\\n\' "$DEF_DAY_HOOK_ESTIMATOR" >> "$W/def_day.final.env"\n',
        ),
        (
            '[ -n "${VCC_DE_MDE_SIGN:-}" ] && ALLOW="$ALLOW,de_mde_sign"\n',
            '[ -n "${VCC_DE_MDE_SIGN:-}" ] && ALLOW="$ALLOW,de_mde_sign"\n'
            '[ -n "${VCC_RC_CALLSET_CPM:-}" ] && ALLOW="$ALLOW,rc_callset_cpm"\n',
        ),
    ]
    return _apply(src, edits, "def_day.sh")


def patch_helpers(src: str) -> str:
    edits = [
        (
            "  def_day_helpers.py riders_match <meta.json> <de_mag variant | -> <clean4:tau | ->   (502; \"-\" = must be absent)\n",
            "  def_day_helpers.py riders_match <meta.json> <de_mag variant | -> <clean4:tau | programme:tau[:tau_c4] |\n"
            "      file:npz:sha:tau[:tau_c4[:leN]] | ->   [rc=<-|N>,stack=<x1|x1.5>,rho=<-|s>,strength=<-|theta>]\n"
            "      (502; 517 adds programme / file; 530's optional 4th arg checks the October meta rows; \"-\" = absent)\n",
        ),
        (
            "def rho_star(sd: float = 0.01) -> float:\n"
            '    """Match relevel_iter.py: q = 0.31 / (0.31 + 1.78), rho* = exp(2 * sd * Phi^-1(q))."""\n'
            "    from scipy.stats import norm\n"
            "    q = 0.31 / (0.31 + 1.78)\n"
            "    return math.exp(2.0 * sd * float(norm.ppf(q)))\n",
            "def rho_star(sd: float = 0.01, scale: float = 1.0) -> float:\n"
            '    """Match relevel_iter.py: q = 0.31 / (0.31 + 1.78), rho* = exp(2 * sd * Phi^-1(q)).\n'
            "    scale 1 (the default, and DEF_DAY_RHO_SCALE unset) returns that value with no multiply.\n"
            '    """\n'
            "    from scipy.stats import norm\n"
            "    q = 0.31 / (0.31 + 1.78)\n"
            "    base = math.exp(2.0 * sd * float(norm.ppf(q)))\n"
            "    return base if scale == 1.0 else base * scale\n",
        ),
        (
            "            rho = _rhos(args[0])\n"
            "            star = rho_star()\n"
            "            ok = bool(rho) and all(\n",
            "            rho = _rhos(args[0])\n"
            "            scale = float(args[1]) if len(args) > 1 else 1.0\n"
            "            star = rho_star(scale=scale)\n"
            "            ok = bool(rho) and all(\n",
        ),
        (
            "        star = rho_star()\n"
            "        best = None\n"
            "        for spec in args:\n",
            "        scale = 1.0\n"
            "        specs = list(args)\n"
            "        if specs and specs[0] == \"--scale\":\n"
            "            scale = float(specs[1])\n"
            "            specs = specs[2:]\n"
            "        star = rho_star(scale=scale)\n"
            "        best = None\n"
            "        for spec in specs:\n",
        ),
        (
            "def main() -> int:\n",
            "X1_G = 0.3443655278649641  # DEF_ENV's VCC_CD4_TARGET_G; the level build records this even when the dump is x1.5\n"
            "\n"
            "\n"
            "def _mde_sign_request(sign: str):\n"
            '    """517: the builder\'s VCC_DE_MDE_SIGN grammar -> the fields its config dict must carry (None = unparseable)."""\n'
            "    p = [x.strip() for x in sign.split(\":\")]\n"
            "    try:\n"
            "        if p[0] == \"clean4\" and len(p) == 2:\n"
            "            return {\"source\": \"clean4\", \"tau\": float(p[1])}\n"
            "        if p[0] == \"programme\" and len(p) in (2, 3):\n"
            "            return {\"source\": \"programme\", \"tau\": float(p[1]), \"tau_c4\": float(p[2]) if len(p) == 3 else 0.05}\n"
            "        if p[0] == \"file\" and len(p) in (4, 5, 6):\n"
            "            out = {\"source\": \"file\", \"sha256\": p[2], \"tau\": float(p[3]), \"tau_c4\": float(p[4]) if len(p) >= 5 else 0.05}\n"
            "            if len(p) == 6:\n"
            "                out[\"le\"] = float(p[5][2:])\n"
            "            return out\n"
            "    except (ValueError, IndexError):\n"
            "        return None\n"
            "    return None\n"
            "\n"
            "\n"
            "def _dump_name530(stack: str, strength: str) -> str:\n"
            '    """The stack-dump filename def_day.sh writes. Unset stack and strength -> today\'s name."""\n'
            "    base = \"dump_h1x2_cd4x15\" if stack == \"x1.5\" else \"dump_h1x2_cd4x1\"\n"
            "    if strength not in (\"-\", \"\"):\n"
            "        base += \"_s\" + strength.replace(\"-\", \"m\")\n"
            "    return base + \".npz\"\n"
            "\n"
            "\n"
            "def _october530_spec(spec: str):\n"
            "    out = {}\n"
            "    for part in spec.split(\",\"):\n"
            "        k, sep, v = part.partition(\"=\")\n"
            "        if sep != \"=\" or k not in {\"rc\", \"stack\", \"rho\", \"strength\"} or v == \"\":\n"
            "            return None\n"
            "        out[k] = v\n"
            "    if set(out) != {\"rc\", \"stack\", \"rho\", \"strength\"}:\n"
            "        return None\n"
            "    return out\n"
            "\n"
            "\n"
            "def _october530_mismatch(cfg: dict, spec: dict) -> str | None:\n"
            "    rc, stack, rho, strength = spec[\"rc\"], spec[\"stack\"], spec[\"rho\"], spec[\"strength\"]\n"
            "    if stack not in (\"x1\", \"x1.5\"):\n"
            "        return f\"stack {stack!r} is not x1 or x1.5\"\n"
            "    want_name = _dump_name530(stack, strength)\n"
            "    if cfg.get(\"bulk_load\") != want_name:\n"
            "        return f\"bulk_load {cfg.get('bulk_load')!r} vs {want_name}\"\n"
            "    g = cfg.get(\"cd4_target_g\")\n"
            "    if g is None or abs(float(g) - X1_G) >= 1e-12:\n"
            "        return f\"cd4_target_g {g!r} vs {X1_G} (the level build keeps DEF_ENV even when the dump is x1.5)\"\n"
            "    got_cpm = cfg.get(\"rc_callset_cpm\")\n"
            "    want_cpm = 2200 if rc == \"-\" else int(rc) if rc.isdigit() else None\n"
            "    if want_cpm is None or got_cpm is None or int(got_cpm) != want_cpm:\n"
            "        return f\"rc_callset_cpm {got_cpm!r} vs {rc}\"\n"
            "    got_rho = cfg.get(\"rho_scale\")\n"
            "    if rho == \"-\":\n"
            "        if got_rho is not None:\n"
            "            return f\"rho_scale present but not requested: {got_rho}\"\n"
            "    else:\n"
            "        s = float(rho)\n"
            "        if not 0.9 <= s <= 1.1:\n"
            "            return f\"rho scale {rho} is outside [0.9, 1.1]\"\n"
            "        # The builder does not record the loop target (530 does not touch it). A value that\n"
            "        # IS recorded must be the one that was requested.\n"
            "        if got_rho is not None and abs(float(got_rho) - s) >= 1e-12:\n"
            "            return f\"rho_scale meta {got_rho} vs {rho}\"\n"
            "    return None\n"
            "\n"
            "\n"
            "def main() -> int:\n",
        ),
        (
            "    if cmd == \"riders_match\":\n"
            "        meta_path, mag, sign = args[0], args[1], args[2]\n"
            "        try:\n"
            "            cfg = json.load(open(meta_path))[\"config\"]\n"
            "        except Exception as exc:\n"
            "            print(f\"riders_match failed: {exc}\", file=sys.stderr)\n"
            "            return 2\n"
            "        got_mag = (cfg.get(\"de_mag\") or {}).get(\"variant\") if isinstance(cfg.get(\"de_mag\"), dict) else None\n"
            "        if (mag == \"-\") != (got_mag is None) or (mag != \"-\" and got_mag != mag):\n"
            "            print(f\"de_mag mismatch: meta {got_mag!r} vs {mag!r}\", file=sys.stderr)\n"
            "            return 1\n"
            "        ms = cfg.get(\"de_mde_sign\")\n"
            "        if sign == \"-\":\n"
            "            if ms is not None:\n"
            "                print(f\"de_mde_sign present but not requested: {ms}\", file=sys.stderr)\n"
            "                return 1\n"
            "        else:\n"
            "            src, _, tau = sign.partition(\":\")\n"
            "            if not isinstance(ms, dict) or ms.get(\"source\") != src or abs(float(ms.get(\"tau\", -1)) - float(tau)) >= 1e-12:\n"
            "                print(f\"de_mde_sign mismatch: meta {ms} vs {sign}\", file=sys.stderr)\n"
            "                return 1\n"
            "        return 0\n",
            "    if cmd == \"riders_match\":\n"
            "        if len(args) not in (3, 4):\n"
            "            print(f\"riders_match: expected 3 or 4 arguments, got {len(args)}\", file=sys.stderr)\n"
            "            return 2\n"
            "        meta_path, mag, sign = args[0], args[1], args[2]\n"
            "        spec = args[3] if len(args) == 4 else None\n"
            "        try:\n"
            "            cfg = json.load(open(meta_path))[\"config\"]\n"
            "        except Exception as exc:\n"
            "            print(f\"riders_match failed: {exc}\", file=sys.stderr)\n"
            "            return 2\n"
            "        got_mag = (cfg.get(\"de_mag\") or {}).get(\"variant\") if isinstance(cfg.get(\"de_mag\"), dict) else None\n"
            "        if (mag == \"-\") != (got_mag is None) or (mag != \"-\" and got_mag != mag):\n"
            "            print(f\"de_mag mismatch: meta {got_mag!r} vs {mag!r}\", file=sys.stderr)\n"
            "            return 1\n"
            "        ms = cfg.get(\"de_mde_sign\")\n"
            "        if sign == \"-\":\n"
            "            if ms is not None:\n"
            "                print(f\"de_mde_sign present but not requested: {ms}\", file=sys.stderr)\n"
            "                return 1\n"
            "        else:\n"
            "            want = _mde_sign_request(sign)\n"
            "            bad = (not isinstance(ms, dict) or want is None\n"
            "                   or any(ms.get(k) != v if isinstance(v, str) else (ms.get(k) is None or abs(float(ms[k]) - v) >= 1e-12)\n"
            "                          for k, v in want.items()))\n"
            "            if bad:\n"
            "                print(f\"de_mde_sign mismatch: meta {ms} vs {sign}\", file=sys.stderr)\n"
            "                return 1\n"
            "        if spec is None:\n"
            "            return 0\n"
            "        parsed = _october530_spec(spec)\n"
            "        if parsed is None:\n"
            "            print(f\"riders_match: bad October spec {spec!r}\", file=sys.stderr)\n"
            "            return 2\n"
            "        msg = _october530_mismatch(cfg, parsed)\n"
            "        if msg:\n"
            "            print(msg, file=sys.stderr)\n"
            "            return 1\n"
            "        return 0\n",
        ),
    ]
    return _apply(src, edits, "def_day_helpers.py")


def patch_relevel(src: str) -> str:
    edits = [
        (
            "usage: uv run python relevel_iter.py <build.ntclevel.json>:<D=f,E=f,F=f> [<...>:<...> ...] [--sd 0.01] [--tol 0.01]\n",
            "usage: uv run python relevel_iter.py <build.ntclevel.json>:<D=f,E=f,F=f> [<...>:<...> ...] [--sd 0.01] [--tol 0.01] [--rho-scale 1]\n"
            "  --rho-scale s (530, default 1, must be in [0.9, 1.1]) targets rho* x s. Omitted, the target is rho*\n"
            "  and the arithmetic is the unscaled formula (no multiply).\n",
        ),
        (
            "ap.add_argument(\"--tol\", type=float, default=0.01, help=\"relative tolerance on rho (0.01 = 0.5 %% in scale)\")\n"
            "a = ap.parse_args()\n"
            "q = 0.31 / (0.31 + 1.78)\n"
            "rho_star = math.exp(2.0 * a.sd * norm.ppf(q))\n",
            "ap.add_argument(\"--tol\", type=float, default=0.01, help=\"relative tolerance on rho (0.01 = 0.5 %% in scale)\")\n"
            "ap.add_argument(\"--rho-scale\", type=float, default=1.0, help=\"530: target rho* x this (default 1)\")\n"
            "a = ap.parse_args()\n"
            "if not 0.9 <= a.rho_scale <= 1.1:\n"
            "    raise SystemExit(f\"--rho-scale must be in [0.9, 1.1] (got {a.rho_scale})\")\n"
            "q = 0.31 / (0.31 + 1.78)\n"
            "base = math.exp(2.0 * a.sd * norm.ppf(q))\n"
            "rho_star = base if a.rho_scale == 1.0 else base * a.rho_scale\n",
        ),
    ]
    return _apply(src, edits, "relevel_iter.py")


def patch_sameline(src: str) -> str:
    """517's --programme-estimator, copied from experiments/517-flex-programme/patch_sameline517.py.
    Default pdown does not enter programme_meanz_ne and adds no provenance key."""
    edits = [
        (
            "def programme(delta, perts, genes, official, panel) -> tuple[np.ndarray, np.ndarray, int, np.ndarray]:\n",
            'ESSENTIAL_LIST = ROOT / "data/donors/k562_essential/perts.txt"   # 517: meanz_ne drops these rows\n\n\n'
            "def programme_meanz_ne(delta, perts, genes, official, panel) -> tuple[np.ndarray, np.ndarray, int]:\n"
            '    """517 (artifact517.py): mean over the non-essential, non-panel rows of delta / sd (sd over all non-panel rows),\n'
            "    as a pseudo-P whose |P - 0.5| >= 0.05 coverage equals the P_down rule's (so the hook's tau keeps its meaning).\"\"\"\n"
            "    ess = {l.strip() for l in open(ESSENTIAL_LIST) if l.strip()}\n"
            "    keep = np.array([p not in panel for p in perts])\n"
            "    X = np.asarray(delta[np.flatnonzero(keep)], dtype=np.float64)\n"
            "    ne = np.array([p not in ess for p in np.asarray(perts)[keep]])\n"
            "    nz = (X != 0).sum(0); dn = (X < 0).sum(0)\n"
            "    with np.errstate(invalid=\"ignore\", divide=\"ignore\"):\n"
            "        pd_all = np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)\n"
            "    pos = {g: j for j, g in enumerate(genes)}\n"
            "    col = np.array([pos.get(g, -1) for g in official])\n"
            "    P_all = np.where(col >= 0, pd_all[np.maximum(col, 0)], np.nan)\n"
            "    sd = X.std(0); sd[sd == 0] = 1.0\n"
            "    zb = (X[ne] / sd).mean(0)\n"
            "    z = np.where(np.isfinite(P_all) & (col >= 0), zb[np.maximum(col, 0)], np.nan)\n"
            "    n_cov = int((np.isfinite(P_all) & (np.abs(P_all - 0.5) >= 0.05)).sum())\n"
            "    z_q = float(np.sort(np.abs(np.nan_to_num(z)))[::-1][n_cov - 1])\n"
            "    pdown = np.where(np.isfinite(z), np.clip(0.5 - 0.05 * z / z_q, 0.001, 0.999), np.nan)\n"
            "    psign = np.where(np.isfinite(pdown) & (pdown <= 0.5), 1.0, -1.0)\n"
            "    return psign, pdown, int(ne.sum())\n\n\n"
            "def programme(delta, perts, genes, official, panel) -> tuple[np.ndarray, np.ndarray, int, np.ndarray]:\n",
        ),
        (
            '    ap.add_argument("--programme-source", default="",\n',
            '    ap.add_argument("--programme-estimator", choices=["pdown", "meanz_ne"], default="pdown",\n'
            '                    help="517: meanz_ne = the hook-source encoding for --programme-source (not for the block prior)")\n'
            '    ap.add_argument("--programme-source", default="",\n',
        ),
        (
            "        spsign, spdown, sn, _ = programme(sdelta, sperts, sgenes, official, panel)\n",
            "        spsign, spdown, sn, _ = programme(sdelta, sperts, sgenes, official, panel)\n"
            "        if a.programme_estimator == \"meanz_ne\":\n"
            "            spsign, spdown, sn_ne = programme_meanz_ne(sdelta, sperts, sgenes, official, panel)\n"
            "            print(f\"[sameline] --programme-estimator meanz_ne: {sn_ne} of {sn} rows (non-essential)\", flush=True)\n",
        ),
        (
            '                "bank": a.programme_source, "bank_delta_sha256": ssha, "programme_rows": sn,\n',
            '                "bank": a.programme_source, "bank_delta_sha256": ssha, "programme_rows": sn,\n'
            '                **({} if a.programme_estimator == "pdown" else {"programme_estimator": a.programme_estimator}),\n',
        ),
    ]
    return _apply(src, edits, "def_sameline.py")


def patch_recipe(src: str) -> str:
    edits = [
        (
            '    "de_mag": "DEF_DAY_DE_MAG",\n'
            '    "mde_sign": "DEF_DAY_MDE_SIGN",\n'
            "}\n",
            '    "de_mag": "DEF_DAY_DE_MAG",\n'
            '    "mde_sign": "DEF_DAY_MDE_SIGN",\n'
            '    "rho_scale": "DEF_DAY_RHO_SCALE",\n'
            '    "rc_callset_cpm": "DEF_DAY_RC_CALLSET_CPM",\n'
            '    "hook_estimator": "DEF_DAY_HOOK_ESTIMATOR",\n'
            '    "lodo_strength": "DEF_DAY_LODO_STRENGTH",\n'
            "}\n",
        ),
        (
            'BUILDER_NEEDS = {"DEF_DAY_DE_MAG": "VCC_DE_MAG", "DEF_DAY_MDE_SIGN": "VCC_DE_MDE_SIGN",\n'
            '                 "DEF_DAY_PRIOR_GATE": "VCC_DE_PRIOR_GATE", "DEF_DAY_SHARPEN_PT": "stercile"}\n',
            'BUILDER_NEEDS = {"DEF_DAY_DE_MAG": "VCC_DE_MAG", "DEF_DAY_MDE_SIGN": "VCC_DE_MDE_SIGN",\n'
            '                 "DEF_DAY_PRIOR_GATE": "VCC_DE_PRIOR_GATE", "DEF_DAY_SHARPEN_PT": "stercile",\n'
            '                 "DEF_DAY_RC_CALLSET_CPM": "VCC_RC_CALLSET_CPM"}\n'
            'SAMELINE = Path(os.environ.get("OCTOBER_SAMELINE", ROOT / "scripts/def_sameline.py"))\n',
        ),
        (
            '    if "DEF_DAY_LODO_HI_PCT" in env and env.get("DEF_DAY_LODO") != "1":\n'
            '        problems.append("lodo_hi_pct needs lodo = 1")\n'
            '    if s.get("cd4_stack", {}).get("value") not in (None, "x1"):\n'
            '        problems.append(f"cd4_stack {s[\'cd4_stack\'][\'value\']!r}: def_day.sh dumps CD4 x1 only; that needs a patch first")\n',
            '    if "DEF_DAY_LODO_HI_PCT" in env and env.get("DEF_DAY_LODO") != "1":\n'
            '        problems.append("lodo_hi_pct needs lodo = 1")\n'
            "    # 530: x1 is today's dump and is not exported, so a DECIDED x1 leaves the command unchanged.\n"
            '    cd4 = s.get("cd4_stack", {})\n'
            '    if cd4.get("status") == "DECIDED":\n'
            '        cd4_val = str(cd4.get("value"))\n'
            '        if cd4_val == "x1.5":\n'
            '            env["DEF_DAY_CD4_STACK"] = "x1.5"\n'
            '        elif cd4_val not in {"x1", "off", "0", "none"}:\n'
            '            problems.append(f"cd4_stack {cd4_val!r} is not x1 or x1.5")\n'
            '    _problems_530(env, problems)\n'
            '    # is_on treats "0" as off, so a DECIDED theta of 0 would vanish instead of being refused.\n'
            '    # 0 still rewrites the LODO npz (522); it is not a way to leave the rider off.\n'
            '    _ls = s.get("lodo_strength", {})\n'
            '    if _ls.get("status") == "DECIDED" and str(_ls.get("value")).strip() == "0":\n'
            '        problems.append("lodo_strength \'0\' must be a non-zero theta in [-1, 1]")\n'
        ),
        (
            "def installed_problems(env: dict, flags: str) -> list:\n"
            "    dd, bu = DEF_DAY.read_text(), BUILDER.read_text()\n"
            "    # The export, not a comment that mentions the knob (519): a substring match let a comment pass.\n"
            "    export_of = {\n"
            '        "DEF_DAY_DE_MAG": \'export VCC_DE_MAG="$DEF_DAY_DE_MAG"\',\n'
            '        "DEF_DAY_MDE_SIGN": \'export VCC_DE_MDE_SIGN="$DEF_DAY_MDE_SIGN"\',\n'
            "    }\n",
            "def _problems_530(env: dict, problems: list) -> None:\n"
            '    """Range and dependency checks for 530\'s knobs. PENDING stays out of env, so this does not run for them."""\n'
            '    if "DEF_DAY_RHO_SCALE" in env:\n'
            "        try:\n"
            '            s = float(env["DEF_DAY_RHO_SCALE"])\n'
            "            ok = 0.9 <= s <= 1.1\n"
            "        except ValueError:\n"
            "            ok = False\n"
            "        if not ok:\n"
            '            problems.append(f"rho_scale {env[\'DEF_DAY_RHO_SCALE\']!r} must be in [0.9, 1.1]")\n'
            '    if "DEF_DAY_RC_CALLSET_CPM" in env:\n'
            '        n = env["DEF_DAY_RC_CALLSET_CPM"]\n'
            "        if not (n.isdigit() and 100 <= int(n) <= 5000):\n"
            '            problems.append(f"rc_callset_cpm {n!r} must be an integer in [100, 5000]")\n'
            '    if "DEF_DAY_HOOK_ESTIMATOR" in env:\n'
            '        if env["DEF_DAY_HOOK_ESTIMATOR"] != "meanz_ne":\n'
            '            problems.append(f"hook_estimator {env[\'DEF_DAY_HOOK_ESTIMATOR\']!r} must be meanz_ne")\n'
            '        if not env.get("DEF_DAY_MDE_SIGN", "").startswith("hookbank:"):\n'
            '            problems.append("hook_estimator needs DEF_DAY_MDE_SIGN=hookbank:<bank>:<tau>...")\n'
            '    if "DEF_DAY_LODO_STRENGTH" in env:\n'
            '        if env.get("DEF_DAY_LODO") != "1":\n'
            '            problems.append("lodo_strength needs lodo = 1")\n'
            '        if "DEF_DAY_LODO_HI_PCT" not in env:\n'
            '            problems.append("lodo_strength needs the percentile cap (lodo_hi_pct)")\n'
            "        import re as _re\n"
            '        th = env["DEF_DAY_LODO_STRENGTH"]\n'
            "        try:\n"
            "            v = float(th)\n"
            '            ok = _re.fullmatch(r"-?(?:0|[1-9]\\d*)(?:\\.\\d+)?", th) is not None and v != 0.0 and abs(v) <= 1.0\n'
            "        except ValueError:\n"
            "            ok = False\n"
            "        if not ok:\n"
            '            problems.append(f"lodo_strength {th!r} must be a non-zero theta in [-1, 1]")\n'
            '    mde = env.get("DEF_DAY_MDE_SIGN", "")\n'
            "    if mde and not (mde.startswith(\"clean4:\") or mde.startswith(\"programme:\") or mde.startswith(\"hookbank:\")):\n"
            '        problems.append(f"mde_sign {mde!r} must be clean4:, programme: or hookbank:")\n'
            "\n"
            "\n"
            "def installed_problems(env: dict, flags: str) -> list:\n"
            "    dd, bu = DEF_DAY.read_text(), BUILDER.read_text()\n"
            "    sameline = SAMELINE.read_text() if SAMELINE.is_file() else \"\"\n"
            "    # The export, not a comment that mentions the knob (519): a substring match let a comment pass.\n"
            "    # 517 rewrites the MDE export to $MDS (clean4 keeps the request; hookbank becomes file:<npz>:<sha>:...).\n"
            "    export_of = {\n"
            '        "DEF_DAY_DE_MAG": \'export VCC_DE_MAG="$DEF_DAY_DE_MAG"\',\n'
            '        "DEF_DAY_MDE_SIGN": \'export VCC_DE_MDE_SIGN="$MDS"\',\n'
            '        "DEF_DAY_RC_CALLSET_CPM": \'export VCC_RC_CALLSET_CPM="$DEF_DAY_RC_CALLSET_CPM"\',\n'
            '        "DEF_DAY_RHO_SCALE": \'--rho-scale "$DEF_DAY_RHO_SCALE"\',\n'
            '        "DEF_DAY_CD4_STACK": \'VCC_CD4_TARGET_G="$CD4G"\',\n'
            '        "DEF_DAY_LODO_STRENGTH": "def_strength_key522.py",\n'
            '        "DEF_DAY_HOOK_ESTIMATOR": \'--programme-estimator "$DEF_DAY_HOOK_ESTIMATOR"\',\n'
            "    }\n",
        ),
        (
            "    out = [f\"scripts/def_day.sh has no {var}\" for var in env if export_of.get(var, var) not in dd]\n"
            "    out += [f\"the builder lacks {need} (needed by {var})\" for var, need in BUILDER_NEEDS.items() if var in env and need not in bu]\n",
            "    out = [f\"scripts/def_day.sh has no {var}\" for var in env if export_of.get(var, var) not in dd]\n"
            "    if \"DEF_DAY_CD4_STACK\" in env and 'VCC_CD4_TARGET_G=\"$CD4G\"' not in dd:\n"
            "        out.append(\"cd4_stack x1.5: scripts/def_day.sh has no DEF_DAY_CD4_STACK dump override (530)\")\n"
            "    out += [f\"the builder lacks {need} (needed by {var})\" for var, need in BUILDER_NEEDS.items() if var in env and need not in bu]\n"
            "    if env.get(\"DEF_DAY_MDE_SIGN\", \"\").startswith(\"hookbank:\") and \"hookbank:\" not in dd:\n"
            "        out.append(\"scripts/def_day.sh has no hookbank form (517)\")\n"
            "    if env.get(\"DEF_DAY_MDE_SIGN\", \"\").startswith(\"hookbank:\") and \"file:<npz>\" not in bu:\n"
            "        out.append(\"the builder lacks file:<npz> (needed by hookbank)\")\n"
            "    if \"DEF_DAY_HOOK_ESTIMATOR\" in env and \"--programme-estimator\" not in sameline:\n"
            "        out.append(\"scripts/def_sameline.py has no --programme-estimator (517)\")\n",
        ),
    ]
    return _apply(src, edits, "october_recipe.py")


def patch_decisions(src: str) -> str:
    """JSON-level edit (coordinator, 2026-09-29): the string anchors broke once october-decisions.json gained 534's
    held_out lines. Parse, edit the two existing switches, insert the four new PENDING switches before cd4_thresholds,
    and write back in the file's own format (indent 1, UTF-8)."""
    d = json.loads(src)
    s = d["switches"]
    ms = s["mde_sign"]
    ms["needs"] = "builder with VCC_DE_MDE_SIGN including the file form (517); 530 accepts clean4, programme and hookbank"
    ms["if_hookbank"] = ("hookbank:<bank>:<tau>[:<tau_c4>[:le<N>]] (the meta records source file; pair with hook_estimator "
                         "meanz_ne)")
    cd = s["cd4_stack"]
    add = (" 530: x1 is today's dump and is not exported; x1.5 exports DEF_DAY_CD4_STACK (457's dump-subshell "
           "VCC_CD4_TARGET_G).")
    if "530:" not in cd.get("evidence", ""):
        cd["evidence"] = cd.get("evidence", "").replace(" def_day.sh dumps H1 x2 / CD4 x1.", "") + add
    new = {
        "rho_scale": {"status": "PENDING", "value": None,
                      "decided_by": "FINAL-PHASE row 2: run rho* x 0.98 / x 1.02 if the members come back. 1 is today's rho*",
                      "if_on": "a scale s with 0.9 <= s <= 1.1", "otherwise": "off",
                      "needs": "530 DEF_DAY_RHO_SCALE (relevel_iter --rho-scale; inband and best_safe take the same scale)",
                      "held_out": "FINAL-PHASE feedback only: rho* is the A/B/C-calibrated NTC-proxy target; no A/B/C board read can move it"},
        "rc_callset_cpm": {"status": "PENDING", "value": None, "decided_by": "523's K pair (10-06/07)",
                           "if_on": "an integer n with 100 <= n <= 5000", "otherwise": "off",
                           "needs": "530 DEF_DAY_RC_CALLSET_CPM on the level builds; builder VCC_RC_CALLSET_CPM (523). Unset = 2200",
                           "held_out": "521: H1 regime analysis (K 1,600 +0.0048 at f 0.35, -0.0155 at f 1); the board's n_real regime is what 523 reads"},
        "hook_estimator": {"status": "PENDING", "value": None,
                           "decided_by": "527's arm 446bj2zmz (10-01 null branch or 10-07) and 496a",
                           "if_on": "meanz_ne", "otherwise": "off",
                           "needs": "DEF_DAY_MDE_SIGN=hookbank:<bank>:<tau>...; 517's def_sameline --programme-estimator (530 folds both patches)",
                           "held_out": "H1 only (the Arc programme's own line; h1-test targets excluded, control halves checked, 517): inherits programme_source's caveat"},
        "lodo_strength": {"status": "PENDING", "value": None, "decided_by": "only if 495 (cap 6) transfers on 10-03",
                          "if_on": "theta, a non-zero plain decimal in [-1, 1] (522's priced arm is 0.25)", "otherwise": "off",
                          "needs": "DEF_DAY_LODO=1 and lodo_hi_pct (the percentile cap); 522's def_strength_key522.py, run between the LODO key and its pin",
                          "held_out": "522: pct x exp(0.25 z_src) +0.7 ... +1.9 % direction on all 9 October-like folds, ~0 on A/B/C (off-board value)"},
    }
    out_sw = {}
    for k, v in s.items():
        if k == "cd4_thresholds":
            out_sw.update(new)
        out_sw[k] = v
    for k, v in new.items():                         # no cd4_thresholds key: append at the end
        out_sw.setdefault(k, v)
    d["switches"] = out_sw
    out = json.dumps(d, indent=1, ensure_ascii=False) + "\n"
    json.loads(out)
    return out


def patch_gate(src: str) -> str:
    edits = [
        (
            '    "DEF_DAY_RC_CALLSET_CPM": {"rc_callset_cpm"},\n'
            "}\n",
            '    "DEF_DAY_RC_CALLSET_CPM": {"rc_callset_cpm"},\n'
            '    "DEF_DAY_RHO_SCALE": {"bulk_ctx_scale"},\n'
            '    "DEF_DAY_CD4_STACK": {"bulk_load", "bulk_load_sha256"},\n'
            '    "DEF_DAY_LODO_STRENGTH": {"bulk_load", "bulk_load_sha256", "pert_scale", "pert_scale_sha256"},\n'
            '    "DEF_DAY_HOOK_ESTIMATOR": {"de_mde_sign"},\n'
            "}\n",
        ),
        (
            "def chk(name: str, cond: bool, detail: str) -> None:\n",
            "def _dump_name530(stack: str, strength: str | None) -> str:\n"
            '    """Same filename rule as def_day.sh / def_day_helpers._dump_name530."""\n'
            '    base = "dump_h1x2_cd4x15" if stack == "x1.5" else "dump_h1x2_cd4x1"\n'
            "    if strength:\n"
            '        base += "_s" + str(strength).replace("-", "m")\n'
            '    return base + ".npz"\n'
            "\n"
            "\n"
            "def _mde_ok530(spec: str, ms) -> bool:\n"
            '    """526 compared source and tau after one partition, which is clean4:<tau> only.\n'
            "    hookbank:<bank>:<tau>[:tau_c4[:leN]] is recorded with source file (the shell rewrites it).\n"
            '    """\n'
            "    if not isinstance(ms, dict):\n"
            "        return False\n"
            '    p = spec.split(":")\n'
            "    try:\n"
            '        if p[0] == "clean4" and len(p) == 2:\n'
            '            return ms.get("source") == "clean4" and abs(float(ms.get("tau", -1)) - float(p[1])) < 1e-12\n'
            '        if p[0] == "programme" and len(p) in (2, 3):\n'
            "            tau_c4 = float(p[2]) if len(p) == 3 else 0.05\n"
            '            return (ms.get("source") == "programme"\n'
            '                    and abs(float(ms.get("tau", -1)) - float(p[1])) < 1e-12\n'
            '                    and abs(float(ms.get("tau_c4", -1)) - tau_c4) < 1e-12)\n'
            '        if p[0] == "hookbank" and len(p) in (3, 4, 5) and p[1]:\n'
            "            tau, tau_c4 = float(p[2]), float(p[3]) if len(p) >= 4 else 0.05\n"
            '            ok = (ms.get("source") == "file" and isinstance(ms.get("sha256"), str) and len(ms["sha256"]) == 64\n'
            '                  and abs(float(ms.get("tau", -1)) - tau) < 1e-12\n'
            '                  and abs(float(ms.get("tau_c4", -1)) - tau_c4) < 1e-12)\n'
            "            if len(p) == 5:\n"
            '                ok = ok and p[4].startswith("le") and ms.get("le") == int(p[4][2:])\n'
            "            else:\n"
            '                ok = ok and "le" not in ms\n'
            "            return bool(ok)\n"
            "    except (TypeError, ValueError, IndexError):\n"
            "        return False\n"
            "    return False\n"
            "\n"
            "\n"
            "def chk(name: str, cond: bool, detail: str) -> None:\n",
        ),
        (
            'if "DEF_DAY_MDE_SIGN" in env:\n'
            "    ms = c.get(\"de_mde_sign\")\n"
            '    src, _, tau = env["DEF_DAY_MDE_SIGN"].partition(":")\n'
            "    try:\n"
            '        tau_ok = isinstance(ms, dict) and ms.get("source") == src and abs(float(ms.get("tau", -1)) - float(tau)) < 1e-12\n'
            "    except (TypeError, ValueError):\n"
            "        tau_ok = False\n"
            '    chk("uncalled-gene sign hook", tau_ok, f"meta de_mde_sign {json.dumps(ms)[:160]} vs {env[\'DEF_DAY_MDE_SIGN\']}")\n',
            'if "DEF_DAY_MDE_SIGN" in env:\n'
            "    ms = c.get(\"de_mde_sign\")\n"
            '    tau_ok = _mde_ok530(env["DEF_DAY_MDE_SIGN"], ms)\n'
            '    chk("uncalled-gene sign hook", tau_ok, f"meta de_mde_sign {json.dumps(ms)[:160]} vs {env[\'DEF_DAY_MDE_SIGN\']}")\n',
        ),
        (
            'if "DEF_DAY_RC_CALLSET_CPM" in env:\n'
            '    chk("call-set CPM block", c.get("rc_callset_cpm") == int(env["DEF_DAY_RC_CALLSET_CPM"]), f"meta {c.get(\'rc_callset_cpm\')}")\n',
            'if "DEF_DAY_RC_CALLSET_CPM" in env:\n'
            '    chk("call-set CPM block", c.get("rc_callset_cpm") == int(env["DEF_DAY_RC_CALLSET_CPM"]), f"meta {c.get(\'rc_callset_cpm\')}")\n'
            "if \"DEF_DAY_CD4_STACK\" in env or \"DEF_DAY_LODO_STRENGTH\" in env:\n"
            '    stack = env.get("DEF_DAY_CD4_STACK", "x1")\n'
            '    theta = env.get("DEF_DAY_LODO_STRENGTH")\n'
            "    want_dump = _dump_name530(stack, theta)\n"
            "    g = c.get(\"cd4_target_g\")\n"
            "    try:\n"
            "        dose_ok = g is not None and abs(float(g) - 0.3443655278649641) < 1e-12\n"
            "    except (TypeError, ValueError):\n"
            "        dose_ok = False\n"
            '    chk("stack dump name", c.get("bulk_load") == want_dump and dose_ok,\n'
            '        f"bulk_load {c.get(\'bulk_load\')} cd4_target_g {g} vs {want_dump} (level g stays DEF_ENV)")\n'
            "else:\n"
            '    chk("CD4 stack x1", c.get("bulk_load") == "dump_h1x2_cd4x1.npz", f"bulk_load {c.get(\'bulk_load\')}")\n'
            'if "DEF_DAY_RHO_SCALE" in env:\n'
            '    scale = env["DEF_DAY_RHO_SCALE"]\n'
            '    chk("rho scale", re.search(rf"rho-scale={re.escape(scale)}(\\s|$)", dd) is not None, f"def_day log names rho-scale={scale}")\n'
            'if "DEF_DAY_LODO_STRENGTH" in env:\n'
            '    side = RUN / "pert_scale_lodo_strength.json"\n'
            "    theta = None\n"
            "    if side.is_file():\n"
            "        try:\n"
            '            theta = json.loads(side.read_text()).get("theta")\n'
            "        except (OSError, json.JSONDecodeError):\n"
            "            theta = None\n"
            "    try:\n"
            '        th_ok = theta is not None and abs(float(theta) - float(env["DEF_DAY_LODO_STRENGTH"])) < 1e-12\n'
            "    except (TypeError, ValueError):\n"
            "        th_ok = False\n"
            '    chk("LODO strength rider", th_ok and f"LODO strength rider: theta {env[\'DEF_DAY_LODO_STRENGTH\']}" in dd,\n'
            '        f"sidecar theta {theta}")\n'
            'if "DEF_DAY_HOOK_ESTIMATOR" in env:\n'
            '    est = env["DEF_DAY_HOOK_ESTIMATOR"]\n'
            '    chk("hook estimator", est == "meanz_ne" and f"estimator {est}" in dd\n'
            '        and str(env.get("DEF_DAY_MDE_SIGN", "")).startswith("hookbank:"),\n'
            '        f"log names estimator {est}; mde {env.get(\'DEF_DAY_MDE_SIGN\')}")\n',
        ),
    ]
    return _apply(src, edits, "gate524.py")


def build(paths: dict[str, Path]) -> dict[str, str]:
    texts = {k: p.read_text() for k, p in paths.items()}
    marks = _markers(texts)
    if all(marks.values()):
        raise SystemExit("already patched")
    if any(marks.values()):
        have = ", ".join(k for k, v in marks.items() if v)
        missing = ", ".join(k for k, v in marks.items() if not v)
        raise SystemExit(f"partial bundle (have {have}; missing {missing}); restore and re-apply")
    out = {
        "def_day": patch_def_day(texts["def_day"]),
        "helpers": patch_helpers(texts["helpers"]),
        "relevel": patch_relevel(texts["relevel"]),
        "sameline": patch_sameline(texts["sameline"]),
        "recipe": patch_recipe(texts["recipe"]),
        "decisions": patch_decisions(texts["decisions"]),
        "gate": patch_gate(texts["gate"]),
    }
    if not all(_markers(out).values()):
        raise SystemExit(f"patch did not land every marker: {_markers(out)}")
    compile(out["helpers"], str(paths["helpers"]), "exec")
    compile(out["relevel"], str(paths["relevel"]), "exec")
    compile(out["sameline"], str(paths["sameline"]), "exec")
    compile(out["recipe"], str(paths["recipe"]), "exec")
    compile(out["gate"], str(paths["gate"]), "exec")
    json.loads(out["decisions"])
    # Syntax-check the shell without touching the target. macOS bash 3.2 accepts the constructs used here
    # (no associative arrays, no empty-array expansion under set -u).
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
        "relevel": Path("experiments/448-def-rehearsal-q175/relevel_iter.py"),
        "sameline": Path("scripts/def_sameline.py"),
        "recipe": Path("scripts/october_recipe.py"),
        "decisions": Path("docs/october-decisions.json"),
        "gate": Path("experiments/524-october-recipe-dr/gate524.py"),
    }
    given = {
        "def_day": a.def_day, "helpers": a.helpers, "relevel": a.relevel, "sameline": a.sameline,
        "recipe": a.recipe, "decisions": a.decisions, "gate": a.gate,
    }
    out = {}
    for k, rel_p in rel.items():
        if given[k] is not None:
            out[k] = Path(given[k])
        elif root is not None:
            out[k] = Path(root) / rel_p
        else:
            raise SystemExit(f"missing --{k.replace('_', '-')} (or pass --root)")
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path)
    ap.add_argument("--def-day", type=Path)
    ap.add_argument("--helpers", type=Path)
    ap.add_argument("--relevel", type=Path)
    ap.add_argument("--sameline", type=Path)
    ap.add_argument("--recipe", type=Path)
    ap.add_argument("--decisions", type=Path)
    ap.add_argument("--gate", type=Path)
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
        tmp = dest.with_name(dest.name + ".530tmp")
        tmp.write_text(text)
        temps.append((tmp, dest))
    for tmp, dest in temps:
        os.replace(tmp, dest)
    print("patched " + ", ".join(str(p) for p in paths.values()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
