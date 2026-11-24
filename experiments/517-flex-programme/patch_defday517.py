"""517: October pass-through for the PROGRAMME-sourced uncalled-gene hook, on top of 502's def_day.sh / def_day_helpers.py
(502 = 497's e29e793a / 60b2a777 + patch_defday502.py -> 610db782 / 20142142). NOT APPLIED. Needs the builder with 496's and
517's patches (`VCC_DE_MDE_SIGN` programme / file forms).

  DEF_DAY_MDE_SIGN=clean4:<tau>                        502's form, unchanged
  DEF_DAY_MDE_SIGN=programme:<tau>[:<tau_c4>]          the hook reads the day's same-line artifact (pdown_<c>); needs
                                                       DEF_DAY_PROGRAMME_SOURCE (the block prior is the programme too);
                                                       follows DEF_DAY_PRIOR_GATE (weak targets -> CLEAN4 at tau_c4)
  DEF_DAY_MDE_SIGN=hookbank:<bank>:<tau>[:<tau_c4>[:le<N>]]
                                                       the hook ALONE (CLEAN4 block): the script builds a second
                                                       def_sameline.py artifact with `--programme-source <bank>` (the panel's
                                                       own rows excluded, 481 F7's row guard) and exports
                                                       VCC_DE_MDE_SIGN=file:<W>/hookprog.npz:<sha>:<tau>[:<tau_c4>[:le<N>]]
  DEF_DAY_HOOK_ESTIMATOR=meanz_ne                      the hookbank artifact's estimator (needs def_sameline.py with
                                                       patch_sameline517.py; default pdown = 396's P_down rule)
riders_match compares the final meta's de_mde_sign with the exported string for every form (source, tau, tau_c4, le, and
the file's sha). Unset / clean4: byte-identical behaviour to 502's script (the added branches do not run).
usage: uv run python experiments/517-flex-programme/patch_defday517.py <def_day.sh with 502> <def_day_helpers.py with 502> [--check]
"""
import os
import sys
from pathlib import Path

dd, hp = Path(sys.argv[1]), Path(sys.argv[2]); s, h = dd.read_text(), hp.read_text()
if "DEF_DAY_MDE_SIGN" not in s or "riders_match" not in h:
    raise SystemExit("apply 502's patch_defday502.py first")
if "hookbank:" in s:
    raise SystemExit("already patched")
E_DD = [
    ("#   DEF_DAY_MDE_SIGN=clean4:<tau> (502, default off; e.g. clean4:0.05): 496's uncalled-gene sign hook (VCC_DE_MDE_SIGN).\n",
     "#   DEF_DAY_MDE_SIGN=clean4:<tau> (502, default off; e.g. clean4:0.05): 496's uncalled-gene sign hook (VCC_DE_MDE_SIGN).\n"
     "#    517: or programme:<tau>[:<tau_c4>] (the hook reads the programme of DEF_DAY_PROGRAMME_SOURCE's artifact and follows\n"
     "#    DEF_DAY_PRIOR_GATE), or hookbank:<bank>:<tau>[:<tau_c4>[:le<N>]] (the hook ALONE from <bank>'s programme, CLEAN4\n"
     "#    block: a second def_sameline.py artifact, exported as VCC_DE_MDE_SIGN=file:<npz>:<sha>:...). Builder >= 517's patch.\n"
     "#    DEF_DAY_HOOK_ESTIMATOR=meanz_ne (517, default pdown): the hook artifact's estimator (def_sameline.py >= 517's patch).\n"),
    ('  case "$DEF_DAY_MDE_SIGN" in clean4:*) ;; *) die "DEF_DAY_MDE_SIGN must be clean4:<tau> (got $DEF_DAY_MDE_SIGN)" ;; esac\n',
     '  case "$DEF_DAY_MDE_SIGN" in\n'
     '    clean4:*) ;;\n'
     '    programme:*) [ -n "${DEF_DAY_PROGRAMME_SOURCE:-}" ] || die "DEF_DAY_MDE_SIGN=$DEF_DAY_MDE_SIGN reads the programme artifact: it needs DEF_DAY_PROGRAMME_SOURCE"\n'
     '                 grep -q "programme:<tau>" scripts/make_submission_pm_cd4.py || die "the builder has no programme hook form (517)" ;;\n'
     '    hookbank:*) grep -q "file:<npz>" scripts/make_submission_pm_cd4.py || die "the builder has no file hook form (517)" ;;\n'
     '    *) die "DEF_DAY_MDE_SIGN must be clean4:<tau>, programme:<tau>[:<tau_c4>] or hookbank:<bank>:<tau>[:<tau_c4>[:le<N>]] (got $DEF_DAY_MDE_SIGN)" ;;\n'
     '  esac\n'),
    ('# 2b. DEF_DAY_LODO=1 (460; default off)',
     '# 2a. 517: DEF_DAY_MDE_SIGN=hookbank:<bank>:... -- the hook\'s own programme artifact (the block prior is untouched)\n'
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
     'esac\n\n'
     '# 2b. DEF_DAY_LODO=1 (460; default off)'),
    ('[ -n "${DEF_DAY_MDE_SIGN:-}" ] && { export VCC_DE_MDE_SIGN="$DEF_DAY_MDE_SIGN"; log "DE rider for the level builds: VCC_DE_MDE_SIGN=$VCC_DE_MDE_SIGN"; }\n',
     '[ -n "${DEF_DAY_MDE_SIGN:-}" ] && { export VCC_DE_MDE_SIGN="$MDS"; log "DE rider for the level builds: VCC_DE_MDE_SIGN=$VCC_DE_MDE_SIGN"; }\n'),
]
E_HP = [
    ('  def_day_helpers.py riders_match <meta.json> <de_mag variant | -> <clean4:tau | ->   (502; "-" = must be absent)\n',
     '  def_day_helpers.py riders_match <meta.json> <de_mag variant | -> <clean4:tau | programme:tau[:tau_c4] |\n'
     '      file:npz:sha:tau[:tau_c4[:leN]] | ->   (502; 517 adds the programme / file forms; "-" = must be absent)\n'),
    ('            src, _, tau = sign.partition(":")\n'
     '            if not isinstance(ms, dict) or ms.get("source") != src or abs(float(ms.get("tau", -1)) - float(tau)) >= 1e-12:\n'
     '                print(f"de_mde_sign mismatch: meta {ms} vs {sign}", file=sys.stderr)\n'
     '                return 1\n',
     '            want = _mde_sign_request(sign)\n'
     '            bad = (not isinstance(ms, dict) or want is None\n'
     '                   or any(ms.get(k) != v if isinstance(v, str) else (ms.get(k) is None or abs(float(ms[k]) - v) >= 1e-12)\n'
     '                          for k, v in want.items()))\n'
     '            if bad:\n'
     '                print(f"de_mde_sign mismatch: meta {ms} vs {sign}", file=sys.stderr)\n'
     '                return 1\n'),
    ('def main() -> int:\n',
     'def _mde_sign_request(sign: str):\n'
     '    """517: the builder\'s VCC_DE_MDE_SIGN grammar -> the fields its config dict must carry (None = unparseable)."""\n'
     '    p = [x.strip() for x in sign.split(":")]\n'
     '    try:\n'
     '        if p[0] == "clean4" and len(p) == 2:\n'
     '            return {"source": "clean4", "tau": float(p[1])}\n'
     '        if p[0] == "programme" and len(p) in (2, 3):\n'
     '            return {"source": "programme", "tau": float(p[1]), "tau_c4": float(p[2]) if len(p) == 3 else 0.05}\n'
     '        if p[0] == "file" and len(p) in (4, 5, 6):\n'
     '            out = {"source": "file", "sha256": p[2], "tau": float(p[3]), "tau_c4": float(p[4]) if len(p) >= 5 else 0.05}\n'
     '            if len(p) == 6:\n'
     '                out["le"] = float(p[5][2:])\n'
     '            return out\n'
     '    except ValueError:\n'
     '        return None\n'
     '    return None\n\n\n'
     'def main() -> int:\n'),
]
for src, edits, name in ((s, E_DD, "def_day.sh"), (h, E_HP, "def_day_helpers.py")):
    for old, new in edits:
        if src.count(old) != 1:
            raise SystemExit(f"{name}: anchor count {src.count(old)} != 1: {old[:80]!r}")
for old, new in E_DD:
    s = s.replace(old, new)
for old, new in E_HP:
    h = h.replace(old, new)
compile(h, str(hp), "exec")
if "--check" in sys.argv:
    print("patches apply cleanly (not written)")
else:
    tmp_dd, tmp_hp = dd.with_name(dd.name + ".517tmp"), hp.with_name(hp.name + ".517tmp")
    tmp_dd.write_text(s); tmp_hp.write_text(h)
    os.replace(tmp_dd, dd); os.replace(tmp_hp, hp)
    print(f"patched {dd} and {hp}")
