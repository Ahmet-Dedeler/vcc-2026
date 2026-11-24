#!/usr/bin/env python3
"""690g — a dump hook for the October path (the PIE day-of step, 700/672). UNAPPLIED; insurance, not the main line.

Why: [bf003b]'s PIE blend (experiments/700-pie-transfer/make_pie_dump700.py) rewrites the stack dump's bulk directions
(norm-preserving) and can write a per-context bulk for 672's VCC_BULK_LOAD_CTX. On D/E/F both need the day's OWN stack dump,
which exists only inside def_day.sh, so the step has to run there, between the dump and the first level build.

`DEF_DAY_DUMP_HOOK=<committed script>`: after the stack dump def_day runs
    bash <script> <bundle> <work> <stack_dump.npz> <work>/identity.json <work>/dumphook
The script must write <work>/dumphook/<the stack dump's file name> (VCC_BULK_LOAD's format: bulk, carriers, perts, genes,
env, with env byte-identical so 598's tercilec still reads the dump key) and MAY write <work>/dumphook/bulk_ctx.npz (672's
format). scripts/def_dumphook_check.py validates both (identity of perts/genes/env, shapes, dtypes, finiteness, contexts)
and logs how far the bulk moved. The hook's dump then becomes VCC_BULK_LOAD; the stack dump stays on disk. bulk_ctx.npz
becomes VCC_BULK_LOAD_CTX (needs 672 in the builder). The final meta must carry that sha; final.env records the hook and
the npz; postbuild allows bulk_load_ctx(_sha256). Since 672's same-line rule (10-15), a context matched to a same-line bank
keeps the same-line rows; the builder logs which, and def_day repeats those lines.

Files: scripts/def_day.sh (6 insertions), scripts/def_dumphook_check.py (new, from this dir), scripts/october_recipe.py (switch
`dump_hook`, PENDING; the value is the script path, which must exist), docs/october-decisions.json (+1 PENDING),
scripts/test_october_recipe.py (fixture). Needs 690a applied first (its anchors). Unset knob = today's command byte for byte.

    uv run python experiments/690-def-readiness-1015/patch690g.py --root <repo> [--check]
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MARK = "690g: dump hook"

DOC_OLD = '''#   DEF_DAY_DE_MAG=<variant> (502, default off; e.g. dn_k250_zs): the zs exact-U pass (VCC_DE_MAG) on the level builds.
'''
DOC_NEW = '''#   DEF_DAY_DUMP_HOOK=<committed script> (690g: dump hook; default off): after the stack dump, `bash <script> <bundle> <work>
#    <dump.npz> <identity.json> <work>/dumphook` writes a replacement dump (same name and format, env unchanged) and maybe
#    bulk_ctx.npz (672's VCC_BULK_LOAD_CTX); scripts/def_dumphook_check.py validates both before any level build.
''' + DOC_OLD
V_OLD = '''# 502: the DE-channel riders (validated here, exported for the level builds after the dump)
'''
V_NEW = '''# 690g: dump hook -- DEF_DAY_DUMP_HOOK=<script> (validated here; it runs right after the stack dump)
if [ -n "${DEF_DAY_DUMP_HOOK:-}" ]; then
  [ -f "$DEF_DAY_DUMP_HOOK" ] || die "DEF_DAY_DUMP_HOOK=$DEF_DAY_DUMP_HOOK is not a file"
  git ls-files --error-unmatch "$DEF_DAY_DUMP_HOOK" >/dev/null 2>&1 || die "DEF_DAY_DUMP_HOOK must be a committed script (got $DEF_DAY_DUMP_HOOK)"
  [ -z "$(git status --porcelain -- "$DEF_DAY_DUMP_HOOK")" ] || die "DEF_DAY_DUMP_HOOK $DEF_DAY_DUMP_HOOK has uncommitted changes"
  [ -f scripts/def_dumphook_check.py ] || die "scripts/def_dumphook_check.py is missing (690g)"
fi
''' + V_OLD
H_OLD = '''DSHA="$(shasum -a 256 "$DUMP" | cut -d' ' -f1)"; log "dump ${DSHA:0:16}"
export VCC_BULK_LOAD="$DUMP:$DSHA" VCC_BULK_SHARPEN_Q="$SQ"
'''
H_NEW = '''DSHA="$(shasum -a 256 "$DUMP" | cut -d' ' -f1)"; log "dump ${DSHA:0:16}"
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
'''
C_OLD = '''if [ -n "${VCC_DE_UPRN:-}" ]; then   # 630dd: rank-neutral UP pass-through -- the final meta must carry the variant and the reviewed module
'''
C_NEW = '''if [ -n "${BLC:-}" ]; then   # 690g: dump hook -- the final meta must carry the per-context bulk it was built with
  python3 - "$FINAL.h5ad.meta.json" "${BLC##*:}" <<'PY' || die "690g: the final meta does not carry bulk_load_ctx_sha256 ${BLC##*:}"
import json, sys
c = json.load(open(sys.argv[1]))["config"]
sys.exit(0 if c.get("bulk_load_ctx_sha256") == sys.argv[2] else 1)
PY
  log "690g: per-context outcome (builder): $(grep -hE '672|bulk.load.ctx' "$W/logs/build.log" 2>/dev/null | tail -4 | tr '\\n' ' ' | cut -c1-400)"
fi
''' + C_OLD
F_OLD = '''[ -n "${VCC_DE_MAG:-}" ] && printf 'VCC_DE_MAG=%s\\n' "$VCC_DE_MAG" >> "$W/def_day.final.env"
'''
F_NEW = '''[ -n "${DEF_DAY_DUMP_HOOK:-}" ] && printf 'DEF_DAY_DUMP_HOOK=%q\\n' "$DEF_DAY_DUMP_HOOK" >> "$W/def_day.final.env"   # 690g: dump hook
[ -n "${BLC:-}" ] && printf 'VCC_BULK_LOAD_CTX=%s\\n' "$BLC" >> "$W/def_day.final.env"   # 690g: dump hook
''' + F_OLD
A_OLD = '''[ -n "${VCC_DE_MAG:-}" ] && ALLOW="$ALLOW,de_mag"
'''
A_NEW = A_OLD + '''[ -n "${BLC:-}" ] && ALLOW="$ALLOW,bulk_load_ctx,bulk_load_ctx_sha256"   # 690g: dump hook
'''
DAY_EDITS = [("DOC", DOC_OLD, DOC_NEW), ("V", V_OLD, V_NEW), ("H", H_OLD, H_NEW), ("C", C_OLD, C_NEW), ("F", F_OLD, F_NEW),
             ("A", A_OLD, A_NEW)]

R_EDITS = [
    ("ENV_OF", '''    "de_dnmag": "DEF_DAY_DE_DNMAG",
}
''', '''    "de_dnmag": "DEF_DAY_DE_DNMAG",
    "dump_hook": "DEF_DAY_DUMP_HOOK",   # 690g: a committed script run after the stack dump (the PIE day-of step)
}
'''),
    ("resolve", '''    _problems_690a(env, problems)
''', '''    _problems_690a(env, problems)
    _problems_690g(env, problems)
'''),
    ("fn", '''def installed_problems(env: dict, flags: str) -> list:
''', '''def _problems_690g(env: dict, problems: list) -> None:
    """690g: the dump hook must name an existing script in the repo (def_day.sh also requires it to be committed)."""
    if "DEF_DAY_DUMP_HOOK" in env and not (ROOT / env["DEF_DAY_DUMP_HOOK"]).is_file():
        problems.append(f"dump_hook {env['DEF_DAY_DUMP_HOOK']!r} is not a file in the repo")


def installed_problems(env: dict, flags: str) -> list:
'''),
]
T_EDITS = [("OFF538", '''          ("de_uprn", "off"), ("mde_sign_keygate", "off"), ("de_dnmag", "off")]
''', '''          ("de_uprn", "off"), ("mde_sign_keygate", "off"), ("de_dnmag", "off"),
          ("dump_hook", "off")]   # 690g
''')]
SWITCH = {"dump_hook": {
    "status": "PENDING", "value": None,
    "decided_by": "a PIE arm's A/B/C board read (700) plus its held-out read; the value is the committed day-of hook script",
    "if_on": "a committed script path, e.g. experiments/700-pie-transfer/<hook>.sh, honouring 690g's interface",
    "otherwise": "off",
    "needs": "def_day 690g; scripts/def_dumphook_check.py; builder 672 if the hook writes bulk_ctx.npz. The PIE assets "
             "(~75 GB, HF_HUB_OFFLINE=1, warm cache) on disk before 10-21",
    "held_out": "700's held-out reads (10-15: bulk blend +0.0008…+0.0017 predicted; PIE's signs lose to ours)",
    "added_by": "690g (2026-10-15)"}}


def apply(text: str, edits: list, label: str) -> tuple[str, list]:
    errs = []
    for name, old, new in edits:
        n = text.count(old)
        if n != 1:
            errs.append(f"{label} {name}: anchor found {n} times (want 1)")
            continue
        text = text.replace(old, new)
    return text, errs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    r = a.root.resolve()
    day, rec, dec, tst = (r / "scripts/def_day.sh", r / "scripts/october_recipe.py", r / "docs/october-decisions.json",
                          r / "scripts/test_october_recipe.py")
    dd, rr, tt = day.read_text(), rec.read_text(), tst.read_text()
    dj = json.loads(dec.read_text())
    if MARK in dd or "_problems_690g" in rr or "dump_hook" in dj["switches"] or (r / "scripts/def_dumphook_check.py").exists():
        print("REFUSED: 690g is already applied", file=sys.stderr)
        return 1
    errs = []
    if "_problems_690a" not in rr:
        errs.append("690a is not applied (690g's recipe anchors come from it)")
    dd2, e1 = apply(dd, DAY_EDITS, "def_day.sh")
    rr2, e2 = apply(rr, R_EDITS, "october_recipe.py")
    tt2, e3 = apply(tt, T_EDITS, "test_october_recipe.py")
    errs += e1 + e2 + e3
    if errs:
        print("REFUSED:\n  " + "\n  ".join(errs), file=sys.stderr)
        return 1
    sw = dj["switches"]
    tail = [k for k in sw if sw[k].get("status") == "DAY-OF"]
    new = {k: v for k, v in sw.items() if k not in tail}
    new.update(SWITCH)
    new.update({k: sw[k] for k in tail})
    dj["switches"] = new
    if a.check:
        print("CHECK OK: 6 def_day.sh anchors, 3 recipe anchors, 1 test anchor, +scripts/def_dumphook_check.py, "
              "+1 PENDING switch; nothing written")
        return 0
    day.write_text(dd2)
    rec.write_text(rr2)
    tst.write_text(tt2)
    dec.write_text(json.dumps(dj, indent=1) + "\n")
    shutil.copy2(HERE / "def_dumphook_check.py", r / "scripts/def_dumphook_check.py")
    print(f"APPLIED 690g to {r}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
