#!/usr/bin/env python3
"""644 -- DEF_DAY_SAMELINE_SKIP=<ctx[,ctx]>: contexts that take the UNMATCHED same-line path whatever their identity call.
UNAPPLIED; two files; applies LAST on 10-06 (606's chain 641 -> 598 -> 612 -> 612dd -> 630 -> 630dd -> 642 -> 604 -> 644;
it anchors on 642's lines, so it refuses a tree without 642). v3 (607's review): names are whitespace-stripped in both
files, and the skip is logged as "treated as unmatched" (a programme source still writes that context at w 0).

Why (644's gap 1): a doubtful identity call in ONE context could only be handled by turning the same-line mode off for every
context, and sameline w 0 is not an off switch (`pt` still applies the bank's per-target signs, builder ~5693).

What it changes:
1. `scripts/def_sameline.py`: `--skip-contexts <ctx[,ctx]>` (default none). A skipped context is recorded as
   "none (644: …)" BEFORE any bank, 642 or 590 logic, so it is not `matched`. Everything downstream then treats it as
   an unmatched context:
   - without a programme source it is absent from the artifact (the builder: `.get` -> None, "not in the artifact --
     unchanged");
   - with `--programme-source` (a bank, or 604's identity prior) it gets the source's programme at w 0 with no rows and
     no signs, exactly like any unmatched context.
   Names outside the call are refused.
2. `scripts/def_day.sh`: the knob's names are checked against the manifest's contexts (else STOP). It is passed to
   def_sameline as `--skip-contexts`, logged at start and recorded in def_day.final.env. Unset: the call is unchanged
   (`${…:+…}` expands to nothing).

    uv run python experiments/644-runbook-amendments/patch_skip644.py --root <repo> --check
    uv run python experiments/644-runbook-amendments/patch_skip644.py --root <repo>
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0] / "642-hek293-contingency"))
import patch_hek293_642 as P642  # noqa: E402  (642's inserted lines are this patch's anchors)

MARK = "644: DEF_DAY_SAMELINE_SKIP"
FILES = (Path("scripts/def_sameline.py"), Path("scripts/def_day.sh"))

# ------------------------------------------------------------------------------------------------ def_sameline.py
SL_ARG_A = '    ap.add_argument("--min-programme-rows", type=int, default=150,\n'
SL_ARG_B = (
    '    ap.add_argument("--skip-contexts", default="",\n'
    '                    help="644: DEF_DAY_SAMELINE_SKIP: comma-separated contexts that take the unmatched path whatever their "\n'
    '                         "identity call (no rows, no pt signs; a --programme-source still gives them its programme at w 0, "\n'
    '                         "like any unmatched context). Default: none")\n'
) + SL_ARG_A
SL_VAL_A = (
    '    if a.contexts:\n'
    '        call = {c: v for c, v in call.items() if c in set(a.contexts.split(","))}\n'
)
SL_VAL_B = SL_VAL_A + (
    '    # 644: DEF_DAY_SAMELINE_SKIP: these contexts are never matched (a doubtful identity call)\n'
    '    _skip644 = {x.strip() for x in a.skip_contexts.split(",") if x.strip()}   # 607: "E, F" works\n'
    '    if _skip644 - set(call):\n'
    '        raise SystemExit(f"--skip-contexts names {sorted(_skip644 - set(call))}, which are not contexts of the call {sorted(call)}")\n'
)
SL_LOOP_A = P642.SL_LOOP_B[len(P642.SL_LOOP_A):].splitlines(keepends=True)[0]     # 642's first inserted line in the loop
SL_LOOP_B = (
    '        # 644: DEF_DAY_SAMELINE_SKIP: the unmatched path whatever the call (before 642 / 590 / any bank logic)\n'
    '        if c in _skip644:\n'
    '            rec["mode"] = "none (644: skipped by --skip-contexts / DEF_DAY_SAMELINE_SKIP; treated as unmatched)"\n'
    '            prov["contexts"][c] = rec\n'
    '            print(f"[sameline] {c}: {sample} skipped (--skip-contexts): treated as unmatched (644)", flush=True)\n'
    '            continue\n'
) + SL_LOOP_A

# ------------------------------------------------------------------------------------------------ def_day.sh
DAY_CALL_A = (
    'uv run python scripts/def_sameline.py --bundle "$B" --identity "$W/identity.json" --out "$W/sameline.npz" \\\n'
    '    ${DEF_DAY_FORCE:+--force "$DEF_DAY_FORCE"} ${DEF_DAY_PROGRAMME_SOURCE:+--programme-source "$DEF_DAY_PROGRAMME_SOURCE"} \\\n'
)
DAY_CALL_B = (
    '# 644: DEF_DAY_SAMELINE_SKIP=<ctx[,ctx]>: contexts that take the unmatched same-line path whatever their identity call\n'
    '# (whitespace is ignored: "E, F" == "E,F"; 607)\n'
    'SLSKIP="$(printf \'%s\' "${DEF_DAY_SAMELINE_SKIP:-}" | tr -d \'[:space:]\')"\n'
    'if [ -n "$SLSKIP" ]; then\n'
    '  for _c in $(printf \'%s\' "$SLSKIP" | tr \',\' \' \'); do\n'
    '    case ",$CTXS," in *",$_c,"*) ;; *) die "644: DEF_DAY_SAMELINE_SKIP names $_c, which is not a manifest context ($CTXS)" ;; esac\n'
    '  done\n'
    '  log "644: same-line skipped for $SLSKIP (def_sameline --skip-contexts): treated as unmatched"\n'
    'fi\n'
    'uv run python scripts/def_sameline.py --bundle "$B" --identity "$W/identity.json" --out "$W/sameline.npz" \\\n'
    '    ${DEF_DAY_FORCE:+--force "$DEF_DAY_FORCE"} ${DEF_DAY_PROGRAMME_SOURCE:+--programme-source "$DEF_DAY_PROGRAMME_SOURCE"} \\\n'
    '    ${SLSKIP:+--skip-contexts "$SLSKIP"} \\\n'
)
DAY_LOG_A = P642.DAY_LOG_B[len(P642.DAY_LOG_A):]
DAY_LOG_B = DAY_LOG_A + ('[ -n "${DEF_DAY_SAMELINE_SKIP:-}" ] && log "644: DEF_DAY_SAMELINE_SKIP=$DEF_DAY_SAMELINE_SKIP: these contexts '
                         'take the unmatched same-line path"\n')
DAY_ENV_A = P642.DAY_ENV_B[len(P642.DAY_ENV_A):]
DAY_ENV_B = DAY_ENV_A + ("# 644: the same-line skip, recorded\n"
                         "[ -n \"${DEF_DAY_SAMELINE_SKIP:-}\" ] && printf 'DEF_DAY_SAMELINE_SKIP=%q\\n' \"$DEF_DAY_SAMELINE_SKIP\" >> \"$W/def_day.final.env\"\n")

EDITS = {FILES[0]: [(SL_ARG_A, SL_ARG_B), (SL_VAL_A, SL_VAL_B), (SL_LOOP_A, SL_LOOP_B)],
         FILES[1]: [(DAY_CALL_A, DAY_CALL_B), (DAY_LOG_A, DAY_LOG_B), (DAY_ENV_A, DAY_ENV_B)]}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def plan(root: Path) -> dict[Path, str]:
    texts = {f: (root / f).read_text() for f in FILES}
    missing = [str(f) for f in FILES if P642.MARK not in texts[f]]
    if missing:
        raise SystemExit(f"642 is not applied to {', '.join(missing)}: 644's skip applies after 642 (last on 10-06)")
    marked = [str(f) for f in FILES if MARK in texts[f]]
    if marked:
        raise SystemExit("already patched: refusing a second apply" if len(marked) == len(FILES)
                         else f"partially patched ({', '.join(marked)}): restore and re-apply")
    new = {}
    for f, pairs in EDITS.items():
        t = texts[f]
        for old, rep in pairs:
            n = t.count(old)
            if n != 1:
                raise SystemExit(f"{f}: anchor found {n} times (need 1): {old.splitlines()[0][:90]!r}")
            t = t.replace(old, rep)
        if f.suffix == ".py":
            compile(t, str(f), "exec")
        new[f] = t
    return new


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    root = a.root.resolve()
    new = plan(root)
    for f in FILES:
        print(f"{f}\n  before {sha(root / f)}")
    if a.check:
        print("patches apply cleanly and compile (not written)")
        return
    for f in FILES:
        dst = root / f
        fd, tmp = tempfile.mkstemp(dir=dst.parent, prefix=f".{dst.name}.644-")
        with os.fdopen(fd, "w") as fh:
            fh.write(new[f])
        os.chmod(tmp, os.stat(dst).st_mode)
        os.replace(tmp, dst)
        print(f"{f}\n  after  {sha(dst)}")
    print("patched " + ", ".join(map(str, FILES)))


if __name__ == "__main__":
    main()
