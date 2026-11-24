#!/usr/bin/env python3
"""554 -- the emitter-dose guard, as an UNAPPLIED patch (the coordinator applies it after the 10-05 bundle, with an inert KA).

Why (553, `experiments/553-r7b-harness/RESULTS.md`, "A production robustness fact for October"): the K-class emitter's
library classes are the context's control-library quartiles. On H1, the block-prior / MDE / uncalled-gene-hook steering
dose collapses to 0 WITHOUT any error somewhere between control-library q75/q25 3.5 and 5, and the placement crashes by
7.24; the MDE / hook dose also erodes quietly when the sub-gate mass per cell is low. A silent dose-0 build would ship the
final entry without its DE design. A/B/C sit at 1.57-1.95.

What it changes (three files; every new knob default-off; silent on PASS, so today's text and exit codes are unchanged on
every input that passes):

1. `scripts/def_preflight.py`: per context, the control-library q25 / q75 / q75/q25 from the library sizes it already
   streams, recorded in the preflight JSON. A FAIL row only when q75/q25 > LIB_SPREAD_MAX (3.0), unless
   `DEF_DAY_ALLOW_WIDE_LIBS=1`, which makes it a passing "ALLOWED" row (recorded in the JSON).
2. `scripts/def_postbuild.py`: `--rc-audit <jsonl>` (the builder's per-block `VCC_RC_AUDIT_DUMP` record, which the meta
   only summarises POOLED over all blocks -- a one-context collapse does not move the pooled median). Per context: median
   `prior_lambda`, `mde_lambda`, `call_set_missed`; FAIL row if any median is past its floor, unless
   `DEF_DAY_ALLOW_LOW_DOSE=1` (passing "ALLOWED" row). The dump must belong to the meta (block count and the pooled
   medians recomputed from it equal the meta's). Without `--rc-audit` the meta's pooled medians are checked (a global
   collapse only). The JSON record gains `dose_guard` (per-context values, floors, margins).
3. `scripts/def_day.sh` (the 538-patched file): every level build writes its own dump
   (`VCC_RC_AUDIT_DUMP=$W/$N-<tag>.rc_audit.jsonl`; the builder's own docs: "diagnostics only, the emission is untouched"),
   the final build's dump is linked as `$W/${N}_final.rc_audit.jsonl` and passed to def_postbuild as `--rc-audit`;
   non-final dumps are deleted with their h5ads; both overrides are logged at start and written to `def_day.final.env`.

    uv run python experiments/554-dose-guard/patch_guard554.py --root <repo> --check     # anchors only, writes nothing
    uv run python experiments/554-dose-guard/patch_guard554.py --root <repo>             # apply (refuses a second time)

It refuses to run on the live repo unless `--root` names it explicitly (the coordinator's apply step).
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys
import tempfile
from pathlib import Path

MARK = "554: the K-class emitter"           # present in every patched file; a second apply is refused
FILES = (Path("scripts/def_preflight.py"), Path("scripts/def_postbuild.py"), Path("scripts/def_day.sh"))

# ------------------------------------------------------------------------------------------------ def_preflight.py
PRE_CONST_A = "MAX_COUNTS_PER_CELL = 1_000_000\n"
PRE_CONST_B = (
    "MAX_COUNTS_PER_CELL = 1_000_000\n"
    "# 554: the K-class emitter's library classes are the control quartiles; above q75/q25 ~3.5 its steering dose (CLEAN4\n"
    "# block prior, MDE, uncalled-gene hook) collapses to 0 WITHOUT an error (553, H1 sweep), by ~7 the placement crashes.\n"
    "# A/B/C 1.57-1.95, H1 1.58, R7b as built 7.24. DEF_DAY_ALLOW_WIDE_LIBS=1 is the explicit override.\n"
    "LIB_SPREAD_MAX = 3.0\n"
)
PRE_SPREAD_A = (
    "            max_library=float(lib.max()), mean_library=round(float(lib.mean()), 1),\n"
    "        )\n"
)
PRE_SPREAD_B = PRE_SPREAD_A + (
    "        # 554: the K-class emitter's depth cliff (LIB_SPREAD_MAX). Always recorded; a row only when it trips, so a passing\n"
    "        # bundle's report text is today's.\n"
    "        _q25, _q75 = (float(v) for v in np.quantile(lib, [0.25, 0.75]))\n"
    "        _spread = _q75 / _q25 if _q25 > 0 else float(\"inf\")\n"
    "        rec.update(library_q25=_q25, library_q75=_q75, library_spread_q75_q25=round(_spread, 4),\n"
    "                   library_spread_max=LIB_SPREAD_MAX)\n"
    "        if _spread > LIB_SPREAD_MAX:\n"
    "            _allow = os.environ.get(\"DEF_DAY_ALLOW_WIDE_LIBS\", \"\") == \"1\"\n"
    "            rec[\"library_spread_allowed\"] = _allow\n"
    "            checks.append((f\"control library q75/q25 <= {LIB_SPREAD_MAX} (554: above ~3.5 the emitter's steering dose \"\n"
    "                           \"collapses silently)\", _allow,\n"
    "                           f\"{_spread:.3f} (q25 {_q25:.0f}, q75 {_q75:.0f})\"\n"
    "                           + (\" -- ALLOWED by DEF_DAY_ALLOW_WIDE_LIBS=1\" if _allow else \"\")))\n"
)

# ------------------------------------------------------------------------------------------------ def_postbuild.py
POST_IMPORT_A = "import json\nfrom pathlib import Path\n"
POST_IMPORT_B = "import json\nimport os  # 554: DEF_DAY_ALLOW_LOW_DOSE\nfrom pathlib import Path\n"
POST_FUNC_A = "def main() -> None:\n"
POST_FUNC_B = '''# 554: the K-class emitter can deliver its CLEAN4 block prior and coarse MDE (with the uncalled-gene hook riding on it) at
# dose 0 WITHOUT an error (553). Floors from 554's calibration (its RESULTS.md): healthy per-context medians are prior 1.0,
# MDE 0.951-1.0 (A / D the lowest), call-set misses 4-21.5; 429a's plain-MDE failure had A's MDE median 0.0.
DOSE_PRIOR_MIN = 0.9
DOSE_MDE_MIN = 0.9
DOSE_MISSED_MAX = 100.0
DOSE_KEYS = ("prior_lambda", "mde_lambda", "call_set_missed")


def dose_guard(meta: dict, rc_audit: Path | None, contexts: list[str]) -> dict:
    """Per-context medians of the emitter's delivered steering dose, from the builder's per-block VCC_RC_AUDIT_DUMP record
    (or, without one, the meta's POOLED medians: a global collapse only). Returns the record and the check rows (a row only
    when a floor trips or the dump does not belong to the meta, so a passing build's text is today's)."""
    rc = meta.get("output", {}).get("rc_callset")
    floors = {"prior_lambda": DOSE_PRIOR_MIN, "mde_lambda": DOSE_MDE_MIN, "call_set_missed": DOSE_MISSED_MAX}
    rec: dict = {"floors": floors, "source": None, "per_context": {}, "rows": []}
    if not isinstance(rc, dict):
        rec["source"] = "none (no rc_callset record: not a K-class build)"
        return rec
    allow = os.environ.get("DEF_DAY_ALLOW_LOW_DOSE", "") == "1"
    problems: list[str] = []
    if rc_audit is not None:
        rows = [json.loads(ln) for ln in rc_audit.read_text().splitlines() if ln.strip()]
        rec["source"] = f"rc_audit {rc_audit.name} ({len(rows)} blocks)"
        rec["rc_audit_sha256"] = sha256_file(rc_audit)
        # the dump must be THIS build's: same block count, and its pooled medians are the meta's (same values, np.median)
        same = len(rows) == rc.get("blocks")
        for k in DOSE_KEYS:
            if k in rc and same:
                v = [float(r[k]) for r in rows if k in r]
                same = bool(v) and float(np.median(v)) == float(rc[k]["median"])
        if not same:
            problems.append(f"--rc-audit {rc_audit.name} does not belong to this meta (block count or pooled medians differ)")
        for c in contexts:
            sub = [r for r in rows if str(r.get("context")) == str(c)]
            rec["per_context"][c] = {"blocks": len(sub), **{k: (float(np.median([float(r[k]) for r in sub if k in r]))
                                                               if any(k in r for r in sub) else None) for k in DOSE_KEYS}}
            if not sub:
                problems.append(f"{c}: no blocks in the dump")
    else:
        rec["source"] = "meta pooled medians (no --rc-audit: a one-context collapse is NOT visible here)"
        rec["per_context"]["pooled"] = {"blocks": rc.get("blocks"),
                                        **{k: (float(rc[k]["median"]) if k in rc else None) for k in DOSE_KEYS}}
    for c, v in rec["per_context"].items():
        m = {}
        for k in DOSE_KEYS:
            if v.get(k) is None:
                continue
            m[k] = (v[k] - floors[k]) if k != "call_set_missed" else (floors[k] - v[k])
            if m[k] < 0:
                problems.append(f"{c}: median {k} {v[k]:.4g} past its floor {floors[k]:g}")
        v["margin"] = m
    rec["allowed"] = allow if problems else None
    if problems:
        rec["rows"].append(("emitter steering dose per context: median prior_lambda >= "
                            f"{DOSE_PRIOR_MIN}, mde_lambda >= {DOSE_MDE_MIN}, call_set_missed <= {DOSE_MISSED_MAX:g} (554)",
                            allow, "; ".join(problems) + (" -- ALLOWED by DEF_DAY_ALLOW_LOW_DOSE=1" if allow else "")))
    return rec


def main() -> None:
'''
POST_ARG_A = '    ap.add_argument("--out", type=Path, default=None)\n    args = ap.parse_args()\n'
POST_ARG_B = (
    '    ap.add_argument("--rc-audit", type=Path, default=None,\n'
    '                    help="554: the build\'s per-block VCC_RC_AUDIT_DUMP jsonl; per-context steering-dose floors. "\n'
    '                         "Absent: the meta\'s pooled medians (a global collapse only); silent on PASS either way")\n'
    '    ap.add_argument("--out", type=Path, default=None)\n    args = ap.parse_args()\n'
)
POST_GUARD_A = (
    "        {c: ctx_out.get(c, {}).get(\"emitted_rows\") for c in sorted(set(ctx_out) | set(manifest[\"contexts\"]))})\n"
    "\n"
    "    print(f\"meta      {args.meta}\")\n"
)
POST_GUARD_B = (
    "        {c: ctx_out.get(c, {}).get(\"emitted_rows\") for c in sorted(set(ctx_out) | set(manifest[\"contexts\"]))})\n"
    "    # 554: the K-class emitter's delivered steering dose (dose_guard). A row only when it trips or is ALLOWED.\n"
    "    dose = dose_guard(meta, args.rc_audit, [str(c) for c in manifest[\"contexts\"]])\n"
    "    for _name, _ok, _value in dose[\"rows\"]:\n"
    "        chk(_name, _ok, _value)\n"
    "\n"
    "    print(f\"meta      {args.meta}\")\n"
)
POST_OUT_A = (
    "            {\"meta\": str(args.meta), \"checks\": [{\"check\": c, \"pass\": p, \"value\": v}\n"
    "                                                for c, p, v in checks]}, indent=2) + \"\\n\")\n"
)
POST_OUT_B = (
    "            {\"meta\": str(args.meta), \"checks\": [{\"check\": c, \"pass\": p, \"value\": v}\n"
    "                                                for c, p, v in checks],\n"
    "             \"dose_guard\": {k: v for k, v in dose.items() if k != \"rows\"}}, indent=2) + \"\\n\")\n"
)

# ------------------------------------------------------------------------------------------------ def_day.sh (538-patched)
DAY_BUILD_A = (
    "  rm -f \"$h5\" \"$h5.meta.json\"\n"
    "  \"${DRV[@]}\" --name \"$N-$tag\" --stages build > \"$W/driver-$tag.log\" 2>&1\n"
)
DAY_BUILD_B = (
    "  rm -f \"$h5\" \"$h5.meta.json\" \"$W/$N-$tag.rc_audit.jsonl\"\n"
    "  # 554: the K-class emitter's per-block steering-dose record (diagnostics only; the emission is untouched), read by\n"
    "  # def_postbuild.py --rc-audit for the final build\n"
    "  VCC_RC_AUDIT_DUMP=\"$W/$N-$tag.rc_audit.jsonl\" \"${DRV[@]}\" --name \"$N-$tag\" --stages build > \"$W/driver-$tag.log\" 2>&1\n"
)
DAY_PRUNE_A = (
    "for P in \"${POINTS[@]}\"; do [ \"${P%%.ntclevel.json:*}\" != \"$FINAL\" ] && rm -f \"${P%%.ntclevel.json:*}.h5ad\" \"${P%%.ntclevel.json:*}.h5ad.meta.json\"; done\n"
    "ln -f \"$FINAL.h5ad\" \"$W/${N}_final.h5ad\" || die \"ln of the final h5ad failed\"\n"
    "ln -f \"$FINAL.h5ad.meta.json\" \"$W/${N}_final.h5ad.meta.json\" || die \"ln of the final meta failed\"\n"
)
DAY_PRUNE_B = (
    "for P in \"${POINTS[@]}\"; do [ \"${P%%.ntclevel.json:*}\" != \"$FINAL\" ] && rm -f \"${P%%.ntclevel.json:*}.h5ad\" \"${P%%.ntclevel.json:*}.h5ad.meta.json\" \"${P%%.ntclevel.json:*}.rc_audit.jsonl\"; done\n"
    "ln -f \"$FINAL.h5ad\" \"$W/${N}_final.h5ad\" || die \"ln of the final h5ad failed\"\n"
    "ln -f \"$FINAL.h5ad.meta.json\" \"$W/${N}_final.h5ad.meta.json\" || die \"ln of the final meta failed\"\n"
    "# 554: the final build's per-block dose record (def_postbuild.py --rc-audit)\n"
    "ln -f \"$FINAL.rc_audit.jsonl\" \"$W/${N}_final.rc_audit.jsonl\" || die \"554: the final build wrote no rc_audit dump\"\n"
)
DAY_ENV_A = "[ -n \"${VCC_RC_CALLSET_BLOCK:-}\" ] && printf 'VCC_RC_CALLSET_BLOCK=%s\\n' \"$VCC_RC_CALLSET_BLOCK\" >> \"$W/def_day.final.env\"\n"
DAY_ENV_B = DAY_ENV_A + (
    "# 554: the dose guard's overrides, recorded\n"
    "[ \"${DEF_DAY_ALLOW_WIDE_LIBS:-}\" = 1 ] && printf 'DEF_DAY_ALLOW_WIDE_LIBS=1\\n' >> \"$W/def_day.final.env\"\n"
    "[ \"${DEF_DAY_ALLOW_LOW_DOSE:-}\" = 1 ] && printf 'DEF_DAY_ALLOW_LOW_DOSE=1\\n' >> \"$W/def_day.final.env\"\n"
)
DAY_POST_A = (
    "uv run python scripts/def_postbuild.py --meta \"$W/${N}_final.h5ad.meta.json\" --bundle \"$B\" --cd4-root \"$W/cd4_compact\" \\\n"
    "    --incumbent \"$INC\" --allow-config-diff \"$ALLOW\" --out \"$W/postbuild.json\" >> \"$LOG\" 2>&1 \\\n"
)
DAY_POST_B = (
    "uv run python scripts/def_postbuild.py --meta \"$W/${N}_final.h5ad.meta.json\" --bundle \"$B\" --cd4-root \"$W/cd4_compact\" \\\n"
    "    --incumbent \"$INC\" --allow-config-diff \"$ALLOW\" --rc-audit \"$W/${N}_final.rc_audit.jsonl\" \\\n"
    "    --out \"$W/postbuild.json\" >> \"$LOG\" 2>&1 \\\n"
)
DAY_LOG_A = "log \"start"
DAY_LOG_ANCHOR = "# 554-start-anchor"   # resolved at run time: the first line that begins with  log "start


def edits() -> dict[Path, list[tuple[str, str]]]:
    return {
        FILES[0]: [(PRE_CONST_A, PRE_CONST_B), (PRE_SPREAD_A, PRE_SPREAD_B)],
        FILES[1]: [(POST_IMPORT_A, POST_IMPORT_B), (POST_FUNC_A, POST_FUNC_B), (POST_ARG_A, POST_ARG_B), (POST_GUARD_A, POST_GUARD_B),
                   (POST_OUT_A, POST_OUT_B)],
        FILES[2]: [(DAY_BUILD_A, DAY_BUILD_B), (DAY_PRUNE_A, DAY_PRUNE_B), (DAY_ENV_A, DAY_ENV_B), (DAY_POST_A, DAY_POST_B)],
    }


def day_start_log(text: str) -> tuple[str, str]:
    """The override log lines go right after def_day's first `log "start` line (its first log call)."""
    lines = text.splitlines(keepends=True)
    for i, ln in enumerate(lines):
        if ln.startswith(DAY_LOG_A):
            add = ("[ \"${DEF_DAY_ALLOW_WIDE_LIBS:-}\" = 1 ] && log \"554: DEF_DAY_ALLOW_WIDE_LIBS=1: preflight's control-library "
                   "spread floor is overridden\"\n"
                   "[ \"${DEF_DAY_ALLOW_LOW_DOSE:-}\" = 1 ] && log \"554: DEF_DAY_ALLOW_LOW_DOSE=1: postbuild's steering-dose "
                   "floors are overridden\"\n")
            return ln, ln + add
    raise SystemExit("def_day.sh: no line starting with  log \"start  (anchor miss)")


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    root = a.root.resolve()
    texts = {f: (root / f).read_text() for f in FILES}
    marked = [str(f) for f in FILES if MARK in texts[f]]
    if marked:
        raise SystemExit(f"already patched ({', '.join(marked)}): refusing a second apply" if len(marked) == len(FILES)
                         else f"partially patched ({', '.join(marked)}): restore the three files from git and re-apply")
    new = {}
    for f, pairs in edits().items():
        t = texts[f]
        for old, rep in pairs:
            n = t.count(old)
            if n != 1:
                raise SystemExit(f"{f}: anchor found {n} times (need 1): {old.splitlines()[0][:90]!r}")
            t = t.replace(old, rep)
        if f == FILES[2]:
            old, rep = day_start_log(t)
            t = t.replace(old, rep, 1)
        new[f] = t
    for f in FILES:
        print(f"{f}\n  before {sha(root / f)}")
    if a.check:
        print("patches apply cleanly (not written)")
        return
    for f in FILES:                               # atomic per file: write a sibling temp, then os.replace
        dst = root / f
        fd, tmp = tempfile.mkstemp(dir=dst.parent, prefix=f".{dst.name}.554-")
        with os.fdopen(fd, "w") as fh:
            fh.write(new[f])
        os.chmod(tmp, os.stat(dst).st_mode)
        os.replace(tmp, dst)
        print(f"{f}\n  after  {sha(dst)}")
    print("patched " + ", ".join(str(f) for f in FILES))


if __name__ == "__main__":
    main()
