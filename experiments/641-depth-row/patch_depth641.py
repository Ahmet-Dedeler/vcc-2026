#!/usr/bin/env python3
"""641 -- the day-of DEPTH row, as an UNAPPLIED patch on top of 554's def_preflight change (the coordinator applies it
after 554, any day from 10-05; inert on every bundle whose controls sit at A/B/C's depth).

Why (Arc, "Behind the Data" 2026-10-01, Fig. 2; experiment 640's notes): the full 2026 data has a median 29,482 UMIs per
cell, but the released A/B/C controls sit at 20,109 / 19,946 / 20,034 (within 0.4 % of each other), so Arc downsampled
every released context to one ~20k median. Every depth-tuned constant of our recipe was chosen on A/B/C at that depth,
above all the DE call-set size (521: the global call-set size is regime-sensitive; K 1,600 is +0.0048 at H1 f 0.35 and
-0.015 at f 1). If D/E/F arrive at another depth (an un-downsampled release would be ~1.47x), the build still runs, every
per-gene quantity is re-fitted from the new controls, but the call-set size no longer sits in the regime it was tuned
in. That is a human decision, not a silent one.

What it changes (two files; the knob is default-off; silent on PASS, so today's text and exit codes are unchanged on every
bundle whose control medians are within 15 % of 20,030):

1. `scripts/def_preflight.py` (554-patched): per context, the control-library median over A/B/C's (LIB_MEDIAN_REF 20,030,
   the mean of A/B/C's three medians), always recorded in the preflight JSON (`library_median_ratio`). A FAIL row only
   when |ratio - 1| > LIB_MEDIAN_TOL (0.15), unless `DEF_DAY_ALLOW_DEPTH=1`, which makes it a passing "ALLOWED" row
   (JSON `library_median_allowed`).
2. `scripts/def_day.sh` (538 + 554-patched): the override is logged at start and written to `def_day.final.env`, next to
   554's two.

Why 15 %: A/B/C's medians agree to 0.4 %, so a correctly downsampled D/E/F context should read within ~1 %; an
un-downsampled one reads +47 %. 15 % is ~37x A/B/C's spread and a third of the raw-data gap. Inside it, the MDE (which
scales as depth^-1/2) moves by <= 7.5 %, inside 554's 10 % dose tolerance; outside it, the DE regime has moved enough to
question the board-tuned call-set size.

    uv run python experiments/641-depth-row/patch_depth641.py --root <repo> --check     # anchors only, writes nothing
    uv run python experiments/641-depth-row/patch_depth641.py --root <repo>             # apply (refuses a second time)

It refuses unless 554 is applied (its anchors are 554's inserted lines) and refuses a second or partial apply.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import tempfile
from pathlib import Path

MARK = "641: the depth row"               # present in every patched file; a second apply is refused
FILES = (Path("scripts/def_preflight.py"), Path("scripts/def_day.sh"))

# ------------------------------------------------------------------------------------------------ def_preflight.py
PRE_CONST_A = (
    "# A/B/C 1.57-1.95, H1 1.58, R7b as built 7.24. DEF_DAY_ALLOW_WIDE_LIBS=1 is the explicit override.\n"
    "LIB_SPREAD_MAX = 3.0\n"
)
PRE_CONST_B = PRE_CONST_A + (
    "# 641: the depth row. Arc downsampled every released context to one median depth (A/B/C control medians 20,109 /\n"
    "# 19,946 / 20,034; the full data's median is 29,482 UMI/cell, Arc's 2026-10-01 post, Fig. 2). The recipe's\n"
    "# depth-tuned constants (the DE call-set size above all) were chosen on A/B/C at that depth. DEF_DAY_ALLOW_DEPTH=1 is\n"
    "# the explicit override (experiments/641-depth-row).\n"
    "LIB_MEDIAN_REF = 20030.0\n"
    "LIB_MEDIAN_TOL = 0.15\n"
)
PRE_ROW_A = "                           + (\" -- ALLOWED by DEF_DAY_ALLOW_WIDE_LIBS=1\" if _allow else \"\")))\n"
PRE_ROW_B = PRE_ROW_A + (
    "        # 641: the depth row (LIB_MEDIAN_REF / LIB_MEDIAN_TOL). Always recorded; a row only when it trips, so a passing\n"
    "        # bundle's report text is today's.\n"
    "        _ratio = float(np.median(lib)) / LIB_MEDIAN_REF\n"
    "        rec.update(library_median_ref=LIB_MEDIAN_REF, library_median_ratio=round(_ratio, 4),\n"
    "                   library_median_tol=LIB_MEDIAN_TOL)\n"
    "        if abs(_ratio - 1.0) > LIB_MEDIAN_TOL:\n"
    "            _allow_d = os.environ.get(\"DEF_DAY_ALLOW_DEPTH\", \"\") == \"1\"\n"
    "            rec[\"library_median_allowed\"] = _allow_d\n"
    "            checks.append((f\"control library median within {LIB_MEDIAN_TOL:.0%} of A/B/C's {LIB_MEDIAN_REF:,.0f} \"\n"
    "                           \"(641: the recipe's depth-tuned constants were chosen at that depth)\", _allow_d,\n"
    "                           f\"{float(np.median(lib)):,.0f} (x{_ratio:.3f})\"\n"
    "                           + (\" -- ALLOWED by DEF_DAY_ALLOW_DEPTH=1\" if _allow_d else \"\")))\n"
)

# ------------------------------------------------------------------------------------------------ def_day.sh (538 + 554)
DAY_LOG_A = ("[ \"${DEF_DAY_ALLOW_LOW_DOSE:-}\" = 1 ] && log \"554: DEF_DAY_ALLOW_LOW_DOSE=1: postbuild's steering-dose "
             "floors are overridden\"\n")
DAY_LOG_B = DAY_LOG_A + (
    "[ \"${DEF_DAY_ALLOW_DEPTH:-}\" = 1 ] && log \"641: the depth row: DEF_DAY_ALLOW_DEPTH=1: preflight's control-library "
    "median check is overridden\"\n"
)
DAY_ENV_A = "[ \"${DEF_DAY_ALLOW_LOW_DOSE:-}\" = 1 ] && printf 'DEF_DAY_ALLOW_LOW_DOSE=1\\n' >> \"$W/def_day.final.env\"\n"
DAY_ENV_B = DAY_ENV_A + (
    "# 641: the depth row's override, recorded\n"
    "[ \"${DEF_DAY_ALLOW_DEPTH:-}\" = 1 ] && printf 'DEF_DAY_ALLOW_DEPTH=1\\n' >> \"$W/def_day.final.env\"\n"
)

NEED_554 = "554: the K-class emitter"


def edits() -> dict[Path, list[tuple[str, str]]]:
    return {FILES[0]: [(PRE_CONST_A, PRE_CONST_B), (PRE_ROW_A, PRE_ROW_B)],
            FILES[1]: [(DAY_LOG_A, DAY_LOG_B), (DAY_ENV_A, DAY_ENV_B)]}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def plan(root: Path) -> dict[Path, str]:
    """The patched texts, or SystemExit with the reason; writes nothing."""
    texts = {f: (root / f).read_text() for f in FILES}
    missing = [str(f) for f in FILES if NEED_554 not in texts[f]]
    if missing:
        raise SystemExit(f"554 is not applied to {', '.join(missing)}: apply patch_guard554.py first (641 anchors on it)")
    marked = [str(f) for f in FILES if MARK in texts[f]]
    if marked:
        raise SystemExit(f"already patched ({', '.join(marked)}): refusing a second apply" if len(marked) == len(FILES)
                         else f"partially patched ({', '.join(marked)}): restore the two files and re-apply")
    new = {}
    for f, pairs in edits().items():
        t = texts[f]
        for old, rep in pairs:
            n = t.count(old)
            if n != 1:
                raise SystemExit(f"{f}: anchor found {n} times (need 1): {old.splitlines()[0][:90]!r}")
            t = t.replace(old, rep)
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
        print("patches apply cleanly (not written)")
        return
    for f in FILES:                               # atomic per file: write a sibling temp, then os.replace
        dst = root / f
        fd, tmp = tempfile.mkstemp(dir=dst.parent, prefix=f".{dst.name}.641-")
        with os.fdopen(fd, "w") as fh:
            fh.write(new[f])
        os.chmod(tmp, os.stat(dst).st_mode)
        os.replace(tmp, dst)
        print(f"{f}\n  after  {sha(dst)}")
    print("patched " + ", ".join(str(f) for f in FILES))


if __name__ == "__main__":
    main()
