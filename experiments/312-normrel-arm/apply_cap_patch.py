#!/usr/bin/env python3
"""Astra's variant: VCC_NORMREL_MODE = plain | cap.

  plain (default, today's shipped arm):  w_j(t) proportional to r^exp
  cap   (the new hypothesis):            w_j(t) proportional to max(1, r)^exp

with r = ||row_j(t)|| / that donor's own panel median. With exp < 0, `cap` SUPPRESSES donor rows that
are oversized relative to their own panel and leaves everything at or below the median untouched,
whereas `plain` additionally UP-WEIGHTS small rows. The two hypotheses are different mechanisms and
today's board point cannot separate them: 143 is the precedent for `pds` rising while signed credit falls.

DO NOT APPLY while a build is running -- the builder is read by the live build process.
"""
import hashlib, pathlib, shutil, sys

B = pathlib.Path("~/Code/virtual-cell-challenge-2026/scripts/make_submission_pm_cd4.py")
src = B.read_text()
if "VCC_NORMREL_MODE" in src:
    print("already patched"); sys.exit(0)

anchor = 'NORMREL_EXP_ENV = os.environ.get("VCC_NORMREL_EXP", "").strip()'
assert anchor in src
src = src.replace(anchor, '''NORMREL_MODE = (os.environ.get("VCC_NORMREL_MODE", "plain").strip() or "plain").lower()
if NORMREL_MODE not in ("plain", "cap"):
    raise SystemExit(f"VCC_NORMREL_MODE must be 'plain' or 'cap'; got {NORMREL_MODE!r}")
''' + anchor, 1)

old = """            if med > 0.0 and nv > 0.0:
                w[name] = float((nv / med) ** exponent) * float(
                    active_weights.get(name, 1.0)
                )"""
new = """            if med > 0.0 and nv > 0.0:
                ratio = nv / med
                if NORMREL_MODE == "cap":
                    # suppress only the OVERSIZED rows; never up-weight a small one
                    ratio = max(1.0, ratio)
                w[name] = float(ratio ** exponent) * float(
                    active_weights.get(name, 1.0)
                )"""
assert old in src and src.count(old) == 1
src = src.replace(old, new, 1)

cfg = '            **({} if NORMREL_EXP is None else {"normrel_exp": NORMREL_EXP}),'
assert cfg in src
src = src.replace(
    cfg,
    cfg + '\n            **({} if NORMREL_EXP is None else {"normrel_mode": NORMREL_MODE}),',
    1,
)

shutil.copyfile(B, "~/Code/virtual-cell-challenge-2026/experiments/312-normrel-arm/work/builder.PRE_CAP.py")
B.write_text(src)
print("patched: VCC_NORMREL_MODE={plain|cap}; sha256", hashlib.sha256(src.encode()).hexdigest()[:16])
