#!/usr/bin/env python3
"""319 -- Astra's variant: VCC_NORMREL_MODE = plain | cap.

  plain (default, today's shipped arm):  w_j(t) proportional to r^exp
  cap   (the new hypothesis):            w_j(t) proportional to max(1, r)^exp

with r = ||row_j(t)|| / that donor's own panel median. With exp < 0, `cap` SUPPRESSES donor rows that
are oversized relative to their own panel and leaves everything at or below the median untouched,
whereas `plain` additionally UP-WEIGHTS small rows. The two hypotheses are different mechanisms and
today's board point cannot separate them: 143 is the precedent for `pds` rising while signed credit falls.

DO NOT APPLY while a build is running -- the builder is read by the live build process.

Provenance: this is `experiments/312-normrel-arm/apply_cap_patch.py` verbatim, with TWO changes,
both required by 319's mandate and neither touching the semantics of the arm:
  (1) the pre-patch byte copy is written into `experiments/319-cap-variant/work/builder.PRE319.py`
      instead of 312's work directory (319 writes only inside its own experiment directory);
  (2) step 4 registers `VCC_NORMREL_MODE` in `def_pipeline.py`'s `ARM_BOARD_DELTA`, in exactly the
      shape 312 used for `VCC_NORMREL_EXP` and 316 for `VCC_NORMREL_PM` (the mandate asks for it and
      312's staged script did not do it).
"""
import hashlib, pathlib, shutil, sys

ROOT = pathlib.Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/319-cap-variant"
B = ROOT / "scripts/make_submission_pm_cd4.py"
D = ROOT / "scripts/def_pipeline.py"

src = B.read_text()
before_sha = hashlib.sha256(src.encode()).hexdigest()
if "VCC_NORMREL_MODE" in src:
    print("already patched")
else:
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

    shutil.copyfile(B, HERE / "work/builder.PRE319.py")
    B.write_text(src)
    print("patched: VCC_NORMREL_MODE={plain|cap}; sha256", hashlib.sha256(src.encode()).hexdigest()[:16])

# ------------------------------------------------------------------ register the arm in the driver
dsrc = D.read_text()
if "VCC_NORMREL_MODE" in dsrc:
    print("def_pipeline.py: already registered")
else:
    danchor = '                   "VCC_CD4_DEFLATE_BASIS": "+0.0042634 (`cd4svd5`) carried to a new panel as the PINNED "'
    assert danchor in dsrc, "ARM_BOARD_DELTA anchor not found"
    dnew = ('''                   "VCC_NORMREL_MODE": "`cap` = the SUPPRESSION-ONLY half of VCC_NORMREL_EXP (319): the "
                                       "ratio is clamped at 1 before the exponent, so a donor row at or "
                                       "below its own panel median keeps weight 1 and only OVERSIZED rows "
                                       "are down-weighted. `plain` (the default) is the shipped two-sided "
                                       "rule, +0.0016245 at exp -1.5 (Q519I2pLWRCW7PAEx5ul) -- `cap` is "
                                       "UNTESTED on the board",
''') + danchor
    shutil.copyfile(D, HERE / "work/def_pipeline.PRE319.py")
    dsrc = dsrc.replace(danchor, dnew, 1)
    D.write_text(dsrc)
    print("def_pipeline.py: VCC_NORMREL_MODE registered as a known arm")

print("builder before sha256:", before_sha)
print("builder after  sha256:", hashlib.sha256(B.read_text().encode()).hexdigest())
print("def_pipeline after sha256:", hashlib.sha256(D.read_text().encode()).hexdigest())
