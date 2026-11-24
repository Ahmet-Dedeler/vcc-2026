"""534 Part C -- the preregistered designs (PREREGISTRATION-C section 2) and steps.

Every B design is the October bulk configuration: composition h2c1 (H1 x2, CD4 x1, D/E/F rules; abc rules on `abc`),
the fold's own LODO key V1 at clip (0.5, 3) re-fitted ON THAT DESIGN'S BASE (as 418 priced its normrel / partial arms),
q 1.75, ballast-realised at KMAX 16, levelled at its own ntclevel factor. Only the base toggles differ:

  B0  normrel cap, -1.5   partial rows g 0.5    production (known answer)
  B1  normrel plain, -1.5 partial rows g 0.5    no cap                 -> the cap step          = B0 - B1  (319)
  B2  normrel exponent 0  partial rows g 0.5    equal donor weights    -> normrel with its cap  = B0 - B2  (312 + 319)
  B3  normrel cap, -1.5   partial rows g 0      no partial rows        -> the partial-rows step = B0 - B3  (319p + 319p2)

B2 uses the production mode with the exponent set to 0 (w = max(1, x)^0 = 1: equal per-target donor weights); its
known answer is bit-identity with lib418.pre_normrel's uniform mean (mode None, the builder with VCC_NORMREL_EXP unset).
The fold's calibration constants are the production base's and are shared by all four designs: the board-emulated pds
beta (414's pattern at the fold's PP|K3 energy), the pricer's X anchors h1c1|K3 / h1c1|none (410b's bulk), and ret_f.
"""
from __future__ import annotations

PROD = dict(NORMREL_MODE="cap", NORMREL_EXP=-1.5, PARTIAL_G=0.5)
BASES = {
    "B0": dict(PROD),
    "B1": dict(NORMREL_MODE="plain", NORMREL_EXP=-1.5, PARTIAL_G=0.5),
    "B2": dict(NORMREL_MODE="cap", NORMREL_EXP=0.0, PARTIAL_G=0.5),
    "B3": dict(NORMREL_MODE="cap", NORMREL_EXP=-1.5, PARTIAL_G=0.0),
}
COMP, KEY, Q = "h2c1", "V1_5_3", 1.75              # the October bulk configuration (the design name h2c1|V1_5_3+sh1.75)
NEAREST_438 = "h2c1|V1_5_3+sh1.5"                  # the nearest design 438 stored on every panel (section 3's fallback KA)
ANCHORS = ("h1c1|none", "h1c1|K3")                 # the pricer's X anchors (production base)
# board reads (FINDINGS section 0's chain: cd4svd5 0.2160944 -> normrel -1.5 0.2177189 -> 319 cap 0.2186597 ->
# 319p 0.2198190 -> 319p2 0.2201159)
BOARD = {"normrel_plain_312": 0.2177189 - 0.2160944, "cap_319": 0.2186597 - 0.2177189,
         "partial_319p": 0.2198190 - 0.2186597, "partial_319p2": 0.2201159 - 0.2198190}
STEPS = (  # (id, arm, ref, what, board read of the same step on A/B/C, its source)
    ("normrel", "B0", "B2", "normrel -1.5 with its cap (per-target donor weights vs equal weights)",
     BOARD["normrel_plain_312"] + BOARD["cap_319"], "312 +0.0016245 (off -> plain -1.5), then 319 +0.0009408 (plain -> cap)"),
    ("cap", "B0", "B1", "the cap form of normrel (cap vs plain at -1.5)", BOARD["cap_319"], "319 +0.0009408"),
    ("partial", "B0", "B3", "the CD4 partial rows at g 0.5 (vs none)", BOARD["partial_319p"] + BOARD["partial_319p2"],
     "319p +0.0011593 (g 0.34), then 319p2 +0.0002969 (g 0.34 -> 0.5)"),
)
PANELS = ("abc", "real", "efp", "efr")
OCT = ("real", "efp", "efr")
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
FOLD_ROUTE = ("efr",)                               # no 3-donor operator: 438's fold-route energies and levels
# 418's stored toggled compositions on A/B/C's folds (rules abc, h1c1, 410 key clip (0.3, 3), unsharpened, design rows)
KA418 = {"B1": "f2:plain-1.5|K3", "B2": "e:exp0|K3", "B3": "d:g0|K3"}
