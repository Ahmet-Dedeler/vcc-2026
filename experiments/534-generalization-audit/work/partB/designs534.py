"""534 Part B -- the preregistered designs (PREREGISTRATION-B section 4) in lib438's naming, and what each step prices.

Milestones are cumulative; each is priced against the one before it (438's LEVELLED route, k = 1):
  M0  h1c1|K3               the 410b-era bulk: per-target 410 key clip (0.3, 3), no sharpening, no stack
  M1  h1c1|K3+sh1.25        + sharpening q 1.25 (417a)
  M2  h2c1|K3+sh1.25        + H1 x2 (417e's H1 half)
  M3  h2c1|K3+sh1.5         q 1.5 (433k)
  M4  h2c1|K3+sh1.75        q 1.75 (446b)
  M5  h2c1|V1_5_3+sh1.75    LODO key, clip (0.5, 3) (446bj2)
  M6  h2c1|V1pct+sh1.75     M5 with the percentile cap P = 96.64697522769474 (522's `pctPROD`): hi = the P-th numpy-linear
                            percentile of key/median over the carried targets, per fold on the fold's own LODO key (522's
                            emulation) and on the full panel's 3-donor LODO key (def_lodo_key.py --hi-pct) for energies
Controls:
  C-a h2c1.5|K3+sh1.75 vs M4                  CD4 x1.5 (the half 438 found A/B/C-only). Predicted <= 0.
  C-b CTRL:h2c1|V1_5_3shuf+sh1.75 vs M5        the LODO key permuted over targets, 438's CTRL construction verbatim
                            (np.random.default_rng(420).permutation of the fold's weight vector, applied after clipping);
                            priced with M5's energies and M5's level factors (420/438's convention: a permutation keeps the
                            key's energy, `en("CTRL:...") = en(<unpermuted design>)`). Predicted < 0.
"h1c1|none" is not a milestone: the pricer's s0 reference (X's relative-cosine anchor), exactly as in 438.
"""
from __future__ import annotations

P_PROD = 96.64697522769474          # 497's unrounded P (== cap 6 on the A/B/C artifact); 522's `pctPROD`
CTRL_SEED = 420                     # 438's CTRL construction (comp438.py: np.random.default_rng(420))

# (name, composition, key, q, realised at KMAX 16)
DESIGNS = (
    ("h1c1|none", "h1c1", None, None, False),
    ("h1c1|K3", "h1c1", "K3", None, True),
    ("h1c1|K3+sh1.25", "h1c1", "K3", 1.25, True),
    ("h2c1|K3+sh1.25", "h2c1", "K3", 1.25, True),
    ("h2c1|K3+sh1.5", "h2c1", "K3", 1.5, True),
    ("h2c1|K3+sh1.75", "h2c1", "K3", 1.75, True),
    ("h2c1|V1_5_3+sh1.75", "h2c1", "V1_5_3", 1.75, True),
    ("h2c1|V1pct+sh1.75", "h2c1", "V1pct", 1.75, True),
    ("h2c1.5|K3+sh1.75", "h2c1.5", "K3", 1.75, True),
    ("CTRL:h2c1|V1_5_3shuf+sh1.75", "h2c1", "V1_5_3shuf", 1.75, True),
)
NAMES = tuple(d[0] for d in DESIGNS)
M = {"M0": "h1c1|K3", "M1": "h1c1|K3+sh1.25", "M2": "h2c1|K3+sh1.25", "M3": "h2c1|K3+sh1.5", "M4": "h2c1|K3+sh1.75",
     "M5": "h2c1|V1_5_3+sh1.75", "M6": "h2c1|V1pct+sh1.75", "C-a": "h2c1.5|K3+sh1.75",
     "C-b": "CTRL:h2c1|V1_5_3shuf+sh1.75"}
CTRL_ENERGY = {M["C-b"]: M["M5"]}   # whose energies / level factors a control uses
# (step id, arm, ref, what it adds, class for the calibrated projection)
STEPS = (
    ("M0->M1", "M1", "M0", "sharpening q 1.25 (417a)", "shape"),
    ("M1->M2", "M2", "M1", "H1 x2 (417e's H1 half)", "composition"),
    ("M2->M3", "M3", "M2", "q 1.25 -> 1.5 (433k)", "shape"),
    ("M3->M4", "M4", "M3", "q 1.5 -> 1.75 (446b)", "shape"),
    ("M4->M5", "M5", "M4", "LODO key clip (0.5, 3) (446bj2)", "key"),
    ("M5->M6", "M6", "M5", "percentile cap P 96.647 (495/522)", "key"),
    ("C-a", "C-a", "M4", "CTRL: CD4 x1.5 on M4", "composition"),
    ("C-b", "C-b", "M5", "CTRL: LODO key permuted over targets on M5", "key"),
)
MAIN_STEPS = STEPS[:6]
OCTOBER_SWITCHES = ("M3->M4", "M4->M5", "M5->M6")      # section 7 rule 3
# the calibrated board-unit transfer (PREREGISTRATION-B section 6): (k for X, k for pds)
K_CAL = {"shape": (0.65, 0.42), "composition": (1.3, 1.0), "key": (1.31, 1.0)}


def energy_name(n: str) -> str:
    return CTRL_ENERGY.get(n, n)
