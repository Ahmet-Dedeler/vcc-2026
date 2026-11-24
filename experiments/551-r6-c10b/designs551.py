"""551 (C10b) -- the designs R6's last access can price, in lib438's naming, and the steps between them.

C10a (534 Part B) priced the October bulk milestones M0 -> M6 on R6. After C10a, the only pending October BULK switches are
`lodo_hi_pct` (M5 = cap 3 vs M6 = the percentile cap; read 10-03) and `sharpen_pt` (per-target sharpening S_P2 or P2g on
top of the base; read 10-02). The October path has no union-source knob (471c stays A/B/C-only unless a pass-through is
written). So the final October bulk configuration F is one of six:

  base   M5  h2c1|V1_5_3+sh1.75      LODO key clip (0.5, 3), q 1.75           (lodo_hi_pct off)
         M6  h2c1|V1pct+sh1.75       LODO key clip (0.5, P 96.647-th pct)     (lodo_hi_pct on)
  per-target q (sharpen_pt), on either base, q = (2.5, 2.0, 1.75) by tercile (LOW, MID, HIGH):
         none   q 1.75 for every target
         S_P2   tercile of response strength: the unkeyed h2c1 row's non-panel norm (493h / 546o)
         P2g    carried-only tercile of the base's key weights; non-carried MID; every floor target LOW (493h
                `P2g_carried` == 546o `p96.6 + P2g`)

Every C10a design is priced again inside the C10b access and must reproduce C10a's stored R6 values bit for bit (KA-REPRO),
so the re-fitted operators and the rebuilt fold caches are proven identical before any new number is read.

Controls: each per-target design with its q vector permuted over targets (np.random.default_rng(493), 493's CTRL seed),
priced at its OWN level factors (a q permutation does not keep the render energy, unlike 534's key permutation).

KA-only designs (efp only, never priced on R6): 546o's rows use the ROUNDED percentile 96.647 (543), not P_PROD, so
`V1p546` reproduces 546o's `p96.6`, `p96.6 + P2g`, `p96.6 + S_P2` exactly; the V1_5_3-based per-target designs reproduce
493h's `P2g_carried` and `S_P2`.
"""
from __future__ import annotations

P_PROD = 96.64697522769474          # 497's unrounded P (== cap 6 on the A/B/C artifact); 522's `pctPROD`; 534's M6
P_546 = 96.647                      # 543 / 546o's rounded P (KA rows only)
CTRL_SEED = 420                     # 534's C-b (438's CTRL construction on the key weights)
PT_CTRL_SEED = 493                  # 493's CTRL seed (the per-target q permuted over targets)

# (name, composition, key, q (scalar or None), pt ("S_P2" | "P2g" | None), perm (bool), realised at KMAX 16)
D534 = (
    ("h1c1|none", "h1c1", None, None, None, False, False),
    ("h1c1|K3", "h1c1", "K3", None, None, False, True),
    ("h1c1|K3+sh1.25", "h1c1", "K3", 1.25, None, False, True),
    ("h2c1|K3+sh1.25", "h2c1", "K3", 1.25, None, False, True),
    ("h2c1|K3+sh1.5", "h2c1", "K3", 1.5, None, False, True),
    ("h2c1|K3+sh1.75", "h2c1", "K3", 1.75, None, False, True),
    ("h2c1|V1_5_3+sh1.75", "h2c1", "V1_5_3", 1.75, None, False, True),
    ("h2c1|V1pct+sh1.75", "h2c1", "V1pct", 1.75, None, False, True),
    ("h2c1.5|K3+sh1.75", "h2c1.5", "K3", 1.75, None, False, True),
    ("CTRL:h2c1|V1_5_3shuf+sh1.75", "h2c1", "V1_5_3shuf", 1.75, None, False, True),
)
NEW = (
    ("h2c1|V1_5_3+ptS_P2", "h2c1", "V1_5_3", None, "S_P2", False, True),
    ("h2c1|V1_5_3+ptP2g", "h2c1", "V1_5_3", None, "P2g", False, True),
    ("h2c1|V1pct+ptS_P2", "h2c1", "V1pct", None, "S_P2", False, True),
    ("h2c1|V1pct+ptP2g", "h2c1", "V1pct", None, "P2g", False, True),
    ("CTRL:h2c1|V1_5_3+ptS_P2perm", "h2c1", "V1_5_3", None, "S_P2", True, True),
    ("CTRL:h2c1|V1_5_3+ptP2gperm", "h2c1", "V1_5_3", None, "P2g", True, True),
    ("CTRL:h2c1|V1pct+ptS_P2perm", "h2c1", "V1pct", None, "S_P2", True, True),
    ("CTRL:h2c1|V1pct+ptP2gperm", "h2c1", "V1pct", None, "P2g", True, True),
)
KA_ONLY = (                          # efp only: the rounded-P rows of 546o
    ("h2c1|V1p546+sh1.75", "h2c1", "V1p546", 1.75, None, False, True),
    ("h2c1|V1p546+ptP2g", "h2c1", "V1p546", None, "P2g", False, True),
    ("h2c1|V1p546+ptS_P2", "h2c1", "V1p546", None, "S_P2", False, True),
)
DESIGNS = D534 + NEW
NAMES = tuple(d[0] for d in DESIGNS)
BY_NAME = {d[0]: d for d in DESIGNS + KA_ONLY}
CTRL_ENERGY = {"CTRL:h2c1|V1_5_3shuf+sh1.75": "h2c1|V1_5_3+sh1.75"}   # 534's convention for the key permutation only

M = {"M0": "h1c1|K3", "M1": "h1c1|K3+sh1.25", "M2": "h2c1|K3+sh1.25", "M3": "h2c1|K3+sh1.5", "M4": "h2c1|K3+sh1.75",
     "M5": "h2c1|V1_5_3+sh1.75", "M6": "h2c1|V1pct+sh1.75", "C-a": "h2c1.5|K3+sh1.75",
     "C-b": "CTRL:h2c1|V1_5_3shuf+sh1.75",
     "M5+S_P2": "h2c1|V1_5_3+ptS_P2", "M5+P2g": "h2c1|V1_5_3+ptP2g",
     "M6+S_P2": "h2c1|V1pct+ptS_P2", "M6+P2g": "h2c1|V1pct+ptP2g",
     "C-S@M5": "CTRL:h2c1|V1_5_3+ptS_P2perm", "C-P@M5": "CTRL:h2c1|V1_5_3+ptP2gperm",
     "C-S@M6": "CTRL:h2c1|V1pct+ptS_P2perm", "C-P@M6": "CTRL:h2c1|V1pct+ptP2gperm"}

# C10a's steps (reproduced inside the access: KA-REPRO) -- (step id, arm, ref, what, class for the calibrated projection)
STEPS_C10A = (
    ("M0->M1", "M1", "M0", "sharpening q 1.25 (417a)", "shape"),
    ("M1->M2", "M2", "M1", "H1 x2 (417e's H1 half)", "composition"),
    ("M2->M3", "M3", "M2", "q 1.25 -> 1.5 (433k)", "shape"),
    ("M3->M4", "M4", "M3", "q 1.5 -> 1.75 (446b)", "shape"),
    ("M4->M5", "M5", "M4", "LODO key clip (0.5, 3) (446bj2)", "key"),
    ("M5->M6", "M6", "M5", "percentile cap P 96.647 (495/522)", "key"),
    ("C-a", "C-a", "M4", "CTRL: CD4 x1.5 on M4", "composition"),
    ("C-b", "C-b", "M5", "CTRL: LODO key permuted over targets on M5", "key"),
)
# C10b's candidate steps (which of them are the claim is fixed by the configuration recorded before the access)
STEPS_C10B = (
    ("S_P2@M6", "M6+S_P2", "M6", "per-target sharpening S_P2 on the percentile-cap base", "shape"),
    ("P2g@M6", "M6+P2g", "M6", "per-target sharpening P2g on the percentile-cap base", "shape"),
    ("S_P2@M5", "M5+S_P2", "M5", "per-target sharpening S_P2 on the cap-3 base", "shape"),
    ("P2g@M5", "M5+P2g", "M5", "per-target sharpening P2g on the cap-3 base", "shape"),
    ("C-S@M6", "C-S@M6", "M6", "CTRL: S_P2's q permuted over targets, on M6", "shape"),
    ("C-P@M6", "C-P@M6", "M6", "CTRL: P2g's q permuted over targets, on M6", "shape"),
    ("C-S@M5", "C-S@M5", "M5", "CTRL: S_P2's q permuted over targets, on M5", "shape"),
    ("C-P@M5", "C-P@M5", "M5", "CTRL: P2g's q permuted over targets, on M5", "shape"),
)
STEPS = STEPS_C10A + STEPS_C10B
# the final configuration F -> (base milestone, per-target step or None)
CONFIGS = {
    "M5": ("M5", None), "M6": ("M6", None),
    "M5+S_P2": ("M5", "S_P2@M5"), "M5+P2g": ("M5", "P2g@M5"),
    "M6+S_P2": ("M6", "S_P2@M6"), "M6+P2g": ("M6", "P2g@M6"),
}
K_CAL = {"shape": (0.65, 0.42), "composition": (1.3, 1.0), "key": (1.31, 1.0)}   # 534-B section 6


def energy_name(n: str) -> str:
    return CTRL_ENERGY.get(n, n)


def designs_for(config: str | None, ka: bool = False):
    """The designs a run prices: every 534 design (KA-REPRO and the cumulative) + F's own per-target design and its
    control. With --ka (efp) every candidate and the KA-only rows, so every code path is checked before the access."""
    if ka:
        return D534 + NEW + KA_ONLY
    assert config in CONFIGS, f"--config must be one of {sorted(CONFIGS)}"
    base, sid = CONFIGS[config]
    if sid is None:
        return D534
    step = {s[0]: s for s in STEPS_C10B}[sid]
    ctrl_id = ("C-S@" if sid.startswith("S_P2") else "C-P@") + base
    ctrl = {s[0]: s for s in STEPS_C10B}[ctrl_id]
    want = {M[step[1]], M[ctrl[1]]}
    return D534 + tuple(d for d in NEW if d[0] in want)
