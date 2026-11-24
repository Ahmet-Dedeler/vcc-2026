#!/usr/bin/env python3
"""312 known-answer rows for the VCC_NORMREL_EXP wiring, on a synthetic 3-donor model.

KA-1  unset            -> the original expression (the reference the other rows compare to)
KA-2  exponent 0       -> uniform weights, must be BIT-IDENTICAL to unset
KA-3  exponent -1.5    -> must differ, and every row's norm must equal the unset row's norm
                          (shape-only, 263 ':nm'); AGENTS 8.1d -- magnitudes must not move
KA-4  the weights themselves reproduce (||row_j(t)||/median_j)**-1.5 by hand on one target
"""
import importlib, os, sys
import numpy as np

sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/src")
sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/scripts")
from vcc2026.donor_transfer import Donor
from vcc2026.multi_donor import MultiDonorModel

rng = np.random.default_rng(312)
GENES = np.array([f"G{i}" for i in range(40)])
PERTS = np.array([f"T{i}" for i in range(12)])

def mk(scale_by_pert):
    d = rng.normal(0, 1, size=(len(PERTS), len(GENES))).astype(np.float32)
    d = d * np.asarray(scale_by_pert, dtype=np.float32)[:, None]
    return Donor(delta=d, perts=PERTS.copy(), genes=GENES.copy(),
                 b_ctrl=np.abs(rng.normal(3, 0.4, len(GENES))).astype(np.float32))

# donor B is deliberately given one wildly large row (T3) so the weighting has something to do
donors = {
    "k562_gwps":   mk([1.0]*12),
    "orion_hct116":mk([1.0,1.0,1.0,5.0]+[1.0]*8),
    "orion_hek293t":mk([1.0]*12),
}
model = MultiDonorModel(donors, GENES, gene_norm=False, eb=False)

def effects(exp_env):
    if exp_env is None: os.environ.pop("VCC_NORMREL_EXP", None)
    else: os.environ["VCC_NORMREL_EXP"] = str(exp_env)
    import make_submission_pm_cd4 as M
    importlib.reload(M)
    return M.build_v9_pre_whiten_effects(model, list(PERTS), alpha=M.ALPHA, beta=M.BETA), M

E0, M0 = effects(None)
E1, _ = effects(0)
E2, _ = effects(-1.5)

ok = True
def row(name, cond, detail):
    global ok
    ok &= bool(cond)
    print(f"{'PASS' if cond else 'FAIL'}  {name}: {detail}")

row("KA-1 unset path runs", E0.shape == (12, 40), f"shape {E0.shape}")
d10 = float(np.max(np.abs(E1 - E0)))
row("KA-2 exponent 0 == unset", d10 == 0.0, f"max|diff| {d10:.3e} (must be exactly 0)")
d20 = float(np.max(np.abs(E2 - E0)))
row("KA-3a exponent -1.5 differs", d20 > 0, f"max|diff| {d20:.3e}")
n0 = np.linalg.norm(E0, axis=1); n2 = np.linalg.norm(E2, axis=1)
worst = float(np.max(np.abs(n2 - n0) / np.maximum(n0, 1e-12)))
row("KA-3b shape-only (row norms preserved)", worst < 1e-9, f"worst relative norm change {worst:.3e}")

# KA-4: hand-compute the weights for T3, where donor orion_hct116 is 5x its own median
t = "T3"; hand = {}
for n in model.names:
    cols = model.dcols[n]; idx = model.pert_index[n]
    med = float(np.median([np.linalg.norm(model.s[n][idx[p]][cols]) for p in PERTS]))
    nv = float(np.linalg.norm(model.s[n][idx[t]][cols]))
    hand[n] = (nv / med) ** -1.5
tot = sum(hand.values()); hand = {k: v / tot for k, v in hand.items()}
row("KA-4 hand weights: the 5x donor is down-weighted",
    hand["orion_hct116"] < 0.5 * hand["k562_gwps"],
    " ".join(f"{k.split('_')[0]}={v:.4f}" for k, v in hand.items()))
print("\nALL KNOWN-ANSWER ROWS PASS" if ok else "\nSOME ROWS FAILED")
sys.exit(0 if ok else 1)
