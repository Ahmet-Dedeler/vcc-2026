"""517 §2: CLEAN4-hook variants for the null/harm branch (analytic, c4 base, 3 seeds): 496's tau .05 (KA -0.00146), the four
banks' unanimous genes (conf4: every bank's P_down on one side of 0.5), unanimous AND tau .05, and tau .02 / .03."""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import price517 as P
RC = P.RC
zr = np.load(RC.REF, allow_pickle=False); genes_full = zr["genes"]; keep = np.asarray(zr["keep"], bool); T = zr["targets"].astype(str)
cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
names = list(cp["bank_names"].astype(str)); PB = cp["P_banks"].astype(float)[[names.index(b) for b in RC.PRIOR_BANKS]]
posc = {g: j for j, g in enumerate(cp["genes"].astype(str))}
colc = np.array([posc.get(str(g), -1) for g in np.asarray(genes_full)[keep]])
Pm = np.where(colc >= 0, PB.mean(0)[np.maximum(colc, 0)], np.nan)
unan = np.where(colc >= 0, np.abs(np.sign(PB - 0.5).sum(0))[np.maximum(colc, 0)] == 4, False)
sg = np.where(Pm > 0.5, -1.0, 1.0)
V = {"tau.05 (496)": np.where(np.isfinite(Pm) & (np.abs(Pm - .5) >= .05), sg, 0.0),
     "unanimous": np.where(unan, sg, 0.0), "unanimous & tau.05": np.where(unan & (np.abs(Pm - .5) >= .05), sg, 0.0),
     "tau.03": np.where(np.isfinite(Pm) & (np.abs(Pm - .5) >= .03), sg, 0.0), "tau.02": np.where(np.isfinite(Pm) & (np.abs(Pm - .5) >= .02), sg, 0.0)}
c4 = RC.prior_signs_full(genes_full)[keep]
res = {k: [] for k in V}
for seed in (0, 1, 2):
    ctx = P.Ctx("c4", seed, lambda i, t: c4); rng = np.random.default_rng(5180 + seed)
    for k, hs in V.items():
        lfc, rec = ctx.hook(lambda i, t, hs=hs: hs); sc, _ = ctx.score(lfc, rng, n_boot=500)
        res[k].append((sc["nmae"], sc["reach"], sc["overall_equiv"], sc["ne_nmae"], ctx.gate_accuracy(lambda i, t, hs=hs: hs, "comp")["acc"], int((hs != 0).sum())))
    del ctx
for k, v in res.items():
    a = np.array(v)
    print(f"{k:20s} genes {int(a[0, 5])}: nmae {a[:, 0].mean():+.5f} reach {a[:, 1].mean():+.5f} ov {a[:, 2].mean():+.5f} ne_nmae {a[:, 3].mean():+.5f} acc {a[0, 4]:.4f}")
