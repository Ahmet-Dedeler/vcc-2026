"""517 §3 follow-up — excess dropout predicts DE CALLS (not signs) at matched expression (dropout517). Would a CPM block
chosen by call propensity capture more of the truth's DE genes than the shipped top-2,200-by-CPM block? 2-fold by target on
h1-test (fit on one half of the targets, evaluate capture on the other), block size fixed at 2,200. Descriptive only (the
block is also where fid's signs are read; 441 PS3: re-choosing it trades fid and jac)."""
import glob, json
from pathlib import Path
import numpy as np
from sklearn.linear_model import LogisticRegression
ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).resolve().parent
import sys; sys.path.insert(0, str(HERE))
from dropout517 import excess_dropout
z = np.load(glob.glob(str(ROOT / "experiments/323-*/work/reference.npz"))[0], allow_pickle=False)
genes, keep, targets = z["genes"].astype(str), np.asarray(z["keep"], bool), list(z["targets"].astype(str))
q, lfc = np.asarray(z["q"], float), np.asarray(z["lfc"], float)
kg = genes[keep]
H = excess_dropout(ROOT / "data/harness/h1-test.h5ad")
x = H["x"][keep]; lc = np.log10(np.maximum(H["mean_cpm"][keep], 1e-3))
own = np.array([[g == t for g in kg] for t in targets])
gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & ~own
tset = set(targets); ok = np.array([g not in tset for g in kg])
rng = np.random.default_rng(517); perm = rng.permutation(len(targets)); halves = (perm[: len(perm) // 2], perm[len(perm) // 2:])
F = np.column_stack([lc, x, lc * x, lc ** 2])
res = {}
for k, (fit, ev) in enumerate((halves, halves[::-1])):
    y = gate[fit].mean(0)                                   # per-gene call rate on the fitting targets
    w = np.full(y.size, len(fit))
    Xr = np.repeat(F, 2, 0); yr = np.tile([1, 0], y.size); wr = np.column_stack([y * w, (1 - y) * w]).ravel()
    m = LogisticRegression(max_iter=2000).fit(Xr, yr, sample_weight=wr)
    prop = m.predict_proba(F)[:, 1]
    for name, score in (("cpm", lc), ("propensity", prop), ("cpm only refit", LogisticRegression(max_iter=2000).fit(np.repeat(F[:, [0, 3]], 2, 0), yr, sample_weight=wr).predict_proba(F[:, [0, 3]])[:, 1])):
        blk = np.zeros(kg.size, bool); blk[np.argsort(-np.where(ok, score, -np.inf))[:2200]] = True
        cap = gate[ev][:, blk].sum() / gate[ev].sum()
        per_t = np.mean([gate[t][blk].sum() / max(gate[t].sum(), 1) for t in ev if gate[t].sum() > 0])
        # sign-accuracy territory: the reference lfc's |value| (sign reliability proxy) on the block, all adjudicable genes
        med_abs = float(np.nanmedian(np.abs(lfc[ev][:, blk])))
        res.setdefault(name, []).append((float(cap), float(per_t), med_abs))
for name, v in res.items():
    print(f"{name:16s} capture of held-out gate genes {np.mean([a for a, _, _ in v]):.4f} (per-target mean {np.mean([b for _, b, _ in v]):.4f}); "
          f"median |ref lfc| on the block {np.mean([c for _, _, c in v]):.4f}")
(HERE / "work/callprop517.json").write_text(json.dumps(res, indent=1))
