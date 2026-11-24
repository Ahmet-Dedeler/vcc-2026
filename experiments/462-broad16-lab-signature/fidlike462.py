"""462f — the fid-like reading (441 §3's last column): sign agreement between a per-gene programme and the reference lfc on
ALL of each h1-test target's block genes (2,200 top control-CPM kept genes, own gene out, lfc != 0), pooled over the 85
targets. fid is precision over a large call set, so this is its best single-programme proxy. Programmes: 458's bank,
458d's full-cell programme (both with h1-test's targets excluded), CLEAN4, and the LOTO oracle (other h1-test targets).
Same line and lab: an upper reference for the cross-line board read. usage: uv run python .../fidlike462.py"""
import glob, json, csv
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).resolve().parent
z = np.load(glob.glob(str(ROOT / "experiments/323-*/work/reference.npz"))[0], allow_pickle=True)
genes, keep, targets, lfc, bctrl = z["genes"].astype(str), z["keep"], list(z["targets"].astype(str)), z["lfc"], z["b_ctrl"]
kg = genes[keep]; e = np.expm1(bctrl[keep]); cpm = 1e6 * e / e.sum()
tset = set(targets)
order = np.argsort(-cpm); block = np.zeros(kg.size, bool); block[[i for i in order if kg[i] not in tset][:2200]] = True
def prog(name):
    d = ROOT / "data/donors" / name; delta = np.load(d / "delta.npy", mmap_mode="r"); perts = [l.strip() for l in open(d / "perts.txt")]
    bg = [l.strip() for l in open(d / "genes.txt")]; rows = [i for i, p in enumerate(perts) if p not in tset]
    M = np.asarray(delta[rows], dtype=np.float64); nz = (M != 0).sum(0); dn = (M < 0).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"): pdn = np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)
    pos = {g: i for i, g in enumerate(bg)}; col = np.array([pos.get(g, -1) for g in kg])
    pk = np.where(col >= 0, pdn[np.maximum(col, 0)], np.nan)
    return np.where(np.isfinite(pk), np.where(pk > 0.5, -1.0, 1.0), np.nan)
c4 = np.load(ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=True)
offi = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]; pos = {g: i for i, g in enumerate(offi)}
c4p = np.asarray(c4["P_down"], float); c4s = np.asarray(c4["psign"], float)
clean4 = np.array([c4s[pos[g]] if g in pos and np.isfinite(c4p[pos[g]]) else np.nan for g in kg])
P = {"bank_h1_2025": prog("h1_2025"), "full_h1_2025_full": prog("h1_2025_full"), "clean4": clean4}
tot = lfc.sum(0)
res = {}
for name in list(P) + ["oracle_loto"]:
    hit = n = 0
    for t in range(len(targets)):
        own = kg == targets[t]
        pred = P[name] if name != "oracle_loto" else np.sign((tot - lfc[t]) / (len(targets) - 1))
        g = block & ~own & np.isfinite(lfc[t]) & (lfc[t] != 0) & np.isfinite(pred) & (pred != 0)
        hit += (np.sign(lfc[t][g]) == pred[g]).sum(); n += g.sum()
    res[name] = dict(acc=round(hit / n, 4), pairs=int(n))
print(json.dumps(res, indent=1)); (HERE / "work/fidlike462.json").write_text(json.dumps(res, indent=1))
