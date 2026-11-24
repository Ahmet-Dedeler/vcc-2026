"""462g — the partial-bin follow-up, priced on H1 before any read: a CLEAN4 + Arc-full blend of per-gene programmes,
z = w * z(P_down_arc - 0.5) + (1 - w) * z(P_down_clean4 - 0.5), sign(-z). Arc alone where CLEAN4 has no value. Readings:
the fid-like row (all block genes) and the gate row (58 targets), as 462e/f. Same line and lab as the Arc source, so
this is an upper reference. usage: uv run python .../blend462.py"""
import glob, json, csv
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).resolve().parent
z = np.load(glob.glob(str(ROOT / "experiments/323-*/work/reference.npz"))[0], allow_pickle=True)
genes, keep, targets, lfc, q, bctrl = z["genes"].astype(str), z["keep"], list(z["targets"].astype(str)), z["lfc"], z["q"], z["b_ctrl"]
kg = genes[keep]; e = np.expm1(bctrl[keep]); cpm = 1e6 * e / e.sum(); tset = set(targets)
order = np.argsort(-cpm); block = np.zeros(kg.size, bool); block[[i for i in order if kg[i] not in tset][:2200]] = True
d = ROOT / "data/donors/h1_2025_full"; delta = np.load(d / "delta.npy", mmap_mode="r"); perts = [l.strip() for l in open(d / "perts.txt")]
bg = [l.strip() for l in open(d / "genes.txt")]; rows = [i for i, p in enumerate(perts) if p not in tset]
M = np.asarray(delta[rows], dtype=np.float64); nz = (M != 0).sum(0); dn = (M < 0).sum(0)
with np.errstate(invalid="ignore", divide="ignore"): pa = np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)
pos = {g: i for i, g in enumerate(bg)}; col = np.array([pos.get(g, -1) for g in kg]); arc = np.where(col >= 0, pa[np.maximum(col, 0)], np.nan)
c4 = np.load(ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=True)
offi = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]; op = {g: i for i, g in enumerate(offi)}
c4p = np.asarray(c4["P_down"], float); cl = np.array([c4p[op[g]] if g in op else np.nan for g in kg])
zs = lambda v: (v - 0.5) / np.nanstd(v[block] - 0.5)
own = np.array([[g == t for g in kg] for t in targets])
gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & ~own; use = np.flatnonzero(gate.sum(1) >= 10)
res = {}
for w in (0.0, 0.25, 0.5, 0.75, 1.0):
    zb = np.where(np.isfinite(cl), w * zs(arc) + (1 - w) * zs(cl), zs(arc)) if w > 0 else zs(cl)
    sign = np.where(np.isfinite(zb), np.where(zb > 0, -1.0, 1.0), np.nan)
    h = n = gh = gn = 0
    for t in range(len(targets)):
        g = block & ~own[t] & np.isfinite(lfc[t]) & (lfc[t] != 0) & np.isfinite(sign)
        h += (np.sign(lfc[t][g]) == sign[g]).sum(); n += g.sum()
        if t in use:
            gg = gate[t] & np.isfinite(sign); gh += (np.sign(lfc[t][gg]) == sign[gg]).sum(); gn += gg.sum()
    res[f"w_arc={w}"] = dict(fidlike=round(h / n, 4), fidlike_pairs=int(n), gate=round(gh / gn, 4), gate_pairs=int(gn))
print(json.dumps(res, indent=1)); (HERE / "work/blend462.json").write_text(json.dumps(res, indent=1))
