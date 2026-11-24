"""462k — Arc's H1 programme from its NON-essential perturbations only? 36 % of the 275 off-panel H1 perturbations are in the
Replogle essential sets; the A/B/C panel has none. On h1-test targets split the same way (essential vs not), compare the
all-row programme (458d) with non-essential-only and essential-only rows (h1-test targets always excluded). Rows: fid-like
(all block genes) and gate (targets with >= 10 gate genes). usage: uv run python .../arcness462.py"""
import glob, json
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).resolve().parent
z = np.load(glob.glob(str(ROOT / "experiments/323-*/work/reference.npz"))[0], allow_pickle=True)
genes, keep, targets, lfc, q, bctrl = z["genes"].astype(str), z["keep"], list(z["targets"].astype(str)), z["lfc"], z["q"], z["b_ctrl"]
kg = genes[keep]; e = np.expm1(bctrl[keep]); cpm = 1e6 * e / e.sum(); tset = set(targets)
order = np.argsort(-cpm); block = np.zeros(kg.size, bool); block[[i for i in order if kg[i] not in tset][:2200]] = True
own = np.array([[g == t for g in kg] for t in targets])
gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & ~own
ess = set(l.strip() for l in open(ROOT / "data/donors/k562_essential/perts.txt")) | set(l.strip() for l in open(ROOT / "data/donors/rpe1/perts.txt"))
d = ROOT / "data/donors/h1_2025_full"; delta = np.load(d / "delta.npy"); perts = [l.strip() for l in open(d / "perts.txt")]
bg = [l.strip() for l in open(d / "genes.txt")]; pos = {g: i for i, g in enumerate(bg)}; col = np.array([pos.get(g, -1) for g in kg])
def prog(f):
    rows = [i for i, p in enumerate(perts) if p not in tset and f(p)]
    M = delta[rows].astype(np.float64); nz = (M != 0).sum(0); dn = (M < 0).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"): pdn = np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)
    pk = np.where(col >= 0, pdn[np.maximum(col, 0)], np.nan)
    return np.where(np.isfinite(pk), np.where(pk > 0.5, -1.0, 1.0), np.nan), len(rows)
def score(sign, tsel):
    h = n = gh = gn = 0
    for t in tsel:
        g = block & ~own[t] & np.isfinite(lfc[t]) & (lfc[t] != 0) & np.isfinite(sign)
        h += (np.sign(lfc[t][g]) == sign[g]).sum(); n += g.sum()
        if gate[t].sum() >= 10:
            gg = gate[t] & np.isfinite(sign); gh += (np.sign(lfc[t][gg]) == sign[gg]).sum(); gn += gg.sum()
    return round(h / n, 4), round(gh / max(gn, 1), 4), int(gn)
T_ne = [t for t in range(len(targets)) if targets[t] not in ess]; T_e = [t for t in range(len(targets)) if targets[t] in ess]
res = {"h1test_targets_nonessential": len(T_ne), "h1test_targets_essential": len(T_e)}
for tag, f in (("all_rows", lambda p: True), ("nonessential_rows", lambda p: p not in ess), ("essential_rows", lambda p: p in ess)):
    s, n = prog(f)
    res[tag] = {"rows": n, "on_nonessential_targets": score(s, T_ne), "on_essential_targets": score(s, T_e), "on_all": score(s, range(len(targets)))}
    print(tag, json.dumps(res[tag]))
print("targets: non-essential", len(T_ne), "essential", len(T_e))
(HERE / "work/arcness462.json").write_text(json.dumps(res, indent=1))
