"""462j — is CLEAN4 (3 of its 4 banks are essential screens) the essential-knockdown programme, and would a programme from
NON-essential perturbations read differently? K562-GWPS split by the K562-essential target list (the Replogle essential
set), and Orion HCT116 likewise, on H1's fid-like row (all 2,200 block genes, 85 targets, 462f) and gate row (58
targets). A cross-line, cross-lab proxy only: H1's truth is an ESC. usage: uv run python .../essnoness462.py"""
import glob, json, csv, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts")); import def_sameline as S  # noqa: E402
z = np.load(glob.glob(str(ROOT / "experiments/323-*/work/reference.npz"))[0], allow_pickle=True)
genes, keep, targets, lfc, q, bctrl = z["genes"].astype(str), z["keep"], list(z["targets"].astype(str)), z["lfc"], z["q"], z["b_ctrl"]
kg = genes[keep]; e = np.expm1(bctrl[keep]); cpm = 1e6 * e / e.sum(); tset = set(targets)
order = np.argsort(-cpm); block = np.zeros(kg.size, bool); block[[i for i in order if kg[i] not in tset][:2200]] = True
own = np.array([[g == t for g in kg] for t in targets])
gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & ~own; use = np.flatnonzero(gate.sum(1) >= 10)
ess = set(l.strip() for l in open(ROOT / "data/donors/k562_essential/perts.txt"))
def prog(bank, keep_row):
    delta, perts, bg, _ = S.load_bank(bank)
    rows = [i for i, p in enumerate(perts) if p not in tset and keep_row(p)]
    M = np.asarray(delta[rows], dtype=np.float64); nz = (M != 0).sum(0); dn = (M < 0).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"): pdn = np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)
    pos = {g: i for i, g in enumerate(bg)}; col = np.array([pos.get(g, -1) for g in kg])
    pk = np.where(col >= 0, pdn[np.maximum(col, 0)], np.nan)
    return np.where(np.isfinite(pk), np.where(pk > 0.5, -1.0, 1.0), np.nan), len(rows)
def score(sign):
    h = n = gh = gn = 0
    for t in range(len(targets)):
        g = block & ~own[t] & np.isfinite(lfc[t]) & (lfc[t] != 0) & np.isfinite(sign)
        h += (np.sign(lfc[t][g]) == sign[g]).sum(); n += g.sum()
        if t in use:
            gg = gate[t] & np.isfinite(sign); gh += (np.sign(lfc[t][gg]) == sign[gg]).sum(); gn += gg.sum()
    return round(h / n, 4), round(gh / gn, 4)
S.BANK_DIRS.update({"orion_hek293t": ROOT / "data/donors/orion_hek293t"})
res = {}
for bank in ("k562_gwps", "orion_hct116"):
    for tag, f in (("all", lambda p: True), ("essential", lambda p: p in ess), ("non_essential", lambda p: p not in ess)):
        sgn, n = prog(bank, f); fl, gt = score(sgn)
        res[f"{bank}:{tag}"] = dict(rows=n, fidlike=fl, gate=gt)
        print(f"{bank:13s} {tag:13s} rows {n:6d}  fid-like {fl:.4f}  gate {gt:.4f}", flush=True)
print("reference: CLEAN4 0.5244 / 0.5733; 458d 0.5545 / 0.6400")
(HERE / "work/essnoness462.json").write_text(json.dumps(res, indent=1))
