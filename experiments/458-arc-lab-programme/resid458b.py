"""458b sizing (analytic, before any 458 read): Arc-H1's programme with Pan-H1's regressed out (z-scored P_down, OLS on
genes both cover). How many block-prior signs does it flip relative to 458's (Arc-H1) programme, per context's CPM block
(2,200 top-CPM genes of its own controls), and how does each agree with CLEAN4? If few signs flip, 458b cannot read
differently from 458 on the board and is not worth a slot.   usage: uv run python .../resid458b.py"""
import csv, json, sys
from pathlib import Path
import h5py, numpy as np
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "scripts"))
import def_sameline as S  # noqa: E402
S.BANK_DIRS.setdefault("h1_2025", ROOT / "data/donors/h1_2025")
official = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
panel = {r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))}
d, p, g, _ = S.load_bank("h1_2025"); arc_sign, arc_pd, _, _ = S.programme(d, p, g, official, panel)
db = ROOT / "data/external/gse295214/pseudo_bulk_perturbation_database.h5ad"
cat = lambda n: n["categories"][:].astype(str)[n["codes"][:]] if isinstance(n, h5py.Group) else n[:].astype(str)
with h5py.File(db, "r") as f:
    ds, pg, vg = cat(f["obs/dataset"]), cat(f["obs/perturbed_gene_name"]), cat(f["var/gene_name"])
    sel = np.sort(np.where((ds == "Pan2025_hESC") & ~np.isin(pg, list(panel)))[0])
    X = np.where(f["layers/measure_mask"][sel, :] > 0, f["X"][sel, :], np.nan)
with np.errstate(invalid="ignore"):
    pp = np.nanmean(np.where(np.isnan(X), np.nan, (X < 0).astype(float)), axis=0)
vi = {x: i for i, x in enumerate(vg)}; col = np.array([vi.get(x, -1) for x in official])
pan_pd = np.where(col >= 0, pp[np.maximum(col, 0)], np.nan)
both = np.isfinite(arc_pd) & np.isfinite(pan_pd)
z = lambda v: (v - np.nanmean(v[both])) / np.nanstd(v[both])
za, zp = z(arc_pd), z(pan_pd)
beta = float(np.sum(za[both] * zp[both]) / np.sum(zp[both] ** 2))
resid = np.where(both, za - beta * zp, np.nan)
# 458b's sign: down (-1) where Arc's residual says "more down than Pan predicts" AND keep Arc's level: threshold the
# residual at Arc's own z of P_down = 0.5, so a gene with no Pan signal keeps exactly 458's sign
thr = float((0.5 - np.nanmean(arc_pd[both])) / np.nanstd(arc_pd[both]))
b_sign = np.where(both, np.where(resid > thr * (1 - beta * 0), -1.0, 1.0), arc_sign)
b_sign = np.where(both, np.where(za - beta * zp > thr - beta * float(np.nanmean(zp[both])), -1.0, 1.0), arc_sign)
c4 = np.load(ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=True)
res = {"beta": beta, "r_arc_pan": float(np.corrcoef(za[both], zp[both])[0, 1]), "genes_both": int(both.sum()),
       "flip_vs_458_all": float(np.mean(b_sign[both] != arc_sign[both])), "contexts": {}}
for c in "ABC":
    b = np.load(ROOT / f"experiments/438-october-stack-repricing/work/abc_frames/bctrl-{c}.npy")
    expr = np.expm1(b); cpm = 1e6 * expr / expr.sum()
    order = np.argsort(-cpm); blk = [i for i in order if official[i] not in panel][:2200]; blk = np.asarray(blk)
    rec = {"flip_vs_458": float(np.mean(b_sign[blk] != arc_sign[blk])),
           "frac_down_458": float(np.mean(arc_sign[blk] < 0)), "frac_down_458b": float(np.mean(b_sign[blk] < 0))}
    for k in c4.files:
        v = np.asarray(c4[k])
        if v.shape == (len(official),) and v.dtype.kind in "fi" and k.lower().startswith(("psign", "sign")) :
            rec[f"agree_458_{k}"] = float(np.mean(np.sign(v[blk]) == arc_sign[blk]))
            rec[f"agree_458b_{k}"] = float(np.mean(np.sign(v[blk]) == b_sign[blk]))
    res["contexts"][c] = rec
print("clean4 keys:", c4.files)
print(json.dumps(res, indent=1))
(Path(__file__).parent / "work/resid458b.json").write_text(json.dumps(res, indent=1))
