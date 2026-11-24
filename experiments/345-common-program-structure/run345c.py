"""345 step 4b-extra + 5: target-matched donor means, common-vector cosines, oracle-best splits."""
from __future__ import annotations
import json
import numpy as np
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DON = REPO / "data/donors"
OUT = Path(__file__).resolve().parent / "work"
st = np.load(OUT / "stage345.npz")
R = st["R"]
Dfull = st["Dfull"]
gate = st["gate"].astype(bool)
elig = st["elig"]
genes = st["genes"].astype(str)
targets = st["targets"].astype(str)
basal = st["basal"].astype(float)
Cmean = st["Cmean"]
G = len(genes)
L = R

def scores_vec(mu):
    pw, pu = [], []
    s = np.sign(mu)
    for ti in elig:
        m_ = gate[ti]
        l = L[ti, m_]
        w = np.abs(l)
        ok = (s[m_] == np.sign(l))
        pw.append(float((w * ok).sum() / w.sum()))
        pu.append(float(ok.mean()))
    return float(np.mean(pw)), float(np.mean(pu)), pw

def load_bank(name):
    d = np.load(DON / name / "delta.npy", mmap_mode="r")
    bg = [x.strip() for x in open(DON / name / "genes.txt")]
    bp = [x.strip() for x in open(DON / name / "perts.txt")]
    lut = {g: i for i, g in enumerate(bg)}
    col = np.array([lut.get(g, -1) for g in genes])
    return d, col, bp

res = json.load(open(OUT / "results_part2.json"))
extra = {}

# ---- target-matched donor means (rows whose target is in the H1 85-panel) ----
tset = set(targets.tolist())
for name in ("k562_gwps", "orion_hct116", "orion_hek293t"):
    d, col, bp = load_bank(name)
    hit = col >= 0
    rows = [i for i, p in enumerate(bp) if p in tset]
    mu_all = np.zeros(G)
    n = d.shape[0]
    for a in range(0, n, 1000):
        mu_all[hit] += d[a:a + 1000][:, col[hit]].astype(np.float64).sum(axis=0)
    mu_all /= n
    mu_m = np.zeros(G)
    if rows:
        blk = d[rows][:, col[hit]].astype(np.float64)
        mu_m[hit] = blk.mean(axis=0)
    w_all, u_all, _ = scores_vec(mu_all)
    out = {"n_rows_matched": len(rows), "all_mean": {"wacc": w_all, "uacc": u_all}}
    if rows:
        w_m, u_m, _ = scores_vec(mu_m)
        out["matched_mean"] = {"wacc": w_m, "uacc": u_m}
        # combination of design with matched mean at same scale rule
        ad = np.median(np.concatenate([np.abs(Dfull[ti, gate[ti]]) for ti in elig]))
        ac = np.median(np.concatenate([np.abs(mu_m[gate[ti]]) for ti in elig]))
        sb = float(ad / ac) if ac > 0 else 0.0
        grid = {}
        for a in (0.25, 0.5, 1, 2):
            pw = []
            for ti in elig:
                m_ = gate[ti]
                l = L[ti, m_]
                w = np.abs(l)
                ok = (np.sign(Dfull[ti, m_] + a * sb * mu_m[m_]) == np.sign(l))
                pw.append(float((w * ok).sum() / w.sum()))
            grid[str(a)] = float(np.mean(pw))
        out["matched_combo_grid"] = grid
        out["s"] = sb
    print(name, "matched rows", len(rows), "all w %.4f u %.4f" % (w_all, u_all),
          ("matched w %.4f u %.4f grid %s" % (out["matched_mean"]["wacc"], out["matched_mean"]["uacc"], out.get("matched_combo_grid"))) if rows else "")
    extra[name] = out

# ---- cosines between common vectors (universe genes, pairwise-complete) ----
vecs = {"H1_LOOmean_full": R.mean(axis=0), "design_own_mean": Dfull.mean(axis=0)}
for name in ("k562_gwps", "orion_hct116", "orion_hek293t"):
    d, col, _ = load_bank(name)
    hit = col >= 0
    mu = np.zeros(G)
    n = d.shape[0]
    for a in range(0, n, 1000):
        mu[hit] += d[a:a + 1000][:, col[hit]].astype(np.float64).sum(axis=0)
    vecs[f"{name}_mean"] = mu / n
names = list(vecs)
cos = {}
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        a, b = vecs[names[i]], vecs[names[j]]
        mk = (a != 0) & (b != 0)
        c = float(a[mk] @ b[mk] / (np.linalg.norm(a[mk]) * np.linalg.norm(b[mk])))
        cos[f"{names[i]}__{names[j]}"] = {"cos": c, "n_genes": int(mk.sum())}
        print(f"{names[i]} x {names[j]}: cos {c:+.4f} on {mk.sum()} genes")
extra["common_vector_cosines"] = cos

# ---- oracle-best splits: UP/DOWN + basal quintiles (mean variant, alpha=1) ----
s_or = res["step3_oracle"]["mean"]["s"]
Csign = np.sign(Dfull + 1.0 * s_or * Cmean)
Dsign = np.sign(Dfull)
spl = {}
for dname, fn in (("ref_UP", lambda l: l > 0), ("ref_DOWN", lambda l: l < 0)):
    wc, wd = [], []
    for ti in elig:
        m_ = gate[ti]
        l = L[ti, m_]
        w = np.abs(l)
        mk = fn(l)
        wc.append(float((w[mk] * (Csign[ti, m_][mk] == np.sign(l[mk]))).sum() / w[mk].sum()))
        wd.append(float((w[mk] * (Dsign[ti, m_][mk] == np.sign(l[mk]))).sum() / w[mk].sum()))
    spl[dname] = {"combo_wacc": float(np.mean(wc)), "design_wacc": float(np.mean(wd))}
    print(dname, "combo %.4f design %.4f" % (np.mean(wc), np.mean(wd)))
logb = np.log10(basal + 1e-3)
qs = np.quantile(np.concatenate([logb[gate[ti]] for ti in elig]), [0, .2, .4, .6, .8, 1.])
qb = []
for i in range(5):
    wc, wd, nn, mass = [], [], [], []
    for ti in elig:
        m_ = gate[ti]
        b = logb[m_]
        mk = (b >= qs[i]) & (b <= qs[i + 1] if i == 4 else b < qs[i + 1])
        if mk.sum() == 0:
            continue
        l = L[ti, m_][mk]
        w = np.abs(l)
        wc.append(float((w * (Csign[ti, m_][mk] == np.sign(l))).sum() / w.sum()))
        wd.append(float((w * (Dsign[ti, m_][mk] == np.sign(l))).sum() / w.sum()))
        nn.append(int(mk.sum())); mass.append(float(w.sum()))
    qb.append({"range": [float(qs[i]), float(qs[i + 1])], "combo_wacc": float(np.mean(wc)),
               "design_wacc": float(np.mean(wd)), "mean_genes": float(np.mean(nn)),
               "mass": float(np.sum(mass))})
tot = sum(q["mass"] for q in qb)
for q in qb:
    q["mass_share"] = q.pop("mass") / tot
    print("basal %s combo %.4f design %.4f mass %.3f" % (q["range"], q["combo_wacc"], q["design_wacc"], q["mass_share"]))
spl["basal_quintiles"] = qb
extra["oracle_best_splits"] = {"variant": "mean", "alpha": 1, "splits": spl}

res["step4b_extra"] = extra
json.dump(res, open(OUT / "results_part2.json", "w"), indent=1)

# unified results.json
p1 = json.load(open(OUT / "results_part1.json"))
p1.update(res)
# results_part2 already includes part1 keys (it was loaded from part1); write canonical copy
json.dump(p1, open(OUT / "results.json", "w"), indent=1)
print("wrote results_part2.json + results.json")
