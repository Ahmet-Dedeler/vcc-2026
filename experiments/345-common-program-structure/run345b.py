"""345 step 3-5: combination grids, transfer vectors (design-own + donor banks), splits.

Reads work/stage345.npz. Donor banks via mmap + chunked means (never densified fully).
"""
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
Lhat = st["Lhat_del"]
genes = st["genes"].astype(str)
targets = st["targets"].astype(str)
basal = st["basal"].astype(float)
Cmean, Cvote, Cwmean = st["Cmean"], st["Cvote"], st["Cwmean"]
G = len(genes)
gi = {g: i for i, g in enumerate(genes)}
L = R  # reference lfc on gate else 0
M = 0.205485  # delivered magnitude ratio, recomputed in part 1 (343: 0.205)

def scores(ti, csign):
    """weighted/unweighted sign acc of a full-G sign vector on gate(ti)."""
    m_ = gate[ti]
    l = L[ti, m_]
    w = np.abs(l)
    ok = (csign[m_] == np.sign(l))
    return float((w * ok).sum() / w.sum()), float(ok.mean())

def exact_nmae(ti, csign):
    """NMAE with delivered magnitudes and counterfactual signs."""
    m_ = gate[ti]
    l = L[ti, m_]
    mag = np.abs(Lhat[ti, m_])
    ln = np.sign(csign[m_]) * mag
    return float(np.abs(ln - l).sum() / np.abs(l).sum())

def grid_eval(D, C, s, label):
    """sign(D + a s C) over the alpha grid; returns dict of rows."""
    alphas = [0, .25, .5, 1, 2, 4, float("inf")]
    rows = []
    for a in alphas:
        comb = C if a == float("inf") else D + a * s * C
        pw, pu, ne = [], [], []
        for ti in elig:
            x, u = scores(ti, np.sign(comb[ti]))
            pw.append(x); pu.append(u); ne.append(exact_nmae(ti, np.sign(comb[ti])))
        pw_m = float(np.mean(pw))
        rows.append({"alpha": a, "wacc": pw_m, "uacc": float(np.mean(pu)),
                     "nmae_exact": float(np.mean(ne)),
                     "nmae_implied": float(1 - M * (2 * pw_m - 1)),
                     "per_target_wacc": pw})
        print(f"{label} a={a}: wacc {pw_m:.4f} uacc {np.mean(pu):.4f} "
              f"nmae_exact {np.mean(ne):.4f} implied {1 - M*(2*pw_m-1):.4f}")
    return rows

def pooled_scale(D, C):
    """s = median|D|/median|C| over pooled eligible gate pairs."""
    ad, ac = [], []
    for ti in elig:
        m_ = gate[ti]
        ad.append(np.abs(D[ti, m_]))
        ac.append(np.abs(C[ti, m_]))
    ad = np.concatenate(ad); ac = np.concatenate(ac)
    return float(np.median(ad) / np.median(ac)), float(np.median(ad)), float(np.median(ac))

res = json.load(open(OUT / "results_part1.json"))
Drow = Dfull  # per-target design rows (full universe)

# ---- step 3: oracle (same-line LOO) combinations ----
step3 = {}
for name, C in (("mean", Cmean), ("vote", Cvote), ("wmean", Cwmean)):
    s, medD, medC = pooled_scale(Drow, C)
    print(f"--- oracle {name}: s={s:.4f} medD={medD:.4f} medC={medC:.4f}")
    step3[name] = {"s": s, "medD": medD, "medC": medC, "rows": grid_eval(Drow, C, s, f"oracle-{name}")}
res["step3_oracle"] = step3

# ---- step 4a: design's own LOO common vector ----
Do = np.zeros_like(Dfull)
for ti in range(len(targets)):
    o = [j for j in range(len(targets)) if j != ti]
    Do[ti] = Dfull[o].mean(axis=0)
pw0, pu0 = [], []
for ti in elig:
    x, u = scores(ti, np.sign(Do[ti]))
    pw0.append(x); pu0.append(u)
print(f"design-own common alone: w {np.mean(pw0):.4f} u {np.mean(pu0):.4f}")
s4a, medD4, medC4 = pooled_scale(Drow, Do)
print(f"--- transfer design-own: s={s4a:.4f} medD={medD4:.4f} medC={medC4:.4f}")
res["step4a"] = {"alone": {"wacc": float(np.mean(pw0)), "uacc": float(np.mean(pu0)),
                           "per_target_wacc": pw0},
                 "s": s4a,
                 "rows": grid_eval(Drow, Do, s4a, "design-own")}

# ---- step 4b: donor-bank mean responses ----
def bank_mean(name):
    d = np.load(DON / name / "delta.npy", mmap_mode="r")
    bg = [x.strip() for x in open(DON / name / "genes.txt")]
    col = np.array([bg.index(g) if g in set(bg) else -1 for g in genes]) \
        if len(bg) < 20000 else np.array([{g: i for i, g in enumerate(bg)}.get(g, -1) for g in genes])
    hit = col >= 0
    print(f"{name}: delta {d.shape}, universe coverage {hit.sum()}/{G}")
    acc = np.zeros(G)
    n = d.shape[0]
    CH = 1000
    for a in range(0, n, CH):
        acc[hit] += d[a:a + CH][:, col[hit]].astype(np.float64).sum(axis=0)
    return acc / n, int(hit.sum())

banks = {}
for name in ("k562_gwps", "orion_hct116", "orion_hek293t"):
    mu, cov = bank_mean(name)
    Cb = np.repeat(mu[None, :], len(targets), axis=0)  # target-independent
    pw, pu = [], []
    for ti in elig:
        x, u = scores(ti, np.sign(mu))
        pw.append(x); pu.append(u)
    sb, _, _ = pooled_scale(Drow, Cb)
    print(f"{name} common alone: w {np.mean(pw):.4f} u {np.mean(pu):.4f}; s={sb:.4f}")
    banks[name] = {"coverage": cov, "alone": {"wacc": float(np.mean(pw)), "uacc": float(np.mean(pu)),
                                              "per_target_wacc": pw},
                   "s": sb, "rows": grid_eval(Drow, Cb, sb, f"bank-{name}")}
res["step4b"] = banks

# ---- step 5: splits for the best transferable combination ----
def pick_best(rows):
    return max(rows, key=lambda r_: r_["wacc"])

cands = [("design-own", res["step4a"]["rows"])] + [(f"bank-{k}", v["rows"]) for k, v in banks.items()]
best_name, best_rows = max(cands, key=lambda kv: pick_best(kv[1])["wacc"])
best = pick_best(best_rows)
print(f"best transferable: {best_name} alpha={best['alpha']} wacc={best['wacc']:.4f}")
# reconstruct its per-target signs: recompute combination
src = {"design-own": Do}.get(best_name, None)
if src is None:
    bname = best_name.replace("bank-", "")
    d = np.load(DON / bname / "delta.npy", mmap_mode="r")
    bg = [x.strip() for x in open(DON / bname / "genes.txt")]
    lut = {g: i for i, g in enumerate(bg)}
    mu = np.zeros(G)
    col = np.array([lut.get(g, -1) for g in genes])
    hit = col >= 0
    n = d.shape[0]
    for a in range(0, n, 1000):
        mu[hit] += d[a:a + 1000][:, col[hit]].astype(np.float64).sum(axis=0)
    mu /= n
    src = np.repeat(mu[None, :], len(targets), axis=0)
    sb = banks[bname]["s"]
else:
    sb = res["step4a"]["s"]
a = best["alpha"]
Csign = np.sign(src) if a == float("inf") else np.sign(Drow + a * sb * src)
Dsign = np.sign(Drow)
splits = {}
for dname, mask_fn in (("ref_UP", lambda l: l > 0), ("ref_DOWN", lambda l: l < 0)):
    ww_c, ww_d, nn = [], [], []
    for ti in elig:
        m_ = gate[ti]
        l = L[ti, m_]
        w = np.abs(l)
        mk = mask_fn(l)
        if mk.sum() == 0:
            continue
        ww_c.append(float((w[mk] * (Csign[ti, m_][mk] == np.sign(l[mk]))).sum() / w[mk].sum()))
        ww_d.append(float((w[mk] * (Dsign[ti, m_][mk] == np.sign(l[mk]))).sum() / w[mk].sum()))
        nn.append(int(mk.sum()))
    splits[dname] = {"combo_wacc": float(np.mean(ww_c)), "design_wacc": float(np.mean(ww_d)),
                     "n_targets": len(nn), "mean_genes": float(np.mean(nn))}
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
               "mass_share": float(np.sum(mass))})
tot_mass = sum(q["mass_share"] for q in qb)
for q in qb:
    q["mass_share"] /= tot_mass
splits["basal_quintiles"] = qb
# oracle-best splits too (ceiling context)
res["step5"] = {"best_transferable": best_name, "best_alpha": best["alpha"],
                "best_wacc": best["wacc"], "splits": splits}
print(json.dumps(splits, indent=1)[:2000])

json.dump(res, open(OUT / "results_part2.json", "w"), indent=1)
print("wrote results_part2.json")
