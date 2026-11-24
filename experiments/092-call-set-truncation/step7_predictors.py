"""Step 7: can we predict n_real (the REFERENCE's significant-gene count) without the reference?

Builds a predictor table from (a) our own prediction's DE statistics -- legitimately available at
submission time, since the scorer's control cells are published -- and (b) donor-bank effect
breadth priors. Saves preds.npz for step8/step9.
"""
import os, json, numpy as np
from scipy import stats

S = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(S))
d = np.load(os.path.join(S, "de_tables.npz"), allow_pickle=True)
genes = d["genes"].astype(str); perts = d["perts"].astype(str)
rq, rl = d["real_padj"], d["real_lfc"]
pp, pq, pl = d["pred_p"], d["pred_padj"], d["pred_lfc"]
A = 0.05; G = len(genes); P = len(perts)
gpos = {g: i for i, g in enumerate(genes)}

# ---- ground truth ---------------------------------------------------------
n_real = np.zeros(P, dtype=float)
n_pred = np.zeros(P, dtype=float)
for i, t in enumerate(perts):
    ok = np.ones(G, bool); j = gpos.get(t)
    if j is not None: ok[j] = False
    n_real[i] = (ok & (rq[i] < A)).sum()
    ind = ok & (rl[i] != 0) & ~np.isnan(rl[i])
    n_pred[i] = (ok & (pq[i] < A) & ind).sum()

# ---- self predictors (no reference used) ----------------------------------
F = {}
lfc = np.nan_to_num(pl, nan=0.0)
absl = np.abs(lfc)
okmask = np.ones((P, G), bool)
for i, t in enumerate(perts):
    j = gpos.get(t)
    if j is not None: okmask[i, j] = False
absl_ok = absl * okmask

F["n_pred_sig"] = n_pred
for a, lab in [(1e-3, "1e3"), (1e-6, "1e6"), (1e-10, "1e10"), (1e-20, "1e20"), (1e-40, "1e40")]:
    F[f"n_sig_{lab}"] = ((pq < a) & okmask).sum(1).astype(float)
F["l2_lfc"] = np.sqrt((absl_ok ** 2).sum(1))
F["l1_lfc"] = absl_ok.sum(1)
F["max_abs_lfc"] = absl_ok.max(1)
for t_ in (0.05, 0.1, 0.25, 0.5, 1.0):
    F[f"n_lfc_gt{t_:g}"] = (absl_ok > t_).sum(1).astype(float)
with np.errstate(divide="ignore", invalid="ignore"):
    nlp = -np.log10(np.clip(np.nan_to_num(pp, nan=1.0), 1e-300, 1.0))
F["sum_neglog10p"] = (nlp * okmask).sum(1)
F["l2_lfc_sig"] = np.sqrt(((absl_ok ** 2) * ((pq < A) & okmask)).sum(1))
sig = (pq < A) & okmask
F["mean_abs_lfc_sig"] = (absl_ok * sig).sum(1) / np.maximum(sig.sum(1), 1)
F["min_padj_nl10"] = -np.log10(np.clip(np.nanmin(np.where(okmask, pq, np.inf), axis=1), 1e-300, 1))
# 90th pct of |lfc| -- a scale-free "how big is this effect" summary
F["q90_abs_lfc"] = np.array([np.percentile(absl_ok[i], 90) for i in range(P)])
F["q99_abs_lfc"] = np.array([np.percentile(absl_ok[i], 99) for i in range(P)])

# ---- donor-bank breadth priors -------------------------------------------
def donor_breadth(name):
    dd = os.path.join(ROOT, "data", "donors", name)
    dperts = np.array(open(os.path.join(dd, "perts.txt")).read().split("\n"))
    dperts = np.array([x for x in dperts if x])
    dgenes = np.array(open(os.path.join(dd, "genes.txt")).read().split("\n"))
    dgenes = np.array([x for x in dgenes if x])
    pos = {p: i for i, p in enumerate(dperts)}
    delta = np.load(os.path.join(dd, "delta.npy"), mmap_mode="r")
    out = {k: np.full(P, np.nan) for k in ("l2", "n_gt", "n_gt2", "max")}
    rows = {}
    for i, t in enumerate(perts):
        j = pos.get(t)
        if j is None: continue
        v = np.asarray(delta[j], dtype=np.float64)
        rows[i] = (v, dgenes)
        av = np.abs(v)
        out["l2"][i] = np.sqrt((av ** 2).sum())
        out["n_gt"][i] = (av > 0.1).sum()
        out["n_gt2"][i] = (av > 0.25).sum()
        out["max"][i] = av.max()
    del delta
    return out, rows, dgenes

donor_rows = {}
for name, tag in [("k562_gwps", "k562"), ("orion_hek293t", "hek"), ("orion_hct116", "hct")]:
    try:
        o, rows, dg = donor_breadth(name)
        donor_rows[tag] = (rows, dg)
        for k, v in o.items():
            F[f"{tag}_{k}"] = v
        cov = np.isfinite(o["l2"]).sum()
        print(f"donor {name:16s}: {cov}/{P} targets covered")
    except Exception as e:
        print(f"donor {name} FAILED: {e}")

# cross-donor cosine agreement (k562 vs hek293t), on shared genes
if "k562" in donor_rows and "hek" in donor_rows:
    (rk, gk), (rh, gh) = donor_rows["k562"], donor_rows["hek"]
    gkp = {g: i for i, g in enumerate(gk)}
    common = [(gkp[g], i) for i, g in enumerate(gh) if g in gkp]
    ik = np.array([a for a, b in common]); ih = np.array([b for a, b in common])
    agree = np.full(P, np.nan)
    for i in range(P):
        if i in rk and i in rh:
            a = rk[i][0][ik]; b = rh[i][0][ih]
            na, nb = np.linalg.norm(a), np.linalg.norm(b)
            if na > 0 and nb > 0: agree[i] = float(a @ b / (na * nb))
    F["xdonor_cos_k562_hek"] = agree
    print(f"cross-donor cosine covered: {np.isfinite(agree).sum()}/{P}")

# ---- correlations + LOO out-of-sample R^2 --------------------------------
y = n_real.copy()
ly = np.log(np.maximum(y, 1))

def loo_r2(x, y, logspace=True):
    """LOO R^2 of a 1-D linear fit. logspace: fit log(y) ~ log(x)."""
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 10: return np.nan, np.nan
    xs = x[m]; ys = y[m]
    X = np.log(np.maximum(xs, 1e-9)) if logspace else xs
    Y = np.log(np.maximum(ys, 1)) if logspace else ys
    n = len(X); pr = np.empty(n)
    for i in range(n):
        k = np.ones(n, bool); k[i] = False
        b, a = np.polyfit(X[k], Y[k], 1)
        pr[i] = a + b * X[i]
    r2_log = 1 - ((Y - pr) ** 2).sum() / ((Y - Y.mean()) ** 2).sum()
    yhat = np.exp(pr) if logspace else pr
    r2_raw = 1 - ((ys - yhat) ** 2).sum() / ((ys - ys.mean()) ** 2).sum()
    return r2_log, r2_raw

rows = []
for k, v in F.items():
    m = np.isfinite(v) & np.isfinite(y)
    if m.sum() < 10:
        rows.append((k, int(m.sum()), np.nan, np.nan, np.nan, np.nan, np.nan)); continue
    sp = stats.spearmanr(v[m], y[m]).statistic
    pe = stats.pearsonr(v[m], y[m]).statistic
    pe_log = stats.pearsonr(np.log(np.maximum(v[m], 1e-9)), ly[m]).statistic
    r2l, r2r = loo_r2(v, y)
    rows.append((k, int(m.sum()), sp, pe, pe_log, r2l, r2r))

rows.sort(key=lambda r: -abs(r[2] if np.isfinite(r[2]) else 0))
print("\n" + "=" * 100)
print("Q1  PREDICTORS OF n_real   (n=132 perturbations, HCT116 harness)")
print("=" * 100)
print(f"{'predictor':<24}{'n':>5}{'spearman':>10}{'pearson':>10}{'pearson_ll':>12}"
      f"{'LOO R2(log)':>13}{'LOO R2(raw)':>13}")
for r in rows:
    print(f"{r[0]:<24}{r[1]:>5}" + "".join(f"{x:>10.3f}" if np.isfinite(x) else f"{'nan':>10}"
                                           for x in r[2:4])
          + "".join(f"{x:>12.3f}" if np.isfinite(x) else f"{'nan':>12}" for x in r[4:5])
          + "".join(f"{x:>13.3f}" if np.isfinite(x) else f"{'nan':>13}" for x in r[5:]))

print(f"\nn_real: min={y.min():.0f} q1={np.percentile(y,25):.0f} med={np.median(y):.0f} "
      f"q3={np.percentile(y,75):.0f} max={y.max():.0f} mean={y.mean():.1f} sd={y.std():.1f}")
print(f"log n_real sd = {ly.std():.3f}  (a predictor must beat this to be useful)")

np.savez(os.path.join(S, "preds.npz"), perts=perts, n_real=n_real, n_pred=n_pred,
         **{k: v for k, v in F.items()})
print(f"\nwrote preds.npz with {len(F)} predictors")
