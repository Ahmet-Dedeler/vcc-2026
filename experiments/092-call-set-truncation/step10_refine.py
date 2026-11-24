"""Step 10: refinements.

  A. fine constant grid + bootstrap CI on the optimum
  B. does ANY predictor beat a constant at matched mean m?  (the real test)
  C. is the fitted predictor just a disguised constant?
  D. practically-implementable variants: a stricter fixed alpha; a fixed fraction of n_pred
  E. robustness of the constant to a shifted n_real distribution
  F. reach robustness if the scorer's secondary sort key is |lfc| desc, not p_adj asc
"""
import os, numpy as np
import trunc_lib as T

S = os.path.dirname(os.path.abspath(__file__))
pre = T.precompute()
P = len(pre)
n_real = np.array([s["n_real"] for s in pre], float)
n_all = np.array([s["n_all"] for s in pre], float)
pr = np.load(os.path.join(S, "preds.npz"), allow_pickle=True)
genes, perts, rp, rq, rl, pp, pq, pl = T.load()

bf, bj, _ = T.fid_jac(pre, [None] * P)
BF, BJ = float(np.nanmean(bf)), float(np.nanmean(bj))
bof, boj = T.scaled(BF, BJ)


def ev(m_vec, w=None):
    f, j, m = T.fid_jac(pre, m_vec)
    if w is None:
        F, J = float(np.nanmean(f)), float(np.nanmean(j))
    else:
        F = float(np.nansum(w * f) / w.sum()); J = float(np.nansum(w * j) / w.sum())
    of, oj = T.scaled(F, J)
    return F, J, ((of - bof) + (oj - boj)) / 6.0, m.mean()


def per_pert(m_vec):
    f, j, _ = T.fid_jac(pre, m_vec)
    return f, j


# ---------------- A. fine constant grid + bootstrap ------------------------
grid = np.arange(300, 2001, 10)
D = np.array([ev([float(c)] * P)[2] for c in grid])
cbest = int(grid[np.argmax(D)])
print("=" * 92)
print("A. GLOBAL CONSTANT -- fine grid + bootstrap")
print("=" * 92)
print(f"argmax c = {cbest}   d_avg6 = {D.max():+.4f}")
FJ = {int(c): per_pert([float(c)] * P) for c in grid}
rng = np.random.default_rng(0)
boot = []
for _ in range(600):
    idx = rng.integers(0, P, P)
    best, bc = -9, None
    for c in grid[::5]:
        f, j = FJ[int(c)]
        of, oj = T.scaled(float(np.nanmean(f[idx])), float(np.nanmean(j[idx])))
        v = (of + oj)
        if v > best: best, bc = v, int(c)
    boot.append(bc)
boot = np.array(boot)
print(f"bootstrap argmax c: median {np.median(boot):.0f}  "
      f"90% CI [{np.percentile(boot,5):.0f}, {np.percentile(boot,95):.0f}]")
# d_avg6 at a few round constants
for c in (500, 700, 800, 900, 1000, 1200):
    print(f"   c={c:<5} d_avg6={ev([float(c)]*P)[2]:+.4f}")

# ---------------- B. matched-mean-m: predictor vs constant -----------------
print("\n" + "=" * 92)
print("B. DOES ANY PREDICTOR BEAT A CONSTANT AT THE SAME MEAN m?")
print("=" * 92)


def loo_fit(x, y):
    X = np.log(np.maximum(np.asarray(x, float), 1e-9)); Y = np.log(np.maximum(y, 1))
    n = len(X); out = np.empty(n); sl = np.empty(n)
    for i in range(n):
        k = np.ones(n, bool); k[i] = False
        b, a = np.polyfit(X[k], Y[k], 1)
        out[i] = np.exp(a + b * X[i]); sl[i] = b
    return out, sl


names = ["l2_lfc_sig", "n_lfc_gt0.05", "n_pred_sig", "k562_n_gt2", "hek_n_gt2",
         "hct_n_gt", "hct_l2"]
print(f"{'predictor':<16}{'mean m':>9}{'d6 pred':>10}{'d6 const@same mean m':>22}{'gain':>10}")
rowsB = []
for nm in names:
    yh, _ = loo_fit(pr[nm], n_real)
    best = (-9,)
    for sf in np.arange(0.5, 8.01, 0.25):
        F, J, d, mm = ev(list(sf * yh))
        if d > best[0]: best = (d, sf, mm)
    d_p, sf, mm = best
    # constant with the SAME mean m
    d_c = ev([float(mm)] * P)[2]
    rowsB.append((nm, mm, d_p, d_c, d_p - d_c))
    print(f"{nm:<16}{mm:>9.1f}{d_p:>+10.4f}{d_c:>+22.4f}{d_p - d_c:>+10.4f}")
# and the ORACLE at matched mean m
best = (-9,)
for sf in np.arange(0.5, 8.01, 0.25):
    F, J, d, mm = ev([sf * s["n_real"] for s in pre])
    if d > best[0]: best = (d, sf, mm)
d_o, sf_o, mm_o = best
print(f"{'ORACLE n_real':<16}{mm_o:>9.1f}{d_o:>+10.4f}{ev([float(mm_o)]*P)[2]:>+22.4f}"
      f"{d_o - ev([float(mm_o)]*P)[2]:>+10.4f}")

# ---------------- C. is the fitted predictor a disguised constant? ---------
print("\n" + "=" * 92)
print("C. IS THE FITTED PREDICTOR JUST A CONSTANT IN DISGUISE?")
print("=" * 92)
print(f"{'predictor':<16}{'LOO slope':>11}{'sd log yhat':>13}{'sd log n_real':>15}"
      f"{'yhat p10/p90':>16}{'n_real p10/p90':>17}")
ln = np.log(np.maximum(n_real, 1))
for nm in names:
    yh, sl = loo_fit(pr[nm], n_real)
    ly = np.log(yh)
    print(f"{nm:<16}{sl.mean():>11.3f}{ly.std():>13.3f}{ln.std():>15.3f}"
          f"{np.percentile(yh,10):>8.0f}/{np.percentile(yh,90):<7.0f}"
          f"{np.percentile(n_real,10):>9.0f}/{np.percentile(n_real,90):<7.0f}")

# ---------------- D. implementable variants --------------------------------
print("\n" + "=" * 92)
print("D. IMPLEMENTABLE VARIANTS")
print("=" * 92)
print("D1. stricter fixed alpha:  m_i = #{genes : our p_adj < alpha'}")
print(f"{'alpha':>10}{'mean m':>10}{'sd m':>9}{'fid':>9}{'jac':>9}{'d avg6':>10}")
G = len(genes)
gpos = {g: i for i, g in enumerate(genes)}
okm = np.ones((P, G), bool)
for i, t in enumerate(perts):
    j = gpos.get(t)
    if j is not None: okm[i, j] = False
for al in [5e-2, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10, 1e-15, 1e-20, 1e-30]:
    mv = ((pq < al) & okm).sum(1).astype(float)
    F, J, d, mm = ev(list(mv))
    print(f"{al:>10.0e}{mm:>10.1f}{mv.std():>9.1f}{F:>9.4f}{J:>9.4f}{d:>+10.4f}")
print("\nD2. fixed fraction of our own call list:  m_i = frac * n_pred_i")
print(f"{'frac':>10}{'mean m':>10}{'sd m':>9}{'fid':>9}{'jac':>9}{'d avg6':>10}")
for fr_ in [0.02, 0.05, 0.08, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0]:
    mv = fr_ * n_all
    F, J, d, mm = ev(list(mv))
    print(f"{fr_:>10g}{mm:>10.1f}{mv.std():>9.1f}{F:>9.4f}{J:>9.4f}{d:>+10.4f}")
print("\nD3. constant with a cap at our own call count (already implied) and a floor:")
for c in (850,):
    for flo in (0, 50, 100, 200):
        mv = np.maximum(np.full(P, float(c)), flo)
        F, J, d, mm = ev(list(mv))
        print(f"   c={c} floor={flo:<5} d_avg6={d:+.4f}")

# ---------------- E. robustness to a shifted n_real distribution -----------
print("\n" + "=" * 92)
print("E. ROBUSTNESS OF THE CONSTANT TO A DIFFERENT n_real DISTRIBUTION")
print("=" * 92)
q = np.quantile(n_real, [0, 1 / 3, 2 / 3, 1.0])
print(f"{'subset':<28}{'n':>5}{'med n_real':>12}{'best c':>9}{'d_avg6@best':>13}"
      f"{'d_avg6@850':>12}{'ret%':>7}")
subsets = [("all 132", np.ones(P, bool)),
           ("tercile 1 (small n_real)", n_real <= q[1]),
           ("tercile 2", (n_real > q[1]) & (n_real <= q[2])),
           ("tercile 3 (large n_real)", n_real > q[2])]
for nm, mask in subsets:
    w = mask.astype(float)
    bb = (-9, None)
    for c in grid[::5]:
        f, j = FJ[int(c)]
        of, oj = T.scaled(float(np.nansum(w * f) / w.sum()), float(np.nansum(w * j) / w.sum()))
        # baseline restricted to the subset
        obf, obj = T.scaled(float(np.nansum(w * bf) / w.sum()), float(np.nansum(w * bj) / w.sum()))
        d = ((of - obf) + (oj - obj)) / 6
        if d > bb[0]: bb = (d, int(c))
    f, j = FJ[850]
    of, oj = T.scaled(float(np.nansum(w * f) / w.sum()), float(np.nansum(w * j) / w.sum()))
    obf, obj = T.scaled(float(np.nansum(w * bf) / w.sum()), float(np.nansum(w * bj) / w.sum()))
    d850 = ((of - obf) + (oj - obj)) / 6
    print(f"{nm:<28}{int(mask.sum()):>5}{np.median(n_real[mask]):>12.0f}{bb[1]:>9}"
          f"{bb[0]:>+13.4f}{d850:>+12.4f}{100*d850/bb[0] if bb[0]>0 else float('nan'):>6.0f}%")

# ---------------- F. reach under an alternative scorer sort key ------------
print("\n" + "=" * 92)
print("F. reach IF THE SCORER'S SECONDARY KEY IS |lfc| desc INSTEAD OF p_adj asc")
print("=" * 92)


def reach_alt(m_vec):
    out = np.empty(P)
    for i, s in enumerate(pre):
        nr = s["n_real"]
        if nr == 0: out[i] = 0.0; continue
        m = s["n_all"] if m_vec[i] is None else int(min(max(round(m_vec[i]), 0), s["n_all"]))
        rr = s["r_rank"]; called = (rr >= 0) & (rr < m)
        order = np.lexsort((s["r_feat"], -s["r_absl"], ~called))
        pur = np.cumsum(s["r_agree"][order]) / np.arange(1, nr + 1)
        good = np.flatnonzero(pur >= 0.9)
        out[i] = (good[-1] + 1) / nr if good.size else 0.0
    return out


b_alt = reach_alt([None] * P).mean()
print(f"{'m policy':<22}{'reach(alt key)':>16}{'delta vs base':>15}{'scaled delta':>14}"
      f"{'d avg6 contrib':>16}")
print(f"{'baseline (no trunc)':<22}{b_alt:>16.6f}{0.0:>15.6f}{0.0:>14.4f}{0.0:>16.4f}")
for lab, mv in [("constant c=850", [850.0] * P), ("constant c=400", [400.0] * P),
                ("constant c=2000", [2000.0] * P),
                ("oracle r=1.2", [1.2 * s["n_real"] for s in pre])]:
    v = reach_alt(mv).mean()
    print(f"{lab:<22}{v:>16.6f}{v - b_alt:>+15.6f}{(v - b_alt) / T.RS:>+14.4f}"
          f"{(v - b_alt) / T.RS / 6:>+16.4f}")

# nmae in scaled units under the two implementations
NB, NR_ = 1.0013, 0.400
def nsc(x): return (x - NB) / (NR_ - NB)
print("\nnmae scaled (lower raw is better):")
for lab, mv in [("baseline", [None] * P), ("c=850 label-mode", [850.0] * P)]:
    raw = T.nmae(pre, mv, "label")
    print(f"  {lab:<20} raw={raw:.6f}  scaled={nsc(raw):+.4f}")
raw_null = T.nmae(pre, [850.0] * P, "null")
print(f"  {'c=850 NULL-mode':<20} raw={raw_null:.6f}  scaled={nsc(raw_null):+.4f}  "
      f"=> d_avg6 {(nsc(raw_null) - nsc(T.nmae(pre, [None]*P))) / 6:+.4f}")
