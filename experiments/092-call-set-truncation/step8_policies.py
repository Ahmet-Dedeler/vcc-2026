"""Step 8 (Q2 + Q3): how much of the truncation gain survives an imperfect estimate of n_real?

Policies:
  (a) ORACLE     m = r * n_real
  (b) PREDICTOR  m = safety * exp(a + b log x),  a,b fit leave-one-out (honest out-of-sample)
  (c) CONSTANT   m = c for every perturbation      <-- the one that matters
"""
import os, numpy as np
import trunc_lib as T

S = os.path.dirname(os.path.abspath(__file__))
pre = T.precompute()
P = len(pre)
n_real = np.array([s["n_real"] for s in pre], float)
n_all = np.array([s["n_all"] for s in pre], float)
pr = np.load(os.path.join(S, "preds.npz"), allow_pickle=True)

# ---- baseline (no truncation) --------------------------------------------
bf, bj, bm = T.fid_jac(pre, [None] * P)
BF, BJ = float(np.nanmean(bf)), float(np.nanmean(bj))
bof, boj = T.scaled(BF, BJ)
print(f"BASELINE (no truncation): fid={BF:.4f} jac={BJ:.4f}  scaled fid={bof:+.4f} jac={boj:+.4f} "
      f"mean m={bm.mean():.0f}")


def ev(m_vec):
    f, j, m = T.fid_jac(pre, m_vec)
    F, J = float(np.nanmean(f)), float(np.nanmean(j))
    of, oj = T.scaled(F, J)
    return F, J, ((of - bof) + (oj - boj)) / 6.0, m.mean()


# ================= (a) ORACLE ==============================================
print("\n" + "=" * 92)
print("(a) ORACLE   m = r * n_real   [uses the reference; upper bound]")
print("=" * 92)
print(f"{'r':>7}{'mean m':>10}{'fid':>10}{'jac':>10}{'d avg6':>10}")
best_o = (-9, None)
for r in [0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 2.0, 3.0]:
    F, J, D, mm = ev([r * s["n_real"] for s in pre])
    if D > best_o[0]: best_o = (D, r, F, J, mm)
    print(f"{r:>7g}{mm:>10.1f}{F:>10.4f}{J:>10.4f}{D:>+10.4f}")
print(f"  best oracle: r={best_o[1]:g}  d_avg6={best_o[0]:+.4f}")

# per-pert argmax (the absolute ceiling)
ceil_f = np.empty(P)
for i, s in enumerate(pre):
    nr, na = s["n_real"], s["n_all"]
    mm = np.arange(1, na + 1)
    fa = s["cum_hit"] / np.maximum(s["cum_adj"], nr)
    ceil_f[i] = fa.max()
_, cj, _ = T.fid_jac(pre, [None] * P)
cof, _ = T.scaled(float(ceil_f.mean()), BJ)
print(f"  per-pert argmax fid ceiling: fid={ceil_f.mean():.4f}  (d_avg6 from fid alone "
      f"{(cof-bof)/6:+.4f})")

# ================= (c) GLOBAL CONSTANT =====================================
print("\n" + "=" * 92)
print("(c) GLOBAL CONSTANT   m = c for every perturbation   [needs no estimator at all]")
print("=" * 92)
grid = list(range(50, 1001, 25)) + list(range(1050, 2001, 50)) + [2250, 2500, 3000, 4000, 5356]
res = []
for c in grid:
    F, J, D, mm = ev([float(c)] * P)
    res.append((c, F, J, D, mm))
best_c = max(res, key=lambda x: x[3])
print(f"{'c':>7}{'mean m':>10}{'fid':>10}{'jac':>10}{'d avg6':>10}")
for c, F, J, D, mm in res:
    if c % 100 == 0 or c in (50, 75, 5356) or c == best_c[0]:
        star = "  <== BEST" if c == best_c[0] else ""
        print(f"{c:>7}{mm:>10.1f}{F:>10.4f}{J:>10.4f}{D:>+10.4f}{star}")
print(f"\n  BEST CONSTANT c={best_c[0]}  fid={best_c[1]:.4f} jac={best_c[2]:.4f} "
      f"d_avg6={best_c[3]:+.4f}   ({100*best_c[3]/best_o[0]:.1f}% of oracle)")

# constant robustness: how flat is the optimum?
ok = [(c, D) for c, F, J, D, mm in res if D >= 0.9 * best_c[3]]
print(f"  within 90% of best over c in [{min(x[0] for x in ok)}, {max(x[0] for x in ok)}]")
ok95 = [(c, D) for c, F, J, D, mm in res if D >= 0.95 * best_c[3]]
print(f"  within 95% of best over c in [{min(x[0] for x in ok95)}, {max(x[0] for x in ok95)}]")

# ================= (b) PREDICTOR-DRIVEN ====================================
def loo_fit_predict(x, y):
    """LOO log-log linear fit; returns out-of-sample yhat (count space)."""
    X = np.log(np.maximum(x, 1e-9)); Y = np.log(np.maximum(y, 1))
    n = len(X); out = np.empty(n)
    for i in range(n):
        k = np.ones(n, bool); k[i] = False
        b, a = np.polyfit(X[k], Y[k], 1)
        out[i] = np.exp(a + b * X[i])
    return out


def loo_fit_predict_multi(Xcols, y):
    X = np.column_stack([np.log(np.maximum(c, 1e-9)) for c in Xcols])
    X = np.column_stack([X, np.ones(len(y))])
    Y = np.log(np.maximum(y, 1)); n = len(Y); out = np.empty(n)
    for i in range(n):
        k = np.ones(n, bool); k[i] = False
        beta, *_ = np.linalg.lstsq(X[k], Y[k], rcond=None)
        out[i] = np.exp(X[i] @ beta)
    return out


cands = {
    "l2_lfc_sig  (self, best clean)": pr["l2_lfc_sig"],
    "n_lfc_gt0.05 (self)": pr["n_lfc_gt0.05"],
    "n_pred_sig  (self)": pr["n_pred_sig"],
    "hek_n_gt2   (cross-line donor)": pr["hek_n_gt2"],
    "k562_n_gt2  (cross-line donor)": pr["k562_n_gt2"],
    "hct_n_gt    (LEAKY same-source donor)": pr["hct_n_gt"],
    "hct_l2      (LEAKY same-source donor)": pr["hct_l2"],
}
multi = {
    "self 3-feature (l2_lfc_sig+n_lfc_gt0.05+max_abs_lfc)":
        [pr["l2_lfc_sig"], pr["n_lfc_gt0.05"], pr["max_abs_lfc"]],
    "self+cross-donor 4-feature":
        [pr["l2_lfc_sig"], pr["n_lfc_gt0.05"], pr["hek_n_gt2"], pr["k562_n_gt2"]],
    "LEAKY hct 2-feature (hct_n_gt+hct_l2)": [pr["hct_n_gt"], pr["hct_l2"]],
}

print("\n" + "=" * 92)
print("(b) PREDICTOR-DRIVEN   m = safety * yhat,  yhat from LOO log-log fit")
print("=" * 92)
print(f"{'predictor':<42}{'safety':>8}{'mean m':>9}{'fid':>9}{'jac':>9}{'d avg6':>9}{'% orc':>8}")
SAFE = [0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
pred_best = {}
for name, x in list(cands.items()):
    yh = loo_fit_predict(np.asarray(x, float), n_real)
    rowbest = (-9,)
    for sf in SAFE:
        F, J, D, mm = ev(list(sf * yh))
        if D > rowbest[0]: rowbest = (D, sf, F, J, mm)
    pred_best[name] = rowbest
    D, sf, F, J, mm = rowbest
    print(f"{name:<42}{sf:>8g}{mm:>9.1f}{F:>9.4f}{J:>9.4f}{D:>+9.4f}{100*D/best_o[0]:>7.0f}%")
for name, cols in multi.items():
    yh = loo_fit_predict_multi([np.asarray(c, float) for c in cols], n_real)
    rowbest = (-9,)
    for sf in SAFE:
        F, J, D, mm = ev(list(sf * yh))
        if D > rowbest[0]: rowbest = (D, sf, F, J, mm)
    pred_best[name] = rowbest
    D, sf, F, J, mm = rowbest
    print(f"{name:<42}{sf:>8g}{mm:>9.1f}{F:>9.4f}{J:>9.4f}{D:>+9.4f}{100*D/best_o[0]:>7.0f}%")

# ================= Q3  SAFETY MULTIPLIER ===================================
print("\n" + "=" * 92)
print("Q3  SAFETY MULTIPLIER SWEEP   (fid punishes undershoot ~linearly, overshoot ~free)")
print("=" * 92)
SW = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0]
best_self = "l2_lfc_sig  (self, best clean)"
yh_self = loo_fit_predict(np.asarray(pr[ "l2_lfc_sig"], float), n_real)
yh_leak = loo_fit_predict(np.asarray(pr["hct_n_gt"], float), n_real)
print(f"{'safety':>8} | {'ORACLE m=s*n_real':>28} | {'SELF-PRED m=s*yhat':>28} | "
      f"{'LEAKY-PRED m=s*yhat':>28}")
print(f"{'':>8} | {'fid':>8}{'jac':>9}{'davg6':>10} | {'fid':>8}{'jac':>9}{'davg6':>10} | "
      f"{'fid':>8}{'jac':>9}{'davg6':>10}")
for s_ in SW:
    a = ev([s_ * s["n_real"] for s in pre])
    b = ev(list(s_ * yh_self))
    c = ev(list(s_ * yh_leak))
    print(f"{s_:>8g} | {a[0]:>8.4f}{a[1]:>9.4f}{a[2]:>+10.4f} | {b[0]:>8.4f}{b[1]:>9.4f}"
          f"{b[2]:>+10.4f} | {c[0]:>8.4f}{c[1]:>9.4f}{c[2]:>+10.4f}")

# expected-fid asymmetry: single-pert view
print("\nasymmetry of fid around the per-pert optimum (pooled over 132 perts):")
print(f"{'m/n_real':>10}{'mean fid':>11}{'loss vs r=1':>13}")
f1 = ev([s["n_real"] for s in pre])[0]
for r in [0.25, 0.5, 0.75, 0.9, 1.0, 1.1, 1.5, 2, 3, 5, 8, 12, 20]:
    F = ev([r * s["n_real"] for s in pre])[0]
    print(f"{r:>10g}{F:>11.4f}{F - f1:>+13.4f}")

np.savez(os.path.join(S, "policy_out.npz"),
         const_grid=np.array([r[0] for r in res]), const_d6=np.array([r[3] for r in res]),
         const_fid=np.array([r[1] for r in res]), const_jac=np.array([r[2] for r in res]),
         best_c=best_c[0], best_oracle_r=best_o[1], yh_self=yh_self, yh_leak=yh_leak,
         n_real=n_real, base=np.array([BF, BJ, bof, boj]))
print("\nwrote policy_out.npz")
