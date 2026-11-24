"""Step 9 (Q4): does truncation move `reach`? does it move `nmae`?

reach orders the REFERENCE's significant pool by (our significant calls first, then p_adj asc),
walks cumulative sign-purity, takes k* = deepest k with purity >= 0.9, R = k*/n_real.

Three models of what "truncation" means in a real submission:
  label : we simply stop *labelling* the tail significant; p_adj and lfc are untouched.
  tail  : the tail is genuinely de-significated (p_adj -> 1) so it sorts to the back;
          lfc/sign survives.
  null  : the tail is de-significated AND its effect is nulled, so it can never sign-agree.
          Worst case -- an implementation that flattens those genes to control.
"""
import os, numpy as np
import trunc_lib as T

S = os.path.dirname(os.path.abspath(__file__))
pre = T.precompute()
P = len(pre)
n_real = np.array([s["n_real"] for s in pre], float)

base_r = {m: T.reach(pre, [None] * P, m) for m in ("label", "tail", "null")}
print(f"BASELINE reach (no truncation) = {base_r['label'].mean():.6f}   "
      f"scaled {(base_r['label'].mean()-T.R0)/T.RS:+.4f}")
assert np.allclose(base_r["label"], base_r["tail"]), "baseline modes must agree"

bf, bj, _ = T.fid_jac(pre, [None] * P)
bof, boj = T.scaled(float(np.nanmean(bf)), float(np.nanmean(bj)))
bor = (base_r["label"].mean() - T.R0) / T.RS
base_nmae = T.nmae(pre, [None] * P)
print(f"BASELINE nmae (raw) = {base_nmae:.6f}")

print("\n" + "=" * 104)
print("reach UNDER TRUNCATION -- oracle m = r * n_real")
print("=" * 104)
print(f"{'r':>7}{'reach(label)':>14}{'reach(tail)':>13}{'reach(null)':>13}"
      f"{'nmae(label)':>13}{'nmae(null)':>12}{'d avg6 label':>14}{'d avg6 null':>13}")
for r in [0.5, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0, 3.0, 5.0]:
    mv = [r * s["n_real"] for s in pre]
    rl_ = T.reach(pre, mv, "label").mean()
    rt_ = T.reach(pre, mv, "tail").mean()
    rn_ = T.reach(pre, mv, "null").mean()
    nl_ = T.nmae(pre, mv, "label"); nn_ = T.nmae(pre, mv, "null")
    f, j, _ = T.fid_jac(pre, mv)
    of, oj = T.scaled(float(np.nanmean(f)), float(np.nanmean(j)))
    d_lab = ((of - bof) + (oj - boj) + ((rl_ - T.R0) / T.RS - bor)) / 6
    d_null = ((of - bof) + (oj - boj) + ((rn_ - T.R0) / T.RS - bor)) / 6
    print(f"{r:>7g}{rl_:>14.6f}{rt_:>13.6f}{rn_:>13.6f}{nl_:>13.6f}{nn_:>12.6f}"
          f"{d_lab:>+14.4f}{d_null:>+13.4f}")

print("\n" + "=" * 104)
print("reach UNDER TRUNCATION -- GLOBAL CONSTANT m = c")
print("=" * 104)
print(f"{'c':>7}{'reach(label)':>14}{'reach(tail)':>13}{'reach(null)':>13}"
      f"{'nmae(label)':>13}{'nmae(null)':>12}{'d avg6 label':>14}{'d avg6 null':>13}")
for c in [200, 400, 600, 700, 850, 1000, 1200, 1500, 2000, 3000]:
    mv = [float(c)] * P
    rl_ = T.reach(pre, mv, "label").mean()
    rt_ = T.reach(pre, mv, "tail").mean()
    rn_ = T.reach(pre, mv, "null").mean()
    nl_ = T.nmae(pre, mv, "label"); nn_ = T.nmae(pre, mv, "null")
    f, j, _ = T.fid_jac(pre, mv)
    of, oj = T.scaled(float(np.nanmean(f)), float(np.nanmean(j)))
    d_lab = ((of - bof) + (oj - boj) + ((rl_ - T.R0) / T.RS - bor)) / 6
    d_null = ((of - bof) + (oj - boj) + ((rn_ - T.R0) / T.RS - bor)) / 6
    print(f"{c:>7}{rl_:>14.6f}{rt_:>13.6f}{rn_:>13.6f}{nl_:>13.6f}{nn_:>12.6f}"
          f"{d_lab:>+14.4f}{d_null:>+13.4f}")

# why label-mode is invariant: the called block IS the p_adj-ascending prefix
print("\n" + "=" * 104)
print("INVARIANCE PROOF CHECK")
print("=" * 104)
mv = [850.0] * P
same = np.allclose(T.reach(pre, mv, "label"), base_r["label"])
print(f"reach(label, c=850) identical to baseline, per-pert : {same}")
mx = max(float(np.max(np.abs(T.reach(pre, [float(c)] * P, 'label') - base_r['label'])))
         for c in (100, 300, 850, 2000))
print(f"max per-pert |delta reach| over c in {{100,300,850,2000}} (label mode) : {mx:.3e}")
print(f"nmae(label) identical to baseline : "
      f"{abs(T.nmae(pre, mv, 'label') - base_nmae) < 1e-12}")
# how much of the reference pool do we even call?
frac_called = np.array([((s['r_rank'] >= 0) & (s['r_rank'] < 850)).mean() for s in pre])
frac_any = np.array([(s['r_rank'] >= 0).mean() for s in pre])
print(f"\nfraction of the reference's DE pool that we call significant:")
print(f"  pre-truncation  : mean {frac_any.mean():.3f}  median {np.median(frac_any):.3f}")
print(f"  at c=850        : mean {frac_called.mean():.3f}  median {np.median(frac_called):.3f}")
