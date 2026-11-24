"""269 -- the target-only lookup, examined: (T1) KOLF2 as the one CLEAN held context (its response is a donor on no other
fold), (T2) which training folds carry the cancer-fold lookup's signal (the held line's own row sits in the training folds'
donor means -- a PLOLO-specific leakage path that does not exist for a challenge line), (T3) shuffled-target controls with
permutation p-values, (T4) the ridge / GBM router with KOLF2 held.  Arithmetic on build269's tables only.

    PYTHONDONTWRITEBYTECODE=1 uv run python lookup269.py
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
CANCER = ["K562", "HCT116", "HEK293T"]
ALL = CANCER + ["KOLF2"]
DONORS = ["K562", "HCT116", "HEK293T"]
B_BOOT = 20_000
N_PERM = 500
SEED = 20260905
t0 = time.time()
LINES: list[str] = []
OUT: dict = {}


def P(s=""):
    print(s, flush=True)
    LINES.append(s)


dfc = pd.read_csv(WORK / "rows_choose269.csv")
dfd = pd.read_csv(WORK / "rows_defl269.csv")
ctx_n = {c: int(dfd[dfd.context == c].t_idx.max() + 1) for c in ALL}
ctx_targets = {c: dfd[dfd.context == c].sort_values("t_idx").target.tolist() for c in ALL}


def candidates(c):
    n = ctx_n[c]
    cand = [[] for _ in range(n)]
    sub = dfc[dfc.context == c]
    for t, d, u in zip(sub.t_idx, sub.donor, sub.u):
        cand[t].append((f"choose {d}", float(u)))
    subd = dfd[dfd.context == c]
    for t, u in zip(subd.t_idx, subd.u_undefl):
        cand[t].append(("undeflate", float(u)))
    return cand


CAND = {c: candidates(c) for c in ALL}
LAB_ORACLE = {c: np.array([max([0.0] + [u for _, u in CAND[c][t]]) for t in range(ctx_n[c])]) for c in ALL}
ORACLE263 = {c: dfd[dfd.context == c].sort_values("t_idx").oracle263.to_numpy() for c in ALL}


def boot(r, rng):
    idx = rng.integers(0, r.size, size=(B_BOOT, r.size))
    m = r[idx].mean(1)
    return dict(point=float(r.mean()), lo=float(np.percentile(m, 2.5)), hi=float(np.percentile(m, 97.5)), sd=float(m.std()), p_le0=float((m <= 0).mean()))


def boot_strat(rs, rng):
    m = np.zeros(B_BOOT)
    for r in rs:
        idx = rng.integers(0, r.size, size=(B_BOOT, r.size))
        m += r[idx].mean(1) / len(rs)
    return dict(point=float(np.mean([r.mean() for r in rs])), lo=float(np.percentile(m, 2.5)), hi=float(np.percentile(m, 97.5)), sd=float(m.std()), p_le0=float((m <= 0).mean()))


def lookup_scores(held, train_ctx, include_undeflate=True):
    """Per target of the held context: the training-derived score of every candidate action (mean label over the training
    contexts where the target has that action)."""
    tr_c = dfc[dfc.context.isin(train_ctx)].groupby(["target", "donor"]).u.mean()
    tr_d = dfd[dfd.context.isin(train_ctx)].groupby("target").u_undefl.mean()
    sc = [[] for _ in range(ctx_n[held])]
    for t in range(ctx_n[held]):
        tg = ctx_targets[held][t]
        for a, u in CAND[held][t]:
            if a == "undeflate":
                if include_undeflate and tg in tr_d.index:
                    sc[t].append((float(tr_d[tg]), a, u))
            else:
                key = (tg, a.split()[1])
                if key in tr_c.index:
                    sc[t].append((float(tr_c[key]), a, u))
    return sc


def decide(sc, tau):
    n = len(sc)
    r = np.zeros(n); a = np.zeros(n, bool); act = ["keep"] * n
    for t in range(n):
        if sc[t]:
            p, aa, u = max(sc[t], key=lambda x: x[0])
            if p > tau:
                r[t] = u; a[t] = True; act[t] = aa
    return r, a, act


def permute_targets(sc, rng):
    """the shuffled-TARGET control: the lookup vectors re-assigned to other targets of the held context (the scores of
    target t' are applied to target t's candidate labels, matched by action name)"""
    n = len(sc)
    perm = rng.permutation(n)
    out = [[] for _ in range(n)]
    for t in range(n):
        src = {aa: p for p, aa, _ in sc[perm[t]]}
        for a, u in CAND[list(ctx_n.keys())[0]][0][:0]:
            pass
        for p, aa, u in sc[t]:
            if aa in src:
                out[t].append((src[aa], aa, u))
    return out


rng = np.random.default_rng(SEED)
res = {}
P("== T1: KOLF2 held (clean: KOLF2's response is a donor on no training fold); training = the three cancer folds ==")
P(f"  KOLF2 ceilings: labelled oracle {LAB_ORACLE['KOLF2'].mean():+.4f}, 263 realised oracle {ORACLE263['KOLF2'].mean():+.4f}; fixed: " +
  ", ".join(f"always choose {d} {np.mean([u for t in range(ctx_n['KOLF2']) for a, u in CAND['KOLF2'][t] if a == f'choose {d}'] + [0.0] * 0):+.4f} (over the carried targets)" for d in DONORS))
fixed_kolf = {}
for d in DONORS + ["undeflate"]:
    key = f"choose {d}" if d != "undeflate" else d
    r = np.zeros(ctx_n["KOLF2"])
    for t in range(ctx_n["KOLF2"]):
        for a, u in CAND["KOLF2"][t]:
            if a == key:
                r[t] = u
    fixed_kolf[key] = boot(r, rng)
P("  fixed actions on KOLF2 (fold mean): " + ", ".join(f"{k} {v['point']:+.4f} [{v['lo']:+.4f}, {v['hi']:+.4f}]" for k, v in fixed_kolf.items()))
res["T1"] = dict(fixed=fixed_kolf, lab_oracle=float(LAB_ORACLE["KOLF2"].mean()), oracle263=float(ORACLE263["KOLF2"].mean()), rows={})
for incl in (True, False):
    sc = lookup_scores("KOLF2", CANCER, include_undeflate=incl)
    for tau in (0.0, 0.01, 0.02, 0.05):
        r, a, act = decide(sc, tau)
        b = boot(r, rng)
        # permutation null: re-assign the lookup vectors across targets
        null = np.array([decide(permute_targets(sc, rng), tau)[0].mean() for _ in range(N_PERM)])
        p_perm = float((1 + (null >= r.mean()).sum()) / (N_PERM + 1))
        acts = pd.Series(act)[a].value_counts().to_dict() if a.any() else {}
        res["T1"]["rows"][f"{'all actions' if incl else 'choose only'} tau {tau}"] = dict(**b, acted=int(a.sum()), acted_frac=float(a.mean()), precision=float((r[a] > 0).mean()) if a.any() else np.nan,
                                                                                         captured_lab=float(r.mean() / LAB_ORACLE["KOLF2"].mean()), null_mean=float(null.mean()), null_sd=float(null.std()), p_perm=p_perm, actions=acts)
        P(f"  lookup {'all actions' if incl else 'choose only ':12s} τ {tau:.2f}: Δpds {b['point']:+.4f} [{b['lo']:+.4f}, {b['hi']:+.4f}] P(<=0) {b['p_le0']:.3f}; acted {int(a.sum())} ({a.mean():.2f}), precision {(r[a] > 0).mean() if a.any() else float('nan'):.2f}; "
          f"captured {r.mean() / LAB_ORACLE['KOLF2'].mean() * 100:+.1f} % of the labelled oracle; shuffled-target null {null.mean():+.4f} ± {null.std():.4f}, p {p_perm:.3f}; actions {acts}")

P("\n== T2: cancer folds held; the lookup by TRAINING SOURCE (all three / the other two cancer folds / KOLF2 only), tau 0.02 ==")
res["T2"] = {}
for src_name, src_fn in [("all three other contexts", lambda c: [x for x in ALL if x != c]), ("the other two cancer folds", lambda c: [x for x in CANCER if x != c]),
                         ("KOLF2 only", lambda c: ["KOLF2"])]:
    for incl in (True, False):
        rs, per, nulls = [], {}, []
        for c in CANCER:
            sc = lookup_scores(c, src_fn(c), include_undeflate=incl)
            r, a, act = decide(sc, 0.02)
            rs.append(r)
            per[c] = dict(**boot(r, rng), acted_frac=float(a.mean()), precision=float((r[a] > 0).mean()) if a.any() else np.nan,
                          actions=pd.Series(act)[a].value_counts().to_dict() if a.any() else {})
            nulls.append(np.array([decide(permute_targets(sc, rng), 0.02)[0].mean() for _ in range(N_PERM)]))
        b = boot_strat(rs, rng)
        null3 = np.mean(nulls, 0)
        p_perm = float((1 + (null3 >= b["point"]).sum()) / (N_PERM + 1))
        key = f"{src_name} | {'all actions' if incl else 'choose only'}"
        res["T2"][key] = dict(mean3=b, per_ctx=per, null3_mean=float(null3.mean()), null3_sd=float(null3.std()), p_perm=p_perm,
                              captured_lab=float(b["point"] / np.mean([LAB_ORACLE[c].mean() for c in CANCER])))
        P(f"  {key:52s} 3-ctx {b['point']:+.4f} [{b['lo']:+.4f}, {b['hi']:+.4f}]  per ctx " + " / ".join(f"{per[c]['point']:+.4f}" for c in CANCER) +
          f"  acted {np.mean([per[c]['acted_frac'] for c in CANCER]):.2f}  captured {b['point'] / np.mean([LAB_ORACLE[c].mean() for c in CANCER]) * 100:+.1f} %  shuffled-target null {null3.mean():+.4f} ± {null3.std():.4f}, p {p_perm:.3f}")

P("\n== T2b: the mechanism -- per target, does the training lookup score track the HELD line's own role in the training folds? ==")
# For held X and action 'choose Z', the training label on fold Y (donors X, Z) is positive when X's row disagrees with Y's truth.
# Measure: Spearman between the training score of 'choose Z' (from fold Y) and the held-fold label of 'choose Z' (X's truth),
# separately for the fold where X is a donor (Y) and for KOLF2 (X is one of three donors).
from scipy.stats import spearmanr
res["T2b"] = {}
for c in CANCER:
    for src_name, src in [("other cancer fold", [x for x in CANCER if x != c]), ("KOLF2", ["KOLF2"])]:
        sc = lookup_scores(c, src, include_undeflate=False)
        xs, ys = [], []
        for t in range(ctx_n[c]):
            for p, a, u in sc[t]:
                xs.append(p); ys.append(u)
        rho = spearmanr(xs, ys).correlation if len(xs) > 10 else np.nan
        sign_agree = float(np.mean((np.array(xs) > 0) == (np.array(ys) > 0))) if xs else np.nan
        res["T2b"][f"{c} <- {src_name}"] = dict(n=len(xs), spearman=float(rho), sign_agreement=sign_agree)
        P(f"  held {c:8s} training source {src_name:18s}: n {len(xs)} Spearman(training score, held label) {rho:+.3f}, sign agreement {sign_agree:.3f}")

P("\n== T4: the ridge / GBM router with KOLF2 held (features as model269; training = the three cancer contexts) ==")
npz = np.load(WORK / "table269.npz", allow_pickle=True)
dm_profile = npz["dm_profile"]
prof = np.where(np.isfinite(dm_profile), dm_profile, np.nanmean(dm_profile, 1, keepdims=True))
prof_c = prof - prof.mean(1, keepdims=True); prof_c /= np.maximum(np.linalg.norm(prof_c, axis=1, keepdims=True), 1e-12)
b_lines = npz["b_ctrl_lines"][:, npz["sh_axis"]]; line_names = npz["b_ctrl_line_names"].tolist()
for d in DONORS:
    dfc[f"d_is_{d}"] = (dfc.donor == d).astype(float)
dfc["s_dm_in_held_missing"] = dfc.s_dm_in_held.isna().astype(float); dfd["s_dm_in_held_missing"] = dfd.s_dm_in_held.isna().astype(float)
dfc["r_dm_in_donor_missing"] = dfc.r_dm_in_donor.isna().astype(float)
SH = ["s_dd_raw", "s_dd_wht", "s_n_whitened", "s_cd4_norm", "s_cd4_energy_removed", "s_union_norm", "s_base_norm", "s_resid_norm",
      "s_resid_cos_base", "s_resid_share", "s_cd4_cos_base", "s_cd4_deflated_cos_base", "s_pred_nn", "s_pred_mean", "s_pred_norm",
      "s_cd4_de_log", "s_target_basal_held", "s_dm_mean", "s_dm_sd", "s_dm_in_held", "s_dm_in_held_missing", "s_collectri_out",
      "s_collectri_in", "s_dorothea_in", "s_is_tf", "s_basal_cc_pc1", "s_basal_cc_pc2", "s_basal_cc_pc3"]
REL = ["r_norm_raw", "r_norm_wht", "r_cd4_cos", "r_cd4_cos_w", "r_cd4d_cos", "r_cd4d_cos_w", "r_selfcons", "r_selfcons_w", "r_basal_sim",
       "r_target_basal_donor", "r_dm_in_donor", "r_dm_in_donor_missing"] + [f"d_is_{d}" for d in DONORS]


def design(df, feats, train_mask, pcs_dm, pcs_b):
    X = df[feats].to_numpy(float)
    X = np.column_stack([X, pcs_dm[df.g_idx.to_numpy()], np.stack([pcs_b[c] for c in df.context])])
    med = np.nanmedian(X[train_mask], 0); X = np.where(np.isfinite(X), X, med[None, :])
    mu, sd = X[train_mask].mean(0), X[train_mask].std(0) + 1e-9
    return (X - mu) / sd


def fit_ridge(X, y, groups):
    best, best_mse = None, np.inf
    for a in [0.3, 1, 3, 10, 30, 100, 300, 1000]:
        mse = sum(float(((Ridge(alpha=a).fit(X[tr], y[tr]).predict(X[va]) - y[va]) ** 2).mean()) for tr, va in GroupKFold(5).split(X, y, groups))
        if mse < best_mse:
            best, best_mse = a, mse
    return Ridge(alpha=best).fit(X, y)


trc = dfc.context.isin(CANCER).to_numpy(); tec = (dfc.context == "KOLF2").to_numpy()
trd = dfd.context.isin(CANCER).to_numpy(); ted = (dfd.context == "KOLF2").to_numpy()
pcs_dm = PCA(n_components=5, random_state=0).fit(prof_c[np.unique(dfc.g_idx.to_numpy()[trc])]).transform(prof_c)
idx = [line_names.index(c) for c in CANCER]; mu_b = b_lines[idx].mean(0)
pb = PCA(n_components=2, random_state=0).fit(b_lines[idx] - mu_b).transform(b_lines - mu_b)
pcs_b = {c: pb[line_names.index(c)] for c in ALL}
Xc = design(dfc, SH + REL, trc, pcs_dm, pcs_b); Xd = design(dfd, SH, trd, pcs_dm, pcs_b)
yc = dfc.u.to_numpy(float); yd = dfd.u_undefl.to_numpy(float)
res["T4"] = {}
for model in ("ridge", "gbm"):
    if model == "ridge":
        pc = fit_ridge(Xc[trc], yc[trc], dfc.target.to_numpy()[trc]).predict(Xc[tec]); pdd = fit_ridge(Xd[trd], yd[trd], dfd.target.to_numpy()[trd]).predict(Xd[ted])
    else:
        g = dict(max_depth=3, max_iter=300, learning_rate=0.05, min_samples_leaf=20, l2_regularization=1.0, random_state=0)
        pc = HistGradientBoostingRegressor(**g).fit(Xc[trc], yc[trc]).predict(Xc[tec]); pdd = HistGradientBoostingRegressor(**g).fit(Xd[trd], yd[trd]).predict(Xd[ted])
    pred_c = {(t, d): p for t, d, p in zip(dfc.t_idx[tec], dfc.donor[tec], pc)}
    pred_d = {t: p for t, p in zip(dfd.t_idx[ted], pdd)}
    sc = [[] for _ in range(ctx_n["KOLF2"])]
    for t in range(ctx_n["KOLF2"]):
        for a, u in CAND["KOLF2"][t]:
            key = t if a == "undeflate" else (t, a.split()[1])
            p = pred_d.get(t) if a == "undeflate" else pred_c.get(key)
            if p is not None and np.isfinite(p):
                sc[t].append((float(p), a, u))
    q = spearmanr(pc, yc[tec]).correlation
    res["T4"][model] = dict(spearman_choose=float(q), rows={})
    for tau in (0.0, 0.02, 0.05, 0.1):
        r, a, act = decide(sc, tau)
        b = boot(r, rng)
        res["T4"][model]["rows"][f"tau {tau}"] = dict(**b, acted_frac=float(a.mean()), precision=float((r[a] > 0).mean()) if a.any() else np.nan, captured_lab=float(r.mean() / LAB_ORACLE["KOLF2"].mean()))
        P(f"  {model:6s} τ {tau:.2f}: Δpds {b['point']:+.4f} [{b['lo']:+.4f}, {b['hi']:+.4f}] acted {a.mean():.2f} precision {(r[a] > 0).mean() if a.any() else float('nan'):.2f} captured {r.mean() / LAB_ORACLE['KOLF2'].mean() * 100:+.1f} %  (Spearman(pred, label) on the choose rows {q:+.3f})")

OUT.update(res, seconds=time.time() - t0)
(WORK / "lookup269.json").write_text(json.dumps(OUT, indent=1, default=float))
(WORK / "lookup269.log").write_text("\n".join(LINES) + "\n")
P(f"done in {time.time() - t0:.0f} s")
