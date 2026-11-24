"""172 — the run: LODO model selection, then ONE h1-test evaluation per model class.

Order of operations is deliberate and is not to be rearranged:

  A  known-answer rows (rule 10) -- abort if the instrument has drifted
  B  LODO grid: every model class x hyperparameter x top-heavy weighting, trained on two
     donor lines and scored on the third.  THIS IS THE ONLY PLACE A CHOICE IS MADE.
  C  h1-test, each class refitted on all three folds at its LODO-selected setting,
     evaluated ONCE.  Every class is reported, so nothing is selected on the eval surface.
  D  permutation null: labels shuffled WITHIN (fold, target), everything else identical
  E  ablation: which feature groups carry the signal, on BOTH surfaces
  F  save the key

AGENTS.md §11.1 rule 9 applies to the fold `reach` numbers in section B: their POOL SIZE is a
proxy (68-101 genes/target against h1-test's 89 median n_real) so their MAGNITUDE is not
comparable to h1-test's.  They are used only to ORDER model classes, which is the one thing an
n_real-denominated statistic can still do when every row shares the same denominator.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

import numpy as np
from scipy.stats import rankdata as rd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
WORK = os.path.join(HERE, "work")

from feats172 import FEATURES                                          # noqa: E402
from fit172 import (PROD, Model, fold_front, fold_reach, h1_features,  # noqa: E402
                    load_folds, topheavy_weights)
from harness172 import H1Test, check                                   # noqa: E402

MAX_TRAIN_ROWS = 300_000

# class -> [(label, model kwargs, top-heavy weighting), ...]
GRID = {
    "M1 logistic (22 params)": [
        (f"C={c} w={w}", dict(kind="logistic", C=c), w)
        for c in (0.1, 1.0, 10.0) for w in ("none", "invrank", "top25")],
    "M2 step-additive (~400)": [
        (f"C={c} w={w}", dict(kind="step", C=c), w)
        for c in (0.01, 0.05, 0.25) for w in ("none", "invrank", "exp10")],
    "M4 gbt (~5x10^3)": [
        (f"it={i} w={w}", dict(kind="gbt", max_iter=i), w)
        for i in (100, 300) for w in ("none", "invrank")],
    "M5 mlp": [
        (f"h={h} w=none", dict(kind="mlp", hidden=h), "none")
        for h in ((16,), (32, 16))],
}
GROUPS = {
    "consensus magnitude": ["rank_cons", "log_cons_rel", "rank_geo", "log_geo_rel"],
    "cross-source agreement": ["concord", "rel_dis", "ratio_mm", "rank_agree", "n_don"],
    "donor-panel responsiveness": ["rank_z", "log_z_rel", "log_marg", "log_sd", "resp"],
    "expression": ["log_cpm_rec", "log_cpm_don"],
    "saturation": ["sat"],
    "target-level": ["t_cos", "t_lognorm", "t_concfrac", "t_ndon"],
}


def subsample(X, y, t, w, cap=MAX_TRAIN_ROWS, seed=0):
    """Keep whole targets, so the top-heavy weights keep their per-target meaning."""
    if X.shape[0] <= cap:
        return X, y, t, w
    rng = np.random.default_rng(seed)
    u = np.unique(t)
    keep = set(u[rng.permutation(u.size)[:max(1, int(u.size * cap / X.shape[0]))]].tolist())
    m = np.isin(t, list(keep))
    return X[m], y[m], t[m], (None if w is None else w[m])


def train(cfg, weighting, Xs, ys, ts, cols=None, seed=0):
    """Fit one config.  Top-heavy weights come from a PRIOR model's ordering, never labels."""
    if cols is not None:
        Xs = Xs[:, cols]
    w = None
    if weighting != "none":
        a, b, c, _ = subsample(Xs, ys, ts, None, seed=seed)
        base = Model(kind="logistic").fit(a, b)
        w = topheavy_weights(weighting, Xs, ts, base.predict(Xs))
    Xr, yr, tr, wr = subsample(Xs, ys, ts, w, seed=seed)
    return Model(**cfg).fit(Xr, yr, wr, tr)


def lodo_score(cfg, weighting, F, cols=None):
    rs, fs, npar = [], [], 0
    for held in PROD:
        tr = [h for h in PROD if h != held]
        Xs = np.concatenate([F[h]["X"] for h in tr])
        ys = np.concatenate([F[h]["y"] for h in tr])
        ts = np.concatenate([F[h]["t"] + 100000 * i for i, h in enumerate(tr)])
        m = train(cfg, weighting, Xs, ys, ts, cols=cols)
        npar = m.n_params
        Xe = F[held]["X"] if cols is None else F[held]["X"][:, cols]
        s = m.predict(Xe)
        rs.append(fold_reach(s, F[held]["t"], F[held]["y"]))
        fs.append(fold_front(s, F[held]["t"], F[held]["y"]))
    return rs, fs, npar


def prod_key(T):
    """debench's production order, written as an explicit key so `deciles` can read it."""
    q = np.nan_to_num(T.S.pred_padj, nan=1.0)
    return np.where(T.S.pred_sig, 1e12 - q, -q)


def main() -> None:
    T = H1Test()
    print("=" * 104)
    print("A. KNOWN-ANSWER ROWS (AGENTS.md §11.1 rule 10)")
    check(T)
    F = load_folds()
    ci = FEATURES.index("rank_cons")
    for h in PROD:
        print(f"  fold {h:<16} {F[h]['X'].shape[0]:>7} rows  "
              f"base sign accuracy {F[h]['y'].mean():.4f}  reach(|cons|) "
              f"{fold_reach(F[h]['X'][:, ci], F[h]['t'], F[h]['y']):.5f}")

    # ---------------------------------------------------------------- B. LODO grid
    print("\n" + "=" * 104)
    print("B. LODO GRID — train on two donor lines, score the third.  Magnitudes are NOT")
    print("   comparable to h1-test (rule 9: the pool denominator differs); the ORDER is.")
    hdr = (f"{'model / setting':<34}" + "".join(f"{h[:12]:>14}" for h in PROD) +
           f"{'mean':>10}{'front@5':>10}{'params':>9}")
    print(hdr)
    out = {"lodo": {}, "baselines": {}, "models": {}}
    for nm, col in (("|cons| (un-whitened)", ci),
                    ("geo x concord", FEATURES.index("rank_agree")),
                    ("expression (log cpm)", FEATURES.index("log_cpm_rec"))):
        r = [fold_reach(F[h]["X"][:, col], F[h]["t"], F[h]["y"]) for h in PROD]
        f5 = [fold_front(F[h]["X"][:, col], F[h]["t"], F[h]["y"]) for h in PROD]
        out["lodo"][nm] = dict(fold_reach=r, front5=f5, n_params=0)
        print(f"{nm:<34}" + "".join(f"{v:>14.5f}" for v in r) +
              f"{np.mean(r):>10.5f}{np.mean(f5):>10.4f}{'--':>9}")
    best_of_class = {}
    for cls, settings in GRID.items():
        rows = []
        for label, cfg, w in settings:
            rs, fs, npar = lodo_score(cfg, w, F)
            rows.append((label, cfg, w, rs, fs, npar))
            out["lodo"][f"{cls} | {label}"] = dict(fold_reach=rs, front5=fs, n_params=npar)
            print(f"{cls.split(' ')[0] + ' ' + label:<34}" + "".join(f"{v:>14.5f}" for v in rs) +
                  f"{np.mean(rs):>10.5f}{np.mean(fs):>10.4f}{npar:>9}")
        b = max(rows, key=lambda r: np.mean(r[3]))
        best_of_class[cls] = b
        print(f"  -> {cls}: LODO-selected `{b[0]}`, mean fold reach {np.mean(b[3]):.5f}\n")
    winner_cls = max(best_of_class, key=lambda c: np.mean(best_of_class[c][3]))
    wlabel, wcfg, wweight, wrs, _, _ = best_of_class[winner_cls]
    print(f"  LODO WINNER: {winner_cls} @ `{wlabel}`  (mean fold reach {np.mean(wrs):.5f})")
    out["winner"] = f"{winner_cls} | {wlabel}"

    # ---------------------------------------------------------------- C. h1-test
    print("\n" + "=" * 104)
    print("C. h1-test — each class refitted on ALL THREE folds at its LODO-selected setting,")
    print("   evaluated ONCE.  `topx` is pooled; the macro variant is in work/run172.json.")
    Xh, _ = h1_features(T)
    shape = T.H.shape
    Xall = np.concatenate([F[h]["X"] for h in PROD])
    yall = np.concatenate([F[h]["y"] for h in PROD])
    tall = np.concatenate([F[h]["t"] + 100000 * i for i, h in enumerate(PROD)])
    print(f"{'key':<34}{'reach':>10}{'d1':>7}{'d2':>7}{'d3':>7}"
          f"{'top1%':>8}{'top2%':>8}{'top5%':>8}{'top10%':>8}{'rho':>8}")

    def row(label, key, store):
        k = np.asarray(key, dtype=np.float64)
        d = T.deciles(k)
        tx = T.topx(k)
        v = k[T.pool]
        rho = float("nan") if v.std() == 0 else float(
            spearmanr(v, T.correct[T.pool].astype(float)).statistic)
        r = T.reach(k)
        print(f"{label:<34}{r:>10.5f}{d[0]:>7.3f}{d[1]:>7.3f}{d[2]:>7.3f}"
              f"{tx[0]:>8.3f}{tx[1]:>8.3f}{tx[2]:>8.3f}{tx[3]:>8.3f}{rho:>8.4f}")
        store[label] = dict(reach=r, deciles=d, topx=tx, topx_macro=T.topx_macro(k), rho=rho)
        return r

    hf = np.load(os.path.join(WORK, "h1test-feats.npz"))
    B = out["baselines"]
    allkeys = {}
    row("constant key (gene-name order)", np.zeros(shape), B)
    row("production p-values", prod_key(T), B)
    allkeys["production p-values"] = prod_key(T)
    row("|cons| (un-whitened consensus)", np.abs(hf["cons"]), B)
    allkeys["|cons| (un-whitened)"] = np.abs(hf["cons"])
    row("geo x concord (013/095)", hf["agree"], B)
    allkeys["geo x concord"] = hf["agree"]
    row("|H| (whitened, THE BENCHMARK)", np.abs(T.H), B)
    allkeys["|H| (whitened)"] = np.abs(T.H)
    keys = {}
    for cls in GRID:
        label, cfg, w, _, _, _ = best_of_class[cls]
        m = train(cfg, w, Xall, yall, tall)
        c = m.predict(Xh).reshape(shape)
        keys[cls] = c
        allkeys[f"{cls.split(' (')[0]} @ {label}"] = c
        row(f"{cls} @ {label}", c, out["models"])
        out["models"][f"{cls} @ {label}"]["n_params"] = m.n_params
    best = keys[winner_cls]
    bl = np.zeros(shape)
    for p in range(shape[0]):
        bl[p] = 0.5 * rd(best[p]) + 0.5 * rd(np.abs(T.H[p]))
    row("50/50 rank blend with |H| (UNTUNED)", bl, out["models"])
    allkeys["50/50 blend with |H|"] = bl
    row("ORACLE (our own signs)", T.key_oracle(), B)
    allkeys["ORACLE"] = T.key_oracle()

    # ---------------------------------------------------------------- D. permutation null
    print("\n" + "=" * 104)
    print("D. PERMUTATION NULLS — two of them, because the briefed one LEAKS.")
    print("   N1 (briefed): labels shuffled WITHIN (fold, target).  This preserves each")
    print("      target's own base sign accuracy, and `rank_cons` is NOT uniform inside a")
    print("      target's DE pool (the pool is a top-K subset of the live genes), so the")
    print("      BETWEEN-target association survives at corr(y, rank_cons) = +0.037 against")
    print("      +0.124 unshuffled -- and a monotone fit converts that into a WITHIN-target")
    print("      ordering.  N1 therefore reads roughly `|cons|`, not the constant-key floor.")
    print("   N2 (strict): labels shuffled across ALL rows, destroying fold and target")
    print("      structure too.  Nothing survives, so the coefficient signs are coin flips")
    print("      and the seed spread is the honest measure of how much `reach` a key can")
    print("      collect from noise alone.")
    nulls = {"N1 within-target": [], "N2 global": []}
    nullfold = []
    for seed in range(6):
        rng = np.random.default_rng(seed)
        y1 = yall.copy()
        for p in np.unique(tall):
            m = tall == p
            y1[m] = rng.permutation(yall[m])
        y2 = rng.permutation(yall)
        for nm, yp in (("N1 within-target", y1), ("N2 global", y2)):
            mm = train(wcfg, wweight, Xall, yp, tall, seed=seed)
            c = mm.predict(Xh).reshape(shape)
            nulls[nm].append(T.reach(c))
            if nm == "N1 within-target":
                allkeys[f"null N1 seed{seed}"] = c
        print(f"   seed {seed}:  N1 {nulls['N1 within-target'][-1]:.5f}   "
              f"N2 {nulls['N2 global'][-1]:.5f}")
    for nm, v in nulls.items():
        print(f"   {nm:<18} mean {np.mean(v):.5f} +/- {np.std(v):.5f}   "
              f"(min {np.min(v):.5f}, max {np.max(v):.5f})")
    print(f"   for reference: fitted {out['models'][f'{winner_cls} @ {wlabel}']['reach']:.5f}, "
          f"|cons| {B['|cons| (un-whitened consensus)']['reach']:.5f}, constant-key floor "
          f"{B['constant key (gene-name order)']['reach']:.5f}")
    out["null"] = {k: [float(x) for x in v] for k, v in nulls.items()}
    out["null"]["corr_y_rankcons"] = dict(
        real=float(np.corrcoef(yall, Xall[:, ci])[0, 1]),
        shuffled_within_target=float(np.corrcoef(y1, Xall[:, ci])[0, 1]),
        shuffled_global=float(np.corrcoef(y2, Xall[:, ci])[0, 1]))

    # ---------------------------------------------------------------- E. ablation
    print("\n" + "=" * 104)
    print(f"E. ABLATION on the LODO winner ({winner_cls} @ {wlabel}) — drop one feature group")
    full_h1 = out["models"][f"{winner_cls} @ {wlabel}"]["reach"]
    full_fold = np.mean(wrs)
    print(f"{'dropped group':<30}{'h1 reach':>10}{'delta':>10}{'mean fold':>12}{'delta':>10}"
          f"{'h1 top1%':>10}")
    abl = {}
    for gname, feats in [("(nothing: full model)", [])] + list(GROUPS.items()):
        cols = [i for i, f in enumerate(FEATURES) if f not in feats]
        m = train(wcfg, wweight, Xall, yall, tall, cols=cols)
        c = m.predict(Xh[:, cols]).reshape(shape)
        r = T.reach(c)
        fr = np.mean(lodo_score(wcfg, wweight, F, cols=cols)[0])
        t1 = T.topx(c)[0]
        abl[gname] = dict(h1=r, fold=float(fr), top1=float(t1))
        print(f"{gname:<30}{r:>10.5f}{r-full_h1:>+10.5f}{fr:>12.5f}{fr-full_fold:>+10.5f}"
              f"{t1:>10.3f}")
    out["ablation"] = abl

    # ---------------------------------------------------------------- F. save
    path = os.path.join(HERE, "key172.npz")
    np.savez_compressed(path, key=best.astype(np.float32), genes=T.genes, targets=T.perts,
                        model=np.asarray(f"{winner_cls} @ {wlabel}"),
                        features=np.asarray(FEATURES), reach=np.float64(full_h1),
                        blend_with_absH=bl.astype(np.float32))
    sha = hashlib.sha256(open(path, "rb").read()).hexdigest()
    out["key"] = dict(path="experiments/172-sign-confidence-key/key172.npz", sha256=sha,
                      model=f"{winner_cls} @ {wlabel}", reach=full_h1)
    np.savez_compressed(os.path.join(WORK, "allkeys.npz"), genes=T.genes, targets=T.perts,
                        **{k: v.astype(np.float32) for k, v in allkeys.items()})
    print("\n" + "=" * 104)
    print(f"F. wrote {path}\n   sha256 {sha}")
    print(f"   wrote work/allkeys.npz ({len(allkeys)} keys, for err172.py / power172.py)")
    with open(os.path.join(WORK, "run172.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)


if __name__ == "__main__":
    main()
