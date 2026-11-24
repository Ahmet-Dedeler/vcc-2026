"""516 step 4: supervised stacking of same-target sources into a per-(target, gene) sign with a calibrated P(up).

Features per (target, gene), identical on every surface (a source the surface may not use reads 0):
  s_i        standardised uncentred score clip(x / (mad_g tau_p), +-3) of source i (K562-GWPS, Orion HCT116, Orion
             HEK293T, VIPerturb K562, CD4 Rest / Stim8hr / Stim48hr), 0 where absent
  det_i      the source's detection of the gene: log10 control CPM (CD4: log10 baseMean + 1), centred per source
  rec        the recipient's own log10 control CPM (centred), block indicator (top 2,200 by recipient control CPM)
Models (no intercept anywhere: the sign must come from the sources, not from a learnt down-lean):
  eq         equal-weight Stouffer (no fit; the vote516 construction)
  lin        logistic on s_i                                   (7 weights)
  reg        logistic on s_i x {block, rest}                   (14)
  det        reg + s_i x det_i + s_i x rec                     (28)
  gbm        HistGradientBoosting on s_i, det_i, rec, block, n_sources (monotone +1 in every s_i), on +-symmetrised
             data (each cell duplicated with every s_i negated and the label flipped) so P(up | s) = 1 - P(up | -s)
Validation: H1 with 5-fold target-grouped CV x 3 repeats (out-of-fold P); transfer H1 -> each fold and folds -> H1
(leave-surface-out). Candidates are priced as hybrids: model sign where |2P - 1| >= thr on non-block genes (rest) or on
genes production does not call (H1), production / its fold analog elsewhere.
    PYTHONDONTWRITEBYTECODE=1 uv run python stack516.py
"""
from __future__ import annotations

import json

import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression

import fold516 as FD
from common516 import WORK, known_answers, load_h1, price, summary_line

SRCS = ["k562_gwps", "orion_hct116", "orion_hek293t", "viperturb_k562", "cd4_Rest", "cd4_Stim8hr", "cd4_Stim48hr"]
SHORT = ["K", "H", "E", "V", "R", "S8", "S48"]


def features(z, ctrl_cpm, use, cap=3.0):
    """Returns dict of (targets x genes) arrays: S (n_src x T x G), DET (n_src x T x G), rec (G,), blk (G,)."""
    T = z["targets"].shape[0]; Gn = z["genes"].shape[0]
    S = np.zeros((len(SRCS), T, Gn), np.float32); DET = np.zeros((len(SRCS), T, Gn), np.float32)
    for i, s in enumerate(SRCS):
        if s not in use:
            continue
        x = (z[f"{s}__log_fc"] if s.startswith("cd4_") else z[f"{s}__x"]).astype(float)
        v = np.clip(x / (z[f"{s}__mad"][None, :] * z[f"{s}__tau"][:, None]), -cap, cap)
        S[i] = np.nan_to_num(v)
        if s.startswith("cd4_"):
            det = np.log10(np.nan_to_num(z[f"{s}__baseMean"].astype(float)) + 1.0)
        else:
            det = np.tile(np.log10(np.nan_to_num(z[f"{s}__cpm"].astype(float)) + 0.1), (T, 1))
        det = det - np.nanmedian(det[np.isfinite(v)]) if np.isfinite(v).any() else det * 0
        DET[i] = np.where(np.isfinite(v), det, 0.0)
    rec = np.log10(np.asarray(ctrl_cpm, float) + 0.1); rec = rec - np.median(rec)
    blk = np.zeros(Gn, bool); blk[np.argsort(-np.asarray(ctrl_cpm))[:2200]] = True
    return dict(S=S, DET=DET, rec=rec, blk=blk)


def design(F, kind, cells):
    """Design matrix rows for flat cell indices `cells` (into T x G)."""
    t, g = np.unravel_index(cells, F["S"].shape[1:])
    s = F["S"][:, t, g].T                                     # n x 7
    if kind == "lin":
        return s
    b = F["blk"][g][:, None].astype(float)
    X = [s * b, s * (1 - b)]
    if kind == "det":
        X += [s * F["DET"][:, t, g].T, s * F["rec"][g][:, None]]
    return np.hstack(X)


def gbm_design(F, cells):
    t, g = np.unravel_index(cells, F["S"].shape[1:])
    s = F["S"][:, t, g].T
    det = F["DET"][:, t, g].T
    nsrc = (s != 0).sum(1, keepdims=True)
    return np.hstack([s, det, F["rec"][g][:, None], F["blk"][g][:, None].astype(float), nsrc])


def fit(kind, Xy_list):
    X = np.vstack([a for a, _ in Xy_list]); y = np.concatenate([b for _, b in Xy_list])
    if kind == "gbm":
        ns = len(SRCS)
        Xs = X.copy(); Xs[:, :ns] *= -1                        # symmetrise: negate scores, flip labels
        X2 = np.vstack([X, Xs]); y2 = np.concatenate([y, 1 - y])
        cst = [1] * ns + [0] * (X.shape[1] - ns)
        m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=200,
                                           l2_regularization=1.0, monotonic_cst=cst, random_state=516)
        m.fit(X2, y2)
        return m
    m = LogisticRegression(fit_intercept=False, C=1.0, max_iter=2000)
    m.fit(X, y)
    return m


def predict(kind, m, F, cells):
    X = gbm_design(F, cells) if kind == "gbm" else design(F, kind, cells)
    return m.predict_proba(X)[:, 1]


def full_P(kind, m, F):
    T, G = F["S"].shape[1:]
    out = np.empty(T * G)
    step = 200000
    for a in range(0, T * G, step):
        c = np.arange(a, min(a + step, T * G))
        out[c] = predict(kind, m, F, c)
    return out.reshape(T, G)


def eq_P(F):
    """Equal-weight Stouffer mapped to a pseudo-probability (for thresholds only; not calibrated)."""
    s = F["S"]; n = (s != 0).sum(0)
    Z = np.where(n > 0, s.sum(0) / np.sqrt(np.maximum(n, 1)), 0.0)
    return 1 / (1 + np.exp(-Z))


def training_cells(pool, truth, targets_mask):
    m = pool & targets_mask[:, None]
    cells = np.flatnonzero(m.ravel())
    y = (truth.ravel()[cells] > 0).astype(int)
    return cells, y


def calib_table(P, pool, truth, label):
    m = pool & np.isfinite(P)
    conf = np.abs(2 * P[m] - 1); right = (np.sign(P[m] - 0.5) == truth[m])
    edges = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.01]
    s = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = (conf >= lo) & (conf < hi)
        if k.sum():
            pred_acc = 0.5 + 0.5 * conf[k].mean()
            s.append(f"[{lo},{hi}) n{k.sum()} acc {right[k].mean():.3f} (pred {pred_acc:.3f})")
    print(f"   calibration {label}: " + " | ".join(s), flush=True)


def main():
    d = load_h1()
    got, ok, _ = known_answers(d)
    print(f"KNOWN ANSWER {'PASS' if ok else 'FAIL'}: fid {got['fid']!r} jac {got['jac']!r} nmae {got['nmae']!r} "
          f"reach {got['reach']!r}", flush=True)
    assert ok
    zh = np.load(WORK / "src-h1.npz", allow_pickle=True)
    ess = np.zeros(len(d["targets"]), bool)
    for bn in ("k562_essential", "rpe1", "nadig_jurkat"):
        ess |= zh[f"{bn}__has"]
    SUBS = {"all85": np.ones(len(d["targets"]), bool), "no_ess": ~ess}
    FH = features(zh, d["ctrl"], set(SRCS))
    # fold surfaces (cross-study sources only; the analog baseline base_mcd4+c4q80)
    folds = {}
    for name, spec in FD.FOLDS.items():
        Fo = FD.load_fold(name)
        donors = [x for x in spec["prod"] if x != spec["sib"]]
        use = set(donors + FD.CD4 + spec["extra"])
        A, blk = FD.analogs(Fo, spec, donors)
        folds[name] = dict(Fo=Fo, F=features(Fo["z"], Fo["ctrl"], use), base=A["base_mcd4+c4q80"],
                           base_mean=A["base_mean"])
    rows = []
    T = len(d["targets"])
    rng = np.random.default_rng(516)
    for kind in ("eq", "lin", "reg", "det", "gbm"):
        # ---- H1 out-of-fold P (5-fold target CV x 3 repeats)
        if kind == "eq":
            P = eq_P(FH)
        else:
            acc = np.zeros((T, len(d["genes"]))); cnt = np.zeros(T)
            for rep in range(3):
                perm = rng.permutation(T); parts = np.array_split(perm, 5)
                for te in parts:
                    trm = np.ones(T, bool); trm[te] = False
                    cells, y = training_cells(d["pool"], d["truth"], trm)
                    X = gbm_design(FH, cells) if kind == "gbm" else design(FH, kind, cells)
                    m = fit(kind, [(X, y)])
                    for k in te:
                        c = np.arange(k * len(d["genes"]), (k + 1) * len(d["genes"]))
                        acc[k] += predict(kind, m, FH, c)
                    cnt[te] += 1
            P = acc / cnt[:, None]
        calib_table(P, d["pool"], d["truth"], f"{kind} H1 out-of-fold")
        S = np.sign(P - 0.5)
        for thr in (0.0, 0.1, 0.2):
            conf = np.abs(2 * P - 1) >= thr
            for region, rm in (("rest", ~FH["blk"][None, :]), ("uncalled", ~d["calls"]), ("all", np.ones_like(d["calls"]))):
                C = np.where(rm & conf, S, 0.0)
                r = price(d, C, f"{kind} H1-OOF |2P-1|>={thr} {region}", subsets=SUBS, verbose=False)
                for sn in SUBS:
                    print(summary_line(r[sn]), flush=True)
                    rows.append(dict(model=kind, surface="h1", thr=thr, region=region, subset=sn, acc=r[sn]["acc"],
                                     members=r[sn]["members"], tmean=r[sn]["acc_tmean"]))
        # ---- transfer: fit on all H1 -> each fold; fit on the folds -> H1
        if kind != "eq":
            cells, y = training_cells(d["pool"], d["truth"], np.ones(T, bool))
            X = gbm_design(FH, cells) if kind == "gbm" else design(FH, kind, cells)
            mH = fit(kind, [(X, y)])
            if kind in ("lin", "reg", "det"):
                w = mH.coef_.ravel()
                print(f"   {kind} weights fit on H1: " + " ".join(f"{w_:+.3f}" for w_ in w), flush=True)
        for name, fo in folds.items():
            Fo, Ff = fo["Fo"], fo["F"]
            Pf = eq_P(Ff) if kind == "eq" else full_P(kind, mH, Ff)
            calib_table(Pf, Fo["pool"], Fo["truth"], f"{kind} H1->{name}")
            for thr in (0.0, 0.1, 0.2):
                C = np.where(~Ff["blk"][None, :] & (np.abs(2 * Pf - 1) >= thr), np.sign(Pf - 0.5), 0.0)
                FD.price_fold(Fo, fo["base"], C, f"{kind} H1->fold |2P-1|>={thr} rest vs mcd4+c4q80", rows)
                FD.price_fold(Fo, fo["base_mean"], C, f"{kind} H1->fold |2P-1|>={thr} rest vs mean", rows)
        if kind != "eq":
            # folds -> H1 (all three folds pooled), and leave-one-fold-out among folds + H1
            tr = []
            for name, fo in folds.items():
                cells, y = training_cells(fo["Fo"]["pool"], fo["Fo"]["truth"], np.ones(len(fo["Fo"]["T"]), bool))
                tr.append(((gbm_design(fo["F"], cells) if kind == "gbm" else design(fo["F"], kind, cells)), y))
            mF = fit(kind, tr)
            if kind in ("lin", "reg", "det"):
                print(f"   {kind} weights fit on the folds: " + " ".join(f"{w_:+.3f}" for w_ in mF.coef_.ravel()), flush=True)
            Pf = full_P(kind, mF, FH)
            calib_table(Pf, d["pool"], d["truth"], f"{kind} folds->H1")
            for thr in (0.0, 0.1):
                for region, rm in (("rest", ~FH["blk"][None, :]), ("uncalled", ~d["calls"])):
                    C = np.where(rm & (np.abs(2 * Pf - 1) >= thr), np.sign(Pf - 0.5), 0.0)
                    r = price(d, C, f"{kind} folds->H1 |2P-1|>={thr} {region}", subsets=SUBS, verbose=False)
                    for sn in SUBS:
                        print(summary_line(r[sn]), flush=True)
                        rows.append(dict(model=kind, surface="folds->h1", thr=thr, region=region, subset=sn,
                                         acc=r[sn]["acc"], members=r[sn]["members"], tmean=r[sn]["acc_tmean"]))
        np.save(WORK / f"stackP-h1-{kind}.npy", P.astype(np.float32))
    json.dump(rows, open(WORK / "stack516.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
