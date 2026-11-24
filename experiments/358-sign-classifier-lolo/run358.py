"""358 sign-classifier LOLO. Actions: cache | run | report. Resumable per block."""
import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
SEED = 3580909
BANKS = ("data/donors/k562_gwps", "data/donors/orion_hct116",
         "data/donors/orion_hek293t", "data/donors/rpe1",
         "data/donors/nadig_jurkat",
         "experiments/166-donor-acquisition/donors/nadig_hepg2")
LINES = ("k562_gwps", "orion_hct116", "orion_hek293t", "rpe1",
         "nadig_jurkat", "nadig_hepg2")
FEATS = ("cons", "cons_sign", "abs_cons", "agree_frac", "basal_g", "basal_t",
         "miss_t", "donor_mean_g", "diff_g", "is_tf_t", "edge_tg_any",
         "edge_tg_sign", "shared_log", "edge_gt_any", "edge_gt_sign")
DONOR_ONLY = [0, 1, 2, 3]
DONOR_PRIOR = [0, 1, 2, 3, 9, 10, 11, 12, 13, 14]
ARMS = ("diag", "lr_full", "hgb_full", "lr_shufrecip", "hgb_shufrecip",
        "lr_donoronly", "lr_donorprior", "lr_shuflabel")


def order(names, salt):
    return sorted(range(len(names)),
                  key=lambda i: hashlib.sha256(
                      (salt + str(names[i])).encode()).hexdigest())


def derangement(n, seed):
    rng = np.random.default_rng(seed)
    cycle = rng.permutation(n)
    out = np.empty(n, dtype=np.int64)
    out[cycle] = np.roll(cycle, 1)
    assert np.all(out != np.arange(n))
    return out


def load_direction(targets, genes):
    mats = []
    for b in BANKS:
        p = ROOT / b
        pi = {t: i for i, t in
              enumerate(p.joinpath("perts.txt").read_text().splitlines())}
        gi = {g: i for i, g in
              enumerate(p.joinpath("genes.txt").read_text().splitlines())}
        ti, ix = [pi[t] for t in targets], [gi[g] for g in genes]
        raw = np.load(p / "delta.npy", mmap_mode="r")
        mat = np.array(raw[np.ix_(ti, ix)], dtype=np.float64)
        basg = np.array(np.load(p / "b_ctrl.npy")[ix], dtype=np.float64)
        # target basal (NaN where the target symbol is absent from genes.txt)
        bfull = np.load(p / "b_ctrl.npy")
        hit = np.array([gi[t] if t in gi else -1 for t in targets])
        ok = hit >= 0
        bt = np.full(len(targets), np.nan)
        bt[ok] = bfull[hit[ok]]
        miss = ~ok
        bt[miss] = float(np.median(bt[ok]))
        mats.append((mat, basg, bt, miss))
    direction = np.array([m[0] for m in mats])
    rms = np.sqrt(np.mean(direction * direction, axis=-1))
    direction = direction / np.maximum(rms[..., None], 1e-12)
    return (direction,
            np.array([m[1] for m in mats]),
            np.array([m[2] for m in mats]),
            np.array([m[3] for m in mats]))


def build_prior(targets, genes):
    d = pd.read_csv(ROOT / "data/external/omnipath-dorothea/dorothea-abc.tsv",
                    sep="\t", usecols=["source_genesymbol", "target_genesymbol",
                                        "is_stimulation", "is_inhibition"])
    tfs = sorted(set(d.source_genesymbol))
    tf_ix = {t: i for i, t in enumerate(tfs)}
    syms = list(dict.fromkeys(list(targets) + list(genes)))
    s_ix = {s: i for i, s in enumerate(syms)}
    t_ix = {t: i for i, t in enumerate(targets)}
    g_ix = {g: i for i, g in enumerate(genes)}
    edge = {}
    for tf, tg, st, inh in zip(d.source_genesymbol, d.target_genesymbol,
                               d.is_stimulation, d.is_inhibition):
        cur = edge.get((tf, tg), [False, False])
        edge[(tf, tg)] = [cur[0] or bool(st), cur[1] or bool(inh)]
    nt, ng = len(targets), len(genes)
    edge_tg_any = np.zeros((nt, ng), np.int8)
    edge_tg_sign = np.zeros((nt, ng), np.int8)
    edge_gt_any = np.zeros((nt, ng), np.int8)
    edge_gt_sign = np.zeros((nt, ng), np.int8)
    for (a, b), (st, inh) in edge.items():
        s = (1 if st and not inh else (-1 if inh and not st else 0))
        if a in t_ix and b in g_ix:
            edge_tg_any[t_ix[a], g_ix[b]] = 1
            edge_tg_sign[t_ix[a], g_ix[b]] = np.int8(s)
        if a in g_ix and b in t_ix:
            edge_gt_any[t_ix[b], g_ix[a]] = 1
            edge_gt_sign[t_ix[b], g_ix[a]] = np.int8(s)
    A = np.zeros((len(syms), len(tfs)), bool)
    for tf, tg in zip(d.source_genesymbol, d.target_genesymbol):
        if tg in s_ix and tf in tf_ix:
            A[s_ix[tg], tf_ix[tf]] = True
    At = A[[s_ix[t] for t in targets]].astype(np.float32)
    Ag = A[[s_ix[g] for g in genes]].astype(np.float32)
    shared = np.log1p(At @ Ag.T).astype(np.float32)
    return {"is_tf_t": np.array([t in set(tfs) for t in targets], np.float32),
            "edge_tg_any": edge_tg_any, "edge_tg_sign": edge_tg_sign,
            "edge_gt_any": edge_gt_any, "edge_gt_sign": edge_gt_sign,
            "shared_log": shared,
            "n_tf": len(tfs), "n_edges": len(edge)}


def action_cache():
    t0 = time.time()
    sel = json.loads((ROOT / "experiments/334-text-conditioned-transfer"
                      / "selection.json").read_text())
    targets, genes = sel["targets"], sel["genes"]
    direction, basal_g, basal_t, miss = load_direction(targets, genes)
    prior = build_prior(targets, genes)
    np.savez_compressed(WORK / "cache358.npz", direction=direction,
                        basal_g=basal_g, basal_t=basal_t, miss_t=miss,
                        targets=np.array(targets), genes=np.array(genes),
                        **{k: v for k, v in prior.items()
                           if k not in ("n_tf", "n_edges")})
    meta = {"targets": len(targets), "genes": len(genes),
            "miss_t_per_bank": {ln: int(miss[bi].sum())
                                for bi, ln in enumerate(LINES)},
            "n_tf": prior["n_tf"], "n_edges": prior["n_edges"],
            "edge_tg_any": int((prior["edge_tg_any"] != 0).sum()),
            "edge_tg_signed": int((prior["edge_tg_sign"] != 0).sum()),
            "edge_gt_any": int((prior["edge_gt_any"] != 0).sum()),
            "edge_gt_signed": int((prior["edge_gt_sign"] != 0).sum()),
            "shared_nonzero": int((prior["shared_log"] > 0).sum()),
            "label_frac_pos": float((direction > 0).mean()),
            "seconds": time.time() - t0}
    (WORK / "cache-meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))


def block_features(direction, basal_g, basal_t, miss_t, pr, R_idx, donors,
                   tids, gperm=None):
    """Feature matrix (len(tids), 512, 13) for recipient rows R_idx."""
    cons = direction[donors][:, tids].mean(0)
    sgn = np.sign(cons)
    agree = np.mean(np.sign(direction[donors][:, tids]) ==
                    np.where(sgn == 0, 1, sgn)[None], axis=0)
    ng = direction.shape[2]
    bg = basal_g[R_idx]
    if gperm is not None:
        bg = bg[gperm]
    dmg = basal_g[donors].mean(0)
    F = np.empty((len(tids), ng, len(FEATS)), np.float32)
    F[..., 0] = cons
    F[..., 1] = np.where(sgn == 0, 1, sgn)
    F[..., 2] = np.abs(cons)
    F[..., 3] = agree
    F[..., 4] = np.broadcast_to(bg[None, :], (len(tids), ng))
    F[..., 5] = np.broadcast_to(basal_t[R_idx][tids][:, None],
                                (len(tids), ng))
    F[..., 6] = np.broadcast_to(
        miss_t[R_idx][tids][:, None].astype(np.float32),
        (len(tids), ng))
    F[..., 7] = np.broadcast_to(dmg[None, :], (len(tids), ng))
    diff = basal_g[R_idx] - dmg
    if gperm is not None:
        diff = diff[gperm]
    F[..., 8] = np.broadcast_to(diff[None, :], (len(tids), ng))
    F[..., 9] = np.broadcast_to(pr["is_tf_t"][tids][:, None],
                                (len(tids), ng))
    F[..., 10] = pr["edge_tg_any"][tids]
    F[..., 11] = pr["edge_tg_sign"][tids]
    F[..., 12] = pr["shared_log"][tids]
    F[..., 13] = pr["edge_gt_any"][tids]
    F[..., 14] = pr["edge_gt_sign"][tids]
    assert np.isfinite(F).all()
    return F, cons


def fit_predict(kind, Xtr, ytr, wtr, Xte, cols):
    if kind == "lr":
        from sklearn.preprocessing import StandardScaler
        from sklearn.linear_model import LogisticRegression
        sc = StandardScaler().fit(Xtr[:, cols])
        clf = LogisticRegression(max_iter=2000, C=1.0)
        clf.fit(sc.transform(Xtr[:, cols]), ytr, sample_weight=wtr)
        return clf.predict_proba(sc.transform(Xte[:, cols]))[:, 1], clf
    from sklearn.ensemble import HistGradientBoostingClassifier
    clf = HistGradientBoostingClassifier(
        loss="log_loss", max_iter=150, learning_rate=0.05,
        max_leaf_nodes=15, max_depth=4, min_samples_leaf=64,
        l2_regularization=10.0, early_stopping=False,
        random_state=SEED)
    clf.fit(Xtr[:, cols], ytr, sample_weight=wtr)
    return clf.predict_proba(Xte[:, cols])[:, 1], clf


def action_run():
    t0 = time.time()
    z = np.load(WORK / "cache358.npz")
    direction = z["direction"]
    basal_g, basal_t = z["basal_g"], z["basal_t"]
    targets = list(z["targets"])
    pr = {k: z[k] for k in ("is_tf_t", "edge_tg_any", "edge_tg_sign",
                             "shared_log", "edge_gt_any", "edge_gt_sign")}
    miss_t = z["miss_t"]
    perm = order(targets, "334-split|")
    blocks = [np.array(perm[::2]), np.array(perm[1::2])]
    ALL = list(range(len(FEATS)))
    for outer in range(6):
        donors = [i for i in range(6) if i != outer]
        for block in (0, 1):
            out = WORK / f"pred_outer{outer}_block{block}.npz"
            if out.exists():
                print(f"skip outer {outer} block {block} (exists)", flush=True)
                continue
            b0 = time.time()
            train, test = blocks[1 - block], blocks[block]
            gperm = derangement(512, SEED + 100 + outer * 2 + block)
            labperm = derangement(len(train), SEED + 200 + outer * 2 + block)
            # --- train assembly over 5 inner recipients ---
            Xs, ys, ws = [], [], []
            for D in donors:
                others = [i for i in donors if i != D]
                F, _ = block_features(direction, basal_g, basal_t, miss_t,
                                      pr, D, others, train)
                # control-A variant uses D's own basal under the same gperm
                Fa, _ = block_features(direction, basal_g, basal_t, miss_t,
                                       pr, D, others, train, gperm=gperm)
                y = (direction[D][train] > 0).astype(np.int8)
                w = np.abs(direction[D][train]).astype(np.float64)
                Xs.append((F, Fa))
                ys.append(y)
                ws.append(w)
            Ftr = np.concatenate([f for f, _ in Xs]).reshape(-1, len(FEATS))
            Ftra = np.concatenate([f for _, f in Xs]).reshape(-1, len(FEATS))
            ytr = np.concatenate([y.ravel() for y in ys])
            ysh = np.concatenate([y[labperm].ravel() for y in ys])
            wtr = np.concatenate([w.ravel() for w in ws])
            del Xs, ys, ws
            # --- test features (recipient = held line) ---
            Fte, cons5 = block_features(direction, basal_g, basal_t, miss_t,
                                        pr, outer, donors, test)
            Ftea, _ = block_features(direction, basal_g, basal_t, miss_t,
                                     pr, outer, donors, test, gperm=gperm)
            Xte = Fte.reshape(-1, len(FEATS))
            Xtea = Ftea.reshape(-1, len(FEATS))
            truth = direction[outer][test]
            # diagonal gain on train (KA-verified construction)
            num = np.zeros(512)
            den = np.zeros(512)
            for D in donors:
                others = [i for i in donors if i != D]
                x = direction[others][:, train].mean(0)
                y = direction[D][train]
                num += (x * y).sum(0)
                den += (x * x).sum(0)
            gain = (num + np.median(den)) / np.maximum(np.median(den) + den,
                                                      1e-12)
            saved = {"test_indices": test, "truth": truth, "cons5": cons5,
                     "gain": gain,
                     "diag": cons5 * gain}
            jobs = [("lr_full", "lr", Ftr, ALL, ytr),
                    ("hgb_full", "hgb", Ftr, ALL, ytr),
                    ("lr_shufrecip", "lr", Ftra, ALL, ytr),
                    ("hgb_shufrecip", "hgb", Ftra, ALL, ytr),
                    ("lr_donoronly", "lr", Ftr, DONOR_ONLY, ytr),
                    ("lr_donorprior", "lr", Ftr, DONOR_PRIOR, ytr),
                    ("lr_shuflabel", "lr", Ftr, ALL, ysh)]
            for name, kind, Xm, cols, yy in jobs:
                Xt = Xtea if name.endswith("shufrecip") else Xte
                p, _ = fit_predict(kind, Xm, yy, wtr, Xt, cols)
                saved[name] = p.reshape(len(test), 512).astype(np.float32)
                print(f"outer {outer} block {block} {name} done "
                      f"({time.time()-b0:.0f}s)", flush=True)
            np.savez_compressed(out, **saved)
            del Ftr, Ftra, Xte, Xtea, saved
    print("RUN seconds:", round(time.time() - t0, 1))


def wsign(pred_sign, truth):
    w = np.abs(truth)
    return ((w * (np.sign(pred_sign) == np.sign(truth))).sum(1) /
            w.sum(1))


def top64(pred_sign, truth):
    w = np.abs(truth)
    top = np.argsort(-w, axis=1, kind="stable")[:, :64]
    sw = np.take_along_axis(w, top, axis=1)
    sc = np.take_along_axis(np.sign(pred_sign) == np.sign(truth), top, axis=1)
    return (sw * sc).sum(1) / sw.sum(1)


def retrieval(pred, truth):
    pn = pred / np.maximum(np.linalg.norm(pred, axis=1, keepdims=True), 1e-12)
    tn = truth / np.maximum(np.linalg.norm(truth, axis=1, keepdims=True),
                            1e-12)
    s = pn @ tn.T
    own = np.diag(s).copy()
    return ((s < own[:, None]).sum(1) +
            0.5 * ((s == own[:, None]).sum(1) - 1)) / (len(truth) - 1)


def interval(delta, draws=10000):
    rng = np.random.default_rng(SEED)
    L, T = delta.shape
    tgt, cro = [], []
    for s in range(0, draws, 100):
        n = min(100, draws - s)
        ti = rng.integers(T, size=(n, T))
        li = rng.integers(L, size=(n, L))
        by = delta[:, ti].mean(-1).T
        tgt.extend(by.mean(1))
        cro.extend(by[np.arange(n)[:, None], li].mean(1))
    m = delta.mean(1)
    return {"mean": float(m.mean()),
            "target_ci95": np.quantile(tgt, [.025, .975]).tolist(),
            "crossed_ci95": np.quantile(cro, [.025, .975]).tolist(),
            "line_means": m.tolist(),
            "positive_lines": int((m > 0).sum())}


def action_report():
    z = np.load(WORK / "cache358.npz")
    targets = list(z["targets"])
    perm = order(targets, "334-split|")
    blocks = [np.array(perm[::2]), np.array(perm[1::2])]
    per = {a: {"w": [], "t64": [], "ret": []} for a in ARMS}
    for outer in range(6):
        asm = {}
        for block in (0, 1):
            d = np.load(WORK / f"pred_outer{outer}_block{block}.npz")
            for a in ARMS:
                if a == "diag":
                    pred = d["diag"]
                else:
                    p = d[a]
                    pred = np.where(p >= 0.5, 1.0, -1.0) * np.abs(d["cons5"])
                asm.setdefault(a, {})[block] = (pred, d["truth"])
        for a in ARMS:
            pred = np.concatenate([asm[a][0][0], asm[a][1][0]])
            # restore target order (block0=test-half-0 etc.); both halves cover all
            ti = np.concatenate([blocks[0], blocks[1]])
            back = np.empty_like(ti)
            back[ti] = np.arange(len(ti))
            pred, truth = pred[back], np.concatenate(
                [asm[a][0][1], asm[a][1][1]])[back]
            assert truth.shape == (1024, 512)
            per[a]["w"].append(wsign(np.sign(pred), truth))
            per[a]["t64"].append(top64(np.sign(pred), truth))
            per[a]["ret"].append(retrieval(pred, truth))
    for k in per:
        for m in per[k]:
            per[k][m] = np.array(per[k][m])
    means = {a: {m: float(per[a][m].mean()) for m in per[a]} for a in ARMS}
    print(f"{'arm':16s} {'w_sign':>10s} {'top64':>10s} {'retr':>10s}")
    for a in ARMS:
        print(f"{a:16s} {means[a]['w']:10.6f} {means[a]['t64']:10.6f} "
              f"{means[a]['ret']:10.6f}")
    comps = {}
    pairs = [("lr_full", "diag"), ("hgb_full", "diag"),
             ("lr_shufrecip", "lr_full"), ("hgb_shufrecip", "hgb_full"),
             ("lr_donoronly", "diag"), ("lr_donorprior", "lr_donoronly"),
             ("lr_shuflabel", "diag")]
    for prim, base in pairs:
        comps[f"{prim}_minus_{base}"] = {}
        for m in ("w", "t64", "ret"):
            comps[f"{prim}_minus_{base}"][m] = interval(per[prim][m] -
                                                       per[base][m])
    for k, v in comps.items():
        for m in ("w", "t64", "ret"):
            c = v[m]
            print(f"{k}/{m}: mean {c['mean']:+.6f} "
                  f"tgt95 [{c['target_ci95'][0]:+.6f},{c['target_ci95'][1]:+.6f}] "
                  f"cro95 [{c['crossed_ci95'][0]:+.6f},{c['crossed_ci95'][1]:+.6f}] "
                  f"pos {c['positive_lines']}/6")
    rec = {"means": means, "comparisons": comps,
           "line_names": list(LINES),
           "per_line": {a: {m: per[a][m].mean(1).tolist() for m in per[a]}
                        for a in ARMS}}
    (WORK / "results358.json").write_text(json.dumps(rec, indent=2) + "\n")
    rows = []
    for li, ln in enumerate(LINES):
        for ti, t in enumerate(targets):
            rows.append({"line": ln, "target": t, **{
                f"{a}__{m}": float(per[a][m][li, ti])
                for a in ARMS for m in ("w", "t64", "ret")}})
    import csv
    with open(WORK / "per-target358.csv", "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("wrote results358.json + per-target358.csv")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=("cache", "run", "report"))
    args = ap.parse_args()
    with threadpool_limits(limits=4):
        globals()[f"action_{args.action}"]()
