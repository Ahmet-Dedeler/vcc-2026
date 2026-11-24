#!/usr/bin/env python3
"""377 score stage: subspaces, corrections, endpoints, paired bootstrap."""
import os
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys, json, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
ROOT = os.path.dirname(os.path.dirname(HERE))
KS = [10, 20, 50, 100]
GAMMAS = [0.25, 0.5]
NBOOT = 10000
RNG = np.random.default_rng(3770909)

def endpoints(pred, truth):
    """wsign (|truth|-weighted), usign, cosine, retrieval (334 formula)."""
    w = np.abs(truth)
    corr = (np.sign(pred) == np.sign(truth))
    qw = (w * corr).sum(1) / w.sum(1)
    qu = corr.mean(1)
    pn = pred / np.maximum(np.linalg.norm(pred, axis=1, keepdims=True), 1e-12)
    yn = truth / np.maximum(np.linalg.norm(truth, axis=1, keepdims=True), 1e-12)
    cos = (pn * yn).sum(1)
    s = pn @ yn.T
    own = np.diag(s)
    beaten = (s < own[:, None]).sum(1)
    ties = (s == own[:, None]).sum(1) - 1
    ret = (beaten + 0.5 * ties) / (len(truth) - 1)
    return qw, qu, cos, ret

def boot_delta(a, b, n=NBOOT, seed=1):
    """Paired bootstrap of mean(a-b); a,b per-target arrays."""
    rng = np.random.default_rng(seed)
    d = a - b
    T = len(d)
    means = np.empty(n)
    for s in range(0, n, 500):
        e = min(s + 500, n)
        ix = rng.integers(T, size=(e - s, T))
        means[s:e] = d[ix].mean(1)
    return float(d.mean()), [float(np.quantile(means, .025)), float(np.quantile(means, .975))]

def subspace_U(X, k_max=100):
    mu = X.mean(0)
    Xc = X - mu
    C = (Xc.T @ Xc) / len(Xc)
    ev, U = np.linalg.eigh(C)
    ix = np.argsort(ev)[::-1]
    ev, U = ev[ix], U[:, ix]
    varexp = ev / ev.sum()
    return mu, U, varexp

def correct(P, mu, U, k, mode, gamma=0.0):
    Uk = U[:, :k]
    Pc = P - mu[None, :]
    proj = (Pc @ Uk) @ Uk.T
    if mode == "proj":
        return proj
    if mode == "shrink":
        return P + gamma * (proj - Pc)
    if mode == "signonly":
        out = P.copy()
        med = np.median(np.abs(proj), axis=1, keepdims=True)
        use = np.abs(proj) > med
        out[use] = np.sign(proj[use]) * np.abs(out[use])
        # avoid exact zeros in output sign ties: keep magnitude, sign from proj
        return out
    raise ValueError(mode)

def main():
    t0 = time.time()
    base = np.load(os.path.join(WORK, "base377.npz"))
    l1000 = {}
    for c in ["HEPG2", "MCF7", "JURKAT"]:
        d = np.load(os.path.join(WORK, f"l1000_{c}.npz"))
        l1000[c] = (d["X"].astype(np.float64), [str(x) for x in d["genes"]],
                     np.array([str(x) for x in d["shrna_target"]]))
    line2cell = {"hepg2": "HEPG2", "jurkat": "JURKAT"}
    report = {"ks": KS, "gammas": GAMMAS, "lines": {}}
    for L in ["hepg2", "jurkat"]:
        cell = line2cell[L]
        axis = [str(x) for x in base[f"{L}_axis"]]
        targets = [str(x) for x in base[f"{L}_targets"]]
        Pmean = base[f"{L}_Pmean"].astype(np.float64)
        Pdiag = base[f"{L}_Pdiag"].astype(np.float64)
        Y = base[f"{L}_Y"].astype(np.float64)
        # restrict L1000 to axis
        Xc_full, syms, tgt = l1000[cell]
        pos = np.array([syms.index(g) for g in axis])
        S_same = Xc_full[:, pos]
        Xm, sm, tm = l1000["MCF7"]
        posm = np.array([sm.index(g) for g in axis])
        S_wrong = Xm[:, posm]
        # n-matched wrong-line subsample
        ixm = RNG.choice(len(S_wrong), size=min(len(S_same), len(S_wrong)), replace=False)
        S_wrong_nm = S_wrong[np.sort(ixm)]
        # held-target subspace (exclude shRNA rows hitting evaluated CRISPRi targets)
        tset = set(targets)
        keep = np.array([g not in tset for g in tgt])
        S_held = S_same[keep]
        print(f"{L}: S_same={S_same.shape} S_wrong={S_wrong.shape} S_held={S_held.shape} "
              f"(kept {keep.mean():.3f})", flush=True)
        subspaces = {}
        mu, U, ve = subspace_U(S_same)
        subspaces["same"] = (mu, U)
        report.setdefault("varexp", {})[L] = {
            "same": {str(k): float(ve[:k].sum()) for k in KS},
            "n_same": int(len(S_same))}
        muw, Uw, vew = subspace_U(S_wrong)
        subspaces["wrong"] = (muw, Uw)
        report["varexp"][L]["wrong"] = {str(k): float(vew[:k].sum()) for k in KS}
        muwn, Uwn, _ = subspace_U(S_wrong_nm)
        subspaces["wrong_nm"] = (muwn, Uwn)
        muh, Uh, _ = subspace_U(S_held)
        subspaces["held"] = (muh, Uh)
        # gene-shuffled same-line
        Ush = {}
        for s in range(3):
            r = np.random.default_rng(1000 + s)
            Xs = S_same[:, r.permutation(S_same.shape[1])]
            _, U2, _ = subspace_U(Xs)
            Ush[f"shuffled_{s}"] = U2
        # random subspaces (share the same-line mean for centering)
        Urnd = {}
        for s in range(3):
            r = np.random.default_rng(2000 + s)
            Q, _ = np.linalg.qr(r.standard_normal((len(axis), max(KS))))
            Urnd[f"random_{s}"] = Q
        # base endpoints (with retrieval)
        res = {}
        for bname, P in [("mean", Pmean), ("diag", Pdiag)]:
            qw, qu, cos, ret = endpoints(P, Y)
            res[bname + "_uncorrected"] = {
                "wsign": float(qw.mean()), "usign": float(qu.mean()),
                "cos": float(cos.mean()), "retrieval": float(ret.mean())}
        lr = {"base": res, "arms": {}, "varexp_same": report["varexp"][L]["same"]}
        # sweep corrections on Pmean (full), Pdiag (endpoints w/o retrieval)
        for bname, P in [("mean", Pmean), ("diag", Pdiag)]:
            qw0, qu0, cos0, ret0 = endpoints(P, Y)
            for sname, (muS, US) in subspaces.items():
                for k in KS:
                    for mode, gamma in ([("proj", 0.0)] + [("shrink", g) for g in GAMMAS]
                                         + [("signonly", 0.0)]):
                        Pp = correct(P, muS, US, k, mode, gamma)
                        w = np.abs(Y)
                        corr = (np.sign(Pp) == np.sign(Y))
                        qw = (w * corr).sum(1) / w.sum(1)
                        qu = corr.mean(1)
                        pn = Pp / np.maximum(np.linalg.norm(Pp, axis=1, keepdims=True), 1e-12)
                        yn = Y / np.maximum(np.linalg.norm(Y, axis=1, keepdims=True), 1e-12)
                        cos = (pn * yn).sum(1)
                        key = f"{bname}|{sname}|{mode}|k{k}|g{gamma}"
                        lr["arms"][key] = {"wsign": float(qw.mean()),
                                           "usign": float(qu.mean()),
                                           "cos": float(cos.mean())}
            # shuffled + random use same-line centering mean
            for sname, UU in list(Ush.items()) + list(Urnd.items()):
                for k in KS:
                    for mode, gamma in ([("proj", 0.0)] + [("shrink", g) for g in GAMMAS]
                                         + [("signonly", 0.0)]):
                        Pp = correct(P, mu, UU, k, mode, gamma)
                        w = np.abs(Y)
                        corr = (np.sign(Pp) == np.sign(Y))
                        qw = (w * corr).sum(1) / w.sum(1)
                        lr["arms"][f"{bname}|{sname}|{mode}|k{k}|g{gamma}"] = {
                            "wsign": float(qw.mean()), "usign": float(corr.mean(1).mean()),
                            "cos": float(((Pp / np.maximum(np.linalg.norm(Pp, axis=1, keepdims=True), 1e-12)) * (Y / np.maximum(np.linalg.norm(Y, axis=1, keepdims=True), 1e-12))).sum(1).mean())}
        # retrieval for shortlist: uncorrected + best-wsign arm per (subspace, mode) on mean base
        qw0 = endpoints(Pmean, Y)[0]
        short = {}
        cands = {}
        for key, v in lr["arms"].items():
            b, s, m, kk, gg = key.split("|")
            if b != "mean":
                continue
            cands.setdefault((s, m), []).append((v["wsign"], key))
        for (s, m), lst in cands.items():
            lst.sort(reverse=True)
            short[f"{s}|{m}"] = lst[0][1]
        retr = {}
        for name, key in list(short.items()) + [("uncorrected", None)]:
            if key is None:
                Pp = Pmean
            else:
                b, s, m, kk, gg = key.split("|")
                k = int(kk[1:]); gamma = float(gg[1:])
                if s in subspaces:
                    muS, US = subspaces[s]
                elif s.startswith("shuffled"):
                    US = Ush[s]; muS = mu
                else:
                    US = Urnd[s]; muS = mu
                Pp = correct(Pmean, muS, US, k, m, gamma)
            _, _, _, ret = endpoints(Pp, Y)
            retr[name] = {"key": key, "retrieval": float(ret.mean()),
                          "dvector": (endpoints(Pp, Y)[0] - qw0)}
        # paired bootstraps: best same-line arm vs uncorrected (wsign), + retrieval delta
        # pick best same-line wsign arm on mean base
        same_keys = [k for k in lr["arms"] if k.startswith("mean|same|")]
        best_same = max(same_keys, key=lambda k: lr["arms"][k]["wsign"])
        b, s, m, kk, gg = best_same.split("|")
        Pp_best = correct(Pmean, mu, U, int(kk[1:]), m, float(gg[1:]))
        qw_best, _, _, ret_best = endpoints(Pp_best, Y)
        qw_u, _, _, ret_u = endpoints(Pmean, Y)
        md, ci = boot_delta(qw_best, qw_u, seed=11)
        mdr, cir = boot_delta(ret_best, ret_u, seed=12)
        # wrong-line best (same mode class for comparability + overall best)
        wrong_keys = [k for k in lr["arms"] if k.startswith("mean|wrong|")]
        best_wrong = max(wrong_keys, key=lambda k: lr["arms"][k]["wsign"])
        lr["best"] = {"same": best_same, "same_wsign": lr["arms"][best_same]["wsign"],
                      "wrong": best_wrong, "wrong_wsign": lr["arms"][best_wrong]["wsign"],
                      "uncorrected_wsign": float(qw_u.mean()),
                      "delta_wsign": md, "delta_wsign_ci95": ci,
                      "delta_retrieval": mdr, "delta_retrieval_ci95": cir}
        lr["shortlist_retrieval"] = {k: {"key": v["key"], "retrieval": v["retrieval"]} for k, v in retr.items()}
        # store per-target vectors needed for pooled bootstrap
        np.savez_compressed(os.path.join(WORK, f"vec_{L}.npz"),
                            qw_best=qw_best, qw_u=qw_u, ret_best=ret_best, ret_u=ret_u)
        report["lines"][L] = lr
        print(f"{L}: uncorr wsign={qw_u.mean():.4f} best_same={best_same} "
              f"{qw_best.mean():.4f} d={md:+.4f} {ci} ret_d={mdr:+.4f} {cir}", flush=True)
    # pooled inverse-variance over lines
    ds, vs = [], []
    for L in ["hepg2", "jurkat"]:
        b = report["lines"][L]["best"]
        lo, hi = b["delta_wsign_ci95"]
        se = (hi - lo) / (2 * 1.96)
        ds.append(b["delta_wsign"]); vs.append(se ** 2)
    w = [1 / v for v in vs]
    pooled = sum(a * c for a, c in zip(ds, w)) / sum(w)
    se_p = sum(w) ** -0.5
    report["pooled"] = {"delta_wsign": pooled,
                        "ci95": [pooled - 1.96 * se_p, pooled + 1.96 * se_p],
                        "line_deltas": ds}
    json.dump(report, open(os.path.join(WORK, "scores377.json"), "w"), indent=1)
    print("pooled:", pooled, [pooled - 1.96 * se_p, pooled + 1.96 * se_p], flush=True)
    print(f"done {time.time()-t0:.1f}s", flush=True)

main()
