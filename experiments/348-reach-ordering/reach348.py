#!/usr/bin/env python3
"""348 — reach ordering on saved cap15 arrays. Read-only inputs from 323/work.

Modes: known (scorer + independent walk vs recorded csv), order (all
alternative orderings), profile (P(k) table + k* diagnostics for current).
Writes only into experiments/348-reach-ordering/work/.
"""
from __future__ import annotations
import argparse
import json
import os
import sys
import time
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
HERE = __import__("pathlib").Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/323-current-family-validation"))
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import polars as pl
import run323 as R
import evaluate323 as E
from cell_eval2.config import EvalConfig
from cell_eval2.de import prepare_de
from cell_eval2.run import dispatch_de_metrics

P0 = 0.9
ALPHA = 0.05
BASE348 = ROOT / "experiments/323-current-family-validation/work"
WORK = HERE / "work"
REACH = "de_wilcoxon_direction_reach_raw"


def load():
    d = {}
    z = np.load(BASE348 / "stats-cap15.npz", allow_pickle=False)
    d["targets"] = z["targets"].astype(str)
    d["pp"], d["qp"], d["lp"] = z["p"], z["q"], z["lfc"]
    z.close()
    r = np.load(BASE348 / "reference.npz", allow_pickle=False)
    d["genes"] = r["genes"].astype(str)
    d["keep"] = r["keep"]
    d["rgenes"] = r["genes"].astype(str)[r["keep"]]
    d["rtargets"] = r["targets"].astype(str)
    d["rqr"] = r["q"]
    d["rlr"] = r["lfc"]
    r.close()
    assert np.array_equal(d["targets"], d["rtargets"]), "target order drift"
    real = pl.read_parquet(BASE348 / "reference_de.parquet")
    d["real"] = real
    # parquet truth vs npz truth must agree exactly
    piv = real.pivot(on="feature", index="target",
                     values=["p_value", "p_adj", "log2_fold_change"])
    piv = piv.sort("target")
    cols = [c for c in piv.columns if c != "target"]
    feats = [c.split("_", 1)[1] if False else c for c in cols]
    # pivot col names are like "p_value_<gene>"? check shape instead via join below
    d["piv"] = piv
    rec = pl.read_csv(BASE348 / "per-target-cap15.csv").filter(
        pl.col("metric") == REACH).sort("perturbation")
    d["rec"] = rec
    ze = np.load(BASE348 / "effects.npz", allow_pickle=False)
    d["e_genes"] = ze["genes"].astype(str)
    d["e_cap15"] = ze["cap15"]
    d["e_targets"] = ze["targets"].astype(str)
    ze.close()
    return d


def scorer_reach(d):
    pred = E.de_frame(d["targets"], d["rgenes"], d["pp"], d["qp"], d["lp"])
    prepared = prepare_de(pred, d["real"], control=E.CTRL,
                          sort_by="abs_log2_fold_change",
                          p_adj_threshold=ALPHA, nan_lfc_policy="mask")
    cfg = EvalConfig.from_preset("vcc2026")
    rows = dispatch_de_metrics([REACH], prepared, cfg)
    return {r["perturbation"]: r["value"] for r in rows}


def build_long(d):
    """One row per (target, universe gene) with pred + truth fields."""
    n, g = len(d["targets"]), len(d["rgenes"])
    # verify parquet agrees with npz truth
    rp = d["real"].sort(["target", "feature"])
    assert rp.height == n * g, (rp.height, n * g)
    t_idx = {t: i for i, t in enumerate(d["targets"])}
    gi = {f: j for j, f in enumerate(d["rgenes"])}
    rt = rp["target"].to_numpy()
    rf = rp["feature"].to_numpy()
    ii = np.array([t_idx[t] for t in rt])
    jj = np.array([gi[f] for f in rf])
    rq = rp["p_adj"].to_numpy().reshape(n, g)  # sorted target-major?
    # rows are sorted by (target, feature); rebuild per-target matrices by index
    order = np.lexsort((jj, ii))
    # simpler: fill via fancy indexing
    qrm = np.empty((n, g)); lrm = np.empty((n, g)); prm = np.empty((n, g))
    qflat = rp["p_adj"].to_numpy()
    lflat = rp["log2_fold_change"].to_numpy()
    pflat = rp["p_value"].to_numpy()
    qrm[ii, jj] = qflat; lrm[ii, jj] = lflat; prm[ii, jj] = pflat
    dq = np.max(np.abs(qrm - d["rqr"]))
    dl = np.nanmax(np.abs(lrm - d["rlr"]))
    print(f"parquet-vs-npz truth: max|dq|={dq:.3e} max|dlfc|={dl:.3e}", flush=True)
    assert dq == 0.0 and (dl == 0.0 or np.isnan(dl)), "truth source drift"
    return qrm, lrm


def design_matrix(d):
    pos = {gg: i for i, gg in enumerate(d["e_genes"])}
    assert np.array_equal(np.sort(d["e_targets"]), np.sort(d["targets"]))
    epos = {t: i for i, t in enumerate(d["e_targets"])}
    D = np.zeros((len(d["targets"]), len(d["rgenes"])))
    missing = []
    for j, gg in enumerate(d["rgenes"]):
        if gg in pos:
            D[:, j] = d["e_cap15"][[epos[t] for t in d["targets"]], pos[gg]]
        else:
            missing.append(gg)
    print(f"design coverage of universe: {len(d['rgenes'])-len(missing)}/{len(d['rgenes'])} "
          f"missing={len(missing)}", flush=True)
    return D


def walk_reach(order_idx, match, adjud, n_conf):
    """order_idx: per-target ascending rank positions over pool genes."""
    n = match.shape[0]
    R = np.full(n, np.nan)
    Ks = np.zeros(n, dtype=int)
    for i in range(n):
        if n_conf[i] == 0:
            continue
        m = match[i][order_idx[i]]
        a = adjud[i][order_idx[i]]
        ma = m[a]
        denom = np.arange(1, ma.size + 1)
        if ma.size == 0:
            R[i] = 0.0
            continue
        pur = np.cumsum(ma) / denom
        ok = np.flatnonzero(pur >= P0)
        k = int(ok.max() + 1) if ok.size else 0
        Ks[i] = k
        R[i] = k / n_conf[i]
    return R, Ks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["known", "order", "profile"])
    a = ap.parse_args()
    t0 = time.monotonic()
    d = load()
    n, g = len(d["targets"]), len(d["rgenes"])
    qrm, lrm = build_long(d)
    pp, qp, lp = d["pp"], d["qp"], d["lp"]

    # shared per-(target,gene) fields
    ref_sig = qrm < ALPHA
    own = np.equal.outer(d["targets"], d["rgenes"])
    pool = ref_sig & ~own
    n_conf = pool.sum(1)
    lr = lrm
    adjud = np.isfinite(lr) & (lr != 0)
    # pred commits iff finite nonzero
    commit = np.isfinite(lp) & (lp != 0)
    match = adjud & commit & (np.sign(lp) == np.sign(lr))
    print(f"targets={n} universe={g} pool-empty={(n_conf==0).sum()} "
          f"unadjudicable-in-pool={(pool & ~adjud).sum()}", flush=True)

    feat = np.tile(d["rgenes"], (n, 1))
    qm = np.where(np.isfinite(qp), qp, np.inf)
    pm = np.where(np.isfinite(pp), pp, np.inf)
    ala = np.abs(lp)
    ala_m = np.where(np.isfinite(ala), -ala, np.inf)  # asc: biggest first, nan last
    sig = (np.isfinite(qp) & (qp < ALPHA))

    def lex(*keys_last_primary):
        # keys given primary-last; feature always least significant
        idx = np.empty((n, g), dtype=int)
        for i in range(n):
            idx[i] = np.lexsort((feat[i],) + tuple(k[i] for k in keys_last_primary))
        return idx

    current_key = (qm, pm, ala_m, ~sig)  # primary-last order: ~sig primary? no:
    # lexsort: LAST key is primary. want primary ~sig, then qm, pm, ala_m, feature.
    cur = lex(ala_m, pm, qm, ~sig)

    if a.mode == "known":
        s = scorer_reach(d)
        rec = {t: v for t, v in zip(d["rec"]["perturbation"], d["rec"]["value"])}
        Rmine, Ks = walk_reach(cur, match & pool, adjud & pool, n_conf)
        diffs = []
        for i, t in enumerate(d["targets"]):
            sv, rv = s[t], rec[t]
            mv = Rmine[i]
            ok_s = (np.isnan(sv) and np.isnan(rv)) or abs(sv - rv) <= 1e-12
            bothnan = np.isnan(mv) and (rv is None or (isinstance(rv, float) and np.isnan(rv)))
            dm = 0.0 if bothnan else abs(mv - rv)
            diffs.append((t, sv, rv, mv, ok_s, dm))
        maxd_s = max(abs(sv - rv) for _, sv, rv, _, _, _ in diffs
                     if not (np.isnan(sv) and np.isnan(rv)))
        maxd_m = max(dm for *_, dm in diffs)
        print(f"scorer-vs-recorded max abs diff: {maxd_s:.3e}", flush=True)
        print(f"mywalk-vs-recorded max abs diff: {maxd_m:.3e}", flush=True)
        bad = [(t, sv, rv, mv) for t, sv, rv, mv, _, dm in diffs if dm > 1e-9]
        print(f"targets off by >1e-9 in my walk: {len(bad)}", flush=True)
        for b in bad[:10]:
            print("  MISS", b, flush=True)
        out = {"scorer_max_abs_diff": maxd_s, "walk_max_abs_diff": maxd_m,
               "n_bad": len(bad),
               "nan_targets": [t for i, t in enumerate(d["targets"]) if n_conf[i] == 0],
               "mean_recorded": float(np.nanmean([rv for _, _, rv, _, _, _ in diffs])),
               "mean_walk": float(np.nanmean(Rmine)),
               "seconds": time.monotonic() - t0}
        WORK.mkdir(parents=True, exist_ok=True)
        (WORK / "known-answer.json").write_text(json.dumps(out, indent=2))
        print(json.dumps(out, indent=2), flush=True)

    elif a.mode == "order":
        D = design_matrix(d)
        ad = np.abs(D)
        ad_m = np.where(np.isfinite(ad), -ad, np.inf)
        common = D.mean(0)
        agree = (D * common) > 0
        lp_neg = lp < 0
        o3 = lex(ala_m)
        o4 = lex(ad_m)
        o5 = lex(ad_m, ~agree)
        o7b = lex(ad_m, agree)
        o7a = lex(ala_m, ~lp_neg)
        o7c = lex(ala_m, ~sig)
        # oracle: correct first, ties by current key
        correct = (match & pool)
        o6 = lex(ala_m, pm, qm, ~correct)
        schemes = {"O1_current": cur, "O3_abs_lfc": o3, "O4_abs_design": o4,
                   "O5_agree_first": o5, "O6_oracle": o6, "O7a_down_first": o7a,
                   "O7b_disagree_first": o7b, "O7c_sig_abs_lfc": o7c}
        res = {}
        for name, idx in schemes.items():
            Rr, _ = walk_reach(idx, match & pool, adjud & pool, n_conf)
            assert set(np.flatnonzero(np.isnan(Rr))) == set(np.flatnonzero(n_conf == 0))
            res[name] = {"mean": float(np.nanmean(Rr)),
                         "values": {t: (None if np.isnan(v) else float(v))
                                    for t, v in zip(d["targets"], Rr)}}
            print(f"{name}: mean reach = {res[name]['mean']:.6f}", flush=True)
        # O2 random
        rng_means = []
        for seed in range(100):
            rng = np.random.default_rng(3480000 + seed)
            idx = np.empty((n, g), dtype=int)
            for i in range(n):
                idx[i] = rng.permutation(g)
            Rr, _ = walk_reach(idx, match & pool, adjud & pool, n_conf)
            rng_means.append(float(np.nanmean(Rr)))
        rng_means = np.array(rng_means)
        res["O2_random"] = {"mean": float(rng_means.mean()),
                            "sd": float(rng_means.std()),
                            "min": float(rng_means.min()),
                            "max": float(rng_means.max())}
        print(f"O2_random: mean={rng_means.mean():.6f} sd={rng_means.std():.6f} "
              f"range=[{rng_means.min():.6f},{rng_means.max():.6f}]", flush=True)
        base = res["O1_current"]["mean"]
        for name in schemes:
            print(f"d({name} - O1) = {res[name]['mean']-base:+.6f}", flush=True)
        out = {"panel_mean_current": base, "schemes": {
            k: ({"mean": v["mean"]} if "values" not in v else {"mean": v["mean"]})
            for k, v in res.items()}, "random_means": rng_means.tolist(),
            "seconds": time.monotonic() - t0}
        (WORK / "orderings.json").write_text(json.dumps(out, indent=2))
        pl.DataFrame({"target": d["targets"], **{
            k: [v["values"][t] for t in d["targets"]]
            for k, v in res.items() if "values" in v}}).write_csv(
            WORK / "per-target-orderings.csv")

    elif a.mode == "profile":
        Rr, Ks = walk_reach(cur, match & pool, adjud & pool, n_conf)
        depth = (adjud & pool).sum(1)
        prof = {}
        for k in (5, 10, 20, 50, 100, 200):
            sel = np.flatnonzero(depth >= k)
            pk = []
            for i in sel:
                m = (match & pool)[i][cur[i]][(adjud & pool)[i][cur[i]]][:k]
                pk.append(m.mean())
            prof[str(k)] = {"n": int(len(sel)),
                            "mean_P": float(np.mean(pk)) if pk else None,
                            "frac_ge90": float(np.mean(np.array(pk) >= P0)) if pk else None}
            print(f"k={k}: n={len(sel)} meanP={prof[str(k)]['mean_P']} "
                  f"frac>=0.9={prof[str(k)]['frac_ge90']}", flush=True)
        scored = n_conf > 0
        kstar = Ks[scored]
        f0 = float(np.mean(kstar == 0))
        flt10 = float(np.mean(kstar < 10))
        print(f"k*==0: {f0:.3f}  k*<10: {flt10:.3f}", flush=True)
        # first-position diagnosis
        pos1miss, rec_depth = [], []
        for i in np.flatnonzero(scored):
            m = (match & pool)[i][cur[i]][(adjud & pool)[i][cur[i]]]
            if m.size == 0:
                continue
            if not m[0]:
                pos1miss.append(str(d["targets"][i]))
                pur = np.cumsum(m) / np.arange(1, m.size + 1)
                ok = np.flatnonzero(pur >= P0)
                rec_depth.append(int(ok.max() + 1) if ok.size else -1)
        print(f"pos1-miss targets: {len(pos1miss)}; recovered k*: "
              f"{sorted(set(rec_depth))[:12]}...", flush=True)
        import collections
        hist = collections.Counter(rec_depth)
        out = {"P_profile": prof, "frac_kstar_0": f0, "frac_kstar_lt10": flt10,
               "pos1_miss_targets": pos1miss,
               "pos1_miss_recovery_kstar_hist": {str(k): v for k, v in sorted(hist.items())},
               "seconds": time.monotonic() - t0}
        (WORK / "profile.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
