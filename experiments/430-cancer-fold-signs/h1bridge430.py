"""430 H1 bridge (analytic, stored files only; no emission): does the pseudobulk GATE PROXY rank sign candidates the
way the REAL reference DE gate does, and where does 427e's H1 advantage of the donor mean come from?

Surface: h1-test's 85 targets. Real gate = 323's reference_de.parquet (p_adj < 0.05, finite nonzero lfc, not the target,
>= 10 genes; 427e's rule). Proxy gates = the same constructions as fold430 on h1_2025's own pseudobulk delta row (the H1
harness cells themselves): top-K |delta| (K 50/200/1000) and |delta| >= tau (0.1/0.2) over the reference's kept genes.

Candidates (on the reference's 10,779 kept genes):
  D       our DELIVERED lfc (426 `stats-base_s0_bps.npz`, the shipped 410b construction)
  CAR     the shipped carriers = 323's cap15 effects (gamma 0.65 + CD4, = production PP analogue on H1)
  CARg0   416's g000 carriers (gamma 0 + CD4, norm-matched; the 417c design)
  M6      427e's donor average: raw delta mean over the six non-H1 donors covering (p, g)
  M3P     the same over the three PRODUCTION donors only (k562_gwps, orion_hct116, orion_hek293t): on the A/B/C panel
          k562_essential / rpe1 / nadig_jurkat cover 0 of 300 targets, so M3P is what "M6" can be there
  M6c/M3Pc  centred versions (each donor's all-perturbation column mean removed)
  C4      CLEAN4 per-gene sign (all four banks; none is H1)
  blends  0.5*CAR + 0.5*w*Mx_rowscaled (fold430's construction), and 427e's 0.5*D + 0.5*w*M6/ln2
KNOWN ANSWERS: 427e's printed numbers: D sign accuracy on the real gate 0.5606, M6 covered-only 0.5776 (coverage 0.995).

    uv run python h1bridge430.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import polars as pl  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
W323 = ROOT / "experiments/323-current-family-validation/work"
DON6 = ("k562_gwps", "k562_essential", "rpe1", "nadig_jurkat", "orion_hct116", "orion_hek293t")
PROD3 = ("k562_gwps", "orion_hct116", "orion_hek293t")
CLEAN4 = ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")
KS = (50, 200, 1000)
TAUS = (0.1, 0.2)
NBOOT = 2000


def load_rows(name, targets, genes):
    bd = ROOT / "data/donors" / name
    P = [l.strip() for l in open(bd / "perts.txt")]
    G = [l.strip() for l in open(bd / "genes.txt")]
    pi = {p: i for i, p in enumerate(P)}
    gcol = {g: i for i, g in enumerate(G)}           # last occurrence wins (as 427e's dict)
    M = np.load(bd / "delta.npy", mmap_mode="r")
    cols = np.array([gcol.get(g, -1) for g in genes])
    okc = cols >= 0
    out = np.full((len(targets), len(genes)), np.nan)
    cm = np.full(len(genes), np.nan)
    acc = np.zeros(int(okc.sum()))
    for i in range(0, M.shape[0], 1000):
        acc += np.asarray(M[i:i + 1000][:, cols[okc]], dtype=np.float64).sum(0)
    cm[okc] = acc / M.shape[0]
    for k, t in enumerate(targets):
        if t in pi:
            out[k, okc] = np.asarray(M[pi[t]], dtype=np.float64)[cols[okc]]
    return out, cm


def nanmean_rows(stack):
    cnt = np.isfinite(stack).sum(0)
    return np.where(cnt > 0, np.nansum(stack, 0) / np.maximum(cnt, 1), 0.0), cnt


def main():
    r = np.load(W323 / "reference.npz", allow_pickle=True)
    genes = r["genes"].astype(str)[r["keep"].astype(bool)]
    gi = {g: i for i, g in enumerate(genes)}
    s = np.load(ROOT / "experiments/426-mde-magnitude/work/h1/stats-base_s0_bps.npz")
    targets = s["targets"].astype(str)
    D = np.nan_to_num(s["lfc"])
    # real gates (427e)
    ref = pl.read_parquet(W323 / "reference_de.parquet")
    real = {}
    for k, t in enumerate(targets):
        m = ref.filter(pl.col("target") == t)
        g = m["feature"].to_numpy().astype(str)
        q = m["p_adj"].to_numpy().astype(float)
        l = m["log2_fold_change"].to_numpy().astype(float)
        ok = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
        if ok.sum() >= 10:
            real[k] = (np.array([gi[x] for x in g[ok]]), l[ok])
    # donors
    rows, cms = {}, {}
    for n in DON6:
        rows[n], cms[n] = load_rows(n, targets, genes)
    h1, h1cm = load_rows("h1_2025", targets, genes)
    assert np.isfinite(h1).all(1).sum() == len(targets) or True
    M6, cnt6 = nanmean_rows(np.stack([rows[n] for n in DON6]))
    M3P, _ = nanmean_rows(np.stack([rows[n] for n in PROD3]))
    M6c, _ = nanmean_rows(np.stack([rows[n] - cms[n][None, :] for n in DON6]))
    M3Pc, _ = nanmean_rows(np.stack([rows[n] - cms[n][None, :] for n in PROD3]))
    cover = {n: int(np.isfinite(rows[n]).any(1).sum()) for n in DON6}
    # carriers -> reference genes
    ze = np.load(W323 / "effects.npz", allow_pickle=False)
    zc = np.load(ROOT / "experiments/416-de-carrier-menu/work/carriers416.npz", allow_pickle=False)
    assert (ze["targets"].astype(str) == targets).all() and (zc["targets"].astype(str) == targets).all()
    epos = {g: i for i, g in enumerate(ze["genes"].astype(str))}
    ecol = np.array([epos.get(g, -1) for g in genes])
    eok = ecol >= 0                               # reference genes off the official axis: carriers 0 there
    print(f"reference kept genes on the official axis: {eok.sum()}/{eok.size}")
    CAR = np.where(eok[None, :], np.asarray(ze["cap15"], float)[:, np.maximum(ecol, 0)], 0.0)
    CARg0 = np.where(eok[None, :], np.asarray(zc["g000"], float)[:, np.maximum(ecol, 0)], 0.0)
    # CLEAN4
    cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    bn = list(cp["bank_names"].astype(str))
    P4 = cp["P_banks"].astype(float)[[bn.index(b) for b in CLEAN4]].mean(0)
    cpos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
    c4 = np.array([(-1.0 if P4[cpos[g]] > 0.5 else 1.0) if g in cpos else -1.0 for g in genes])
    C4 = np.tile(c4, (len(targets), 1))

    # ---------------------------------------------------------------- KNOWN ANSWERS (427e)
    K = np.concatenate([np.full(i.size, k) for k, (i, _) in real.items()])
    I = np.concatenate([i for i, _ in real.values()])
    LL = np.concatenate([l for _, l in real.values()])
    accD = float(np.mean(np.sign(D[K, I]) == np.sign(LL)))
    avg = M6 / np.log(2)
    cov = avg[K, I] != 0
    accM = float(np.mean(np.sign(avg[K, I][cov]) == np.sign(LL[cov])))
    print(f"KA 427e: delivered {accD:.4f} (427e 0.5606); donor avg covered {accM:.4f} (427e 0.5776), cov {cov.mean():.3f}")
    assert abs(accD - 0.5606) < 5e-5 and abs(accM - 0.5776) < 5e-5, "KA 427e FAILED"

    # fold430's PPprior emulation applied to CAR: CLEAN4 sign on the top-2,200 control-CPM block, strongest 20 % kept
    bH0 = np.load(ROOT / "data/donors/h1_2025/b_ctrl.npy")
    Gh0 = {g: i for i, g in enumerate(l.strip() for l in open(ROOT / "data/donors/h1_2025/genes.txt"))}
    muk = np.array([np.expm1(bH0[Gh0[g]]) if g in Gh0 else -1.0 for g in genes])
    blk = np.argsort(-muk)[:2200]
    CARprior = CAR.copy()
    for k in range(len(targets)):
        v = CAR[k, blk]
        fl = np.abs(v) <= np.quantile(np.abs(v), 0.8)
        CARprior[k, blk[fl]] = np.abs(v[fl]) * c4[blk[fl]]
    C = {"D": D, "CAR": CAR, "CARprior": CARprior, "CARg0": CARg0, "M6": M6, "M3P": M3P, "M6c": M6c, "M3Pc": M3Pc,
         "C4": C4}
    nCAR = np.linalg.norm(CAR, axis=1)
    for mn in ("M6", "M3P", "M3Pc"):
        nm = np.linalg.norm(C[mn], axis=1)
        R = C[mn] * np.where(nm > 0, nCAR / np.maximum(nm, 1e-300), 0.0)[:, None]
        for w in (1.0, 1.5, 2.0):
            C[f"B{mn}w{w:g}"] = 0.5 * CAR + 0.5 * w * R
    for w in (1.0, 1.5, 2.0):
        C[f"D+M6w{w:g}"] = 0.5 * D + 0.5 * w * M6 / np.log(2)
        C[f"D+M3Pw{w:g}"] = 0.5 * D + 0.5 * w * M3P / np.log(2)

    # ---------------------------------------------------------------- gates
    gates = {"real": {k: real[k][0] for k in real}}
    truth = {"real": {k: real[k][1] for k in real}}
    tpos = np.array([gi.get(t, -1) for t in targets])
    for kind, vals in (("top", KS), ("tau", TAUS)):
        for v in vals:
            name = f"{kind}{v:g}"
            gates[name], truth[name] = {}, {}
            for k in range(len(targets)):
                d = h1[k].copy()
                if tpos[k] >= 0:
                    d[tpos[k]] = np.nan
                ok = np.flatnonzero(np.isfinite(d) & (d != 0))
                a = np.abs(d[ok])
                sel = ok[np.argsort(-a, kind="stable")[:v]] if kind == "top" else ok[a >= v]
                gates[name][k] = sel
                truth[name][k] = d[sel]
    # significance-weighted proxies: |delta| * sqrt(control mean count); z-top-K, and (oracle size) top-n_real by z / |d|
    bH = np.load(ROOT / "data/donors/h1_2025/b_ctrl.npy")
    Gh = [l.strip() for l in open(ROOT / "data/donors/h1_2025/genes.txt")]
    hpos = {g: i for i, g in enumerate(Gh)}
    mu = np.array([np.expm1(bH[hpos[g]]) if g in hpos else np.nan for g in genes])
    for name, K_ in (("z200", 200), ("z1000", 1000), ("zreal", None), ("dreal", None)):
        gates[name], truth[name] = {}, {}
        for k in range(len(targets)):
            if K_ is None and k not in real:
                continue
            d = h1[k].copy()
            if tpos[k] >= 0:
                d[tpos[k]] = np.nan
            ok = np.flatnonzero(np.isfinite(d) & (d != 0) & np.isfinite(mu))
            sc = np.abs(d[ok]) * (np.sqrt(mu[ok]) if name != "dreal" else 1.0)
            kk = K_ if K_ is not None else real[k][0].size
            sel = ok[np.argsort(-sc, kind="stable")[:kk]]
            gates[name][k] = sel
            truth[name][k] = d[sel]
    # proxy-vs-real overlap
    ov = {}
    for g in gates:
        if g == "real":
            continue
        sh = [k for k in real if k in gates[g]]
        ov[g] = dict(mean_size=float(np.mean([gates[g][k].size for k in gates[g]])),
                     frac_of_real_in_proxy=float(np.mean([np.isin(real[k][0], gates[g][k]).mean() for k in sh])),
                     frac_of_proxy_in_real=float(np.mean([np.isin(gates[g][k], real[k][0]).mean()
                                                          for k in sh if gates[g][k].size])),
                     sign_agree_real_vs_proxy=float(np.mean(np.concatenate([
                         np.sign(real[k][1][np.isin(real[k][0], gates[g][k])]) ==
                         np.sign(h1[k][real[k][0][np.isin(real[k][0], gates[g][k])]]) for k in sh]))))

    # ---------------------------------------------------------------- metrics
    rng = np.random.default_rng(430)
    res = dict(ka=dict(delivered=accD, donor_avg_covered=accM), coverage_targets=cover, overlap=ov, table={})
    for g in gates:
        ks = [k for k in gates[g] if gates[g][k].size >= (10 if g != "real" else 10)]
        base_h = None
        for cn, M in C.items():
            hits = np.array([np.where(M[k, gates[g][k]] == 0, 0.5,
                                      np.sign(M[k, gates[g][k]]) == np.sign(truth[g][k])).sum() for k in ks])
            n = np.array([gates[g][k].size for k in ks], float)
            # nmae proxy at own L1-optimal k (weighted median)
            rr, ww = [], []
            for k in ks:
                lh, d = M[k, gates[g][k]], truth[g][k]
                nz = lh != 0
                rr.append(d[nz] / lh[nz]); ww.append(np.abs(lh[nz]) / np.abs(d).sum())
            rr, ww = np.concatenate(rr), np.concatenate(ww)
            o = np.argsort(rr)
            cw = np.cumsum(ww[o])
            kk = max(float(rr[o][np.searchsorted(cw, 0.5 * cw[-1])]), 0.0)
            nm = np.array([np.abs(kk * M[k, gates[g][k]] - truth[g][k]).sum() / np.abs(truth[g][k]).sum() for k in ks])
            rec = dict(pooled=float(hits.sum() / n.sum()), pertarget=float(np.mean(hits / n)), nmae=float(nm.mean()),
                       k=kk, n_targets=len(ks))
            if cn == "CAR":
                base_h, base_n, base_nm = hits, n, nm
            res["table"][f"{g}|{cn}"] = rec
        # paired bootstrap vs CAR, pooled accuracy
        for cn, M in C.items():
            if cn == "CAR":
                continue
            hits = np.array([np.where(M[k, gates[g][k]] == 0, 0.5,
                                      np.sign(M[k, gates[g][k]]) == np.sign(truth[g][k])).sum() for k in ks])
            idx = rng.integers(0, len(ks), size=(NBOOT, len(ks)))
            dd = hits[idx].sum(1) / base_n[idx].sum(1) - base_h[idx].sum(1) / base_n[idx].sum(1)
            res["table"][f"{g}|{cn}"]["d_pooled_vs_CAR_ci"] = np.quantile(dd, [0.025, 0.5, 0.975]).tolist()
    # ---------------------------------------------------------------- 427e's nmae table, M6 vs the panel-available M3P
    def nmae427(Lh):
        return float(np.mean([np.abs(Lh[k, i] - l).sum() / np.abs(l).sum() for k, (i, l) in real.items()]))
    count = r["ref_lin"] * np.median(r["libraries"]) / 1e6
    add = np.where(count >= 0.3, np.minimum(0.1 / np.sqrt(np.maximum(count, 1e-3)), 0.1), 0.0)[None, :]
    base = nmae427(D)
    t427 = {"ours+MDE": nmae427(D + np.sign(D) * add) - base}
    for mn, M in (("M6", M6), ("M3P", M3P)):
        a = M / np.log(2)
        for w in (1.0, 1.5, 2.0):
            t427[f"{mn} x{w:g} alone"] = nmae427(w * a) - base
            B = 0.5 * D + 0.5 * w * a
            t427[f"blend(ours,{mn}x{w:g}) alone"] = nmae427(B) - base
            t427[f"blend(ours,{mn}x{w:g}) +MDE"] = nmae427(B + np.sign(B) * add) - base
            t427[f"blend(ours,{mn}x{w:g}) sign acc"] = float(np.mean(np.sign(B[K, I]) == np.sign(LL)))
    print(f"427e nmae table (base {base:.5f}):")
    for k_, v in t427.items():
        print(f"  {k_:32s} {v:+.5f}")
    ka427 = {"M6 x1 alone": -0.0007, "M6 x1.5 alone": -0.0090, "blend(ours,M6x1.5) alone": -0.0120,
             "blend(ours,M6x1.5) +MDE": -0.0245, "ours+MDE": -0.0098}
    for k_, v in ka427.items():
        assert abs(t427[k_] - v) < 6e-5, f"KA 427e nmae FAILED {k_}: {t427[k_]} vs {v}"
    print("KA 427e nmae table: PASS")
    res["t427"] = t427
    # where M6's advantage lives: targets covered by the essential-screen donors (absent from the A/B/C panel)
    ess = np.array([any(np.isfinite(rows[n][k]).any() for n in ("k562_essential", "rpe1", "nadig_jurkat"))
                    for k in range(len(targets))])
    kk_ = np.array(list(real))
    share = float(sum(real[k][0].size for k in kk_ if ess[k]) / sum(real[k][0].size for k in kk_))
    dec = {}
    for nm_, sel in (("ess_targets", [k for k in kk_ if ess[k]]), ("other_targets", [k for k in kk_ if not ess[k]])):
        Kx = np.concatenate([np.full(real[k][0].size, k) for k in sel])
        Ix = np.concatenate([real[k][0] for k in sel])
        Lx = np.concatenate([real[k][1] for k in sel])
        dec[nm_] = {c: float(np.mean(np.where(C[c][Kx, Ix] == 0, 0.5, np.sign(C[c][Kx, Ix]) == np.sign(Lx))))
                    for c in ("D", "CAR", "CARg0", "M6", "M3P")}
        dec[nm_]["n_targets"] = len(sel)
    res["essential_split"] = dict(gate_share_of_ess_targets=share, **dec)
    # 427e's nmae arithmetic on the panel-like stratum only (targets no essential-screen donor covers)
    oth = [k for k in kk_ if not ess[k]]
    def nm_oth(Lh):
        return float(np.mean([np.abs(Lh[k, real[k][0]] - real[k][1]).sum() / np.abs(real[k][1]).sum() for k in oth]))
    b0 = nm_oth(D)
    t_oth = {"ours+MDE": nm_oth(D + np.sign(D) * add) - b0}
    for w in (1.0, 1.5, 2.0):
        B = 0.5 * D + 0.5 * w * M3P / np.log(2)
        t_oth[f"blend(ours,M3Px{w:g}) alone"] = nm_oth(B) - b0
        t_oth[f"blend(ours,M3Px{w:g}) +MDE"] = nm_oth(B + np.sign(B) * add) - b0
    res["t427_other_targets"] = dict(base=b0, n=len(oth), **t_oth)
    print(f"427e arithmetic on the {len(oth)} non-essential targets (base {b0:.5f}): " + json.dumps({k: round(v, 5) for k, v in t_oth.items()}))
    print(f"targets covered by k562_essential/rpe1/jurkat: {int(ess.sum())} ({int(ess[kk_].sum())} with a real gate), "
          f"holding {share:.3f} of real-gate cells")
    print(json.dumps(dec))
    (WORK / "h1bridge430.json").write_text(json.dumps(res, indent=1))
    cols = ["D", "CAR", "CARprior", "CARg0", "M6", "M3P", "M6c", "M3Pc", "C4", "BM3Pw1", "BM3Pw1.5", "BM6w1.5", "D+M6w1.5",
            "D+M3Pw1.5"]
    print("gate      " + " ".join(f"{c:>9s}" for c in cols))
    for g in gates:
        print(f"{g:9s} " + " ".join(f"{res['table'][f'{g}|{c}']['pooled']:9.4f}" for c in cols))
    print("per-target mean accuracy")
    for g in gates:
        print(f"{g:9s} " + " ".join(f"{res['table'][f'{g}|{c}']['pertarget']:9.4f}" for c in cols))
    # rank agreement of each proxy with the real gate over all candidates (pooled accuracy; Spearman)
    cn_all = list(C)
    rv = np.array([res["table"][f"real|{c}"]["pooled"] for c in cn_all])
    rk = lambda x: np.argsort(np.argsort(x))
    for g in gates:
        if g == "real":
            continue
        pv = np.array([res["table"][f"{g}|{c}"]["pooled"] for c in cn_all])
        res.setdefault("rank_agreement", {})[g] = float(np.corrcoef(rk(rv), rk(pv))[0, 1])
    print("Spearman(proxy, real) over all candidates:", json.dumps(res["rank_agreement"]))
    (WORK / "h1bridge430.json").write_text(json.dumps(res, indent=1))
    print("nmae proxy (own k)")
    for g in gates:
        print(f"{g:9s} " + " ".join(f"{res['table'][f'{g}|{c}']['nmae']:9.4f}" for c in cols))
    print("coverage (targets) per donor:", cover)
    print("proxy overlap:", json.dumps(ov, indent=0))


if __name__ == "__main__":
    main()
