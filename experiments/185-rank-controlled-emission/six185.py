"""185 Parts B/C/D — read the members, the delivery fidelity, the P(k) curve, the call-set
frontier and the a21/a22 retrodiction off the h1-test DE tables written by `realize185.py`.

Known-answer rows (abort if any fails):
  a15  reach 0.18283 / median n_pred 6035        armC  reach 0.20291 / n_pred 3581
  |H| key on a15's table 0.23533                  ORACLE 0.64553
  this file's own walk == debench.reach() on every arm
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
from scipy.stats import rankdata, spearmanr

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
from debench import ANCHOR, Surface  # noqa: E402

ALPHA = 0.05
TAG = os.environ.get("OUTTAG", "de185")
P = os.path.join(WORK, TAG + ".%s.npz")
ARMS = ["a15", "a15_b", "armC_a8_flat", "rc_hist_a", "rc_hist_a_b", "rc_hist_s", "rc_perm",
        "rc_hybrid", "rc_lfc", "rc_k1000", "rc_k2000", "rc_k4000", "a21_like", "a22_like"]
ARMS = [a for a in ARMS if os.path.exists(P % a)]
# board-fitted 1/span (181 board181.py; inside AGENTS.md §4's published spans)
INV = dict(reach=1.12927, fid=3.32482, jac=2.72081, nmae=1.64787)
SHRINK = 0.4277                       # 164 §6 / 171 §7 h1-test -> board reach delta
A21 = dict(nmae=-0.064214, fid=-0.029464, reach=-0.018482, jac=-0.016583)
A22 = dict(nmae=+0.001520, fid=-0.004447, reach=-0.010656, jac=-0.001523)


class PowerSurface(Surface):
    """157's proportional power matching: pool_p = reference top-round(f * n_real_p)."""

    def __init__(self, path, factor=1.0):
        super().__init__(path)
        z = np.load(path, allow_pickle=True)
        rp = np.where(np.isfinite(z["real_p"]), z["real_p"], np.inf)
        base = (z["real_padj"] < ALPHA) & np.isfinite(z["real_padj"])
        sig = np.zeros(rp.shape, bool)
        for i in range(rp.shape[0]):
            k = int(round(factor * int(base[i].sum())))
            k = max(0, min(k, int(np.isfinite(z["real_p"][i]).sum())))
            if k:
                sig[i, np.argsort(rp[i], kind="stable")[:k]] = True
        self.real_sig = sig & np.isfinite(z["real_p"])
        for p, j in enumerate(self.self_idx):
            if j >= 0:
                self.real_sig[p, j] = False


def walk(s, key=None):
    """My own reach walk (BRIEF §5) + P(k) curve.  Returns reach, P(k) dict, k* list."""
    ks, out = [], []
    pk = {k: [] for k in (1, 2, 3, 5, 10, 20, 50, 100)}
    for p in range(len(s.perts)):
        pool = np.flatnonzero(s.real_sig[p])
        if pool.size == 0:
            continue
        if key is None:
            order = np.lexsort((s.genes[pool], -np.abs(s.pred_lfc[p, pool]), s.pred_p[p, pool],
                               s.pred_padj[p, pool], ~s.pred_sig[p, pool]))
        else:
            order = np.lexsort((s.genes[pool], -key[p][pool]))
        pool = pool[order]
        pool = pool[s.adj[p, pool]]
        if pool.size == 0:
            out.append(0.0)
            ks.append(0)
            continue
        ok = (np.sign(s.pred_lfc[p, pool]) == np.sign(s.real_lfc[p, pool])).astype(float)
        pur = np.cumsum(ok) / np.arange(1, pool.size + 1)
        deep = np.flatnonzero(pur >= 0.9)
        k = int(deep[-1] + 1) if deep.size else 0
        ks.append(k)
        out.append(k / int(s.real_sig[p].sum()))
        for kk in pk:
            if pool.size >= kk:
                pk[kk].append(pur[kk - 1])
    return float(np.mean(out)), {k: float(np.mean(v)) for k, v in pk.items()}, ks


def members(s):
    nr, npd = s.counts()
    prec, rec = [], []
    for p in range(len(s.perts)):
        sel = np.flatnonzero(s.pred_sig[p] & s.adj[p])
        if sel.size:
            prec.append((np.sign(s.pred_lfc[p, sel]) == np.sign(s.real_lfc[p, sel])).mean())
        a_ = s.real_sig[p]
        if a_.sum():
            rec.append(int((a_ & s.pred_sig[p]).sum()) / int(a_.sum()))
    return dict(n_pred=int(np.median(npd)), reach=s.reach(), fid=s.fid(), prec=float(np.mean(prec)),
                jac=s.jac(), nmae=s.nmae(), recall=float(np.mean(rec)))


def spearmans(z, s, absH, expr):
    """171 §12a / 177 §4: median over perturbations of Spearman(-log10 pred_padj, x),
    over the DE universe, inside our call set, and over the reach pool."""
    T = absH.shape[0]
    res = {k: [] for k in ("H_univ", "H_call", "H_pool", "E_univ", "E_call", "E_pool")}
    lp = -np.log10(np.maximum(z["pred_padj"], 1e-300))
    lpr = -np.log10(np.maximum(z["pred_p"], 1e-300))
    for p in range(T):
        call = s.pred_sig[p]
        pool = s.real_sig[p]
        for lab, x in (("H", absH[p]), ("E", expr)):
            res[f"{lab}_univ"].append(spearmanr(lpr[p], x).correlation)
            if call.sum() > 10:
                res[f"{lab}_call"].append(spearmanr(lpr[p][call], x[call]).correlation)
            if pool.sum() > 10:
                res[f"{lab}_pool"].append(spearmanr(lpr[p][pool], x[pool]).correlation)
    return {k: float(np.median(v)) for k, v in res.items()}


def main():
    z15 = np.load(P % "a15", allow_pickle=True)
    absH = z15["absH"].astype(np.float64)
    genes = np.asarray(z15["genes"]).astype(str)
    meta = np.load(os.path.join(WORK, "ctrl_meta.npz"), allow_pickle=True)
    expr = np.log1p(meta["cpm"][meta["keep"]])
    S = {a: Surface(P % a) for a in ARMS}

    # ---- known-answer rows
    s15 = S["a15"]
    ka = dict(a15_reach=s15.reach(), a15_npred=members(s15)["n_pred"])
    ka["H_key"] = s15.reach(key=absH)
    orc = (np.sign(s15.pred_lfc) == np.sign(s15.real_lfc)).astype(float)
    ka["oracle"] = s15.reach(key=orc)
    if "armC_a8_flat" in S:
        ka["armC_reach"] = S["armC_a8_flat"].reach()
        ka["armC_npred"] = members(S["armC_a8_flat"])["n_pred"]
    mywalk = walk(s15)[0]
    print("KNOWN-ANSWER ROWS")
    rows = [("a15 reach (native BH, fresh emission)", ka["a15_reach"], 0.18283, 2e-5),
            ("a15 median n_pred (debench conv.; raw BH 6035)", ka["a15_npred"], 6034, 0.5),
            ("|H| key reach on a15's pool", ka["H_key"], 0.23533, 2e-5),
            ("ORACLE", ka["oracle"], 0.64553, 2e-5),
            ("my walk vs debench.reach (a15)", mywalk, ka["a15_reach"], 1e-9)]
    if "armC_reach" in ka:
        rows += [("a8-flat reach", ka["armC_reach"], 0.20291, 2e-5),
                 ("a8-flat median n_pred", ka["armC_npred"], 3581, 0.5)]
    bad = []
    for lab, got, want, tol in rows:
        ok = abs(got - want) <= tol
        print(f"  {lab:<42} {got:>10.5f}   expect {want:<10.5f} {'PASS' if ok else 'FAIL'}")
        if not ok:
            bad.append(lab)
    if bad and not os.environ.get("NOABORT"):
        sys.exit(f"ABORT (rule 10): {bad}")

    # ---- Part B: members + delivery fidelity + P(k)
    print("\nPART B — members on h1-test (85 perturbations), reach at proportional power f")
    print(f"{'arm':<14}{'n_pred':>7}{'reach f1':>10}{'f2':>8}{'f3.5':>8}{'fid':>8}{'prec':>8}"
          f"{'jac':>8}{'nmae':>8}{'recall':>8}{'P(1)':>7}{'P(3)':>7}{'P(10)':>7}{'k*med':>6}")
    res = {}
    for a in ARMS:
        s = S[a]
        mb = members(s)
        r2 = PowerSurface(P % a, 2.0).reach()
        r35 = PowerSurface(P % a, 3.5).reach()
        rw, pk, ks = walk(s)
        assert abs(rw - mb["reach"]) < 1e-9, a
        mb.update(reach_f2=r2, reach_f35=r35, pk=pk, kstar_med=float(np.median(ks)),
                  kstar0=float(np.mean(np.array(ks) == 0)))
        res[a] = mb
        print(f"{a:<14}{mb['n_pred']:>7}{mb['reach']:>10.5f}{r2:>8.5f}{r35:>8.5f}{mb['fid']:>8.5f}"
              f"{mb['prec']:>8.5f}{mb['jac']:>8.5f}{mb['nmae']:>8.5f}{mb['recall']:>8.4f}"
              f"{pk[1]:>7.4f}{pk[3]:>7.4f}{pk[10]:>7.4f}{mb['kstar_med']:>6.0f}")

    print("\nDELIVERY FIDELITY — median over perturbations of Spearman(-log10 p, x); "
          "x = |H| (the key) or control expression; universe / inside our call set / reach pool")
    print(f"{'arm':<14}{'|H| univ':>10}{'|H| call':>10}{'|H| pool':>10}{'expr univ':>11}"
          f"{'expr call':>11}{'expr pool':>11}")
    for a in ARMS:
        z = np.load(P % a, allow_pickle=True)
        sp = spearmans(z, S[a], absH, expr)
        res[a]["spearman"] = sp
        print(f"{a:<14}{sp['H_univ']:>+10.3f}{sp['H_call']:>+10.3f}{sp['H_pool']:>+10.3f}"
              f"{sp['E_univ']:>+11.3f}{sp['E_call']:>+11.3f}{sp['E_pool']:>+11.3f}")

    # the designed-vs-realised check for the rc arms, and what the |H| key reads on THEIR tables
    print("\nRC ARMS — counterfactual |H| key on each arm's own table (its ceiling) vs realised")
    for a in ARMS:
        if not a.startswith("rc"):
            continue
        s = S[a]
        print(f"  {a:<12} realised {s.reach():.5f}   |H| key on this table {s.reach(key=absH):.5f}"
              f"   oracle {s.reach(key=(np.sign(s.pred_lfc)==np.sign(s.real_lfc)).astype(float)):.5f}")

    # sign agreement between emitted lfc and H (172 §4.2 read 0.945 on the pool for a15)
    z15 = np.load(P % "a15", allow_pickle=True)
    Hs = np.load(f"{ROOT}/experiments/169-sharpening-six-member/work/h1test-effect.npz", allow_pickle=True)
    gi = {g: i for i, g in enumerate(np.asarray(Hs["genes"]).astype(str))}
    ti = {t: i for i, t in enumerate(np.asarray(Hs["targets"]).astype(str))}
    Hm = Hs["H"][np.ix_([ti[t] for t in s15.perts], [gi[g] for g in genes])]
    for a in ("a15", "rc_hist_a", "rc_hist_s"):
        if a in S:
            s = S[a]
            pool = s.real_sig & s.adj
            agree = (np.sign(s.pred_lfc[pool]) == np.sign(Hm[pool])).mean()
            print(f"  P(emitted sign == sign H) over the pool, {a}: {agree:.4f}")

    # ---- Part C: the call-set frontier
    print("\nPART C — fid / jac as a function of n_pred (top-K per perturbation by a key), at "
          "proportional power f (n_real scale); the h1-test n_real distribution is kept, only "
          "scaled")
    Ks = [50, 100, 200, 400, 800, 1600, 3200, 6400]
    keys = {"|H|": absH}
    rc = "rc_hist_a" if "rc_hist_a" in S else "a15"
    for a, lab in (("a15", "a15 p-order"), (rc, rc + " p-order")):
        zz = np.load(P % a, allow_pickle=True)
        keys[lab] = -np.log(np.maximum(zz["pred_p"], 1e-300))
    keys["|pred lfc|"] = np.abs(z15["pred_lfc"])
    keys["random"] = np.random.default_rng(0).random(absH.shape)
    frontier = {}
    for f in (0.5, 1.0, 2.0, 3.5):
        sf = PowerSurface(P % rc, f)
        nr = np.array([int(sf.real_sig[p].sum()) for p in range(len(sf.perts))])
        print(f"\n  f = {f}: median n_real {np.median(nr):.0f}, mean {nr.mean():.0f}")
        print(f"  {'key':<20}" + "".join(f"{'K='+str(K):>14}" for K in Ks))
        for kname, kv in keys.items():
            line_f, line_j = [], []
            for K in Ks:
                call = np.zeros(kv.shape, bool)
                for p in range(kv.shape[0]):
                    call[p, np.argpartition(-kv[p], K - 1)[:K]] = True
                fv, jv = sf.fid(call=call), sf.jac(call=call)
                line_f.append(fv)
                line_j.append(jv)
                frontier[(f, kname, K)] = (fv, jv)
            print(f"  {kname:<20}" + "".join(f"{fv:>7.4f}/{jv:<6.4f}" for fv, jv in zip(line_f, line_j)))
        # realised rc_kK arms on the same surface
        for a in ("rc_k1000", "rc_k2000", "rc_k4000"):
            if a in S:
                sa = PowerSurface(P % a, f)
                print(f"  realised {a:<11}{'':<3} fid {sa.fid():.4f}  jac {sa.jac():.4f}  "
                      f"n_pred {members(S[a])['n_pred']}")
        # optimum in board units: (d fid/span + d jac/span)/6 relative to a15's realised set
        s15f = PowerSurface(P % "a15", f)
        f0, j0 = s15f.fid(), s15f.jac()
        best = max(((K, kname, (frontier[(f, kname, K)][0] - f0) * INV["fid"] / 6
                     + (frontier[(f, kname, K)][1] - j0) * INV["jac"] / 6)
                    for kname in keys for K in Ks), key=lambda x: x[2])
        print(f"  a15 realised set at f={f}: fid {f0:.4f} jac {j0:.4f};  best top-K: K={best[0]} "
              f"by {best[1]} -> d(fid+jac) overall {best[2]:+.5f}")

    # ---- Part D: retrodiction of a21 / a22
    print("\nPART D — board retrodiction of the resample family through THIS instrument")
    print(f"  {'arm':<10}{'member':<7}{'h1 raw d':>10}{'-> scaled':>11}{'x0.4277':>10}"
          f"{'board scaled':>14}{'sign':>7}")
    base = res["a15"]
    retro = {}
    for a, board in (("a21_like", A21), ("a22_like", A22)):
        if a not in res:
            continue
        for mname in ("reach", "fid", "jac", "nmae"):
            d = res[a][mname] - base[mname]
            if mname == "nmae":
                d = -d
            sc = d * INV[mname]
            shr = sc * SHRINK if mname == "reach" else sc
            ok = "OK" if np.sign(sc) == np.sign(board[mname]) else "WRONG"
            retro[(a, mname)] = (d, sc, board[mname], ok)
            print(f"  {a:<10}{mname:<7}{d:>+10.5f}{sc:>+11.5f}{shr:>+10.5f}{board[mname]:>+14.6f}{ok:>7}")
    with open(os.path.join(WORK, "six185.json"), "w") as fh:
        json.dump(dict(known=ka, members=res,
                       frontier={f"{k[0]}|{k[1]}|{k[2]}": v for k, v in frontier.items()},
                       retro={f"{k[0]}|{k[1]}": v for k, v in retro.items()}), fh, default=float)


if __name__ == "__main__":
    main()
