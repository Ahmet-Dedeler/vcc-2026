"""Audit the SHIPPED a12 artifact against the model that predicted it.

Falsifier 2 fired on a12 (predicted raw `expr_mse_unbiased_capped_norm` 0.8882, board
0.9013492). Every term of

    MSE - 1 = ( A - 2X + Kc - rho*claim ) / sumD

except `X` and `sumD` is computable from the shipped file plus `data/challenge/context_*.h5ad`.
This script measures A, B, Kc, C_hat_p, the per-row cap margin and `claim` on the artifact that
produced the board number -- NOT on a re-emission of the policy -- and localises the shortfall.

    uv run python mech_a12_audit.py

Reads `work/{a8,a10,a11,a12}-stats.npz` (produced by `mech_chal_stats.py` from the .vcc files)
and `work/ctrl-stats.npz`. Writes `work/a12_audit.npz`.
"""
from __future__ import annotations

import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from cell_eval2.metrics.delta import _across_pert_budget, _row_weights  # noqa: E402
from cell_eval2.prep import bulk_lognorm_means  # noqa: E402

TS = 50_000.0
CTX = ("A", "B", "C")
B_ANCHOR, SPAN = 0.987862879, 0.950162566
ALPHA_PROD = 1.0796

# board raws, as scored (docs/LOG.md + the arm result files)
BOARD = {"a8": 1.362666960949, "a10": 1.3524148, "a11": 0.9099893, "a12": 0.9013492,
         "a12_overall": 0.16798741}
# nominal VCC_ALPHA each arm was built at
ALPHA = {"a8": 1.0796, "a10": 1.0796, "a11": 0.838, "a12": 1.0108}

CTRL = np.load(f"{ROOT}/experiments/096-mse-floor/work/ctrl-stats.npz", allow_pickle=True)


def cp_proxies():
    """The two `C_p` proxies already on disk, as per-context scalars.

    `ntc`  -- delete-1 jackknife of ONE ntc construct's 400 cells, 16 constructs per context.
              This is the one the published cap 0.246529 was computed from.
    `rand` -- 400 cells drawn at random from the whole 18,400-cell control pool, 12 draws.
    """
    out = {}
    for tag, f in (("ntc", "cp-proxy-ntc.npz"), ("rand", "cp-proxy.npz")):
        z = np.load(f"{ROOT}/experiments/096-mse-floor/work/{f}")
        out[tag] = {c: (float(z[c].mean()), float(z[c].std()), z[c].size) for c in CTX}
    return out


def measure(stats_path, cp):
    """Every computable term of the identity, per context, for one shipped artifact.

    `cp` is `{ctx: C_p}`. Returns per-context dicts plus the per-row jackknife vector.
    """
    z = np.load(stats_path, allow_pickle=True)
    perts = z["perts"].astype(str)
    genes = z["genes"].astype(str)
    assert np.array_equal(genes, CTRL["genes"].astype(str)), "gene axis mismatch"
    G = genes.size
    gpos = {g: i for i, g in enumerate(genes)}
    pm_all = bulk_lognorm_means(z["sums"].astype(np.float64), TS)
    rows = []
    for c in CTX:
        sel = np.flatnonzero(np.char.startswith(perts, f"{c}|"))
        tgt = np.asarray([p.split("|", 1)[1] for p in perts[sel]])
        bctrl = bulk_lognorm_means(CTRL[f"sum_{c}"][None, :].astype(np.float64), TS)[0]
        cols = {k: gpos[t] for k, t in enumerate(tgt) if t in gpos}
        w = _row_weights(sel.size, G, cols)
        pm = pm_all[sel]
        dh = pm - bctrl[None, :]
        a = np.einsum("ij,ij->i", dh, dh)
        for k, col in cols.items():
            a[k] -= dh[k, col] ** 2
        Chat = z["jk"][sel]
        Cp = cp[c]
        capped = np.minimum(Chat, Cp)
        claim = float((capped * w).sum())
        B = _across_pert_budget(pm, cols, w)
        rows.append(dict(
            ctx=c, n=int(sel.size), resolved=len(cols), G=G,
            A=float((a * w).sum()), B=B,
            Kc=float((np.full(sel.size, float(CTRL[f"jk_{c}"])) * w).sum()),
            Cctrl=float(CTRL[f"jk_{c}"]),
            Chat=Chat, Cp=Cp, claim=claim,
            claim_uncapped=float((Chat * w).sum()),
            K=float((np.full(sel.size, Cp) * w).sum()),
            cap_binds=int((Chat > Cp).sum()),
            rho=min(1.0, B / claim) if claim > 0 else 1.0,
        ))
    return rows


def panel(rows, keys=("A", "B", "Kc", "claim", "claim_uncapped", "K")):
    m = {k: float(np.mean([r[k] for r in rows])) for k in keys}
    # ⚠️ #348's budget and #247's cap are formed inside ONE scoring run, and the challenge is
    # scored PER CONTEXT (three anchor sets in `anchor_version`; the We527 two-raws-one-score
    # proof in MSE-MECHANISM §(d)). So `min(claim, B)` is taken per context and averaged --
    # NOT applied to the context means, which is what the 096 model did. The two agree only
    # while every context is on the same side of the cliff.
    m["rho_claim"] = float(np.mean([min(r["claim"], r["B"]) for r in rows]))
    m["rho_claim_meanform"] = min(m["claim"], m["B"])
    m["rho_min"] = float(np.min([min(1.0, r["B"] / r["claim"]) for r in rows]))
    return m


def mse_of(m, X, sD):
    """`(★)` with the metric's own two branches, from measured context-mean terms."""
    return 1.0 + (m["A"] + m["Kc"] - m["rho_claim"] - 2 * X) / sD


def score_mse(raw):
    return min(max((B_ANCHOR - raw) / SPAN, 0.0), 1.0)


def solve_pair(m1, y1, m2, y2):
    """(X, sumD) from two arms' measured terms and their two board raws.

    `X` here is the arm-INDEPENDENT alignment at VCC_ALPHA = 1.0796; each arm enters with its
    own amplitude factor `s = alpha/1.0796`, since `X(s) = s*X`.
    """
    def row(m, y):
        c = m["A"] + m["Kc"] - m["rho_claim"]
        return (-2 * m["s"], -(y - 1.0), -c)
    k1, d1, v1 = row(m1, y1)
    k2, d2, v2 = row(m2, y2)
    X, sD = np.linalg.solve(np.array([[k1, d1], [k2, d2]]), np.array([v1, v2]))
    return float(X), float(sD)


def solve_ls(ms, ys):
    """Least-squares (X, sumD) over >2 arms: (y-1)*sD + 2 s X = A + Kc - min(claim, B)."""
    Amat = np.array([[-2 * m["s"], -(y - 1.0)] for m, y in zip(ms, ys)])
    b = np.array([-(m["A"] + m["Kc"] - m["rho_claim"]) for m in ms])
    sol, *_ = np.linalg.lstsq(Amat, b, rcond=None)
    return float(sol[0]), float(sol[1])


def main():
    cps = cp_proxies()
    print("C_p proxies (per-context scalar, cells with 400 cells / 20k median depth):")
    for tag in ("ntc", "rand"):
        d = cps[tag]
        K = float(np.mean([300 * d[c][0] / 18532 for c in CTX]))
        print(f"  {tag:5s}  " + "  ".join(f"{c}={d[c][0]:7.4f}+-{d[c][1]:.4f}(n={d[c][2]})"
                                          for c in CTX) + f"   -> K = {K:.6f}")

    arms = {}
    for arm in ("a8", "a10", "a11", "a12"):
        p = f"{ROOT}/experiments/096-mse-floor/work/{arm}-stats.npz"
        try:
            rows = measure(p, {c: cps["ntc"][c][0] for c in CTX})
        except FileNotFoundError:
            print(f"  (missing {p}, skipping {arm})")
            continue
        m = panel(rows)
        m["s"] = ALPHA[arm] / ALPHA_PROD
        m["rows"] = rows
        arms[arm] = m

    print("\n===== measured on each SHIPPED artifact, real challenge panel, context-mean =====")
    print(f"{'arm':>5}{'alpha':>8}{'s_nom':>8}{'s_meas':>8}{'A':>10}{'B':>10}{'A-B':>9}"
          f"{'Kc':>9}{'claim':>9}{'uncapd':>9}{'K':>9}{'B/claim':>9}{'rho':>7}{'board':>9}")
    a8 = arms.get("a8")
    for arm, m in arms.items():
        s_meas = np.sqrt(m["A"] / a8["A"]) if a8 else float("nan")
        print(f"{arm:>5}{ALPHA[arm]:>8.4f}{m['s']:>8.4f}{s_meas:>8.4f}{m['A']:>10.6f}"
              f"{m['B']:>10.6f}{m['A']-m['B']:>9.6f}{m['Kc']:>9.6f}{m['claim']:>9.6f}"
              f"{m['claim_uncapped']:>9.4f}{m['K']:>9.6f}{m['B']/m['claim']:>9.4f}"
              f"{min(1.0, m['B']/m['claim']):>7.4f}{BOARD[arm]:>9.4f}")

    print("\n===== the two #348 branches, PER CONTEXT (this is where the model broke) =====")
    print(f"{'arm':>5}{'ctx':>4}{'A_c':>10}{'B_c':>10}{'claim_c':>10}{'Kc_c':>9}"
          f"{'rho_c':>8}{'min(claim,B)':>13}{'num_c':>10}")
    for arm, m in arms.items():
        for r in m["rows"]:
            rc = min(r["claim"], r["B"])
            print(f"{arm:>5}{r['ctx']:>4}{r['A']:>10.6f}{r['B']:>10.6f}{r['claim']:>10.6f}"
                  f"{r['Kc']:>9.6f}{min(1.0, r['B']/r['claim']):>8.4f}{rc:>13.6f}"
                  f"{r['A']+r['Kc']-rc:>10.6f}")
        print(f"{'':>5}{'mean':>4}{m['A']:>10.6f}{m['B']:>10.6f}{m['claim']:>10.6f}"
              f"{m['Kc']:>9.6f}{m['rho_min']:>8.4f}{m['rho_claim']:>13.6f}"
              f"{m['A']+m['Kc']-m['rho_claim']:>10.6f}"
              f"   <- per-context; the 096 model used min of the MEANS = "
              f"{m['rho_claim_meanform']:.6f}"
              f" (diff {m['rho_claim']-m['rho_claim_meanform']:+.6f})")

    print("\n===== a12: the cap, per context and per row =====")
    m12 = arms["a12"]
    for r in m12["rows"]:
        ch = r["Chat"]
        ratio = ch / r["Cp"]
        marg = ch - r["Cp"]
        print(f"  ctx {r['ctx']}  n={r['n']} resolved={r['resolved']}  "
              f"C_p(proxy)={r['Cp']:.4f}")
        print(f"      C_hat_p: mean {ch.mean():8.3f}  median {np.median(ch):8.3f}  "
              f"min {ch.min():8.3f}  p5 {np.percentile(ch,5):8.3f}  "
              f"p95 {np.percentile(ch,95):8.3f}  max {ch.max():9.3f}")
        print(f"      C_hat/C_p: mean {ratio.mean():7.3f}  median {np.median(ratio):7.3f}  "
              f"min {ratio.min():7.3f}  p5 {np.percentile(ratio,5):7.3f}")
        print(f"      cap binds (C_hat > C_p) on {int((marg>0).sum())}/{r['n']} rows; "
              f"rows below cap: {int((marg<=0).sum())}  "
              f"forfeited = sum w*(C_p - C_hat)+ = "
              f"{float((np.maximum(-marg,0)).sum())/ (r['G']-1) :.6f}")
        print(f"      claim {r['claim']:.6f}  vs K (= all rows at cap) {r['K']:.6f}  "
              f"-> {100*r['claim']/r['K']:.2f}% of cap")

    print("\n===== board reconciliation with the a8+a11 solve (X, sumD as shipped) =====")
    bs = np.load(f"{ROOT}/experiments/096-mse-floor/work/board_solve2.npz")
    X0, sD0 = float(bs["X"]), float(bs["sumD"])
    print(f"  a8+a11 re-anchored solve as used for a12's prediction: X={X0:.6f} sumD={sD0:.6f}")
    print(f"{'arm':>5}{'pred raw':>11}{'board raw':>11}{'resid':>10}{'num pred':>11}"
          f"{'num board':>11}{'d num':>10}")
    for arm, m in arms.items():
        mm = dict(m); mm["A"] = m["A"]; mm["Kc"] = m["Kc"]
        pred = 1.0 + (m["A"] + m["Kc"] - m["rho_claim"] - 2 * m["s"] * X0) / sD0
        npred = m["A"] + m["Kc"] - m["rho_claim"] - 2 * m["s"] * X0
        nboard = (BOARD[arm] - 1.0) * sD0
        print(f"{arm:>5}{pred:>11.6f}{BOARD[arm]:>11.6f}{pred-BOARD[arm]:>10.6f}"
              f"{npred:>11.6f}{nboard:>11.6f}{nboard-npred:>10.6f}")

    print("\n===== what each single term would have to be, alone, to close a12's residual =====")
    m = arms["a12"]
    nboard = (BOARD["a12"] - 1.0) * sD0
    base = m["A"] + m["Kc"] - m["rho_claim"]
    print(f"  measured A={m['A']:.6f} B={m['B']:.6f} Kc={m['Kc']:.6f} "
          f"claim={m['claim']:.6f} rho*claim={m['rho_claim']:.6f} s={m['s']:.6f}")
    need_rc = base - min(m["claim"], m["B"]) * 0 - (nboard + 2 * m["s"] * X0)
    # rho*claim that reproduces the board at fixed X, sumD:
    rc = m["A"] + m["Kc"] - 2 * m["s"] * X0 - nboard
    print(f"  rho*claim required        = {rc:.6f}   "
          f"(measured {m['rho_claim']:.6f}, cap K={m['K']:.6f})")
    Xreq = (m["A"] + m["Kc"] - m["rho_claim"] - nboard) / (2 * m["s"])
    print(f"  X required                = {Xreq:.6f}   (solve says {X0:.6f}, "
          f"ratio {Xreq/X0:.4f})")
    sDreq = (m["A"] + m["Kc"] - m["rho_claim"] - 2 * m["s"] * X0) / (BOARD["a12"] - 1.0)
    print(f"  sumD required             = {sDreq:.6f}   (solve says {sD0:.6f}, "
          f"ratio {sDreq/sD0:.4f})")
    Areq = nboard + 2 * m["s"] * X0 - m["Kc"] + m["rho_claim"]
    print(f"  A required (B,claim held) = {Areq:.6f}   (measured {m['A']:.6f}, "
          f"ratio {Areq/m['A']:.4f})")
    _ = need_rc

    print("\n===== re-solve (X, sumD) on the points now available =====")
    combos = [("a8+a11 (as shipped)", ["a8", "a11"]),
              ("a8+a12", ["a8", "a12"]),
              ("a11+a12", ["a11", "a12"]),
              ("a10+a11", ["a10", "a11"]),
              ("a10+a12", ["a10", "a12"])]
    print(f"{'points':>22}{'X':>11}{'sumD':>11}{'cos*':>9}   residual on the others")
    sols = {}
    for name, ks in combos:
        if not all(k in arms for k in ks):
            continue
        X, sD = solve_pair(arms[ks[0]], BOARD[ks[0]], arms[ks[1]], BOARD[ks[1]])
        cos = X / np.sqrt(a8["A"] * sD)
        res = []
        for arm in arms:
            if arm in ks:
                continue
            mm = arms[arm]
            pr = 1.0 + (mm["A"] + mm["Kc"] - mm["rho_claim"] - 2 * mm["s"] * X) / sD
            res.append(f"{arm} {pr-BOARD[arm]:+.4f}")
        sols[name] = (X, sD)
        print(f"{name:>22}{X:>11.6f}{sD:>11.6f}{cos:>9.5f}   " + "  ".join(res))

    for name, ks in (("LS a8+a11+a12 (3 pt)", ["a8", "a11", "a12"]),
                     ("LS all four", ["a8", "a10", "a11", "a12"]),
                     # a8, a10 and a11 all claim their OWN measured C_hat_p -- the #247 cap
                     # never binds on them -- so this solve is free of the C_p proxy. a12 is
                     # the first arm where the proxy is load-bearing, so it is held out.
                     ("LS a8+a10+a11 (C_p-free)", ["a8", "a10", "a11"])):
        if not all(k in arms for k in ks):
            continue
        X, sD = solve_ls([arms[k] for k in ks], [BOARD[k] for k in ks])
        cos = X / np.sqrt(a8["A"] * sD)
        res = []
        for k in ks:
            mm = arms[k]
            pr = 1.0 + (mm["A"] + mm["Kc"] - mm["rho_claim"] - 2 * mm["s"] * X) / sD
            res.append(f"{k} {pr-BOARD[k]:+.4f}")
        sols[name] = (X, sD)
        print(f"{name:>22}{X:>11.6f}{sD:>11.6f}{cos:>9.5f}   " + "  ".join(res))

    print("\n===== the a12 prediction error, attributed step by step =====")
    m = arms["a12"]
    a8r = arms["a8"]["rows"]
    s = m["s"]
    # (a) the prediction of record: assumed terms, min applied to the context MEANS
    A_as = s * s * arms["a8"]["A"]
    B_as = s * s * arms["a8"]["B"]
    cl_as = arms["a8"]["K"]                       # `claim` at the #247 cap
    step = [("as registered (assumed A,B; claim at the cap; min of the MEANS; a8+a11 solve)",
             1.0 + (A_as + m["Kc"] - min(cl_as, B_as) - 2 * s * X0) / sD0)]
    # (b) swap in the terms measured on the artifact, still min-of-means
    step.append(("+ A, B, claim measured on the shipped artifact",
                 1.0 + (m["A"] + m["Kc"] - m["rho_claim_meanform"] - 2 * s * X0) / sD0))
    # (c) evaluate #348's budget PER CONTEXT, as the scorer does
    step.append(("+ #348's budget evaluated PER CONTEXT (the model's error)",
                 1.0 + (m["A"] + m["Kc"] - m["rho_claim"] - 2 * s * X0) / sD0))
    # (d) re-solve (X, sumD) with the same correction applied to every anchor point
    Xn, sDn = sols["LS all four"]
    step.append(("+ (X, sumD) re-solved with the correction applied to the anchors too",
                 1.0 + (m["A"] + m["Kc"] - m["rho_claim"] - 2 * s * Xn) / sDn))
    prev = None
    for name, v in step:
        d = "" if prev is None else f"{v - prev:+10.6f}"
        print(f"  {v:>9.6f}{d:>11}   {name}")
        prev = v
    print(f"  {BOARD['a12']:>9.6f}{BOARD['a12']-prev:>11.6f}   BOARD")

    np.savez(f"{ROOT}/experiments/096-mse-floor/work/a12_audit.npz",
             **{f"{a}_{k}": arms[a][k] for a in arms for k in
                ("A", "B", "Kc", "claim", "claim_uncapped", "K", "s")},
             **{f"sol_{n.replace(' ','_')}": np.array(v) for n, v in sols.items()},
             a12_Chat=np.concatenate([r["Chat"] for r in arms["a12"]["rows"]]),
             a12_Cp=np.array([r["Cp"] for r in arms["a12"]["rows"]]))
    return arms, sols


if __name__ == "__main__":
    main()
