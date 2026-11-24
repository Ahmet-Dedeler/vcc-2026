"""189 Phase 2.2/2.3 — the prior as a CALL-SET filter (`fid`, `jac`) and as a weight on `q_w`.

USE 2 — call-set composition.  `fid_p = k / max(n_pred_p, n_real_p)`.  I own WHICH genes go in
the call set at a given size; 185's INBOX-1 owns the SIZE `n_pred_p` and is not duplicated here.
The transferable quantity is therefore the **sign precision of the top-n call set as a function
of n and of the key** -- `fid`'s numerator per called gene.  Absolute `fid`/`jac` are printed but
carry AGENTS.md §11.1 rule 9: `n_real` is 89 median here against a right-skewed challenge
distribution, so only the ordering of rows is readable.

USE 3 — `q_w`.  187 §5: `nmae = 1 - 2 q_w W_hi - (2 q_lo - 1) K_lo + K_hi`, and with 99% of
genes below the kink this is `nmae = 1 - (2 q_w - 1) K`, `K = sum|p|/sum|t|`.  So

        d nmae = -2 K d q_w        at fixed total predicted mass.

The counterfactual measured here is: if 188 re-shaped per-gene magnitudes by `u(theta)`, at
CONSTANT total |p| over the pool (so K does not move and this stays a direction question, not
an amplitude one), what `q_w` is achievable?  Magnitudes are not changed by this script.
"""
from __future__ import annotations

import json
import os

import numpy as np

from common189 import WORK
from harness189 import Surface189, check

NS = [10, 25, 50, 100, 200, 400, 800, 1600, 3200]
LAMS = [-1.0, -0.5, 0.0, 0.25, 0.5, 1.0, 2.0, 4.0]


def zs(x):
    m = np.isfinite(x)
    o = np.zeros_like(x, dtype=np.float64)
    o[m] = (x[m] - x[m].mean()) / max(x[m].std(), 1e-12)
    return o


def topn_call(T, key, n):
    """Per target, the top-n genes by `key` among the CPM-gated candidates."""
    call = np.zeros(T.H.shape, bool)
    cand = np.flatnonzero(T.gate)
    for p in range(len(T.perts)):
        k = np.asarray(key[p])[cand]
        m = min(n, cand.size)
        call[p, cand[np.argpartition(-k, m - 1)[:m]]] = True
    return call


def precision(T, call):
    """Sign precision over the ADJUDICABLE called genes, pooled and macro."""
    S = T.S
    num = den = 0.0
    macro = []
    for p in range(len(T.perts)):
        idx = np.flatnonzero(call[p] & S.adj[p])
        if idx.size == 0:
            continue
        ok = (np.sign(S.pred_lfc[p, idx]) == np.sign(S.real_lfc[p, idx])).astype(float)
        num += ok.sum()
        den += idx.size
        macro.append(float(ok.mean()))
    return num / max(den, 1), float(np.mean(macro))


def main() -> None:
    T = Surface189()
    check(T)
    th = T.theta
    ze = zs(T.b_ctrl)
    base = np.log(np.maximum(np.abs(T.H), 1e-30))
    rng = np.random.default_rng(20260903)

    keys = {
        "|H| (incumbent)": np.abs(T.H),
        "|H| x exp(+0.50 theta)": np.exp(base + 0.5 * th),
        "|H| x exp(+1.00 theta)": np.exp(base + 1.0 * th),
        "|H| x exp(+0.50 SHUFFLED theta)": np.exp(base + 0.5 * rng.permutation(th)),
        "|H| x exp(+0.50 z_expr)": np.exp(base + 0.5 * ze),
        "|H| x exp(-0.50 z_expr)": np.exp(base - 0.5 * ze),
    }

    print(f"\n{'=' * 112}\nUSE 2 — call-set sign PRECISION of the top-n genes (pooled / macro)")
    print(f"{'key':<34}" + "".join(f"{n:>9}" for n in NS))
    prec = {}
    for name, k in keys.items():
        row = [precision(T, topn_call(T, k, n)) for n in NS]
        prec[name] = [(float(a), float(b)) for a, b in row]
        print(f"{name:<34}" + "".join(f"{a:>9.4f}" for a, _ in row))
    print(f"{'  (macro variant, |H|)':<34}"
          + "".join(f"{b:>9.4f}" for _, b in prec["|H| (incumbent)"]))
    print(f"{'  (macro, |H| x exp(+0.50 theta))':<34}"
          + "".join(f"{b:>9.4f}" for _, b in prec["|H| x exp(+0.50 theta)"]))

    print(f"\n{'=' * 112}\nUSE 2b — `fid` and `jac` at the same call sets "
          f"(rule 9: n_real median {np.median(T.n_real):.0f} here, ORDERING ONLY)")
    print(f"{'key':<34}" + "".join(f"{n:>9}" for n in NS))
    fj = {}
    for name, k in keys.items():
        f = [T.S.fid(call=topn_call(T, k, n)) for n in NS]
        j = [T.S.jac(call=topn_call(T, k, n)) for n in NS]
        fj[name] = dict(fid=list(map(float, f)), jac=list(map(float, j)))
        print(f"{name:<34}" + "".join(f"{x:>9.4f}" for x in f))
    print(f"{'  jac, |H|':<34}" + "".join(f"{x:>9.4f}" for x in fj["|H| (incumbent)"]["jac"]))
    print(f"{'  jac, |H| x exp(+0.50 theta)':<34}"
          + "".join(f"{x:>9.4f}" for x in fj["|H| x exp(+0.50 theta)"]["jac"]))

    # ---- USE 3 -----------------------------------------------------------------------
    K0 = T.mass_ratio()
    q0 = T.qw()
    print(f"\n{'=' * 112}\nUSE 3 — achievable `q_w` under a per-gene weight u(theta), at "
          f"CONSTANT total |p| mass")
    print(f"    baseline: q_w {q0:.5f}   K {K0:.5f}   count accuracy "
          f"{float(T.correct[T.pool].mean()):.5f}")
    print(f"    pricing (187 §5, all-k<1 regime): d nmae_raw = -2 K d q_w ; "
          f"scaled = -1.6532 d raw ; overall = scaled/6")
    print(f"{'weight u':<30}{'q_w':>10}{'d q_w':>10}{'K':>9}{'d nmae raw':>12}"
          f"{'d scaled':>10}{'d overall':>11}")
    qrows = []
    for tag, w in ([(f"exp({l:+.2f} theta)", np.exp(l * th)) for l in LAMS]
                   + [("exp(+0.50 SHUFFLED theta)", np.exp(0.5 * rng.permutation(th))),
                      ("exp(+0.50 z_expr)", np.exp(0.5 * ze)),
                      ("exp(-0.50 z_expr)", np.exp(-0.5 * ze)),
                      ("raw r_g", np.nan_to_num(T.raw, nan=float(np.nanmean(T.raw))))]):
        w = np.asarray(w, dtype=np.float64)
        Kw = T.mass_ratio(w)
        w = w * (K0 / Kw)                       # renormalise: total |p| mass unchanged
        q = T.qw(w)
        dn = -2.0 * K0 * (q - q0)
        qrows.append(dict(weight=tag, qw=q, dq=q - q0, K=float(T.mass_ratio(w)),
                          dnmae=dn, dscaled=-1.6532 * dn, doverall=-1.6532 * dn / 6))
        print(f"{tag:<30}{q:>10.5f}{q - q0:>+10.5f}{T.mass_ratio(w):>9.5f}{dn:>+12.5f}"
              f"{-1.6532 * dn:>+10.5f}{-1.6532 * dn / 6:>+11.6f}")

    # the ceiling: perfect per-gene weighting (oracle on the gene axis)
    acc_g = np.zeros(len(T.genes))
    cnt_g = np.zeros(len(T.genes))
    for p in range(len(T.perts)):
        idx = np.flatnonzero(T.pool[p])
        acc_g[idx] += T.correct[p, idx]
        cnt_g[idx] += 1
    gacc = np.where(cnt_g > 0, acc_g / np.maximum(cnt_g, 1), 0.5)
    for name, w in (("ORACLE gene-level: u = 1[acc>0.5]", (gacc > 0.5).astype(float) + 1e-6),
                    ("ORACLE gene-level: u = exp(8*(acc-0.5))", np.exp(8 * (gacc - 0.5)))):
        Kw = T.mass_ratio(w)
        w2 = w * (K0 / Kw)
        q = T.qw(w2)
        dn = -2.0 * K0 * (q - q0)
        qrows.append(dict(weight=name, qw=q, dq=q - q0, dnmae=dn,
                          dscaled=-1.6532 * dn, doverall=-1.6532 * dn / 6))
        print(f"{name:<30}{q:>10.5f}{q - q0:>+10.5f}{T.mass_ratio(w2):>9.5f}{dn:>+12.5f}"
              f"{-1.6532 * dn:>+10.5f}{-1.6532 * dn / 6:>+11.6f}")

    json.dump(dict(ns=NS, precision=prec, fid_jac=fj, q0=q0, K0=K0, qrows=qrows),
              open(os.path.join(WORK, "members189.json"), "w"), indent=2)
    print("\nwrote work/members189.json")


if __name__ == "__main__":
    main()
