"""521 — the per-target row table: every target's fid / jac / reach / nmae rows under each call-set size K, on the
shipped H1 state (435 production stats), for the two emitter variants, at each power regime.

Because each DE member is a mean of per-target rows and a target's row depends only on its own call set and delivered
lfc, any per-target policy K_p is priced EXACTLY by picking rows from this table (fid's cohort drops a target only when
n_pred = n_real = 0; nmae's and reach's cohorts are reference-only).

Regimes (the truth's pools): f 1.0 / 0.5 / 0.35 (PowerSurface: top round(f n_real) by the reference's p), and 'cmp' — a
compressed map n' = round(c n^beta) that matches BOTH board facts: the upper tail at f 0.425 (192 §4.2) and the weak end
at the board's nmae cohort (13-30 % of targets with a gate < 10 → 22 % mid; H1 has 32 % at f 1).
Content when shrunk: 'headcpm' (min(K, 400) head genes by |key|, then block genes by control CPM; 192 §5.3's own
construction), and 'keyord' (today's calls in realized p order = |key| order) as a check.
Variants: b = the emitter keeps steering dropped block genes (lfc pinned); a = production emitter, dropped block genes
take the unsteered lfc (416 prior_off) + the MDE add (ceiling521.lfc_variant_a).

    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/521-strength-nreal/table521.py [seed]
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

import lib521 as L  # noqa: E402
from ceiling521 import content_rank, lfc_variant_a, shrink_calls  # noqa: E402

KGRID = (0, 10, 25, 50, 100, 200, 400, 800, 1200, 1600, 2000, 2200)     # plus "base" = today's call set
REGIMES = ("1.0", "0.5", "0.35", "cmp")


def cmp_counts(S, frac_lt10=0.22, tail_f=0.425):
    """n' = round(c * n^beta): n' = tail_f * 1000 at n = 1000, and a fraction frac_lt10 of targets with a gate < 10."""
    base = ((S["rq"] < L.ALPHA) & np.isfinite(S["rq"]) & ~S["own"]).sum(1).astype(float)
    best = None
    for beta in np.arange(0.40, 1.301, 0.005):
        c = tail_f * 1000.0 / 1000.0 ** beta
        npr = np.where(base > 0, np.round(c * np.power(base, beta)), 0.0)
        fr = float(np.mean(npr < 10))
        if best is None or abs(fr - frac_lt10) < best[0]:
            best = (abs(fr - frac_lt10), beta, c, npr)
    return best[3].astype(int), float(best[1]), float(best[2])


def pools_k(S, k):
    rp, own = S["rp"], S["own"]
    rpx = np.where(np.isfinite(rp) & ~own, rp, np.inf)
    out = np.zeros(rp.shape, bool)
    for i in range(rp.shape[0]):
        if k[i]:
            out[i, np.argsort(rpx[i], kind="stable")[: int(k[i])]] = True
    return out & ~own


def regime_pool(S, reg):
    if reg == "cmp":
        k, beta, c = cmp_counts(S)
        return pools_k(S, k), dict(beta=beta, c=c)
    return L.pools(S, float(reg)), {}


def rows_for(S, call, lfc, pool):
    adj = L.adjud(S)
    rl = S["rl"]
    fr = L.fid_rows(call, lfc, pool, adj, rl)
    jr = L.jac_rows(call, pool)
    nr = L.nmae_rows(lfc, pool, rl)
    rr = np.full(call.shape[0], np.nan)
    for i in range(call.shape[0]):
        rr[i] = L.reach_row(L.model_order_pos(S, i, call[i]), pool[i], adj[i], lfc[i], rl[i])
    return np.vstack([fr, jr, rr, nr])       # members x targets


def main(seed=0):
    S = L.load_state(seed)
    sig = (S["oq"] < L.ALPHA) & ~S["own"]
    ranks = {kind: [content_rank(S, i, kind) for i in range(85)] for kind in ("headcpm", "keyord")}
    tab = {}
    meta = {}
    for reg in REGIMES:
        pool, info = regime_pool(S, reg)
        meta[reg] = dict(info, nreal=pool.sum(1).tolist())
        tab[f"{reg}|base"] = rows_for(S, sig, S["ol"], pool)
        for kind in ("headcpm", "keyord"):
            for K in KGRID:
                call = shrink_calls(S, np.full(85, float(K)), kind, ranks[kind])
                tab[f"{reg}|{kind}|b|{K}"] = rows_for(S, call, S["ol"], pool)
                tab[f"{reg}|{kind}|a|{K}"] = rows_for(S, call, lfc_variant_a(S, call), pool)
        n = pool.sum(1)
        print(f"regime {reg}: n_real median {np.median(n):.0f} mean {n.mean():.0f}, zero {int((n == 0).sum())}, "
              f"gate<10 {int(((pool & np.isfinite(S['rl'])).sum(1) < 10).sum())} {info}", flush=True)
    np.savez(L.WORK / f"table521_s{seed}.npz", keys=np.array(list(tab.keys())), rows=np.stack(list(tab.values())),
             nreal=np.array([meta[r]["nreal"] for r in REGIMES]), regimes=np.array(REGIMES),
             cmp_beta=meta["cmp"]["beta"], cmp_c=meta["cmp"]["c"], kgrid=np.array(KGRID))
    print(f"wrote table521_s{seed}.npz with {len(tab)} row blocks")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
