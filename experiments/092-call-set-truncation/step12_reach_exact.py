"""Step 12: reach + nmae computed EXACTLY per cell_eval2 0.16.0 source.

reach (direction.py:_purity_curve/_k_star/_reference_stats/de_direction_reach):
  pool     = reference rows with p_adj_real < alpha, target's own gene removed
  order    = [_sig_pred desc, rank_p_adj asc, rank_p_value asc, abs_lfc_pred desc, feature asc]
  in_denom = real lfc not null/nan/zero ; match = in_denom & pred lfc defined & signs equal
  n_denom  = cumsum(in_denom) ; purity = cumsum(match)/n_denom ; k = n_denom
  k*       = MAX k with purity >= 0.9  (deepest, purity is not monotone)
  R        = k* / n_conf,  n_conf = #(p_adj_real < alpha) [significance only]; NaN if 0
nmae (de.py:de_lfc_nmae): mean|lfc_pred-lfc_real| / mean|lfc_real| over {p_adj_real<alpha
  & lfc_real finite & own gene removed}; pred null/non-finite -> 0.0; target dropped if
  n_gate < 10 or den == 0.
"""
import os, numpy as np, trunc_lib as T

S = os.path.dirname(os.path.abspath(__file__))
genes, perts, rp, rq, rl, pp, pq, pl_ = T.load()
A = 0.05; G = len(genes); P = len(perts)
gpos = {g: i for i, g in enumerate(genes)}
fr = np.empty(G, np.int64); fr[np.argsort(genes)] = np.arange(G)

# ---- build the reference pool per target, with every sort key ----
pool = []
for i, t in enumerate(perts):
    ok = np.ones(G, bool); j = gpos.get(t)
    if j is not None: ok[j] = False
    gsel = np.flatnonzero(ok & (rq[i] < A))                 # n_conf rows
    rlv = rl[i][gsel]; plv = pl_[i][gsel]
    in_denom = (~np.isnan(rlv)) & (rlv != 0)
    committed = (~np.isnan(plv)) & (plv != 0)
    match = in_denom & committed & (np.sign(plv) == np.sign(rlv))
    # our own call-list rank (p_adj asc, p asc, |lfc| desc, feature asc) over our sig set
    ps = ok & (pq[i] < A)
    idx = np.flatnonzero(ps)
    key = np.lexsort((fr[idx], -np.abs(pl_[i][idx]), pp[i][idx], pq[i][idx]))
    rank_of = np.full(G, -1, np.int64); rank_of[idx[key]] = np.arange(idx.size)
    pool.append(dict(t=t, n_conf=gsel.size, n_all=int(idx.size),
                     padj=pq[i][gsel], pval=pp[i][gsel],
                     absl=np.abs(np.nan_to_num(pl_[i][gsel], nan=-np.inf)),
                     feat=fr[gsel], in_denom=in_denom, match=match,
                     rank=rank_of[gsel], rlv=rlv, plv=plv))

unadj = sum(int((~s["in_denom"]).sum()) for s in pool)
print(f"unadjudicable rows among {sum(s['n_conf'] for s in pool)} reference-significant "
      f"genes: {unadj}   (source notes 0 on all six official lines)")


def reach(m_vec, mode="tail", seed=None):
    """mode: 'tail' = truncated genes de-significated, p_adj order among them preserved;
             'shuffle' = de-significated AND their relative order randomised;
             'null' = de-significated AND their effect nulled (can never match)."""
    rng = np.random.default_rng(seed)
    out = np.empty(P)
    for i, s in enumerate(pool):
        nc = s["n_conf"]
        if nc == 0: out[i] = np.nan; continue
        m = s["n_all"] if m_vec[i] is None else int(min(max(round(m_vec[i]), 0), s["n_all"]))
        called = (s["rank"] >= 0) & (s["rank"] < m)
        mt = s["match"] & called if mode == "null" else s["match"]
        if mode == "shuffle":
            tie = rng.random(nc)
            order = np.lexsort((s["feat"], tie, ~called))
        else:
            order = np.lexsort((s["feat"], -s["absl"], s["pval"], s["padj"], ~called))
        nd = np.cumsum(s["in_denom"][order]); nm = np.cumsum(mt[order])
        with np.errstate(invalid="ignore", divide="ignore"):
            pur = np.where(nd > 0, nm / np.maximum(nd, 1), np.nan)
        good = np.flatnonzero((nd > 0) & (pur >= 0.9))
        out[i] = (nd[good[-1]] / nc) if good.size else 0.0
    return out


def nmae(m_vec, mode="tail"):
    vals = []
    for i, s in enumerate(pool):
        gate = np.isfinite(s["rlv"])
        if gate.sum() < 10: continue
        m = s["n_all"] if m_vec[i] is None else int(min(max(round(m_vec[i]), 0), s["n_all"]))
        called = (s["rank"] >= 0) & (s["rank"] < m)
        p = np.where(np.isfinite(s["plv"]), s["plv"], 0.0)
        if mode == "null": p = np.where(called, p, 0.0)
        r = s["rlv"][gate]; p = p[gate]
        den = np.abs(r).mean()
        if den == 0 or not np.isfinite(den): continue
        vals.append(np.abs(p - r).mean() / den)
    return float(np.mean(vals))


base = np.nanmean(reach([None] * P))
bn = nmae([None] * P)
NB, NRp = 1.0013, 0.400
print(f"\nBASELINE reach = {base:.6f}  scaled {(base - T.R0)/T.RS:+.4f}")
print(f"BASELINE nmae  = {bn:.6f}  scaled {(bn - NB)/(NRp - NB):+.4f}")

print("\n" + "=" * 96)
print("reach / nmae vs TRUNCATION LEVEL   (exact cell_eval2 semantics)")
print("=" * 96)
print(f"{'m policy':<20}{'reach tail':>12}{'reach shuf':>12}{'reach null':>12}"
      f"{'d6 reach(tail)':>16}{'nmae tail':>11}{'nmae null':>11}{'d6 nmae(null)':>15}")
pols = [("baseline", [None] * P)] + \
       [(f"constant c={c}", [float(c)] * P) for c in (400, 700, 850, 1000, 1500, 2000)] + \
       [("oracle r=1.2", [1.2 * s["n_conf"] for s in pool]),
        ("alpha 1e-8", None)]
mv_a = None
for lab, mv in pols:
    if mv is None: continue
    rt = np.nanmean(reach(mv, "tail"))
    rs_ = np.nanmean([np.nanmean(reach(mv, "shuffle", seed=k)) for k in range(3)])
    rn = np.nanmean(reach(mv, "null"))
    nt = nmae(mv, "tail"); nn = nmae(mv, "null")
    print(f"{lab:<20}{rt:>12.6f}{rs_:>12.6f}{rn:>12.6f}"
          f"{(rt - base)/T.RS/6:>+16.4f}{nt:>11.6f}{nn:>11.6f}"
          f"{((nn - NB)/(NRp - NB) - (bn - NB)/(NRp - NB))/6:>+15.4f}")
