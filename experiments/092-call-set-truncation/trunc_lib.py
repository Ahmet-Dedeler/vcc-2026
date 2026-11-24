"""Shared precompute for the truncation counterfactual. Ranking + official scaling."""
import os, numpy as np

S = os.path.dirname(os.path.abspath(__file__))
A = 0.05
# official anchors (AGENTS.md sec 4, midpoints; identical to step5_fine.py)
F0, FS = 0.5135, 0.3
J0, JS = 0.029, 0.37
R0, RS = 0.072, 0.968      # reach baseline / replicate anchor


def load():
    d = np.load(os.path.join(S, "de_tables.npz"), allow_pickle=True)
    return (d["genes"].astype(str), d["perts"].astype(str),
            d["real_p"], d["real_padj"], d["real_lfc"],
            d["pred_p"], d["pred_padj"], d["pred_lfc"])


def precompute():
    """Per-pert arrays over our significance-ranked call list, plus reach inputs."""
    genes, perts, rp, rq, rl, pp, pq, pl = load()
    G = len(genes); P = len(perts)
    gpos = {g: i for i, g in enumerate(genes)}
    fr = np.empty(G, np.int64); fr[np.argsort(genes)] = np.arange(G)
    out = []
    for i, t in enumerate(perts):
        ok = np.ones(G, bool); j = gpos.get(t)
        if j is not None: ok[j] = False
        rs = ok & (rq[i] < A)                       # reference significant set R
        ps = ok & (pq[i] < A)                       # our significant set (pre-truncation)
        ind = ok & (rl[i] != 0) & ~np.isnan(rl[i])  # adjudicable
        idx = np.flatnonzero(ps)
        # cell_eval2 ranking key: p_adj asc, p asc, |lfc| desc, feature asc
        key = np.lexsort((fr[idx], -np.abs(pl[i][idx]), pp[i][idx], pq[i][idx]))
        r = idx[key]
        adj = ind[r]
        hit = adj & (pl[i][r] != 0) & (np.sign(pl[i][r]) == np.sign(rl[i][r]))
        # rank position of every gene inside our call list (G-length, -1 if never called)
        rank_of = np.full(G, -1, np.int64); rank_of[r] = np.arange(r.size)

        # ---- reach inputs: the REFERENCE's significant pool ----
        gsel = np.flatnonzero(rs)
        r_padj = pq[i][gsel]
        r_agree = ((pl[i][gsel] != 0) & ~np.isnan(pl[i][gsel]) &
                   (np.sign(pl[i][gsel]) == np.sign(rl[i][gsel])))
        r_rank = rank_of[gsel]                       # -1 => not in our call list at all
        r_absl = np.abs(np.nan_to_num(pl[i][gsel]))
        r_feat = fr[gsel]

        out.append(dict(
            t=t, n_real=int(rs.sum()), n_all=int(r.size),
            cum_adj=np.cumsum(adj), cum_hit=np.cumsum(hit), cum_in=np.cumsum(rs[r]),
            r_padj=r_padj, r_agree=r_agree, r_rank=r_rank, r_absl=r_absl, r_feat=r_feat,
            real_lfc_R=rl[i][gsel], pred_lfc_R=np.nan_to_num(pl[i][gsel]),
        ))
    return out


def fid_jac(pre, m_vec):
    """fid & jac per pert for a per-pert call-set size m_vec (None entry => no truncation)."""
    F = np.empty(len(pre)); J = np.empty(len(pre)); M = np.empty(len(pre))
    for i, s in enumerate(pre):
        nr, na = s["n_real"], s["n_all"]
        m = na if m_vec[i] is None else int(min(max(round(m_vec[i]), 0), na))
        M[i] = m
        if m <= 0:
            F[i] = 0.0; J[i] = 1.0 if nr == 0 else 0.0; continue
        npd = int(s["cum_adj"][m - 1]); k = int(s["cum_hit"][m - 1]); it = int(s["cum_in"][m - 1])
        F[i] = k / max(npd, nr) if max(npd, nr) else np.nan
        u = nr + m - it
        J[i] = 1.0 if u == 0 else it / u
    return F, J, M


def scaled(f, j):
    return (f - F0) / FS, (j - J0) / JS


def reach(pre, m_vec, mode="label"):
    """reach per pert.

    mode='label'  : truncation only relabels significance; p_adj/lfc untouched.
                    order = (called-first, then p_adj asc).
    mode='tail'   : truncated genes are de-significated (p_adj -> 1) so they sort to the
                    back; their lfc / sign is preserved.
    mode='null'   : truncated genes are de-significated AND their effect is nulled, so they
                    can never sign-agree. Worst-case implementation.
    """
    out = np.empty(len(pre))
    for i, s in enumerate(pre):
        nr = s["n_real"]
        if nr == 0: out[i] = 0.0; continue
        na = s["n_all"]
        m = na if m_vec[i] is None else int(min(max(round(m_vec[i]), 0), na))
        rr = s["r_rank"]
        called = (rr >= 0) & (rr < m)
        agree = s["r_agree"].copy()
        if mode == "null":
            agree = agree & called
        # sort: called-block first, then by p_adj asc (feature asc as final tiebreak)
        order = np.lexsort((s["r_feat"], s["r_padj"], ~called))
        pur = np.cumsum(agree[order]) / np.arange(1, nr + 1)
        good = np.flatnonzero(pur >= 0.9)
        out[i] = (good[-1] + 1) / nr if good.size else 0.0
    return out


def nmae(pre, m_vec, mode="label"):
    """Normalized MAE of log2FC over the reference's DE genes.
    raw = mean_p [ mean_{g in R_p} |pred_lfc - real_lfc| / mean_{g in R_p} |real_lfc| ]"""
    out = []
    for i, s in enumerate(pre):
        nr = s["n_real"]
        if nr == 0: continue
        na = s["n_all"]
        m = na if m_vec[i] is None else int(min(max(round(m_vec[i]), 0), na))
        rr = s["r_rank"]; called = (rr >= 0) & (rr < m)
        p = s["pred_lfc_R"].copy()
        if mode == "null": p = np.where(called, p, 0.0)
        den = np.abs(s["real_lfc_R"]).mean()
        out.append(np.abs(p - s["real_lfc_R"]).mean() / den if den > 0 else np.nan)
    return float(np.nanmean(out))
