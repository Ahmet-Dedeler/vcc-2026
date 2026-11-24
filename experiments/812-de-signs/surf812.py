"""812 — shared loader for the two DE-H1 surfaces (h1-test, h1-trainval) on the production emitter's realized stats.

State per surface and seed (all on the scorer's 5-CPM universe, the reference's `keep` axis):
  rq, rp, rl     the truth's BH q, p, lfc (323's / 618's reference.npz)
  oq, op, ol     our emitted q, p, lfc (production emitter rc_emit_k, CLEAN4 q80 block prior, coarse MDE, zs;
                 h1test = 435's `..._dnk250zs`, h1tv = 618's g1 `..._zr-5_..._dnk250zs`, i.e. the hook-free zs arm)
  key            the cap15 design key gated at control CPM > 10 (lib521's reconstruction of run435's `eff_de`)
  head, block    the emitter's prescribed call set: top 400 |key|, then top 2,200 emit_ref_lin among the rest
  cpm            control CPM per gene (from b_ctrl), c4      CLEAN4 P_down per gene (NaN where uncovered)
Read-only on every input; nothing here writes.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
ALPHA = 0.05
N_HEAD, N_CPM = 400, 2200
SURF = {
    "h1test": dict(
        ref=ROOT / "experiments/323-current-family-validation/work/reference.npz",
        eff=ROOT / "experiments/323-current-family-validation/work/effects.npz",
        st=ROOT / "experiments/435-zs-on-coarse/work/h1/stats-base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz"),
    "h1tv": dict(
        ref=ROOT / "experiments/618-h1-trainval-surface/work/ref_tv/reference.npz",
        eff=ROOT / "experiments/618-h1-trainval-surface/work/tv/effects.npz",
        st=ROOT / "experiments/618-h1-trainval-surface/work/g1/stats-base_s{seed}_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz"),
}
C4 = ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz"


def sha(p, n=1 << 20):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while b := f.read(n):
            h.update(b)
    return h.hexdigest()


def load(surface: str, seed: int = 0, stats_path: str | None = None):
    P = SURF[surface]
    compact = HERE / f"work/stage/ref_{surface}.npz"   # the 8 small arrays of reference.npz (EC2 runs; same values)
    ref = np.load(compact if compact.exists() else P["ref"])
    genes_full = ref["genes"].astype(str)
    keep = ref["keep"].astype(bool)
    targets = ref["targets"].astype(str)
    genes = genes_full[keep]
    rq, rp, rl = (np.asarray(ref[k], float) for k in ("q", "p", "lfc"))
    ref_lin = np.asarray(ref["emit_ref_lin"], float)
    b_ctrl = np.asarray(ref["b_ctrl"], float)
    stp = Path(stats_path) if stats_path else Path(str(P["st"]).format(seed=seed))
    st = np.load(stp)
    assert list(st["targets"].astype(str)) == list(targets), "stats targets differ from the reference's"
    oq, op, ol = (np.asarray(st[k], float) for k in ("q", "p", "lfc"))
    e = np.load(P["eff"])
    assert list(e["targets"].astype(str)) == list(targets)
    pos = {g: i for i, g in enumerate(e["genes"].astype(str))}
    cols = np.array([pos.get(g, -1) for g in genes_full])
    eff = np.zeros((len(targets), len(genes_full)))
    eff[:, cols >= 0] = np.asarray(e["cap15"], float)[:, cols[cols >= 0]]
    ex = np.expm1(b_ctrl)
    cpm_full = 1e6 * ex / ex.sum()
    gate10 = (cpm_full > 10.0) & (b_ctrl > 0)
    key = np.where(gate10[None, :], eff, 0.0)[:, keep]
    T, G = key.shape
    own = np.equal.outer(targets, genes)
    head = np.zeros((T, G), bool)
    block = np.zeros((T, G), bool)
    for i in range(T):
        head[i, np.argsort(-np.abs(key[i]), kind="stable")[:N_HEAD]] = True
        block[i, np.argsort(-np.where(head[i], -np.inf, ref_lin), kind="stable")[:N_CPM]] = True
    c4 = np.load(C4, allow_pickle=False)
    cp = {g: i for i, g in enumerate(map(str, c4["genes"]))}
    c4p = np.array([c4["P_down"][cp[g]] if g in cp else np.nan for g in genes], float)
    return dict(surface=surface, seed=seed, stats=str(stp), targets=targets, genes=genes, rq=rq, rp=rp, rl=rl,
                oq=oq, op=op, ol=ol, key=key, head=head, block=block, own=own, cpm=cpm_full[keep], ref_lin=ref_lin,
                c4=c4p)


def masks(S):
    adj = np.isfinite(S["rl"]) & (S["rl"] != 0) & ~S["own"]
    sig = (S["rq"] < ALPHA) & np.isfinite(S["rq"]) & ~S["own"]
    call = (S["oq"] < ALPHA) & np.isfinite(S["oq"]) & ~S["own"]
    return adj, sig, call


def fid_rows(call, sgn, pool, adj, rl):
    """fid per target (cell-eval2 de_wilcoxon_direction_fidelity_yield_raw at f 1): k / max(n_pred, n_real)."""
    sel = call & adj
    k = ((np.sign(sgn) == np.sign(rl)) & sel).sum(1)
    npred = sel.sum(1)
    nreal = pool.sum(1)
    den = np.maximum(npred, nreal)
    return np.where(den > 0, k / np.maximum(den, 1), np.nan)


def realized_order_pos(oq, op, ol, genes):
    """the scorer's ranking of the prediction: sig first, p_adj, p, |lfc| desc, name (cell-eval2 direction.py)."""
    sig = oq < ALPHA
    o = np.lexsort((genes, -np.abs(np.nan_to_num(ol)), np.nan_to_num(op, nan=1.0), np.nan_to_num(oq, nan=1.0), ~sig))
    pos = np.empty(o.size, np.int64)
    pos[o] = np.arange(o.size)
    return pos


def reach_row(order_pos, pool_i, adj_i, lfc_i, rl_i, p0=0.9):
    """de_wilcoxon_direction_reach_raw for one target (lib521.reach_row)."""
    idx = np.flatnonzero(pool_i)
    n_real = idx.size
    if n_real == 0:
        return np.nan
    idx = idx[np.argsort(order_pos[idx], kind="stable")]
    idx = idx[adj_i[idx]]
    if idx.size == 0:
        return 0.0
    ok = (np.sign(lfc_i[idx]) == np.sign(rl_i[idx])) & (lfc_i[idx] != 0)
    pur = np.cumsum(ok) / np.arange(1, idx.size + 1)
    deep = np.flatnonzero(pur >= p0)
    return (int(deep[-1]) + 1 if deep.size else 0) / n_real
