"""521 — shared machinery: the shipped H1 DE state (435's production-emitter stats) with its call set reconstructed,
and the four DE members as exact functions of (call set, delivered lfc, order) under proportional power f.

The shipped construction on H1 is 435's `base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse` (production emitter
`rc_emit_k`, CallSet(400, 2200, 1.9), clean4 q80 block prior, coarse MDE; `docs/FINDINGS.md` §1). Its stats hold our
per-gene p, q, lfc on the scorer's 5-CPM universe (10,779 genes); 323's reference.npz holds the truth's p, q, lfc.

The emitter's call set and order are reconstructed from their definitions (src/vcc2026/rc_emit_k.py):
  key_p     = where(cpm_ctrl > 10 & b_ctrl > 0, cap15_effect_p, 0)          (run435's `np.where(gate, eff_de[i], 0.)`)
  head_p    = top 400 of the gated universe by |key_p| (stable argsort)
  block_p   = top 2,200 by emit_ref_lin among the rest (stable)             -> C_p = head_p | block_p
  order     = |key| + jitter(seed 1000 + i, 1e-12) + 1e6 * [g in C_p]; the complement's |z| is one common factor
              below z_cap, so the p order is |key| order inside C_p and inside its complement.
Nothing here writes outside experiments/521-strength-nreal/.
"""
from __future__ import annotations

import sys
from functools import lru_cache
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
REF = ROOT / "experiments/323-current-family-validation/work/reference.npz"
EFF = ROOT / "experiments/323-current-family-validation/work/effects.npz"
ST435 = ROOT / "experiments/435-zs-on-coarse/work/h1/stats-base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse.npz"
MET435 = ROOT / "experiments/435-zs-on-coarse/work/h1/metrics-base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse.json"
DIAG435 = ROOT / "experiments/435-zs-on-coarse/work/h1/diag-base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse.json"
OFF416 = ROOT / "experiments/416-de-carrier-menu/work/h1/stats-prior_off_s{seed}.npz"
ALPHA, P0 = 0.05, 0.9
N_HEAD, N_CPM = 400, 2200
JITTER = 1.0e-12
OFFSET = 1.0e6
# board-fitted 1/span (192's INV; inside the BRIEF's published spans). pds/mse do not move in this family.
INV = dict(reach=1.12927, fid=3.32482, jac=2.72081, nmae=1.64787)


@lru_cache(maxsize=4)
def load_state(seed: int = 0):
    ref = np.load(REF)
    genes_full = ref["genes"].astype(str)
    keep = ref["keep"].astype(bool)
    targets = ref["targets"].astype(str)
    genes = genes_full[keep]
    rq = np.asarray(ref["q"], float)
    rp = np.asarray(ref["p"], float)
    rl = np.asarray(ref["lfc"], float)
    ref_lin = np.asarray(ref["emit_ref_lin"], float)
    b_ctrl = np.asarray(ref["b_ctrl"], float)
    st = np.load(str(ST435).format(seed=seed))
    assert list(st["targets"].astype(str)) == list(targets)
    oq, op, ol = (np.asarray(st[k], float) for k in ("q", "p", "lfc"))
    # the key (run435: eff_de = map_axis(cap15, genes, r["genes"]); key = where(gate, eff_de, 0))
    e = np.load(EFF)
    eg = e["genes"].astype(str)
    assert list(e["targets"].astype(str)) == list(targets)
    pos = {g: i for i, g in enumerate(eg)}
    cols = np.array([pos.get(g, -1) for g in genes_full])
    eff = np.zeros((len(targets), len(genes_full)))
    eff[:, cols >= 0] = np.asarray(e["cap15"], float)[:, cols[cols >= 0]]
    ex = np.expm1(b_ctrl)
    cpm = 1e6 * ex / ex.sum()
    gate10 = (cpm > 10.0) & (b_ctrl > 0)
    key = np.where(gate10[None, :], eff, 0.0)[:, keep]
    T, G = key.shape
    own = np.equal.outer(targets, genes)
    head = np.zeros((T, G), bool)
    block = np.zeros((T, G), bool)
    jit = np.zeros((T, G))
    for i in range(T):
        ka = np.abs(key[i])
        head[i, np.argsort(-ka, kind="stable")[:N_HEAD]] = True
        block[i, np.argsort(-np.where(head[i], -np.inf, ref_lin), kind="stable")[:N_CPM]] = True
        jit[i] = np.random.default_rng(1000 + i).random(G) * JITTER
    return dict(targets=targets, genes=genes, rq=rq, rp=rp, rl=rl, oq=oq, op=op, ol=ol, ref_lin=ref_lin,
                key=key, absk=np.abs(key), head=head, block=block, C=head | block, jit=jit, own=own, seed=seed)


def pools(S, f: float = 1.0):
    """PowerSurface (157/192): pool_p = the reference's top round(f * n_real_p) genes by its own p (stable), own target
    removed. f = 1 is the BH set itself (n_real_p = #{q < 0.05} without the own target)."""
    rq, rp, own = S["rq"], S["rp"], S["own"]
    base = (rq < ALPHA) & np.isfinite(rq) & ~own
    if f == 1.0:
        return base
    out = np.zeros_like(base)
    rpx = np.where(np.isfinite(rp) & ~own, rp, np.inf)
    for i in range(base.shape[0]):
        k = int(round(f * int(base[i].sum())))
        if k:
            out[i, np.argsort(rpx[i], kind="stable")[:k]] = True
    return out & ~own


def adjud(S):
    return np.isfinite(S["rl"]) & (S["rl"] != 0) & ~S["own"]


def fid_rows(call, lfc, pool, adj, rl):
    """BRIEF 4, per target (NaN where n_pred = n_real = 0: dropped from the mean)."""
    out = np.full(call.shape[0], np.nan)
    for i in range(call.shape[0]):
        n_real = int(pool[i].sum())
        sel = call[i] & adj[i]
        n_pred = int(sel.sum())
        if n_pred == 0 and n_real == 0:
            continue
        k = int((np.sign(lfc[i, sel]) == np.sign(rl[i, sel])).sum())
        out[i] = k / max(n_pred, n_real)
    return out


def jac_rows(call, pool):
    inter = (call & pool).sum(1)
    uni = (call | pool).sum(1)
    return np.where(uni == 0, 1.0, inter / np.maximum(uni, 1))


def nmae_rows(lfc, pool, rl):
    out = np.full(lfc.shape[0], np.nan)
    for i in range(lfc.shape[0]):
        s = pool[i] & np.isfinite(rl[i])
        if s.sum() < 10:
            continue
        den = np.abs(rl[i, s]).sum()
        if not np.isfinite(den) or den == 0:
            continue
        h = np.nan_to_num(lfc[i, s], nan=0.0, posinf=0.0, neginf=0.0)
        out[i] = np.abs(h - rl[i, s]).sum() / den
    return out


def reach_row(order_pos, pool_i, adj_i, lfc_i, rl_i):
    """BRIEF 5 for one target; order_pos = each gene's position in our confidence order (smaller = earlier)."""
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
    deep = np.flatnonzero(pur >= P0)
    return (int(deep[-1]) + 1 if deep.size else 0) / n_real


def realized_order_pos(S, i):
    """The scorer's order from the realized stats: sig first, adj p, p, |lfc| desc, name (BRIEF 5)."""
    oq, op, ol, g = S["oq"][i], S["op"][i], S["ol"][i], S["genes"]
    sig = oq < ALPHA
    o = np.lexsort((g, -np.abs(np.nan_to_num(ol)), np.nan_to_num(op, nan=1.0), np.nan_to_num(oq, nan=1.0), ~sig))
    pos = np.empty(o.size, np.int64)
    pos[o] = np.arange(o.size)
    return pos


def model_order_pos(S, i, call_i):
    """The emitter's designed order for a delivered call set `call_i` that is a SUBSET of today's realized calls:
    the called genes first in their REALIZED (q, p) order (inside the prescribed set that is |key| order, and a subset
    keeps it), then every other gene by |key| + jitter, because the complement's target |z| is one common monotone
    factor of |key| order (rc_emit_k: `z_target[complement] *= z_cap / z_worst`)."""
    oq, op, g = S["oq"][i], S["op"][i], S["genes"]
    assert not np.any(call_i & ~(oq < ALPHA)), "model order handles subsets of today's calls only"
    front = np.flatnonzero(call_i)
    front = front[np.lexsort((g[front], np.nan_to_num(op[front], nan=1.0), np.nan_to_num(oq[front], nan=1.0)))]
    tail = np.flatnonzero(~call_i)
    tail = tail[np.lexsort((g[tail], -(S["absk"][i][tail] + S["jit"][i][tail])))]
    ordr = np.concatenate([front, tail])
    pos = np.empty(ordr.size, np.int64)
    pos[ordr] = np.arange(ordr.size)
    return pos


def members(S, call, lfc, f=1.0, order="model", rows=False):
    """fid, jac, reach, nmae for a prescribed call set (bool T x G) and delivered lfc, at proportional power f."""
    pool = pools(S, f)
    adj = adjud(S)
    rl = S["rl"]
    fr = fid_rows(call, lfc, pool, adj, rl)
    jr = jac_rows(call, pool)
    nr = nmae_rows(lfc, pool, rl)
    rr = np.full(call.shape[0], np.nan)
    for i in range(call.shape[0]):
        pos = realized_order_pos(S, i) if order == "realized" else model_order_pos(S, i, call[i])
        rr[i] = reach_row(pos, pool[i], adj[i], lfc[i], rl[i])
    out = dict(fid=float(np.nanmean(fr)), jac=float(np.mean(jr)), reach=float(np.nanmean(rr)), nmae=float(np.nanmean(nr)))
    if rows:
        out["rows"] = dict(fid=fr, jac=jr, reach=rr, nmae=nr)
    return out


def overall_delta(m1, m0):
    """board-overall delta of the four DE members (192's INV; nmae lower is better)."""
    return ((m1["fid"] - m0["fid"]) * INV["fid"] + (m1["jac"] - m0["jac"]) * INV["jac"]
            + (m1["reach"] - m0["reach"]) * INV["reach"] - (m1["nmae"] - m0["nmae"]) * INV["nmae"]) / 6.0
