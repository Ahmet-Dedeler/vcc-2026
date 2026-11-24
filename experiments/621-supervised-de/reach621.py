#!/usr/bin/env python3
"""621c — analytic ceiling of a DIRECTION-CONFIDENCE REACH ORDER on h1 (sweep candidate #2), order-only.

Production's delivered lfc (signs) and calls are PINNED (521's state = 435's production emitter on h1-test); only the
scorer's confidence order changes. Orders compared, with the scorer's structure (our significant calls first) kept:
  realized        today's order (sig first, then q, p, |lfc|)
  zabs            inside the called set by |Zbal| desc (then realized); uncalled kept in realized order
  zagree          inside the called set by agreement-confidence a = |Zbal| if sign(Zbal) == sign(our lfc) else -|Zbal|
                  (promote calls the consensus confidently backs, demote the ones it confidently contradicts)
  zagree_all      as zagree, applied to the uncalled tail too
  oracle_called   inside the called set, the genes whose delivered sign is right first (order-only ceiling)
Reach per BRIEF 5 via lib521.reach_row, at f = 1 and f = 0.35 (pools), all and non-essential targets, with a
2,000-draw target bootstrap of each order's delta vs realized. Zbal from 621's h1 features (no H1 bank).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
from inst636 import ESS, Inst  # noqa: E402
import lib521 as L  # noqa: E402

RNG = np.random.default_rng(6212)


def main() -> None:
    I = Inst(); D = I.S["h1"]
    Zb636 = np.load(HERE / "work/pred-h1.npz", allow_pickle=True)["Zbal"].astype(float)
    S = L.load_state(0)
    T, G = S["targets"], S["genes"]
    ti = {t: i for i, t in enumerate(D["T"])}; gi = {g: j for j, g in enumerate(D["G"])}
    rt = np.array([ti[t] for t in T]); cg = np.array([gi.get(g, -1) for g in G]); okg = cg >= 0
    Zb = np.full((len(T), len(G)), np.nan)
    Zb[:, okg] = Zb636[rt][:, cg[okg]]
    ol, oq, op = S["ol"], S["oq"], S["op"]
    adj = L.adjud(S)
    sig = oq < L.ALPHA
    noness = np.array([t not in ESS for t in T])
    absz = np.nan_to_num(np.abs(Zb), nan=0.0)
    agree = np.where(np.sign(np.nan_to_num(Zb)) == np.sign(np.nan_to_num(ol)), absz, -absz)

    def order_pos(i, key_called, key_tail=None):
        """sig first; inside sig by key desc then realized; tail by key_tail desc (or realized)."""
        rpos = L.realized_order_pos(S, i)
        front = np.flatnonzero(sig[i]); tail = np.flatnonzero(~sig[i])
        front = front[np.lexsort((rpos[front], -key_called[front]))] if key_called is not None else front[np.argsort(rpos[front])]
        tail = tail[np.lexsort((rpos[tail], -key_tail[tail]))] if key_tail is not None else tail[np.argsort(rpos[tail])]
        ordr = np.concatenate([front, tail]); pos = np.empty(ordr.size, np.int64); pos[ordr] = np.arange(ordr.size)
        return pos

    out = {}
    for f in (1.0, 0.35):
        pool = L.pools(S, f)
        rows = {k: np.full(len(T), np.nan) for k in ("realized", "zabs", "zagree", "zagree_all", "oracle_called")}
        for i in range(len(T)):
            right = ((np.sign(ol[i]) == np.sign(S["rl"][i])) & (ol[i] != 0)).astype(float)
            orders = {
                "realized": L.realized_order_pos(S, i),
                "zabs": order_pos(i, absz[i]),
                "zagree": order_pos(i, agree[i]),
                "zagree_all": order_pos(i, agree[i], agree[i]),
                "oracle_called": order_pos(i, right),
            }
            for k, pos in orders.items():
                rows[k][i] = L.reach_row(pos, pool[i], adj[i], ol[i], S["rl"][i])
        res = {}
        for subset, tm in (("all", np.ones(len(T), bool)), ("noness", noness)):
            base = rows["realized"]
            for k in rows:
                v = rows[k]
                ok = tm & np.isfinite(v) & np.isfinite(base)
                d = v[ok] - base[ok]
                idx = np.arange(d.size)
                bs = [d[RNG.choice(idx, idx.size)].mean() for _ in range(2000)]
                res[f"{subset}:{k}"] = dict(reach=float(np.nanmean(v[tm])), delta=float(d.mean()),
                                            ci=[float(np.quantile(bs, .025)), float(np.quantile(bs, .975))],
                                            overall_equiv=float(d.mean() * L.INV["reach"] / 6))
        out[f"f{f}"] = res
        print(f"\n== f = {f}")
        for k, r in res.items():
            print(f"  {k:24s} reach {r['reach']:.4f}  Δ vs realized {r['delta']:+.4f} [{r['ci'][0]:+.4f}, {r['ci'][1]:+.4f}]  overall-equiv {r['overall_equiv']:+.5f}")
    (HERE / "work/reach621.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
