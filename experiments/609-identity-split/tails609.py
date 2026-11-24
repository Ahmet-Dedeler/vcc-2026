#!/usr/bin/env python3
"""609b ([4d0207], 2026-10-02 ~14:20Z): does 622n's "low-identity tail goes UP" (found on UNCALLED genes) also hold in
the CPM BLOCK and the head? If so, a targeted override (low-I block genes -> UP, CLEAN4 elsewhere) would be a cheap,
exactly-readable prior-only arm; if not, the low-I signal is uncalled-only and 630's UP pass stays its only route.
Same surfaces, masks and identity as split609 (639d's machinery; KA'd there). Per region and identity decile:
gate pairs, truth-UP fraction, CLEAN4's UP-call rate, and the accuracy of CLEAN4 vs "force UP" on the bottom-I decile
(and bottom quintile). Expression-only control: the same tail cut on E = log1p(own control CPM) instead of I.
usage: uv run python experiments/609-identity-split/tails609.py   -> work/tails609.json"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/639-identity-axis-lolo"))
import decomp639 as D  # noqa: E402


def main() -> None:
    I = D.Inst()
    prof = D.flex_profiles()
    out = {}
    for s in I.S:
        X = I.S[s]
        T, G, lfc, gate = X["T"], X["G"], X["lfc"], X["gate"]
        if s == "h1":
            own = prof["h1"]
            refs = [prof[k] for k in ("ctxA", "ctxB", "ctxC")]
        else:
            zz = np.load(D.G480 / D.SURF[s], allow_pickle=True)
            own = dict(zip(zz["genes"].astype(str), zz["ctrl_cpm"].astype(float)))
            refs = [prof[k] for k in ("ctxA", "ctxB", "ctxC", "h1")]
        Iv = D.identity(own, refs, list(G))
        Ev = np.array([np.log1p(own[g]) if (g in own and np.isfinite(Iv[j])) else np.nan for j, g in enumerate(G)])
        C4, _ = D.clean4_lolo(I, s)
        c4up = (C4 < 0.5)                                   # CLEAN4 says UP
        head, block = D.callset_for(I, s)
        unc = ~head & ~block
        up = lfc > 0
        rec = {}
        for rname, R in (("head", head), ("block", block), ("uncalled", unc)):
            rr = {}
            for axis_name, ax in (("I", Iv), ("E", Ev)):
                fin = np.isfinite(ax) & np.isfinite(C4)
                for qname, q in (("bottom10", 0.10), ("bottom20", 0.20), ("top10", None)):
                    if q is None:
                        cut = np.nanquantile(ax[fin], 0.90)
                        tail = fin & (ax >= cut)
                    else:
                        cut = np.nanquantile(ax[fin], q)
                        tail = fin & (ax <= cut)
                    m = gate & R & tail[None, :]
                    n = int(m.sum())
                    if n == 0:
                        continue
                    tu = float(up[m].mean())
                    acc_c4 = float(((c4up[None, :] == up) & m).sum() / n)
                    acc_force = tu if q is not None else 1.0 - tu       # bottom tails: force UP; top tail: force DOWN
                    rr[f"{axis_name}_{qname}"] = dict(n=n, truth_up=tu, c4_up_rate=float(c4up[None, :].repeat(len(T), 0)[m].mean()),
                                                      acc_clean4=acc_c4, acc_force=acc_force, gain=acc_force - acc_c4)
            rec[rname] = rr
        out[s] = rec
        print(f"\n== {s}")
        for rname in ("head", "block", "uncalled"):
            for key in ("I_bottom10", "I_bottom20", "E_bottom10", "I_top10"):
                v = rec[rname].get(key)
                if v:
                    print(f"  {rname:8s} {key:11s} n {v['n']:6d} truth_up {v['truth_up']:.3f} C4 says up {v['c4_up_rate']:.3f} "
                          f"acc C4 {v['acc_clean4']:.3f} force {v['acc_force']:.3f} gain {v['gain']:+.3f}", flush=True)
    (HERE / "work/tails609.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
