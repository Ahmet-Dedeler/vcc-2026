#!/usr/bin/env python3
"""639d — WHERE production would deliver the identity axis's sign gain: head / CPM block / uncalled (coordinator's ask via
[21dd4c], 2026-10-02 ~13:05Z).

Call set per surface:
  h1     : lib521.load_state(0)["head"], ["block"] (the real design key's head 400 + the CPM block 2,200);
  folds  : ds612c.callset(cons0, cnt, 400, 2200) (the 3-donor consensus key's head + the block by control expression),
           cons0 / cnt exactly as ds612.surface_forms builds them.
Sign sources (636's instrument, gate genes, LOLO): CLEAN4 (P_down; forbidden banks dropped) vs score = z(P_C4 - .5) +
lam * z(I), I from a FLEX-ONLY reference set (A/B/C + h1, never the surface itself; 622: 3' references hurt), lam 0.5 / 1.
Per region: gate pairs, paired sign accuracy (both sources defined), the paired gain, the coverage extension (identity
only where CLEAN4 is undefined vs no prior), |lfc|-weighted versions, and each region's SHARE of the total net gain
(net = (right - wrong) count, and |lfc|-weighted). Analysis only, < 4 GB.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
sys.path.insert(0, str(ROOT / "experiments/612-key-gated-hook"))
from inst636 import G480, SURF, Inst, clean4_lolo  # noqa: E402
import ident639 as M  # noqa: E402
import lib521 as L  # noqa: E402

FOLD = {"hct116": "orion_hct116", "hek293t": "orion_hek293t", "k562": "k562_gwps"}


def flex_profiles():
    prof = M.abc_profiles()
    z = np.load(G480 / SURF["h1"], allow_pickle=True)
    prof["h1"] = dict(zip(z["genes"].astype(str), z["ctrl_cpm"].astype(float)))
    return prof


def identity(own_prof, refs, genes):
    common = set(own_prof).intersection(*(set(r) for r in refs))
    cg = sorted(common)
    I = np.log1p(np.array([own_prof[g] for g in cg])) - np.mean([np.log1p(np.array([r[g] for g in cg])) for r in refs], 0)
    pos = {g: i for i, g in enumerate(genes)}
    out = np.full(len(genes), np.nan)
    for g, v in zip(cg, I):
        j = pos.get(g)
        if j is not None:
            out[j] = v
    return out


def callset_for(I, s):
    D = I.S[s]
    T, G = D["T"], D["G"]
    if s == "h1":
        st = L.load_state(0)
        assert list(st["targets"]) == list(T)
        pos = {g: i for i, g in enumerate(st["genes"])}
        cols = np.array([pos.get(g, -1) for g in G])
        head = np.zeros((len(T), len(G)), bool)
        block = np.zeros((len(T), len(G)), bool)
        ok = cols >= 0
        head[:, ok] = st["head"][:, cols[ok]]
        block[:, ok] = st["block"][:, cols[ok]]
        return head, block
    import ds612 as DS
    import ds612c as DC
    Sx = DS.F.load_surface(FOLD[s])
    genes_f = [str(g) for g in Sx["genes"]]
    targets_f = list(Sx["targets"])
    cons0 = np.nan_to_num(np.nanmean(np.stack([DS.load_donor_rows(d, targets_f, np.array(genes_f), pc_cpm=1.0)[0]
                                               for d in Sx["donors"]]), 0))
    Cf = DC.callset(cons0, Sx["cnt"], 400, 0)              # head only
    Call = DC.callset(cons0, Sx["cnt"], 400, 2200)         # head + block
    # map onto the instrument's axes (targets, genes) by name
    tpos = {t: i for i, t in enumerate(targets_f)}
    gpos = {g: i for i, g in enumerate(genes_f)}
    ti = np.array([tpos.get(t, -1) for t in T])
    gi = np.array([gpos.get(g, -1) for g in G])
    head = np.zeros((len(T), len(G)), bool)
    block = np.zeros((len(T), len(G)), bool)
    for a, t in enumerate(ti):
        if t < 0:
            continue
        okg = gi >= 0
        head[a, okg] = Cf[t, gi[okg]]
        block[a, okg] = Call[t, gi[okg]] & ~Cf[t, gi[okg]]
    return head, block


def main():
    I = Inst()
    prof = flex_profiles()
    out = {}
    for s in I.S:
        D = I.S[s]
        T, G, lfc, gate = D["T"], D["G"], D["lfc"], D["gate"]
        if s == "h1":
            own = prof["h1"]
            refs = [prof[k] for k in ("ctxA", "ctxB", "ctxC")]
        else:
            z = np.load(G480 / SURF[s], allow_pickle=True)
            own = dict(zip(z["genes"].astype(str), z["ctrl_cpm"].astype(float)))
            refs = [prof[k] for k in ("ctxA", "ctxB", "ctxC", "h1")]
        Iv = identity(own, refs, list(G))
        zi = (Iv - np.nanmean(Iv)) / np.nanstd(Iv)
        C4, _ = clean4_lolo(I, s)
        zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)
        head, block = callset_for(I, s)
        unc = ~head & ~block
        dn = lfc < 0
        w = np.abs(np.nan_to_num(lfc))
        rec = {"n_gate": int(gate.sum()), "regions": {}}
        for lam in (0.5, 1.0):
            sc = zc + lam * zi
            both = np.isfinite(zc) & np.isfinite(zi)
            ext = ~np.isfinite(zc) & np.isfinite(zi)
            totals = {"paired_net": 0.0, "paired_netw": 0.0, "ext_net": 0.0, "ext_netw": 0.0}
            regs = {}
            for rname, R in (("head", head), ("block", block), ("uncalled", unc)):
                m = gate & R & both[None, :]
                okc = ((zc > 0)[None, :] == dn) & m
                oks = ((sc > 0)[None, :] == dn) & m
                n = int(m.sum())
                d_net = float((2 * oks.sum() - n) - (2 * okc.sum() - n))           # change in (right - wrong)
                d_netw = float(((2 * oks - 1) * w * m).sum() - ((2 * okc - 1) * w * m).sum())
                me = gate & R & ext[None, :]
                oke = ((zi > 0)[None, :] == dn) & me
                ext_net = float(2 * oke.sum() - me.sum())                          # vs no prior (0)
                ext_netw = float(((2 * oke - 1) * w * me).sum())
                regs[rname] = dict(gate_pairs=int((gate & R).sum()), paired_n=n,
                                   acc_C4=float(okc.sum() / max(n, 1)), acc_comb=float(oks.sum() / max(n, 1)),
                                   gain=float((oks.sum() - okc.sum()) / max(n, 1)), d_net=d_net, d_netw=d_netw,
                                   ext_n=int(me.sum()), ext_acc=float(oke.sum() / max(me.sum(), 1)), ext_net=ext_net,
                                   ext_netw=ext_netw, gate_mass_share=float((w * gate * R).sum() / max((w * gate).sum(), 1e-12)))
                totals["paired_net"] += d_net
                totals["paired_netw"] += d_netw
                totals["ext_net"] += ext_net
                totals["ext_netw"] += ext_netw
            for rname in regs:
                r = regs[rname]
                r["share_paired_net"] = r["d_net"] / totals["paired_net"] if totals["paired_net"] else float("nan")
                r["share_paired_netw"] = r["d_netw"] / totals["paired_netw"] if totals["paired_netw"] else float("nan")
                r["share_all_netw"] = ((r["d_netw"] + r["ext_netw"]) / (totals["paired_netw"] + totals["ext_netw"])
                                       if (totals["paired_netw"] + totals["ext_netw"]) else float("nan"))
            rec["regions"][f"lam{lam:g}"] = regs
        out[s] = rec
        print(f"\n== {s} (gate pairs {rec['n_gate']})", flush=True)
        for lam in ("lam0.5", "lam1"):
            print(f"  {lam}: {'region':9s} {'gate pairs':>10s} {'mass':>5s} {'acc C4':>7s} {'acc comb':>8s} {'gain':>7s} "
                  f"{'share net':>9s} {'share |l|-w':>11s} | {'ext n':>6s} {'ext acc':>7s} {'share all |l|-w':>15s}")
            for rname, r in rec["regions"][lam].items():
                print(f"        {rname:9s} {r['gate_pairs']:10d} {r['gate_mass_share']:5.2f} {r['acc_C4']:7.4f} {r['acc_comb']:8.4f} "
                      f"{r['gain']:+7.4f} {r['share_paired_net']:9.2f} {r['share_paired_netw']:11.2f} | {r['ext_n']:6d} "
                      f"{r['ext_acc']:7.4f} {r['share_all_netw']:15.2f}")
    (HERE / "work/decomp639.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    sys.exit(main())
