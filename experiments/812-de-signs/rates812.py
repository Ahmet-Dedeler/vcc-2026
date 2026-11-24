"""812 — exchange rates: what +x of sign accuracy on each gene class buys in the four DE members (H1, f 1, analytic on the
realized stats; calls, order and magnitudes fixed, only signs move). Answers "where could the field's +0.01 come from"
in board units, and prices any future sign source (or 811's model) by class before an emitter run.

Classes (per target, own gene excluded):
  a0      called, truth non-significant          (fid's filler, 721's lever)
  a1      called, truth significant               (fid's TP part, reach, nmae on called genes)
  a1head  called head (top-400 |key|), truth sig   (the front of the order: reach)
  u1      uncalled, truth significant             (nmae's bulk: the hook's genes)
Move: a random share of the WRONG-signed cells in the class is flipped to the truth's sign so that class accuracy rises by
+x (x 0.01 / 0.05); 3 random draws. Flipping an uncalled gene flips its delivered lfc sign (nmae); flipping a called gene
flips its lfc sign (fid, reach purity, nmae).
KA: the unflipped members equal the arm's recorded raw fid / reach / nmae (1e-12).

usage: uv run python experiments/812-de-signs/rates812.py   -> work/rates812.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import surf812 as SF  # noqa: E402

MET = {"h1test": SF.ROOT / "experiments/435-zs-on-coarse/work/h1/metrics-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.json",
       "h1tv": SF.ROOT / "experiments/618-h1-trainval-surface/work/g1/metrics-base_s{s}_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.json"}
INV = dict(fid=3.32482, reach=1.12927, nmae=1.64787, jac=2.72081)


def nmae_rows(lfc, pool, rl):
    out = np.full(lfc.shape[0], np.nan)
    for i in range(lfc.shape[0]):
        s = pool[i] & np.isfinite(rl[i])
        if s.sum() < 10:
            continue
        den = np.abs(rl[i, s]).sum()
        if den == 0:
            continue
        h = np.nan_to_num(lfc[i, s], nan=0.0, posinf=0.0, neginf=0.0)
        out[i] = np.abs(h - rl[i, s]).sum() / den
    return out


def members(S, sig, adj, call, ol, pos):
    rl = S["rl"]
    T = rl.shape[0]
    fid = float(np.nanmean(SF.fid_rows(call, ol, sig, adj, rl)))
    reach = float(np.nanmean([SF.reach_row(pos[i], sig[i], adj[i], ol[i], rl[i]) for i in range(T)]))
    nm = float(np.nanmean(nmae_rows(ol, sig, rl)))
    return dict(fid=fid, reach=reach, nmae=nm)


out = {}
for surf in ("h1test", "h1tv"):
    S = SF.load(surf, 0)
    adj, sig, call = SF.masks(S)
    rl, ol, oq, op, genes = S["rl"], S["ol"], S["oq"], S["op"], S["genes"]
    T = rl.shape[0]
    pos = [SF.realized_order_pos(oq[i], op[i], ol[i], genes) for i in range(T)]   # order is |lfc|-tied only; sign flips keep it
    m0 = members(S, sig, adj, call, ol, pos)
    rec = json.load(open(str(MET[surf]).format(s=0)))["raw"]
    assert all(abs(m0[k] - rec[k]) < 1e-9 for k in ("fid", "reach", "nmae")), f"KA FAIL {surf}: {m0} vs {rec}"
    right = (np.sign(ol) == np.sign(rl)) & (ol != 0)
    classes = {"a0": call & adj & ~sig, "a1": call & sig, "a1head": call & sig & S["head"], "u1": ~call & sig & adj}
    R = {"base": m0, "class_acc": {}, "rates": {}}
    for cn, cm in classes.items():
        n = int(cm.sum())
        acc = float(right[cm].mean())
        R["class_acc"][cn] = dict(n=n, acc=round(acc, 4))
        for x in (0.01, 0.05):
            wrong = np.flatnonzero((cm & ~right).ravel())
            k = int(round(x * n))
            if k > wrong.size:
                continue
            ds = []
            for d in range(3):
                rng = np.random.default_rng(812000 + d)
                fl = rng.choice(wrong, k, replace=False)
                ol1 = ol.copy().ravel()
                ol1[fl] = np.abs(ol1[fl]) * np.sign(rl.ravel()[fl])
                ol1[fl] = np.where(ol1[fl] == 0, 1e-3 * np.sign(rl.ravel()[fl]), ol1[fl])
                m1 = members(S, sig, adj, call, ol1.reshape(ol.shape), pos)
                ds.append({k2: m1[k2] - m0[k2] for k2 in m0})
            mean = {k2: float(np.mean([d_[k2] for d_ in ds])) for k2 in m0}
            ov = (mean["fid"] * INV["fid"] + mean["reach"] * INV["reach"] - mean["nmae"] * INV["nmae"]) / 6
            R["rates"][f"{cn}+{x}"] = dict(cells=k, dfid=round(mean["fid"], 6), dreach=round(mean["reach"], 6),
                                           dnmae=round(mean["nmae"], 6), d_overall=round(ov, 6),
                                           draws_overall=[round((d_["fid"] * INV["fid"] + d_["reach"] * INV["reach"] - d_["nmae"] * INV["nmae"]) / 6, 6) for d_ in ds])
            print(surf, cn, x, R["rates"][f"{cn}+{x}"], flush=True)
    out[surf] = R
json.dump(out, open(SF.HERE / "work/rates812.json", "w"), indent=1)
