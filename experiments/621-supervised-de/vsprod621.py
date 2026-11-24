#!/usr/bin/env python3
"""621 — the switch vs what we SHIP: production's delivered lfc signs on h1-test (435's production emitter, via 521's
lib) against CLEAN4 and the 621 consensus switch, on identical gate pairs, all targets and the non-essential subset,
with a 2,000-draw target bootstrap of the paired differences. The switch sign is the family-balanced donor consensus
where |Zbal| >= c and CLEAN4 elsewhere. h1 donors exclude every H1 bank (636's rule). Only h1 has a delivered
production surface; the other three surfaces have no production emitter run, so they read vs CLEAN4 only (switch621).
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
from inst636 import ESS, Inst, clean4_lolo  # noqa: E402
import lib521 as L  # noqa: E402

RNG = np.random.default_rng(6211)


def main() -> None:
    I = Inst()
    D = I.S["h1"]
    T, G = D["T"], D["G"]
    z = np.load(HERE / "work/pred-h1.npz", allow_pickle=True)
    Zb = z["Zbal"].astype(float)
    C4 = np.broadcast_to(clean4_lolo(I, "h1")[0], D["lfc"].shape).astype(float)
    S = L.load_state(0)
    ti = {t: i for i, t in enumerate(S["targets"])}
    gi = {g: j for j, g in enumerate(S["genes"])}
    rows_t = np.array([ti.get(t, -1) for t in T]); cols_g = np.array([gi.get(g, -1) for g in G])
    assert (rows_t >= 0).all(), "h1 targets missing from 521 state"
    OL = np.full(D["lfc"].shape, np.nan)
    okc = cols_g >= 0
    OL[:, okc] = S["ol"][rows_t][:, cols_g[okc]]
    m = D["gate"] & np.isfinite(C4) & np.isfinite(Zb) & np.isfinite(OL) & (OL != 0)
    truth = np.sign(D["lfc"])
    sgn_c4 = np.where(C4 > 0.5, -1.0, 1.0)
    sgn_cons = np.sign(Zb)  # Zbal < 0 = donors say down = sign -1
    sgn_prod = np.sign(OL)
    noness = np.array([t not in ESS for t in T])
    out = {"n_gate_pairs": int(m.sum()), "rows": []}

    def per_target_correct(sgn, mask):
        ok = mask & (sgn != 0)
        c = ((sgn == truth) & ok).sum(1).astype(float)
        n = ok.sum(1).astype(float)
        return c, n

    for cset_name, tmask in (("all", np.ones(len(T), bool)), ("noness", noness)):
        mm = m & tmask[:, None]
        cp, npp = per_target_correct(sgn_prod, mm)
        res = {"targets": cset_name, "n": int(mm.sum()), "prod": float(cp.sum() / npp.sum())}
        c4c, c4n = per_target_correct(sgn_c4, mm)
        res["clean4"] = float(c4c.sum() / c4n.sum())
        for c in (0, 1, 2, 3):
            sw = np.where(np.abs(Zb) >= c, sgn_cons, sgn_c4)
            sc, sn = per_target_correct(sw, mm)
            acc = float(sc.sum() / sn.sum())
            # paired bootstrap over targets of (switch - prod)
            idx = np.flatnonzero(tmask & (npp > 0))
            bs = []
            for _ in range(2000):
                b = RNG.choice(idx, size=len(idx), replace=True)
                bs.append(sc[b].sum() / sn[b].sum() - cp[b].sum() / npp[b].sum())
            res[f"switch_c{c}"] = acc
            res[f"switch_c{c}_minus_prod_CI"] = [float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))]
            res[f"switch_c{c}_share_cons"] = float(((np.abs(Zb) >= c) & mm).sum() / mm.sum())
        out["rows"].append(res)
        print(json.dumps(res, indent=1))
    (HERE / "work/vsprod621.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
