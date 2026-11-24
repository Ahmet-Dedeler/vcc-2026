"""493 (third session) — does 492's per-target sharpening hold on October-like panels? 438's panel-agnostic PLOLO
machinery (lib438), the October stack h2c1 (H1 x2, CD4 x1), K3 key, ballast-realised at KMAX 16.

Rows are sharpened per target (169's sharpen(mode="row") on each frame's gate with the row's own q), then rendered and
realised exactly as lib438.realised_rows does. Arms vs the q 1.75 base: P2 (2.5, 2.0, 1.75 by the fold's K3 tercile),
L30 (3.0, 1.75, 1.75), S_P2 (P2's q by response-strength tercile), CTRL (P2's q permuted over targets).
KNOWN ANSWER (abort): the per-target path at q 1.5 for every row == 438's stored `h2c1|K3+sh1.5` row on this panel/fold
(design and realised 16; cos_agg, pds, pds lam 0.5/1, render A/B), |d| < 1e-12.

    uv run python experiments/493-pertarget-q-october/comp493.py <panel> <fold>   # panel in efp efr real abc; ~6 GB
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/438-october-stack-repricing"))
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import lib438 as L  # noqa: E402
from lib438 import ALPHA, GAMMA, LAM, SLOTS, E062  # noqa: E402
import arms414  # noqa: E402
from plolo import TS, profile_from_effect  # noqa: E402

W438 = ROOT / "experiments/438-october-stack-repricing/work"
TOL = 1e-12
SPECS = {"P2": (2.5, 2.0, 1.75), "L30": (3.0, 1.75, 1.75)}


def terc(v, spec, floor=None):
    t = np.digitize(v, np.quantile(v, [1 / 3, 2 / 3]))
    if floor is not None:                    # the October guard: every target at the key floor is LOW
        t = np.where(v <= floor + 1e-12, 0, t)
    return np.asarray([spec[k] for k in t], float)


def main(panel, held):
    t0 = time.time()
    log = lambda *a: print(f"[{panel}/{held}]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    ctx = L.bind(panel)
    fold = L.build_fold438(ctx, held, "held")
    scored, pp_pos, wcols = fold["scored"], fold["panel_positions"], fold["wcols"]
    truth = fold["truth"][:, scored]
    targets, sources, op = list(fold["targets"]), list(fold["sources"]), fold["op"]
    nonpanel = ctx.nonpanel
    dn = L.Donors438(targets, sources)
    cd4, partial, h1, _ = L.residuals(ctx, targets)
    pre = L.pre_normrel(dn, sources, "cap", -1.5, fallback_rows=fold["pre"])
    whit = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    FC = L.Finals(ctx, whit, cd4, partial, h1, wcols, pp_pos)
    key410, carried410, _ = L.shared_energy_key(targets, sources, op, wcols, cols=None)
    W3 = L.weights(key410, carried410, 0.3, 3.0)[0]
    P_DE = L.carrier_profiles(FC.final("h1c1"))
    bulk_k3 = L.key_scale(FC.final("h1c1"), W3)[0]
    dk3 = L.realised_rows(bulk_k3, P_DE, q=None, kmaxes=(), zero_cols=None)
    bet = L.betas(ctx, held, nonpanel, A_pp_k3={c: L.energy(dk3[c]["design"], nonpanel)["A"] for c in SLOTS})
    fin = FC.final("h2c1")
    bulk = L.key_scale(fin, W3)[0]
    holder = {"bulk": bulk}
    strength = np.linalg.norm(fin[:, nonpanel], axis=1)
    res = dict(panel=panel, fold=held, n_targets=len(targets), ka={}, arms={})
    arrays = {}

    def rows_q(qv):
        out = {}
        for c, (b, gate) in arms414.frames().items():
            B = holder["bulk"]
            F = np.empty_like(B)
            for q in np.unique(qv):
                sel = qv == q
                F[sel] = L.sharpen(B[sel], float(q), gate, mode="row")
            des, rea, K = np.empty_like(F), np.empty_like(F), np.empty(F.shape[0], int)
            for i in range(F.shape[0]):
                p_bk = profile_from_effect(F[i], b, gate, "minnorm")
                des[i] = np.log1p(TS * p_bk) - b
                p_re, _, K[i], _, _ = L.realise_k(P_DE[c][i], p_bk, 16)
                rea[i] = np.log1p(TS * p_re) - b
            out[c] = dict(design=des, realised=rea, K=K)
        return out

    def rec(name, qv):
        dr = rows_q(qv)
        sd, ad = L.row_score({c: dr[c]["design"] for c in SLOTS}, scored, truth, bet, nonpanel)
        sr, ar = L.row_score({c: dr[c]["realised"] for c in SLOTS}, scored, truth, bet, nonpanel)
        K = np.concatenate([dr[c]["K"] for c in SLOTS])
        sr["frac_capped"], sr["K_mean"] = float(np.mean(K >= 16)), float(np.mean(K))
        res["arms"][name] = dict(design=sd, realised={"16": sr})
        arrays[f"{name}@R16"] = ar
        log(f"{name:18s} R16 cos {sr['cos_agg']:.6f} pds {sr['pds']:.6f} pds@.5 {sr['pds_lam0.5']:.6f} K {sr['K_mean']:.2f}")
        return res["arms"][name]

    n = len(targets)
    e = rec("h2c1|K3+sh1.5", np.full(n, 1.5))
    ref = json.loads((W438 / panel / f"comp438-{held}-held.json").read_text())["arms"]["h2c1|K3+sh1.5"]
    dD, dR = L.cmp(e["design"], ref["design"]), L.cmp(e["realised"]["16"], ref["realised"]["16"])
    log(f"KA h2c1|K3+sh1.5 vs 438: design |d| {dD:.1e}, realised |d| {dR:.1e}")
    assert dD < TOL and dR < TOL, f"KA FAILED {panel}/{held}: {dD} {dR}"
    res["ka"]["h2c1_sh1.5"] = dict(design=dD, realised=dR)
    rec("base q1.75", np.full(n, 1.75))
    for name, spec in SPECS.items():
        rec(name, terc(W3, spec))
    rec("S_P2", terc(strength, SPECS["P2"]))
    for name, spec in SPECS.items():            # 493b: the knob with the floor guard (identical where the floor < 1/3)
        rec(f"{name}g", terc(W3, spec, floor=0.3))
    res["frac_floor"] = float(np.mean(W3 <= 0.3 + 1e-12))
    if "--stack" in sys.argv:                   # 493c: floor 0.5 stacked with the per-target q (tercile from K3, guarded)
        W5 = L.weights(key410, carried410, 0.5, 3.0)[0]
        bulk5 = L.key_scale(fin, W5)[0]
        holder["bulk"] = bulk5
        rec("F5 base q1.75", np.full(n, 1.75))
        rec("F5 P2g", terc(W3, SPECS["P2"], floor=0.3))
        rec("F5 S_P2", terc(strength, SPECS["P2"]))
        holder["bulk"] = bulk
    rec("CTRL P2 permuted", terc(W3, SPECS["P2"])[np.random.default_rng(493).permutation(n)])
    res["peak_rss_gib"] = E062.peak_rss_gib()
    (HERE / "work" / (f"comp493c-{panel}-{held}.json" if "--stack" in sys.argv else f"comp493-{panel}-{held}.json")).write_text(json.dumps(res, indent=1, default=float))
    np.savez_compressed(HERE / "work" / (f"arrays493c-{panel}-{held}.npz" if "--stack" in sys.argv else f"arrays493-{panel}-{held}.npz"), names=np.asarray(list(arrays)),
                        d=np.asarray((truth * truth).sum(1)), **{f"{k}__{kk}": v for k, a in arrays.items() for kk, v in a.items()})
    log(f"done; peak RSS {res['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
