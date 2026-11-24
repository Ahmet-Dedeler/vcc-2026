"""543o (third session) — the PERCENTILE cap on 438's October-like panels (543's rule: hi = the p-th percentile of the
carried key/median), vs fixed caps 3 and 6 (540o's rows, known answers).

(540o's docstring follows.)
540o (third session) — the LODO key's CAP on 438's October-like panels (540 found cap 6 best on the A/B/C folds).
comp493g.py's path (lib438; October stack h2c1; the fold's LODO-V1 fit exactly as comp438.py), at q 1.75 for every row,
with the key re-fitted at clip (0.5, hi) for hi in 3 (base = the October default), 4, 6, 8 and uncapped, and cap 6 with
the per-target arms (P2g with carried-only cut points and the floor guard; S_P2).
KNOWN ANSWER (abort): clip (0.5, 3) at q 1.5 == 438's stored `h2c1|V1_5_3+sh1.5` (design and realised 16), |d| < 1e-12.

(comp493g.py's docstring follows.)
493g (third session) — 492's per-target sharpening on October-like panels under the October key: the LODO key at clip
(0.5, 3) (`DEF_DAY_LODO=1`, scripts/def_lodo_key.py), instead of 493's K3 key.

comp493.py's path (438's panel-agnostic lib438; the October stack h2c1 = H1 x2, CD4 x1; ballast-realised at KMAX 16) with
the bulk keyed by the fold's LODO-V1 weights at (0.5, 3), computed exactly as comp438.py does (lodo_bases + lodo_key on
the h1c1 composition; weights(kv1, carried>=2, 0.5, 3)). On these panels the LODO floor holds 20-45 % of targets
(comp438's `keys`), so the floor guard decides the knob wherever it holds more than a third.

Arms vs the q 1.75 base: P2g (2.5, 2.0, 1.75 by the fold's LODO tercile, every floor target LOW), P2 without the guard,
S_P2 (P2's q by response-strength tercile: the unkeyed h2c1 row's non-panel norm), L30g (3.0, 1.75, 1.75, guarded),
CTRL (P2g's q permuted, seed 493).
KNOWN ANSWER (abort): the per-target path at q 1.5 for every row == 438's stored `h2c1|V1_5_3+sh1.5` row on this
panel/fold (design and realised 16; cos_agg, pds, pds lam 0.5/1, render A/B), |d| < 1e-12.

    uv run python experiments/543-quantile-cap/comp543o.py <panel> <fold>   # panel in efp efr real; ~6 GB
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
P2, L30 = (2.5, 2.0, 1.75), (3.0, 1.75, 1.75)
FLOOR = 0.5


def terc(v, spec, floor=None):
    t = np.digitize(v, np.quantile(v, [1 / 3, 2 / 3]))
    if floor is not None:                    # the October guard: every target at the key floor is LOW
        t = np.where(v <= floor + 1e-12, 0, t)
    return np.asarray([spec[k] for k in t], float), t


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
    # the fold's LODO-V1 at (0.5, 3), comp438.py's arithmetic
    d_rows = {s: L.donor_rows(s, targets)[0] for s in sources}
    A_ = {s: dn.A[s] for s in sources}
    carried_all = sum(A_[s].astype(int) for s in sources) >= 2
    bases = L.lodo_bases(op, dn, sources, "cap", -1.5, pre, wcols)
    kv1 = L.lodo_key(FC.final("h1c1"), FC.parts("h1c1"), bases, d_rows, A_, carried_all, sources, pp_pos)
    kvv = np.nan_to_num(kv1, nan=0.0)
    WV5 = L.weights(kvv, carried_all, 0.5, 3.0)[0]
    del bases, d_rows
    P_DE = L.carrier_profiles(FC.final("h1c1"))
    bulk_k3 = L.key_scale(FC.final("h1c1"), W3)[0]
    dk3 = L.realised_rows(bulk_k3, P_DE, q=None, kmaxes=(), zero_cols=None)
    bet = L.betas(ctx, held, nonpanel, A_pp_k3={c: L.energy(dk3[c]["design"], nonpanel)["A"] for c in SLOTS})
    fin = FC.final("h2c1")
    bulk = L.key_scale(fin, WV5)[0]
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
    e = rec("h2c1|V1_5_3+sh1.5", np.full(n, 1.5))
    ref = json.loads((W438 / panel / f"comp438-{held}-held.json").read_text())["arms"]["h2c1|V1_5_3+sh1.5"]
    dD, dR = L.cmp(e["design"], ref["design"]), L.cmp(e["realised"]["16"], ref["realised"]["16"])
    log(f"KA h2c1|V1_5_3+sh1.5 vs 438: design |d| {dD:.1e}, realised |d| {dR:.1e}")
    assert dD < TOL and dR < TOL, f"KA FAILED {panel}/{held}: {dD} {dR}"
    res["ka"]["h2c1_V1_5_3_sh1.5"] = dict(design=dD, realised=dR)
    rec("base q1.75", np.full(n, 1.75))
    ref540 = json.loads((ROOT / f"experiments/540-lodo-clip-scan/work/comp540o-{panel}-{held}.json").read_text())["arms"]
    kn = kvv[carried_all] / float(np.median(kvv[carried_all]))
    PCAPS = {f"p{pc:g}": float(np.percentile(kn, pc)) for pc in (90, 95, 96, 97.5, 99)}
    res["pcaps"] = PCAPS
    res["keys"] = {}
    for name, hi in [("fixed 6", 6.0)] + list(PCAPS.items()):
        W = L.weights(kvv, carried_all, 0.5, max(hi, 0.5))[0]
        holder["bulk"] = L.key_scale(fin, W)[0]
        e = rec(name, np.full(n, 1.75))
        res["keys"][name] = dict(hi=hi, frac_cap=float(np.mean(carried_all & (W >= hi - 1e-12))))
        if name == "fixed 6":
            dd = L.cmp(e["realised"]["16"], ref540["V1_0.5_6"]["realised"]["16"])
            log(f"KA fixed 6 vs 540o: realised |d| {dd:.1e}")
            assert dd < TOL, f"KA FAILED {panel}/{held}: {dd}"
    holder["bulk"] = bulk
    floor = WV5 <= FLOOR + 1e-12
    res["frac_floor"] = float(np.mean(floor))
    res["peak_rss_gib"] = E062.peak_rss_gib()
    log(f"floor share {res['frac_floor']:.3f}; keys {res['keys']}")
    (ROOT / "experiments/543-quantile-cap/work" / f"comp543o-{panel}-{held}.json").write_text(json.dumps(res, indent=1, default=float))
    np.savez_compressed(ROOT / "experiments/543-quantile-cap/work" / f"arrays543o-{panel}-{held}.npz", names=np.asarray(list(arrays)),
                        d=np.asarray((truth * truth).sum(1)), **{f"{k}__{kk}": v for k, a in arrays.items() for kk, v in a.items()})
    log(f"done; peak RSS {res['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
