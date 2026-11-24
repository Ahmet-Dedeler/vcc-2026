"""522 — one heavy pass per (October-like panel, PLOLO fold): dump everything a per-target KEY arm needs, so every later
arm (caps, percentile caps, the strength key, combinations, sub-panels) runs from the dump without the whitening
operator (~1-2 GB instead of ~6 GB).

The path is comp546o.py's (546, itself 540o/543o's) verbatim: lib438 bound to the panel; the fold's held-out operator;
Donors438; the panel's CD4 banks under the D/E/F rules; the October stack h2c1 (H1 x2, CD4 x1); the fold's own LODO-V1
fit on the h1c1 composition (lodo_bases + lodo_key, weights(kv, carried>=2, lo, hi)); q 1.75 row sharpening on each
frame's gate; render (minnorm); ballast realisation at KMAX 16 against PP's unkeyed h1c1 carriers; L.row_score.

Saved (work/dump522-<panel>-<fold>.npz): targets, sources, fin (h2c1 final, panel zeroed), fPP (h1c1 final, the carriers),
kv (LODO key), carried, avail (per source), V (per-donor LODO terms v_j = <s(b^-j + radd), d_j> / |F|^2), their parts
(|d_j| non-panel, cos(P^-j, d_j), |P^-j|), cos(d_1, d_2), |F|, scored, panel positions, truth (float32, as the fold cache
holds it), betas, frames.

KNOWN ANSWERS (abort unless |d| < 1e-12, L.cmp's fields: cos_agg, pds, pds lam 0.5/1, render A/B), BEFORE saving:
  438's stored `h2c1|V1_5_3+sh1.5`;  546o's `base q1.75`, `fixed 6` and `p96.6` rows;  540o's `V1_0.5_6` row;
  the per-donor terms' mean == lodo_key (1e-12 relative).

    uv run python experiments/522-october-bulk-key/dump522.py <panel> <fold>    # panel in efp efr real; ~6 GB peak
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
W540 = ROOT / "experiments/540-lodo-clip-scan/work"
W546 = ROOT / "experiments/546-october-cap-pertarget/work"
WORK = HERE / "work"
TOL = 1e-12


def main(panel: str, held: str) -> None:
    t0 = time.time()
    log = lambda *a: print(f"[{panel}/{held}]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    ctx = L.bind(panel)
    fold = L.build_fold438(ctx, held, "held")
    scored, pp_pos, wcols = fold["scored"], fold["panel_positions"], fold["wcols"]
    truth_full = fold["truth"]
    truth = truth_full[:, scored]
    targets, sources, op = list(fold["targets"]), list(fold["sources"]), fold["op"]
    nonpanel = ctx.nonpanel
    dn = L.Donors438(targets, sources)
    cd4, partial, h1, _ = L.residuals(ctx, targets)
    pre = L.pre_normrel(dn, sources, "cap", -1.5, fallback_rows=fold["pre"])
    whit = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    FC = L.Finals(ctx, whit, cd4, partial, h1, wcols, pp_pos)
    key410, carried410, _ = L.shared_energy_key(targets, sources, op, wcols, cols=None)
    W3 = L.weights(key410, carried410, 0.3, 3.0)[0]
    d_rows = {s: L.donor_rows(s, targets)[0] for s in sources}
    A_ = {s: dn.A[s] for s in sources}
    carried_all = sum(A_[s].astype(int) for s in sources) >= 2
    bases = L.lodo_bases(op, dn, sources, "cap", -1.5, pre, wcols)
    fPP, pPP = FC.final("h1c1"), FC.parts("h1c1")
    kv1 = L.lodo_key(fPP, pPP, bases, d_rows, A_, carried_all, sources, pp_pos)
    kvv = np.nan_to_num(kv1, nan=0.0)
    # ---- per-donor LODO terms (lodo_key's arithmetic, cols=None), and their parts
    pp = np.einsum("ij,ij->i", fPP, fPP)
    ns = len(sources)
    V = np.full((ns, len(targets)), np.nan)
    dnorm, pnorm, cosPd = (np.full((ns, len(targets)), np.nan) for _ in range(3))
    for a, j in enumerate(sources):
        Pj = pPP["s"] * (bases[j] + pPP["radd"])
        Pj[:, pp_pos] = 0.0
        num = np.einsum("ij,ij->i", Pj, d_rows[j])
        V[a] = np.where(A_[j], num / pp, np.nan)
        dj = d_rows[j].copy()
        dj[:, pp_pos] = 0.0
        nd, npj = np.linalg.norm(dj, axis=1), np.linalg.norm(Pj, axis=1)
        dnorm[a] = np.where(A_[j], nd, np.nan)
        pnorm[a] = np.where(A_[j], npj, np.nan)
        cosPd[a] = np.where(A_[j], num / np.maximum(nd * npj, 1e-300), np.nan)
    V[:, ~carried_all] = np.nan
    with np.errstate(all="ignore"):
        kmean = np.nan_to_num(np.nanmean(V, axis=0), nan=0.0)
    dmax = float(np.max(np.abs(kmean[carried_all] - kvv[carried_all]) / np.maximum(np.abs(kvv[carried_all]), 1e-300)))
    log(f"KA per-donor mean == lodo_key: max rel |d| {dmax:.1e}")
    assert dmax < 1e-12, "KA per-donor terms"
    d0, d1 = d_rows[sources[0]].copy(), d_rows[sources[1]].copy()
    d0[:, pp_pos] = 0.0
    d1[:, pp_pos] = 0.0
    both = A_[sources[0]] & A_[sources[1]]
    cos_dd = np.where(both, np.einsum("ij,ij->i", d0, d1)
                      / np.maximum(np.linalg.norm(d0, axis=1) * np.linalg.norm(d1, axis=1), 1e-300), np.nan)
    del bases, d_rows, d0, d1
    P_DE = L.carrier_profiles(fPP)
    bulk_k3 = L.key_scale(fPP, W3)[0]
    dk3 = L.realised_rows(bulk_k3, P_DE, q=None, kmaxes=(), zero_cols=None)
    bet = L.betas(ctx, held, nonpanel, A_pp_k3={c: L.energy(dk3[c]["design"], nonpanel)["A"] for c in SLOTS})
    fin = FC.final("h2c1")
    frames = {c: (np.asarray(b), np.asarray(g)) for c, (b, g) in arms414.frames().items()}
    holder = {}

    def rows_q(q):
        out = {}
        for c, (b, gate) in frames.items():
            F = L.sharpen(holder["bulk"], float(q), gate, mode="row")
            des, rea = np.empty_like(F), np.empty_like(F)
            for i in range(F.shape[0]):
                p_bk = profile_from_effect(F[i], b, gate, "minnorm")
                des[i] = np.log1p(TS * p_bk) - b
                p_re, _, _, _, _ = L.realise_k(P_DE[c][i], p_bk, 16)
                rea[i] = np.log1p(TS * p_re) - b
            out[c] = dict(design=des, realised=rea)
        return out

    def rec(w, q):
        holder["bulk"] = L.key_scale(fin, w)[0]
        dr = rows_q(q)
        sd, _ = L.row_score({c: dr[c]["design"] for c in SLOTS}, scored, truth, bet, nonpanel)
        sr, _ = L.row_score({c: dr[c]["realised"] for c in SLOTS}, scored, truth, bet, nonpanel)
        return sd, sr

    ka = {}
    ref438 = json.loads((W438 / panel / f"comp438-{held}-held.json").read_text())["arms"]["h2c1|V1_5_3+sh1.5"]
    ref546 = json.loads((W546 / f"comp546o-{panel}-{held}.json").read_text())["arms"]
    ref540 = json.loads((W540 / f"comp540o-{panel}-{held}.json").read_text())["arms"]
    kn = kvv[carried_all] / float(np.median(kvv[carried_all]))
    hi546 = float(np.percentile(kn, 96.647))
    checks = [("438 h2c1|V1_5_3+sh1.5", L.weights(kvv, carried_all, 0.5, 3.0)[0], 1.5, ref438, True),
              ("546 base q1.75", L.weights(kvv, carried_all, 0.5, 3.0)[0], 1.75, ref546["base q1.75"], False),
              ("546 fixed 6", L.weights(kvv, carried_all, 0.5, 6.0)[0], 1.75, ref546["fixed 6"], False),
              ("540o V1_0.5_6", L.weights(kvv, carried_all, 0.5, 6.0)[0], 1.75, ref540["V1_0.5_6"], False),
              ("546 p96.6", L.weights(kvv, carried_all, 0.5, hi546)[0], 1.75, ref546["p96.6"], False)]
    for name, w, q, ref, with_design in checks:
        sd, sr = rec(w, q)
        dR = L.cmp(sr, ref["realised"]["16"])
        dD = L.cmp(sd, ref["design"]) if with_design else 0.0
        ka[name] = dict(realised=dR, design=dD, cos_agg=sr["cos_agg"], pds=sr["pds"])
        log(f"KA {name:24s} cos {sr['cos_agg']:.9f} pds {sr['pds']:.6f}: realised |d| {dR:.1e} design |d| {dD:.1e}")
        assert dR < TOL and dD < TOL, f"KA FAILED {panel}/{held} {name}: {dR} {dD}"
    out = WORK / f"dump522-{panel}-{held}.npz"
    np.savez_compressed(out, targets=np.asarray(targets), sources=np.asarray(sources), fin=fin, fPP=fPP, kv=kvv,
                        carried=carried_all, avail=np.stack([np.asarray(A_[s], bool) for s in sources]), V=V,
                        dnorm=dnorm, pnorm=pnorm, cosPd=cosPd, cos_dd=cos_dd, Fnorm=np.sqrt(pp), scored=np.asarray(scored),
                        panel_positions=np.asarray(pp_pos), truth=np.asarray(truth_full, np.float32),
                        **{f"beta_{c}": bet[c] for c in SLOTS}, **{f"b_{c}": frames[c][0] for c in SLOTS},
                        **{f"gate_{c}": frames[c][1] for c in SLOTS})
    rec_ = dict(panel=panel, fold=held, n=len(targets), sources=sources, ka=ka, hi546=hi546,
                n_carried=int(carried_all.sum()), peak_rss_gib=E062.peak_rss_gib(), seconds=time.time() - t0)
    (WORK / f"dump522-{panel}-{held}.json").write_text(json.dumps(rec_, indent=1, default=float))
    log(f"dumped {out.name}; peak RSS {rec_['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
