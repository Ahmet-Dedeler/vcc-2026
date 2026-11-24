"""430 — DE-gene SIGNS on a held-out line: which composition gets them right? One PLOLO fold per process.

Truth = the held-out line L's own donor delta row (pseudobulk, pert minus control; primary = RAW, i.e. what a DE test
against controls reads; secondary = CENTRED on L's all-perturbation column mean, i.e. PLOLO's `truth`).
Gate proxies per target (universe = PLOLO's scored columns: L measures the gene, panel columns out; plus L's control
CPM > 5, the scorer's expression floor; kolf2's bank has no control profile, so no floor there):
  top-K by |delta_L|, K in {50, 200, 1000};  |delta_L| >= tau, tau in {0.1, 0.2} (log1p-CP50k units; not for kolf2's
  z-scores).
  zK (K 200, 1000): top-K by |delta_L| * sqrt(L's control CPM), a significance proxy; on H1 it ranks 23 sign candidates
  like the real reference DE gate (Spearman 0.94 for z1000, 24 candidates), where top-50 |delta| does not (0.22) (h1bridge430.py).
Per candidate x gate: per-target sign hits (a zero prediction scores 1/2) and gate size, and the L1 nmae proxy
sum|k*lhat - delta| / sum|delta| averaged over targets, at the candidate's own L1-optimal global k (exact: the weighted
median of delta/lhat with weights |lhat| / S_p).

Candidates, every one built WITHOUT L (sources = FOLDS[L]["sources"], the fold's held operator):
  PP       the production composition (418/420 "PP"): normrel cap -1.5 + PM50 donor base, whitened at gamma 0.65,
           + CD4 (deflated k5) + H1 union + CD4 partial residual, one post-add preserve, panel columns zeroed.
  PPg0     the same with the whitening skipped (gamma 0): the residual doses re-solved against the unwhitened base,
           each row norm-matched to PP's row (416's g000 construction).
  Wonly    PP's whitened donor base with no residual.       pre   the unwhitened donor base (normrel + PM), no residual.
  PPxs     (Orion folds) PP with the donor base restricted to the cross-study source (k562_gwps): tests whether the
           same-study partner line inflates PP.
  M3L      equal-weight mean of all other donor tables (data/donors/*: k562_gwps, k562_essential, rpe1, nadig_jurkat,
           orion_hct116, orion_hek293t minus L; never h1_2025), RAW delta, nanmean over donors covering (p, g).
  M3S      the same excluding L's whole study (Replogle 2022 = k562_gwps/k562_essential/rpe1; Orion = hct116/hek293t).
  M3Lc/M3Sc  the same on CENTRED donor rows (each donor's all-perturbation column mean removed, as production does).
  C4L/C4S  the CLEAN4 per-gene prior sign (371's P_banks; mean P_down over K562-GWPS, K562-ess, Nadig Jurkat, Nadig
           HepG2 minus L's bank / minus L's study), -1 where uncovered (the builder's rule), constant across targets.
  PPprior  PP with the CLEAN4 (LOLO) sign on L's CPM block (top 2,200 expressed genes by L's control CPM) except the
           strongest 20 % |PP| per target (an emulation of the delivered clean4q80 sign structure; H1 bridge checks it).
  PPxsg0   (Orion folds) PPxs with the whitening skipped (gamma 0), rows norm-matched to PPxs.
  BXw{w}   (Orion folds) 0.5*PPxs + 0.5*w*M3S_rowscaled (a cross-study-only blend).
  B{3}w{w} 0.5*PP + 0.5*w*M3x_rowscaled, M3x each row scaled to PP's row norm (non-panel columns), w in {1, 1.5, 2}.

KNOWN ANSWERS (abort on failure, before any candidate number):
  KA-PP     compose_arm(defaults) == sens415.compose_pp bit-for-bit, and the effect-level aggregate cosine of PP with
            PLOLO's truth on the live scored rows == sens415.json[L]["base_agg"] (1e-12).
  KA-DLV    PP delivered through the a23 read-back in frames A/B/C reproduces 418's stored per-target arrays
            `PP|none__x/__a` (P.T and P.P per target per frame) and 420's `PP|K3__x/__a` (the 410 key applied), |d| 1e-12.
  KA-TRUTH  my raw delta_L - colmean_L == PLOLO's centred truth on the scored columns (float32 tolerance).

    uv run python fold430.py <fold>     # k562_gwps | orion_hct116 | orion_hek293t | kolf2
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True           # importing 414/415/418 modules must not drop __pycache__ into their dirs
import json  # noqa: E402
import os  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/418-bulk-composition"))
import lib418  # noqa: E402
from lib418 import (ALPHA, E062, FOLDS, GAMMA, LAM, PROD_OPS, UNION, W415, Donors, Operator, agg, agg_parts,  # noqa: E402
                    betas_for, build_fold, cd4_frozen_and_partial, compose_arm, compose_pp, deflate, key_scale,
                    pre_normrel, score, source_residual, union_rows)
from plolo import OFFICIAL, PANEL, TS  # noqa: E402

W418 = ROOT / "experiments/418-bulk-composition/work"
W420 = ROOT / "experiments/420-bulk-stack/work"
DON = ("k562_gwps", "k562_essential", "rpe1", "nadig_jurkat", "orion_hct116", "orion_hek293t")
STUDY = {"k562_gwps": "replogle22", "k562_essential": "replogle22", "rpe1": "replogle22", "nadig_jurkat": "nadig25",
         "nadig_hepg2": "nadig25", "orion_hct116": "orion", "orion_hek293t": "orion", "kolf2": "kolf2"}
CLEAN4 = ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")
KS = (50, 200, 1000)
ZKS = (200, 1000)          # top-K by |delta_L| * sqrt(control mean expression): Spearman 0.82 / 0.94 with the real gate on H1
TAUS = (0.1, 0.2)
BLOCK = 2200
PROTECT_Q = 0.8
WS = (1.0, 1.5, 2.0)
TOL = 1e-12


def l1_opt(pairs):
    """Exact L1-optimal k >= 0 for mean_p sum_g |k lhat - d| / S_p. pairs = [(lhat, d)] per target (gate values)."""
    r, w = [], []
    for lh, d in pairs:
        S = np.abs(d).sum()
        nz = lh != 0
        r.append(d[nz] / lh[nz])
        w.append(np.abs(lh[nz]) / S)
    r, w = np.concatenate(r), np.concatenate(w)
    if r.size == 0:
        return 0.0
    o = np.argsort(r)
    cw = np.cumsum(w[o])
    k = float(r[o][np.searchsorted(cw, 0.5 * cw[-1])])
    return max(k, 0.0)


def main(held: str) -> None:
    t0 = time.time()
    log = lambda *a: print(f"[{held}]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    fold = build_fold(held)
    scored = np.asarray(fold["scored"])
    pp_pos = np.asarray(fold["panel_positions"])
    wcols = np.asarray(fold["wcols"])
    truth_c = np.asarray(fold["truth"], dtype=np.float64)          # centred, full official axis, panel zeroed
    targets = list(fold["targets"])
    sources = list(FOLDS[held]["sources"])
    op = Operator(PROD_OPS / "whitener_union_pruned.npz" if FOLDS[held]["op_kind"] == "production"
                  else FOLDS[held]["op"])
    nonpanel = np.ones(len(OFFICIAL), bool)
    nonpanel[pp_pos] = False
    res = dict(fold=held, n_targets=len(targets), sources=sources, ka={}, gates={}, cands={})

    # ------------------------------------------------------------------ the PP composition (418's KA-415 path)
    dn = Donors(targets, sources)
    frozen, partial = cd4_frozen_and_partial(targets)
    cd4_5, _ = deflate(frozen, 5)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    pre_cached = np.asarray(fold["pre"], dtype=np.float64)
    pre228 = dn.pre(sources, False)
    pre228[:, pp_pos] = 0.0
    nocov = sum(dn.A[s].astype(int) for s in sources) == 0
    pre228[nocov] = pre_cached[nocov]
    fb = pre228.copy()
    pre = pre_normrel(dn, sources, "cap", -1.5, fallback_rows=fb)
    fpp, _, c415, radd415, _ = compose_pp(op, pre, cd4_5, partial, h1, wcols, pp_pos)
    whit = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    fpp2, parts = compose_arm(whit, cd4_5, partial, h1, wcols, pp_pos)
    assert np.abs(fpp2 - fpp).max() == 0.0 and parts["s"] == c415["s"], "KA-PP FAILED: compose_arm != compose_pp"
    T = truth_c[:, scored]
    P = fpp[:, scored]
    live = (np.abs(T).sum(1) > 0) & (np.abs(P).sum(1) > 0)
    pt, pp_, tt = agg_parts(P[live], T[live])
    base_agg = agg(np.ones(int(live.sum())), pt, pp_, tt)
    s415 = json.loads((W415 / "sens415.json").read_text())[held]
    d = abs(base_agg - s415["base_agg"])
    log(f"KA-PP base_agg {base_agg:.12f} vs sens415 {s415['base_agg']:.12f} |d| {d:.1e}")
    assert d < TOL, "KA-PP FAILED: base_agg"
    res["ka"]["base_agg"] = [base_agg, s415["base_agg"]]
    pcos = pt / np.sqrt(pp_ * (T[live] ** 2).sum(1))
    res["ka"]["pp_percos_median"] = float(np.median(pcos))

    # ------------------------------------------------------------------ KA-DLV: delivered per-target arrays
    betas = betas_for(held, nonpanel)
    a418 = np.load(W418 / f"arrays418-{held}.npz")
    a420 = np.load(W420 / f"arrays420-{held}.npz")
    worst = 0.0
    _, arr = score(fpp, fold, T, betas, None)
    for k in ("x", "a"):
        worst = max(worst, float(np.abs(arr[k] - a418[f"PP|none__{k}"]).max() / np.abs(a418[f"PP|none__{k}"]).max()))
    W3 = np.asarray(a418["W3"])
    _, arr3 = score(key_scale(fpp, W3)[0], fold, T, betas, None)
    for k in ("x", "a"):
        worst = max(worst, float(np.abs(arr3[k] - a420[f"PP|K3__{k}"]).max() / np.abs(a420[f"PP|K3__{k}"]).max()))
    dcos = arr["x"] / np.sqrt(arr["a"] * a418["d"][None, :])
    log(f"KA-DLV PP|none (418) and PP|K3 (420) per-target P.T, P.P in A/B/C: worst rel |d| {worst:.1e}; "
        f"delivered per-target cos median A/B/C {np.round(np.median(dcos, 1), 5).tolist()}")
    assert worst < TOL, "KA-DLV FAILED"
    res["ka"]["dlv_worst_rel"] = worst
    res["ka"]["dlv_percos_median_ABC"] = np.median(dcos, 1).tolist()
    del a418, a420, arr, arr3

    # ------------------------------------------------------------------ truth (raw) and KA-TRUTH
    tix = np.asarray([PANEL.index(t) for t in targets])
    Z = {n: np.load(WORK / f"donor430-{n}.npz") for n in (*DON, "kolf2")}
    zL = Z[held]
    raw_L = np.asarray(zL["raw"][tix], dtype=np.float64)
    cm_L = np.asarray(zL["colmean"])
    assert zL["avail"][tix].all()
    cen_L = raw_L - cm_L[None, :]
    dt = float(np.abs(cen_L[:, scored] - truth_c[:, scored]).max())
    log(f"KA-TRUTH raw_L - colmean_L vs PLOLO truth on scored: max |d| {dt:.2e} (truth max {np.abs(T).max():.2f})")
    assert dt < 1e-4, "KA-TRUTH FAILED"
    res["ka"]["truth_max_abs_d"] = dt
    assert np.isfinite(raw_L[:, scored]).all()

    # universe: scored, expressed in L (CPM > 5) where L has a control profile
    if held == "kolf2":
        expr = np.ones(len(OFFICIAL), bool)
        cpmL = None
    else:
        bL = np.asarray(fold["b_ctrl"], dtype=np.float64)
        lin = np.expm1(bL)
        cpmL = 1e6 * lin / lin.sum()
        expr = cpmL > 5.0
    uni = np.zeros(len(OFFICIAL), bool)
    uni[scored] = True
    uni &= expr
    U = np.flatnonzero(uni)
    res["universe"] = dict(n_scored=int(scored.size), n_expressed=int(U.size))
    # L's CPM block (for the PPprior emulation): top-2200 expressed scored genes by L's control CPM
    block = U[np.argsort(-cpmL[U])[:BLOCK]] if cpmL is not None else U[:0]

    # ------------------------------------------------------------------ candidates
    C = {}
    C["PP"] = fpp
    whit0 = op.apply(pre * ALPHA, wcols, 0.0, LAM)
    g0, parts0 = compose_arm(whit0, cd4_5, partial, h1, wcols, pp_pos)
    nr = np.linalg.norm(fpp, axis=1) / np.maximum(np.linalg.norm(g0, axis=1), 1e-300)
    C["PPg0"] = g0 * nr[:, None]
    res["g0"] = dict(fractions=parts0["fractions"], fractions_pp=parts["fractions"])
    w_only = whit.copy()
    w_only[:, pp_pos] = 0.0
    C["Wonly"] = w_only
    p_only = pre * ALPHA
    p_only[:, pp_pos] = 0.0
    C["pre"] = p_only
    if held in ("orion_hct116", "orion_hek293t"):
        pre_xs = pre_normrel(dn, ["k562_gwps"], "cap", -1.5, fallback_rows=None)   # no partner-line fallback
        C["PPxs"], _ = compose_arm(op.apply(pre_xs * ALPHA, wcols, GAMMA, LAM), cd4_5, partial, h1, wcols, pp_pos)
        res["ppxs_uncovered"] = int((~dn.A["k562_gwps"]).sum())
        gx, _ = compose_arm(op.apply(pre_xs * ALPHA, wcols, 0.0, LAM), cd4_5, partial, h1, wcols, pp_pos)
        nx = np.linalg.norm(C["PPxs"], axis=1) / np.maximum(np.linalg.norm(gx, axis=1), 1e-300)
        C["PPxsg0"] = gx * nx[:, None]          # the cross-study gamma-0 composition, rows norm-matched to PPxs

    def dmean(names, centred):
        acc = np.zeros((len(targets), len(OFFICIAL)))
        cnt = np.zeros((len(targets), len(OFFICIAL)))
        for n in names:
            r = np.asarray(Z[n]["raw"][tix], dtype=np.float64)
            if centred:
                r = r - np.asarray(Z[n]["colmean"])[None, :]
            ok = np.isfinite(r)
            acc += np.where(ok, r, 0.0)
            cnt += ok
        m = np.where(cnt > 0, acc / np.maximum(cnt, 1), 0.0)
        m[:, pp_pos] = 0.0
        return m, cnt

    others = [n for n in DON if n != held]
    xstudy = [n for n in DON if STUDY[n] != STUDY[held]]
    res["m3"] = dict(literal=others, xstudy=xstudy,
                     covering_literal={n: int(Z[n]["avail"][tix].sum()) for n in others})
    C["M3L"], cntL = dmean(others, False)
    C["M3S"], _ = dmean(xstudy, False)
    C["M3Lc"], _ = dmean(others, True)
    C["M3Sc"], _ = dmean(xstudy, True)
    res["m3"]["targets_uncovered_literal"] = int((cntL.max(1) == 0).sum())

    cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    bnames = list(cp["bank_names"].astype(str))
    PB = cp["P_banks"].astype(float)
    gpos = {g: j for j, g in enumerate(OFFICIAL)}
    ccol = np.asarray([gpos.get(str(g), -1) for g in cp["genes"].astype(str)])

    def clean4(banks):
        P = PB[[bnames.index(b) for b in banks]].mean(0)
        s = np.full(len(OFFICIAL), -1.0)
        ok = ccol >= 0
        s[ccol[ok]] = np.where(P[ok] > 0.5, -1.0, 1.0)
        return s

    c4L_banks = [b for b in CLEAN4 if b != held]
    c4S_banks = [b for b in CLEAN4 if STUDY[b] != STUDY[held]]
    res["clean4"] = dict(literal=c4L_banks, xstudy=c4S_banks)
    for tag, banks in (("C4L", c4L_banks), ("C4S", c4S_banks)):
        s = clean4(banks)
        M = np.tile(s, (len(targets), 1))
        M[:, pp_pos] = 0.0
        C[tag] = M
    # PPprior: CLEAN4 (literal LOLO) sign on L's block, strongest 20 % |PP| per target protected
    s4 = clean4(c4L_banks)
    pr = fpp.copy()
    if block.size:
        for i in range(len(targets)):
            v = fpp[i, block]
            thr = np.quantile(np.abs(v), PROTECT_Q)
            flip = np.abs(v) <= thr
            pr[i, block[flip]] = np.abs(v[flip]) * s4[block[flip]]
    C["PPprior"] = pr

    npp = np.linalg.norm(fpp[:, nonpanel], axis=1)
    for m3 in ("M3L", "M3S", "M3Sc"):
        n3 = np.linalg.norm(C[m3][:, nonpanel], axis=1)
        sc = np.where(n3 > 0, npp / np.maximum(n3, 1e-300), 0.0)
        R3 = C[m3] * sc[:, None]
        for w in WS:
            C[f"B{m3[2:]}w{w:g}"] = 0.5 * fpp + 0.5 * w * R3

    if "PPxs" in C:                             # cross-study blend on the cross-study composition
        nx = np.linalg.norm(C["PPxs"][:, nonpanel], axis=1)
        n3 = np.linalg.norm(C["M3S"][:, nonpanel], axis=1)
        R3 = C["M3S"] * np.where(n3 > 0, nx / np.maximum(n3, 1e-300), 0.0)[:, None]
        for w in WS:
            C[f"BXw{w:g}"] = 0.5 * C["PPxs"] + 0.5 * w * R3

    # ------------------------------------------------------------------ evaluation
    rows = {"raw": raw_L, "cen": cen_L}
    covered = np.asarray([any(dn.A[s][i] for s in sources) for i in range(len(targets))])
    res["targets_covered_by_sources"] = int(covered.sum())
    evalt = np.flatnonzero(covered)
    out = {}
    for tname, TR in rows.items():
        gates = {f"top{k}": ("top", k) for k in KS}
        if held != "kolf2":
            gates.update({f"tau{t:g}": ("tau", t) for t in TAUS})
            gates.update({f"z{k}": ("z", k) for k in ZKS})     # H1-validated significance proxy (h1bridge430)
        gidx = {}
        for g, (kind, v) in gates.items():
            gl = []
            for i in evalt:
                dv = TR[i, U]
                nz = dv != 0
                cand = U[nz]
                a = np.abs(TR[i, cand])
                if kind == "top":
                    sel = cand[np.argsort(-a, kind="stable")[:v]]
                elif kind == "z":
                    sel = cand[np.argsort(-(a * np.sqrt(cpmL[cand])), kind="stable")[:v]]
                else:
                    sel = cand[a >= v]
                gl.append(sel)
            gidx[g] = gl
            res["gates"][f"{tname}:{g}"] = dict(mean_size=float(np.mean([x.size for x in gl])),
                                               targets_ge10=int(sum(x.size >= 10 for x in gl)))
        for cname, M in C.items():
            for g, gl in gidx.items():
                hits, ns, pairs = [], [], []
                for j, i in enumerate(evalt):
                    sel = gl[j]
                    lh, dv = M[i, sel], TR[i, sel]
                    h = np.where(lh == 0, 0.5, (np.sign(lh) == np.sign(dv)).astype(float))
                    hits.append(h.sum())
                    ns.append(sel.size)
                    pairs.append((lh, dv))
                k = l1_opt([p for p in pairs if p[1].size])
                nm = np.asarray([np.abs(k * lh - dv).sum() / np.abs(dv).sum() if dv.size else np.nan
                                 for lh, dv in pairs])
                out[f"{tname}|{cname}|{g}"] = dict(hits=np.asarray(hits), n=np.asarray(ns), nmae=nm, k=k)
        log(f"evaluated truth={tname}")
    # summary print (pooled accuracy; raw truth)
    for g in ("z1000", "z200", "top50", "top200", "top1000", "tau0.1", "tau0.2"):
        line = []
        for cname in C:
            key = f"raw|{cname}|{g}"
            if key not in out:
                continue
            o = out[key]
            line.append(f"{cname} {o['hits'].sum() / max(o['n'].sum(), 1):.4f}/{np.nanmean(o['nmae']):.4f}")
        if line:
            log(f"raw {g}: " + "  ".join(line))
    res["cands"] = {k: dict(pooled=float(v["hits"].sum() / max(v["n"].sum(), 1)),
                            pertarget=float(np.mean((v["hits"] / np.maximum(v["n"], 1))[v["n"] >= 10]))
                            if (v["n"] >= 10).any() else None,
                            nmae=float(np.nanmean(v["nmae"])), k=v["k"], n_targets=int(v["n"].size),
                            n_ge10=int((v["n"] >= 10).sum()))
                    for k, v in out.items()}
    res["peak_rss_gib"] = E062.peak_rss_gib()
    (WORK / f"fold430-{held}.json").write_text(json.dumps(res, indent=1))
    cov_xs = np.asarray([any(Z[n]["avail"][tix[i]] for n in xstudy) for i in range(len(targets))])
    np.savez_compressed(WORK / f"fold430-{held}.npz", evalt=evalt, targets=np.asarray(targets), cov_xs=cov_xs[evalt],
                        **{f"{k}__{f}": v[f] for k, v in out.items() for f in ("hits", "n", "nmae")},
                        **{f"{k}__k": np.float64(v["k"]) for k, v in out.items()})
    log(f"done; peak RSS {res['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main(sys.argv[1])
