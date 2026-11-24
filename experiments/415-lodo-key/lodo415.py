"""415 — leave-one-donor-out per-target key for the PRODUCTION-COMPOSITION bulk row, validated on 162's PLOLO folds
against 410's key evaluated on the same prediction. Preregistration: PREREGISTRATION.md (fixed before this ran).

    key_p = mean_j <P^{-j}_p, d_{j,p}> / |P_p|^2          (mean over source donors j carrying p; carried = >= 2 donors)

P      = 228's build_final(m=None): c*W(ALPHA*pre) on the operator axis (ALPHA*pre off it) + frac_cd4*CD4 + frac_u*H1,
         times the post-add preserve s, panel zeroed.
P^{-j} = the same with pre -> pre^{-j} (mean of the OTHER carrying sources' prepared rows), constants c, frac_cd4,
         frac_u, s held at P's values (linear LOO; V1/V2), or everything recomputed through the pipeline (V3).
d_j    = V1: donor j's raw bank row (410's donor_rows); V2/V3: donor j's prepared row (the row `pre` averages).

Known answers (abort on failure): KA-1 check410.json, KA-2 cached pre, KA-3 228/223 base222 cos+pds, KA-4 algebra.
Output: work/lodo415.json (+ printed tables).

    uv run python experiments/415-lodo-key/lodo415.py
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
for p in ("experiments/410-pertscale-revisit", "experiments/162-retrodicting-validator",
          "experiments/218-independent-donors"):
    sys.path.insert(0, str(ROOT / p))
import plolo  # noqa: E402

from resid218 import ALPHA, GAMMA, LAM, score_effects, source_residual  # noqa: E402  (sets plolo.WORK; reset below)
from compose218 import merged, solve_fraction_carrying, union_rows  # noqa: E402

plolo.WORK = WORK            # challenge b_ctrl caches were copied here; nothing outside 415 is written
from plolo import E062, FOLDS, OFFICIAL, PANEL, cd4_for, solve_fraction  # noqa: E402
from key410 import raw_whiten, shared_energy_key, weights  # noqa: E402
from size410 import donor_rows  # noqa: E402

CD4_F, UNION_G = 0.34, 0.28                       # arms228 constants (not imported: arms228 repoints plolo.WORK)
UNION = ("gse295214_h1", "h1_2025")
KA228 = {  # arms228.KA, 223 focus228 base222 (a23 cos, a23 pds)
    "k562_gwps": (0.03943565390646537, 0.8053968077826482),
    "orion_hct116": (0.032080093772601896, 0.7922125442469077),
    "orion_hek293t": (0.026156872210626163, 0.7636228552854588),
    "kolf2": (0.012520242539233312, 0.6361551030668678),
}
RULES = ((0.3, 3.0), (0.3, 1e9))
NBOOT, NPERM = 1000, 200
PREP_NAME = {"k562_gwps": "k562_gwps_pm50", "orion_hct116": "orion_hct116", "orion_hek293t": "orion_hek293t"}
CHECK410 = json.loads((ROOT / "experiments/410-pertscale-revisit/work/check410.json").read_text())


def ranks(x):
    return np.argsort(np.argsort(x))


def spearman(a, b):
    return float(np.corrcoef(ranks(a), ranks(b))[0, 1])


def agg_parts(P, T):
    return (P * T).sum(1), (P * P).sum(1), float((T * T).sum())


def agg(w, pt, pp, tt, idx=None):
    if idx is None:
        return float((w * pt).sum() / np.sqrt((w * w * pp).sum() * tt))
    return float((w[idx] * pt[idx]).sum() / np.sqrt((w[idx] ** 2 * pp[idx]).sum() * tt))


def whiten_split(op, M, wcols):
    """Production Operator.apply decomposed: returns (c, W M on the axis without preserve). apply == M off-axis,
    c * WM on it. Asserted against op.apply by the caller."""
    WM = raw_whiten(op, M, wcols)
    before = float(np.median(np.linalg.norm(np.asarray(M, float)[:, wcols], axis=1)))
    after = float(np.median(np.linalg.norm(WM[:, wcols], axis=1)))
    return before / after, WM


def compose(op, pre, cd4, h1, wcols, panel_positions):
    """228 build_final(m=None, a=1) with every constant returned."""
    whitened = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    frac_cd4, _ = solve_fraction(whitened, cd4, CD4_F)
    frac_u, _ = solve_fraction_carrying(whitened, h1, UNION_G)
    res = merged(cd4, frac_cd4, h1, frac_u)
    final, info = E062.add_post_whiten_residual(whitened, res, 1.0, wcols)
    final[:, panel_positions] = 0.0
    radd = np.where(np.asarray(res.gene_votes) > 0, np.asarray(res.mean_effect), 0.0)
    return final, dict(frac_cd4=float(frac_cd4), frac_u=float(frac_u),
                       s=float(info["post_add_global_preserve_scale"])), radd, whitened


def boot(keys_w, pt, pp, tt, rng, names):
    n = len(pt)
    one = np.ones(n)
    out = {k: [] for k in names}
    diffs = {k: [] for k in names if k != "K410"}
    for _ in range(NBOOT):
        idx = rng.integers(0, n, n)
        b = agg(one, pt, pp, tt, idx)
        g = {k: agg(keys_w[k], pt, pp, tt, idx) / b - 1 for k in names}
        for k in names:
            out[k].append(g[k])
            if k != "K410":
                diffs[k].append(g[k] - g["K410"])
    ci = lambda v: [float(np.quantile(v, .025)), float(np.quantile(v, .975))]
    return {k: ci(v) for k, v in out.items()}, {k: ci(v) for k, v in diffs.items()}, \
        {k: float(np.mean(np.asarray(v) > 0)) for k, v in diffs.items()}


def main() -> None:
    t_all = time.time()
    rng = np.random.default_rng(415)
    pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([pos[g] for g in PANEL if g in pos], dtype=np.int64)
    prep = {}
    for src, nm in PREP_NAME.items():
        z = np.load(WORK / f"prep-{nm}.npz", allow_pickle=True)
        prep[src] = dict(eff=z["effects"], avail=z["available"], mean_on_axis=z["mean_on_axis"],
                         tidx={str(t): i for i, t in enumerate(z["targets"])})
    res = {}
    for held, fold in FOLDS.items():
        t0 = time.time()
        z = np.load(ROOT / f"experiments/162-retrodicting-validator/work/fold-{held}.npz", allow_pickle=True)
        pre = np.asarray(z["pre"], float)
        truth = np.asarray(z["truth"], float)
        scored = np.asarray(z["scored"])
        wcols = np.asarray(z["wcols"])
        targets = [str(t) for t in z["targets"]]
        sources = list(fold["sources"])
        op = plolo.Operator(fold["op"])
        out = dict(held=held, sources=sources)

        # ---------------- KA-1: 410's key on the ORIGINAL fold prediction reproduces check410.json
        P0 = op.apply(pre, wcols, GAMMA, 0.0)[:, scored]
        T0 = truth[:, scored]
        live0 = (np.abs(T0).sum(1) > 0) & (np.abs(P0).sum(1) > 0)
        tg0 = [t for t, k in zip(targets, live0) if k]
        pt0, pp0, tt0 = agg_parts(P0[live0], T0[live0])
        k0, car0, _ = shared_energy_key(tg0, sources, op, wcols, cols=scored)
        rho0 = spearman(k0[car0], (pt0 / pp0)[car0])
        base0 = agg(np.ones(len(tg0)), pt0, pp0, tt0)
        ka1 = dict(rho=rho0, rho_410=CHECK410[held]["spearman_sstar"], gains={})
        assert abs(rho0 - ka1["rho_410"]) < 1e-9, f"KA-1 rho {held}"
        for lo, hi in RULES:
            w, _ = weights(k0, car0, lo, hi)
            g = agg(w, pt0, pp0, tt0) / base0 - 1
            g410 = CHECK410[held]["rules"][f"lo{lo}_hi{hi if hi < 1e8 else 1000000000.0}"]["gain"]
            ka1["gains"][f"lo{lo}_hi{hi:g}"] = [g, g410]
            assert abs(g - g410) < 1e-9, f"KA-1 gain {held} {lo} {hi}: {g} vs {g410}"
        out["KA1"] = ka1
        print(f"[{held}] KA-1 PASS rho {rho0:+.6f}  gains " +
              " ".join(f"{k} {v[0]:+.6f}" for k, v in ka1["gains"].items()), flush=True)

        # ---------------- KA-2: re-prepared source rows reproduce the cached `pre`
        E, A = {}, {}
        for s in sources:
            ii = np.asarray([prep[s]["tidx"][t] for t in targets])
            E[s] = prep[s]["eff"][ii]
            A[s] = prep[s]["avail"][ii]
        votes = sum(A[s].astype(int) for s in sources)
        pre_re = np.zeros_like(pre)
        for s in sources:
            pre_re[A[s]] += E[s][A[s]]
        cov = votes > 0
        pre_re[cov] /= votes[cov, None]
        if (~cov).any():
            mr = sum(prep[s]["mean_on_axis"] for s in sources) / len(sources)
            pre_re[~cov] = E062.EFFECT_SCALE * mr / (float(np.linalg.norm(mr)) or 1.0)
        rel = float(np.abs(pre_re - pre).max() / np.abs(pre).max())
        assert rel < 1e-5, f"KA-2 {held}: rel max diff {rel}"
        out["KA2_rel_maxdiff"] = rel
        print(f"[{held}] KA-2 PASS pre rebuilt, rel max |diff| {rel:.2e} (float32 cache)", flush=True)
        pre = pre_re                                   # float64 exact rows from here on

        # ---------------- production composition + KA-3
        cd4 = cd4_for(targets)
        h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
        final, const, radd, whitened = compose(op, np.asarray(z["pre"], float), cd4, h1, wcols, panel_positions)
        e = score_effects(final, dict(scored=scored), truth[:, scored])
        dc, dp = abs(e["cos"] - KA228[held][0]), abs(e["pds"] - KA228[held][1])
        assert dc < 1e-12 and dp < 1e-12, f"KA-3 {held}: dcos {dc} dpds {dp}"
        out["KA3"] = dict(cos=e["cos"], pds=e["pds"], dcos=dc, dpds=dp)
        print(f"[{held}] KA-3 PASS 228 base222 cos {e['cos']:.15f} (d {dc:.1e}) pds {e['pds']:.15f} (d {dp:.1e})",
              flush=True)
        # the analysis composition on the float64 pre (constants re-derived, reported)
        final, const, radd, whitened = compose(op, pre, cd4, h1, wcols, panel_positions)
        c, WM = whiten_split(op, pre * ALPHA, wcols)
        chk = pre * ALPHA
        chk[:, wcols] = c * WM[:, wcols]
        assert np.allclose(chk, whitened, rtol=1e-10, atol=1e-12), "whiten_split != Operator.apply"
        s_add = const["s"]
        const["c"] = c
        # residual share of the production row on the scored axis (reported, mandate says ~35 % in production)
        Pb = s_add * whitened
        Pb[:, panel_positions] = 0.0
        Rr = s_add * radd
        Rr[:, panel_positions] = 0.0

        P = final[:, scored]
        T = truth[:, scored]
        live = (np.abs(T).sum(1) > 0) & (np.abs(P).sum(1) > 0)
        tg = [t for t, k in zip(targets, live) if k]
        n = len(tg)
        P, T = P[live], T[live]
        pt, pp, tt = agg_parts(P, T)
        sstar = pt / pp
        base = agg(np.ones(n), pt, pp, tt)
        oracle = float(np.sqrt((np.maximum(pt, 0) ** 2 / pp).sum() / tt) / base - 1)
        out.update(n=n, base_agg=base, oracle_gain=oracle, constants=const,
                   residual_norm_share_median=float(np.median(
                       np.linalg.norm(Rr[live][:, scored], axis=1) / np.linalg.norm(P, axis=1))),
                   base_orig_agg=float(base0))

        keys = {}
        # K410 on the production-composition live set
        k410, car410, nd = shared_energy_key(tg, sources, op, wcols, cols=scored)
        keys["K410"] = (k410, car410)

        # LOO rows (live targets only), linear construction
        Al = {s: A[s][live] for s in sources}
        El = {s: E[s][live] for s in sources}
        nsrc = sum(Al[s].astype(int) for s in sources)
        carried = nsrc >= 2
        assert np.array_equal(carried, car410), "carried sets differ between K410 and LODO"
        raw = {s: donor_rows(s, tg) for s in sources}
        for s in sources:
            assert np.array_equal(raw[s][1], Al[s]), f"availability mismatch {s}"
        Pm = {}          # P^{-j}, linear, scored columns
        Pm_base = {}     # its base-only part
        Pm_full = {}     # V3: full pipeline rebuild
        for j in sources:
            prej = np.zeros((n, pre.shape[1]))
            cnt = np.zeros(n, int)
            for i in sources:
                if i == j:
                    continue
                prej[Al[i]] += El[i][Al[i]]
                cnt += Al[i].astype(int)
            ok = cnt > 0
            prej[ok] /= cnt[ok, None]
            bj = prej * ALPHA
            WMj = raw_whiten(op, bj, wcols)
            bj[:, wcols] = c * WMj[:, wcols]
            basej = s_add * bj
            basej[:, panel_positions] = 0.0
            Pm_base[j] = basej[:, scored]
            full = s_add * (bj + radd[live])
            full[:, panel_positions] = 0.0
            Pm[j] = full[:, scored]
            # V3: production rebuild without donor j over ALL fold targets (build_baseline semantics incl. fallback)
            others = [i for i in sources if i != j]
            pre_v3 = np.zeros_like(pre)
            v3 = sum(A[i].astype(int) for i in others)
            for i in others:
                pre_v3[A[i]] += E[i][A[i]]
            cv = v3 > 0
            pre_v3[cv] /= v3[cv, None]
            if (~cv).any():
                mr = sum(prep[i]["mean_on_axis"] for i in others) / len(others)
                pre_v3[~cv] = E062.EFFECT_SCALE * mr / (float(np.linalg.norm(mr)) or 1.0)
            fin_j, cj, _, _ = compose(op, pre_v3, cd4, h1, wcols, panel_positions)
            Pm_full[j] = fin_j[live][:, scored]
            out.setdefault("v3_constants", {})[j] = cj

        def lodo(num_rows, d_rows, denom):
            key = np.full(n, np.nan)
            for p in np.flatnonzero(carried):
                ips = [float(num_rows[j][p] @ d_rows[j][p]) for j in sources if Al[j][p]]
                key[p] = np.mean(ips) / denom[p]
            return key

        d_raw = {s: raw[s][0][:, scored] for s in sources}
        d_prep = {s: El[s][:, scored] for s in sources}
        keys["LODO_V1"] = (lodo(Pm, d_raw, pp), carried)
        keys["LODO_V2"] = (lodo(Pm, d_prep, pp), carried)
        keys["LODO_V3"] = (lodo(Pm_full, d_prep, pp), carried)
        keys["D2_V2_baseonly_num"] = (lodo(Pm_base, d_prep, pp), carried)
        # D1: 410's formula on the prepared rows (base-only numerator AND denominator), + KA-4 algebra check
        Pb_live = (s_add * whitened)[live]
        Pb_live[:, panel_positions] = 0.0
        ppb = (Pb_live[:, scored] ** 2).sum(1)
        keys["D1_410formula_prep"] = (lodo(Pm_base, d_prep, ppb), carried)
        # KA-4: with 2 sources and every target carried by both, lodo(base) == pairwise formula on prepared rows
        if len(sources) == 2:
            a, b = sources
            both = Al[a] & Al[b]
            Wa = s_add * c * ALPHA * raw_whiten(op, El[a], wcols)
            Wb = s_add * c * ALPHA * raw_whiten(op, El[b], wcols)
            Wa[:, panel_positions] = 0.0
            Wb[:, panel_positions] = 0.0
            # off-axis columns of P^{-j} carry the unwhitened other row; include them exactly as P does
            offax = np.setdiff1d(np.arange(pre.shape[1]), wcols)
            Wa[:, offax] = s_add * ALPHA * El[a][:, offax]
            Wb[:, offax] = s_add * ALPHA * El[b][:, offax]
            Wa[:, panel_positions] = 0.0
            Wb[:, panel_positions] = 0.0
            Wa, Wb = Wa[:, scored], Wb[:, scored]
            pair = 0.5 * ((Wa * d_prep[b]).sum(1) + (Wb * d_prep[a]).sum(1)) / ppb
            kd1 = keys["D1_410formula_prep"][0]
            relka4 = float(np.nanmax(np.abs(pair[both] - kd1[both]) / np.abs(kd1[both]).max()))
            assert relka4 < 1e-9, f"KA-4 {held}: {relka4}"
            out["KA4_rel"] = relka4
            print(f"[{held}] KA-4 PASS lodo(base, 2 sources) == pairwise 410 formula, rel {relka4:.1e}", flush=True)

        # ---------------- endpoints
        names = list(keys)
        out["keys"] = {}
        for rule_lo, rule_hi in RULES:
            tag = f"lo{rule_lo}_hi{rule_hi:g}"
            W_ = {}
            for k in names:
                key, car = keys[k]
                w, med = weights(np.nan_to_num(key, nan=0.0), car, rule_lo, rule_hi)
                W_[k] = w
            gains = {k: agg(W_[k], pt, pp, tt) / base - 1 for k in names}
            ci, dci, pwin = boot(W_, pt, pp, tt, rng, names)
            # permutation control on V1 weights
            wv1 = W_["LODO_V1"]
            cidx = np.flatnonzero(carried)
            perm = []
            for _ in range(NPERM):
                wp = wv1.copy()
                wp[cidx] = rng.permutation(wv1[cidx])
                perm.append(agg(wp, pt, pp, tt) / base - 1)
            for k in names:
                d = out["keys"].setdefault(k, dict(spearman_sstar=spearman(keys[k][0][carried], sstar[carried]),
                                                   frac_key_le0=float((keys[k][0][carried] <= 0).mean())))
                d[tag] = dict(gain=gains[k], ci=ci[k])
                if k != "K410":
                    d[tag].update(diff_vs_410=gains[k] - gains["K410"], diff_ci=dci[k], p_boot_gt=pwin[k])
            out.setdefault("perm_V1", {})[tag] = dict(mean=float(np.mean(perm)),
                                                      q025=float(np.quantile(perm, .025)),
                                                      q975=float(np.quantile(perm, .975)))
        out["spearman_key_vs_410"] = {k: spearman(keys[k][0][carried], k410[carried]) for k in names if k != "K410"}
        out["seconds"] = time.time() - t0
        out["peak_rss_gib"] = E062.peak_rss_gib()
        res[held] = out
        print(f"[{held}] n {n} carried {int(carried.sum())} base agg {base:.5f} (orig {base0:.5f}) oracle {oracle:+.3f} "
              f"const {json.dumps({k: round(v, 4) for k, v in const.items()})} resid share "
              f"{out['residual_norm_share_median']:.3f}", flush=True)
        for k in names:
            d = out["keys"][k]
            line = f"   {k:22s} rho {d['spearman_sstar']:+.3f}"
            for rule_lo, rule_hi in RULES:
                tag = f"lo{rule_lo}_hi{rule_hi:g}"
                g = d[tag]
                line += f" | {tag} {g['gain']:+.3f} [{g['ci'][0]:+.3f},{g['ci'][1]:+.3f}]"
                if k != "K410":
                    line += f" d {g['diff_vs_410']:+.4f} [{g['diff_ci'][0]:+.4f},{g['diff_ci'][1]:+.4f}]"
            print(line, flush=True)
        print(f"   perm V1: {json.dumps(out['perm_V1'])}  ({out['seconds']:.0f}s, peak {out['peak_rss_gib']:.2f} GiB)",
              flush=True)
        (WORK / "lodo415.json").write_text(json.dumps(res, indent=1))
        del op, final, whitened, Pm, Pm_full, Pm_base
    print(f"done in {time.time() - t_all:.0f}s -> {WORK / 'lodo415.json'}")


if __name__ == "__main__":
    main()
