"""420 — bulk STACKS on PLOLO, one fold per process: composition x key x sharpening q, scored as DESIGNED (418's read-back)
and as REALISED by the ballast at KMAX 16 / 24 / 32 (414's emulation; the carriers keep PP unkeyed, as in 417a).

KNOWN ANSWERS FIRST (abort on failure, before any arm number is computed):
  KA-414   228's composition: `Kinf+sh@1.25` == arms414 + board414 and `Kinf+sh@1.5` == extra414 (cos_agg, pds,
           board-emulated pds lam 0.5 / 1, render A/B per frame), 1e-12.
  KA-PP    compose_arm(defaults) == sens415.compose_pp bit-for-bit (418's KA-415b).
  KA-V1    my LODO-V1 weights on PP (clip 0.3,3) == 418's saved `V1__PP` bit-for-bit.
  KA-418   `PP|K3`, `PP|K3+sh`, `b:h1x2|K3+sh`, `a:cd4x1.5|K3+sh` == comp418-<fold>.json (every field, |d| = 0).
  KA-REAL  scoring the DESIGNED rows through my row-scorer == lib418.score (1e-12): the realised numbers use one path.

GRID (every design): compositions {PP, h1x2, cd4x1.5, stack, stack_h1x1.5, stack_cd4x1.25} x keys {K3 (410 clip .3,3),
K5_3 (410 clip .5,3), V1 (415 LODO-V1 on PP, clip .3,3, held fixed), V1_5_3 (same key, clip .5,3)} x q {none, 1.25, 1.35,
1.5, 1.75, 2.0}. Secondary: V1own (LODO-V1 recomputed on the arm's composition) for the non-PP compositions at q 1.25/1.5.

    uv run python comp420.py <fold>            # k562_gwps | orion_hct116 | orion_hek293t | kolf2
    uv run python comp420.py <fold> --level    # per-context scaled variants of chosen designs (after analyse420)
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np

from lib420 import (ALPHA, COMP418, COMPS, E062, FOLDS, GAMMA, KEYS, KMAXES, LAM, PROD_OPS, QS, SHARP_Q, UNION, W414,
                    W418, WORK, Donors, Operator, betas_for, build_fold, cd4_for, cd4_frozen_and_partial, compose_arm,
                    compose_pp, de_profiles, deflate, design_and_realised, donor_rows, endpoints, frames, key_scale,
                    lodo_bases, lodo_key, pre_normrel, score, shared_energy_key, source_residual, union_rows, weights)
import lodo415

TOL = 1e-12
LAMS = (0.5, 1.0)


def score_rows(dls, scored, truth, betas, nonpanel):
    """lib418.score's read-back on given delivered rows {ctx: rows} (identical arithmetic; KA-REAL checks it)."""
    out, X, Aa, pdsp, pdsl, AB = {}, [], [], [], [], {}
    for ctx in "ABC":
        dl = dls[ctx]
        P = dl[:, scored]
        e = endpoints(P, truth)
        e["cos_agg"] = float(np.sqrt(max(0.0, 1.0 - e["sse"])))
        for lam in LAMS:
            e[f"pds_lam{lam:g}"] = endpoints((dl + lam * betas[ctx][None, :])[:, scored], truth)["pds"]
        out[ctx] = e
        X.append((P * truth).sum(1))
        Aa.append((P * P).sum(1))
        Tu = truth / np.maximum(np.linalg.norm(truth, axis=1, keepdims=True), 1e-300)
        for Pm, acc in ((P, pdsp), ((dl + 0.5 * betas[ctx][None, :])[:, scored], pdsl)):
            Pu = Pm / np.maximum(np.linalg.norm(Pm, axis=1, keepdims=True), 1e-300)
            S = Pu @ Tu.T
            own = np.diag(S).copy()
            n = S.shape[0]
            acc.append(((S < own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1)) / (n - 1))
        H = dl[:, nonpanel]
        a_ = float((H * H).sum(1).mean())
        Hc = H - H.mean(0, keepdims=True)
        b_ = float((Hc * Hc).sum(1).mean())
        AB[ctx] = dict(A=a_, B=b_, BoverA=b_ / a_)
    keys = ("pds", "cos", "sse", "cos_agg", "sgn_w", "pds_lam0.5", "pds_lam1")
    summ = {k: float(np.mean([out[c][k] for c in "ABC"])) for k in keys}
    summ["per_frame"] = {c: {k: out[c][k] for k in ("pds", "cos_agg", "pds_lam0.5")} for c in "ABC"}
    summ["AB"] = AB
    arrays = dict(x=np.asarray(X), a=np.asarray(Aa), pds=np.asarray(pdsp), pdsl=np.asarray(pdsl))
    return summ, arrays


def cmp(s, r, fields=("cos_agg", "pds", "pds_lam0.5", "pds_lam1")):
    d = max(abs(s[k] - r[k]) for k in fields)
    return max(d, *[abs(s["AB"][c][k] / r["AB"][c][k] - 1) for c in "ABC" for k in ("A", "B")])


def run(held: str, level_specs=None) -> None:
    t0 = time.time()
    log = lambda *a: print(f"[{held}]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    fold = build_fold(held)
    fold["scored"] = np.asarray(fold["scored"])
    scored = fold["scored"]
    pp_pos = fold["panel_positions"]
    wcols = np.asarray(fold["wcols"])
    truth_full = np.asarray(fold["truth"], dtype=np.float64)
    truth = truth_full[:, scored]
    pre_cached = np.asarray(fold["pre"], dtype=np.float64)
    targets = list(fold["targets"])
    sources = list(FOLDS[held]["sources"])
    op = Operator(PROD_OPS / "whitener_union_pruned.npz" if FOLDS[held]["op_kind"] == "production"
                  else FOLDS[held]["op"])
    nonpanel = np.ones(truth_full.shape[1], bool)
    nonpanel[pp_pos] = False
    betas = betas_for(held, nonpanel)
    res = dict(fold=held, n_targets=len(targets), sources=sources, ka={}, arms={}, comp={})
    arrays = {}

    # ============================================================== KA-414 (228's composition, Kinf, sharpening 1.25 / 1.5)
    cd4_old = cd4_for(targets)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    base228, _, _, _ = lodo415.compose(op, pre_cached, cd4_old, h1, wcols, pp_pos)
    key410, carried410, _ = shared_energy_key(targets, sources, op, wcols, cols=None)
    Kinf = key_scale(base228, weights(key410, carried410, 0.3, 1e9)[0])[0]
    a414 = json.loads((W414 / f"arms414-{held}.json").read_text())["arms"]
    b414 = json.loads((W414 / f"board414-{held}.json").read_text())["arms"]
    x414 = json.loads((W414 / f"extra414-{held}.json").read_text())["arms"]
    s125, _ = score(Kinf, fold, truth, betas, 1.25)
    r = dict(a414["Kinf+sh@1.25"])
    r["pds_lam0.5"], r["pds_lam1"] = b414["Kinf+sh@1.25"]["lam0.5"]["pds"], b414["Kinf+sh@1.25"]["lam1"]["pds"]
    d125 = cmp(s125, r)
    s150, _ = score(Kinf, fold, truth, betas, 1.5)
    d150 = cmp(s150, x414["Kinf+sh@1.5"])
    log(f"KA-414 Kinf+sh@1.25 worst |d| {d125:.1e}; Kinf+sh@1.5 (extra414) worst |d| {d150:.1e}")
    assert d125 < TOL and d150 < TOL, f"KA-414 FAILED {held}: {d125} {d150}"
    res["ka"]["ka414"] = dict(sh125=d125, sh150=d150)
    del base228, Kinf

    # ============================================================== PP composition (418's path) + KA-PP
    dn = Donors(targets, sources)
    frozen, partial = cd4_frozen_and_partial(targets)
    cd4, _ = deflate(frozen, 5)
    pre228 = dn.pre(sources, False)
    pre228[:, pp_pos] = 0.0
    nocov = sum(dn.A[s].astype(int) for s in sources) == 0
    pre228[nocov] = pre_cached[nocov]
    fb = pre228.copy()
    pre = pre_normrel(dn, sources, "cap", -1.5, fallback_rows=fb)
    assert np.array_equal(pre, dn.pre(sources, True, fallback_rows=fb)), "KA-PP pre"
    f415, _, c415, radd415, _ = compose_pp(op, pre, cd4, partial, h1, wcols, pp_pos)
    whit = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    finals, parts = {}, {}
    finals["PP"], parts["PP"] = compose_arm(whit, cd4, partial, h1, wcols, pp_pos)
    assert np.array_equal(finals["PP"], f415) and np.array_equal(parts["PP"]["radd"], radd415), "KA-PP compose"
    for name, (m, g) in COMPS.items():
        if name != "PP":
            finals[name], parts[name] = compose_arm(whit, cd4, partial, h1, wcols, pp_pos, m_cd4=m, g_h1=g)
        res["comp"][name] = dict(m_cd4=m, g_h1=g, fractions=parts[name]["fractions"], s=parts[name]["s"])
    log("KA-PP compose_arm(defaults) == compose_pp: PASS")

    # ============================================================== keys + KA-V1
    d_rows = {s: donor_rows(s, targets)[0] for s in sources}
    A_ = {s: dn.A[s] for s in sources}
    carried_all = sum(A_[s].astype(int) for s in sources) >= 2
    bases = lodo_bases(op, dn, sources, "cap", -1.5, pre, wcols)
    kv1 = lodo_key(finals["PP"], parts["PP"], bases, d_rows, A_, carried_all, sources, pp_pos)
    W = {"K3": weights(key410, carried410, 0.3, 3.0)[0], "K5_3": weights(key410, carried410, 0.5, 3.0)[0],
         "V1": weights(np.nan_to_num(kv1, nan=0.0), carried_all, 0.3, 3.0)[0],
         "V1_5_3": weights(np.nan_to_num(kv1, nan=0.0), carried_all, 0.5, 3.0)[0]}
    z418 = np.load(W418 / f"arrays418-{held}.npz")
    assert np.array_equal(W["V1"], z418["V1__PP"]) and np.array_equal(W["K3"], z418["W3"]), "KA-V1 / K3 weights"
    log("KA-V1: LODO-V1 (0.3,3) and K3 weights == 418's saved arrays: PASS")
    res["keys"] = {k: dict(q=np.quantile(w, [0, .1, .25, .5, .75, .9, 1]).tolist(), frac_floor=float(np.mean(w <= KEYS[k][1] + 1e-12)))
                   for k, w in W.items()}
    res["keys"]["spearman_K3_V1"] = float(np.corrcoef(np.argsort(np.argsort(W["K3"])), np.argsort(np.argsort(W["V1"])))[0, 1])

    P_DE = de_profiles(finals["PP"])

    ncheck = [0]

    def rec(name, comp, w, q, realised=True, scale=None, note="", check=False):
        """check=True also runs lib418.score on the design and asserts the row-scorer reproduces it (KA-REAL); it is on
        for every known-answer row and every 7th grid row."""
        bulk = key_scale(finals[comp], w)[0]
        ncheck[0] += 1
        if scale is None and (check or not realised or ncheck[0] % 7 == 0):
            s, a = score(bulk, fold, truth, betas, q)
        else:
            s = None
        entry = dict(comp=comp, q=q, note=note)
        if realised or scale is not None:
            dr = design_and_realised(bulk, P_DE, q=q, scale=scale, kmaxes=KMAXES)
            sd, ad = score_rows({c: dr[c]["design"] for c in "ABC"}, scored, truth, betas, nonpanel)
            if s is not None:
                dka = cmp(sd, s, ("cos_agg", "pds", "pds_lam0.5", "pds_lam1"))
                assert dka < TOL, f"KA-REAL FAILED {name}: {dka}"
                res["ka"]["real_worst"] = max(res["ka"].get("real_worst", 0.0), dka)
                res["ka"]["real_n"] = res["ka"].get("real_n", 0) + 1
            else:
                s, a = sd, ad
            entry["realised"] = {}
            for k in KMAXES:
                sr, ar = score_rows({c: dr[c]["realised"][k] for c in "ABC"}, scored, truth, betas, nonpanel)
                K = np.concatenate([dr[c]["stats"][k]["K"] for c in "ABC"])
                sr["frac_capped"] = float(np.mean(K >= k))
                sr["K_mean"] = float(np.mean(K))
                entry["realised"][str(k)] = sr
                arrays[f"{name}@R{k}"] = ar
        entry["design"] = s
        res["arms"][name] = entry
        arrays[name] = a
        rr = entry.get("realised", {}).get("16")
        log(f"{name:34s} cos {s['cos_agg']:.6f} pds@.5 {s['pds_lam0.5']:.6f} A {np.mean([s['AB'][c]['A'] for c in 'ABC']):.3f}"
            + (f" | R16 cos {rr['cos_agg']:.6f} pds@.5 {rr['pds_lam0.5']:.6f} capped {100 * rr['frac_capped']:.0f}%" if rr else ""))
        return entry

    if level_specs:
        for spec in level_specs:     # "comp|key+shq:A=..,B=..,C=.."
            nm, sc = spec.split(":")
            comp, rest = nm.split("|")
            key, q = (rest.split("+sh") + [None])[:2]
            scale = {kv.split("=")[0]: float(kv.split("=")[1]) for kv in sc.split(",")}
            rec(f"{nm}:{sc}", comp, W[key], float(q) if q else None, scale=scale, note="per-context bulk scale")
            if f"{nm}" not in res["arms"]:
                e = rec(nm, comp, W[key], float(q) if q else None)
                main = json.loads((WORK / f"comp420-{held}.json").read_text())["arms"][nm]
                d = max(cmp(e["realised"][k], main["realised"][k]) for k in e["realised"])
                assert d == 0.0, f"KA-LEVEL: unscaled {nm} != comp420 main grid: {d}"
                res["ka"]["level_unscaled_vs_main"] = max(res["ka"].get("level_unscaled_vs_main", 0.0), d)
        res["peak_rss_gib"] = E062.peak_rss_gib()
        (WORK / f"comp420L-{held}.json").write_text(json.dumps(res, indent=1, default=float))
        np.savez_compressed(WORK / f"arrays420L-{held}.npz", names=np.asarray(list(arrays)), d=np.asarray((truth * truth).sum(1)),
                            **{f"{n}__{k}": v for n, a in arrays.items() for k, v in a.items()})
        log(f"level done; peak RSS {res['peak_rss_gib']:.2f} GiB")
        return

    # ============================================================== KA-418 (before the grid)
    c418 = json.loads((W418 / f"comp418-{held}.json").read_text())["arms"]
    worst = 0.0
    for comp, q, ref in (("PP", None, "PP|K3"), ("PP", SHARP_Q, "PP|K3+sh"), ("h1x2", SHARP_Q, "b:h1x2|K3+sh"),
                         ("cd4x1.5", SHARP_Q, "a:cd4x1.5|K3+sh")):
        e = rec(f"{comp}|K3" + ("" if q is None else f"+sh{q:g}"), comp, W["K3"], q, check=True)
        d = cmp(e["design"], c418[ref])
        worst = max(worst, d)
        log(f"KA-418 {ref}: |d| {d:.1e}")
        assert d == 0.0, f"KA-418 FAILED {held} {ref}: {d}"
    res["ka"]["ka418_worst"] = worst
    log("ALL KNOWN ANSWERS PASS")

    # ============================================================== the grid
    for comp in COMPS:
        for key in KEYS:
            for q in (None, *QS):
                name = f"{comp}|{key}" + ("" if q is None else f"+sh{q:g}")
                if name not in res["arms"]:
                    rec(name, comp, W[key], q)
        (WORK / f"comp420-{held}.json").write_text(json.dumps(res, indent=1, default=float))
    # secondary: the arm's own LODO-V1 key, recomputed on its composition
    for comp in COMPS:
        if comp == "PP":
            continue
        kv = lodo_key(finals[comp], parts[comp], bases, d_rows, A_, carried_all, sources, pp_pos)
        wv = weights(np.nan_to_num(kv, nan=0.0), carried_all, 0.3, 3.0)[0]
        for q in (SHARP_Q, 1.5):
            rec(f"{comp}|V1own+sh{q:g}", comp, wv, q, note="own LODO-V1 recomputed on this composition")
    # control: the V1 key permuted over targets on the sharpened PP bulk (must lose)
    rng = np.random.default_rng(420)
    rec("CTRL:PP|V1shuf+sh1.25", "PP", W["V1"][rng.permutation(len(W["V1"]))], SHARP_Q, note="CONTROL: V1 permuted")
    res["peak_rss_gib"] = E062.peak_rss_gib()
    (WORK / f"comp420-{held}.json").write_text(json.dumps(res, indent=1, default=float))
    np.savez_compressed(WORK / f"arrays420-{held}.npz", names=np.asarray(list(arrays)), d=np.asarray((truth * truth).sum(1)),
                        **{f"{n}__{k}": v for n, a in arrays.items() for k, v in a.items()})
    log(f"done; {len(res['arms'])} arms; peak RSS {res['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    specs = [a for a in sys.argv[2:] if a != "--level"]
    run(sys.argv[1], level_specs=specs if "--level" in sys.argv[2:] else None)
