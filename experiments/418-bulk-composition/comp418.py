"""418 — one-variable changes of the BULK composition, priced on PLOLO with 414's instruments, one fold per process.

KNOWN ANSWERS FIRST (abort on failure, before any arm number is computed):
  KA-414  228's composition: `base`, `Kinf`, `Kinf+sh@1.25` reproduce arms414-<fold>.json (cos_agg, pds, render A/B)
          and board414-<fold>.json (board-emulated pds lam 0.5 / 1) to 1e-12.
  KA-415  the PP composition: my pre (normrel cap -1.5) == Donors.pre bit-for-bit; my compose_arm at the defaults ==
          sens415.compose_pp bit-for-bit (final row, max |d| = 0); the composition constants == sens415.json
          constants_pp; effect-level base_agg, K410 and V1_PP gains at clip(0.3,3)/(0.3,inf) == sens415.json (1e-12).

ARMS (all on PP; key = the fold's 410-formula key at clip(0.3,3) held fixed, i.e. 410b's bulk; each also with 417a's
sharpening q 1.25; secondary: the arm's own LODO-V1 key recomputed on its composition):
  (a) CD4 target f x {0.5, 0.75, 1.25, 1.5, 2.0}   (b) H1-union g x {0, 0.5, 1.5, 2.0}   (c) CD4 deflation k {0, 3, 8}
  (d) CD4 partial g {0, 0.25, 0.75, 1.0}           (e) normrel cap exponent {0 (off), -0.75, -3}
  (f1) the key scales the donor base only; the CD4/H1/partial residual is added unkeyed
  (f2) normrel PLAIN mode at -1.5 (the cap was chosen on the welded board by its DE members)
RETRO (no key, no sharpening): the historical one-variable composition steps that have board pds/mse readings.

    uv run python comp418.py <fold>          # k562_gwps | orion_hct116 | orion_hek293t | kolf2
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np

from lib418 import (ALPHA, E062, FOLDS, GAMMA, LAM, PROD_OPS, SHARP_Q, UNION, W414, W415, WORK, Donors, Operator,
                    agg, agg_parts, betas_for, build_fold, cd4_for, cd4_frozen_and_partial, compose_arm, compose_pp,
                    deflate, donor_rows, key_on_base_only, key_scale, lodo_bases, lodo_key, pre_normrel, score,
                    shared_energy_key, source_residual, union_rows, weights)
import lodo415

NR = {  # normrel settings: (mode, exponent)
    "cap-1.5": ("cap", -1.5), "off": (None, None), "cap-0.75": ("cap", -0.75), "cap-3": ("cap", -3.0),
    "plain-1.5": ("plain", -1.5), "plain-3": ("plain", -3.0),
}
# name -> (normrel, deflate k, m_cd4, g_h1, g_part)
ARMS = {
    "PP": ("cap-1.5", 5, 1.0, 0.28, 0.5),
    **{f"a:cd4x{m:g}": ("cap-1.5", 5, m, 0.28, 0.5) for m in (0.5, 0.75, 1.25, 1.5, 2.0)},
    **{f"b:h1x{m:g}": ("cap-1.5", 5, 1.0, 0.28 * m, 0.5) for m in (0.0, 0.5, 1.5, 2.0)},
    **{f"c:k{k}": ("cap-1.5", k, 1.0, 0.28, 0.5) for k in (0, 3, 8)},
    **{f"d:g{g:g}": ("cap-1.5", 5, 1.0, 0.28, g) for g in (0.0, 0.25, 0.75, 1.0)},
    "e:exp0": ("off", 5, 1.0, 0.28, 0.5),
    "e:exp-0.75": ("cap-0.75", 5, 1.0, 0.28, 0.5),
    "e:exp-3": ("cap-3", 5, 1.0, 0.28, 0.5),
    "f2:plain-1.5": ("plain-1.5", 5, 1.0, 0.28, 0.5),
}
RETRO = {  # historical board steps (pre-key, pre-ballast; bulk members are bulk functionals either way)
    "R:rc_callset(f0.10)": ("off", 0, 0.10 / 0.34, 0.0, 0.0),
    "R:cd4f034": ("off", 0, 1.0, 0.0, 0.0),
    "R:h1union": ("off", 0, 1.0, 0.28, 0.0),
    "R:cd4svd5": ("off", 5, 1.0, 0.28, 0.0),
    "R:312(plain-1.5)": ("plain-1.5", 5, 1.0, 0.28, 0.0),
    "R:312b(plain-3)": ("plain-3", 5, 1.0, 0.28, 0.0),
    "R:319(cap-1.5)": ("cap-1.5", 5, 1.0, 0.28, 0.0),
    "R:319d3(cap-3)": ("cap-3", 5, 1.0, 0.28, 0.0),
    "R:254(svd5+p0.34)": ("off", 5, 1.0, 0.28, 0.34),
    "R:319p(p0.34)": ("cap-1.5", 5, 1.0, 0.28, 0.34),
    "R:319p2(p0.5)": ("cap-1.5", 5, 1.0, 0.28, 0.5),
}
TOL = 1e-12
# dose continuation (run with --extra after the main grid): the H1/CD4 doses were still rising at the grid's edge.
EXTRA = {
    "PP": ARMS["PP"],
    **{f"b:h1x{m:g}": ("cap-1.5", 5, 1.0, 0.28 * m, 0.5) for m in (2.5, 3.0, 4.0)},
    "a:cd4x1.75": ("cap-1.5", 5, 1.75, 0.28, 0.5),
}


def run(held: str, extra: bool = False) -> None:
    global ARMS, RETRO
    if extra:
        ARMS, RETRO = EXTRA, {}
    t0 = time.time()
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
    log = lambda *a: print(f"[{held}]", *a, f"({time.time() - t0:.0f}s)", flush=True)

    def rec(name, bulk, sharp_q=None, note=""):
        s, a = score(bulk, fold, truth, betas, sharp_q)
        s["note"] = note
        res["arms"][name] = s
        arrays[name] = a
        log(f"{name:34s} cos_agg {s['cos_agg']:.6f} pds {s['pds']:.6f} pds@0.5 {s['pds_lam0.5']:.6f} "
            f"A {np.mean([s['AB'][c]['A'] for c in 'ABC']):.4f}")
        return s

    # ============================================================== KA-414 (228's composition)
    cd4_old = cd4_for(targets)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    base228, _, _, _ = lodo415.compose(op, pre_cached, cd4_old, h1, wcols, pp_pos)
    key, carried, nd = shared_energy_key(targets, sources, op, wcols, cols=None)
    W3 = weights(key, carried, 0.3, 3.0)[0]
    Wi = weights(key, carried, 0.3, 1e9)[0]
    Kinf = key_scale(base228, Wi)[0]
    a414 = json.loads((W414 / f"arms414-{held}.json").read_text())["arms"]
    b414 = json.loads((W414 / f"board414-{held}.json").read_text())["arms"]
    worst = 0.0
    for name, bulk, q in (("base", base228, None), ("Kinf", Kinf, None), ("Kinf+sh@1.25", Kinf, 1.25)):
        s = rec(f"KA414:{name}", bulk, q, "228 composition, KA-414")
        r, rb = a414[name], b414[name]
        d = max(abs(s["cos_agg"] - r["cos_agg"]), abs(s["pds"] - r["pds"]),
                abs(s["pds_lam0.5"] - rb["lam0.5"]["pds"]), abs(s["pds_lam1"] - rb["lam1"]["pds"]),
                *[abs(s["AB"][c][k] / r["AB"][c][k] - 1) for c in "ABC" for k in ("A", "B")])
        worst = max(worst, d)
        log(f"KA-414 {name}: worst |d| {d:.1e}")
        assert d < TOL, f"KA-414 FAILED {held} {name}: {d}"
    res["ka"]["ka414_worst"] = worst
    res["key"] = dict(carried=int(carried.sum()), w3_q=np.quantile(W3, [0, .1, .25, .5, .75, .9, 1]).tolist())

    # ============================================================== KA-415 (the PP composition)
    dn = Donors(targets, sources)
    frozen, partial = cd4_frozen_and_partial(targets)
    assert np.array_equal(frozen.mean_effect, cd4_old.mean_effect) and np.array_equal(frozen.gene_votes, cd4_old.gene_votes)
    cd4 = {0: frozen}
    cd4[5], daud = deflate(frozen, 5)
    pre228 = dn.pre(sources, False)
    pre228[:, pp_pos] = 0.0
    nocov = sum(dn.A[s].astype(int) for s in sources) == 0
    pre228[nocov] = pre_cached[nocov]
    assert np.abs(pre228 - pre_cached).max() / np.abs(pre_cached).max() < 1e-5, "228 pre rebuild"
    fb = pre228.copy()
    pre415 = dn.pre(sources, True, fallback_rows=fb)
    pre = {"cap-1.5": pre_normrel(dn, sources, "cap", -1.5, fallback_rows=fb)}
    assert np.array_equal(pre415, pre["cap-1.5"]), "KA-415 FAILED: pre_normrel(cap,-1.5) != Donors.pre"
    f415, _, c415, radd415, cc415 = compose_pp(op, pre415, cd4[5], partial, h1, wcols, pp_pos)
    whit = {"cap-1.5": op.apply(pre["cap-1.5"] * ALPHA, wcols, GAMMA, LAM)}
    fpp, ppp = compose_arm(whit["cap-1.5"], cd4[5], partial, h1, wcols, pp_pos)
    dfin = float(np.abs(fpp - f415).max())
    log(f"KA-415 compose_arm(defaults) vs sens415.compose_pp: max |d| {dfin:.1e}; radd |d| "
        f"{np.abs(ppp['radd'] - radd415).max():.1e}; s {ppp['s']:.12f} vs {c415['s']:.12f}")
    assert dfin == 0.0 and np.array_equal(ppp["radd"], radd415) and ppp["s"] == c415["s"], "KA-415 FAILED: compose"
    s415 = json.loads((W415 / "sens415.json").read_text())[held]
    dc = max(abs(c415[k] - v) for k, v in s415["constants_pp"].items())
    log(f"KA-415 constants_pp vs sens415.json: max |d| {dc:.1e}")
    assert dc < TOL, f"KA-415 FAILED constants {dc}"
    # effect-level reads on the live scored rows (sens415.plolo_fold's convention)
    P, T = fpp[:, scored], truth
    live = (np.abs(T).sum(1) > 0) & (np.abs(P).sum(1) > 0)
    pt, pp_, tt = agg_parts(P[live], T[live])
    base_agg = agg(np.ones(int(live.sum())), pt, pp_, tt)
    tg = [t for t, k in zip(targets, live) if k]
    k410s, car_s, _ = shared_energy_key(tg, sources, op, wcols, cols=scored)
    d_rows = {s: donor_rows(s, targets)[0] for s in sources}
    A_ = {s: dn.A[s] for s in sources}
    carried_all = sum(A_[s].astype(int) for s in sources) >= 2
    bases = {"cap-1.5": lodo_bases(op, dn, sources, "cap", -1.5, pre["cap-1.5"], wcols)}
    kv1s = lodo_key(fpp, ppp, bases["cap-1.5"], d_rows, A_, carried_all, sources, pp_pos, cols=scored)[live]
    ka_rows = [("base_agg", base_agg, s415["base_agg"])]
    for lo, hi in ((0.3, 3.0), (0.3, 1e9)):
        tag = f"lo{lo}_hi{hi:g}"
        for kname, kk in (("K410", k410s), ("V1_PP", kv1s)):
            g = agg(weights(np.nan_to_num(kk, nan=0.0), car_s, lo, hi)[0], pt, pp_, tt) / base_agg - 1
            ka_rows.append((f"{kname} {tag}", g, s415["keys"][kname][tag]["gain"]))
    for lab, got, exp in ka_rows:
        log(f"KA-415 {lab:22s} {got:+.9f} (sens415 {exp:+.9f}, |d| {abs(got - exp):.1e})")
        assert abs(got - exp) < TOL, f"KA-415 FAILED {lab}"
    res["ka"]["ka415"] = {lab: [got, exp] for lab, got, exp in ka_rows}
    res["ka"]["deflation"] = daud
    log("ALL KNOWN ANSWERS PASS")

    # ============================================================== bases for every normrel setting
    need_nr = sorted({v[0] for v in ARMS.values()} | {v[0] for v in RETRO.values()})
    for nr in need_nr:
        if nr not in pre:
            mode, exp = NR[nr]
            pre[nr] = pre_normrel(dn, sources, mode, exp, fallback_rows=fb)
            whit[nr] = op.apply(pre[nr] * ALPHA, wcols, GAMMA, LAM)
    for k in {v[1] for v in ARMS.values()} | {v[1] for v in RETRO.values()}:
        if k not in cd4:
            cd4[k], _ = deflate(frozen, k)

    def comp(cfg):
        nr, k, m, gh, gp = cfg
        return compose_arm(whit[nr], cd4[k], partial, h1, wcols, pp_pos, m_cd4=m, g_h1=gh, g_part=gp)

    # reference rows: PP unkeyed (= 395's bulk), PP + K3 (= 410b), PP + K3 + sh (= 417a)
    rec("PP|none", fpp, None, "PP unkeyed = 395/319p2 bulk")
    V1 = {}
    for name, cfg in ARMS.items():
        final, parts = (fpp, ppp) if name == "PP" else comp(cfg)
        if name != "PP":
            assert not np.array_equal(final, fpp), f"{name} is not a change"
        res["comp"][name] = dict(cfg=list(cfg), fractions=parts["fractions"], s=parts["s"],
                                 row_norm_q=np.quantile(np.linalg.norm(final, axis=1), [.1, .5, .9]).tolist(),
                                 frob=float(np.linalg.norm(final)),
                                 cos_with_PP_median=float(np.median(np.einsum("ij,ij->i", final, fpp) / np.maximum(
                                     np.linalg.norm(final, axis=1) * np.linalg.norm(fpp, axis=1), 1e-300))))
        bulk = key_scale(final, W3)[0]
        rec(f"{name}|K3", bulk, None, "fixed 410-formula key clip(0.3,3)")
        rec(f"{name}|K3+sh", bulk, SHARP_Q, "fixed key + 417a sharpening")
        nr = cfg[0]
        if nr not in bases:
            mode, exp = NR[nr]
            bases[nr] = lodo_bases(op, dn, sources, mode, exp, pre[nr], wcols)
        kv = lodo_key(final, parts, bases[nr], d_rows, A_, carried_all, sources, pp_pos)
        wv = weights(np.nan_to_num(kv, nan=0.0), carried_all, 0.3, 3.0)[0]
        V1[name] = wv
        bv = key_scale(final, wv)[0]
        rec(f"{name}|V1", bv, None, "arm's own LODO-V1 key clip(0.3,3)")
        rec(f"{name}|V1+sh", bv, SHARP_Q, "own LODO-V1 key + 417a sharpening")
        if not extra:
            (WORK / f"comp418-{held}.json").write_text(json.dumps(res, indent=1))
    # f1: key on the donor base only (PP composition, same w)
    b1 = key_on_base_only(whit["cap-1.5"], ppp, W3, pp_pos, fpp)
    res["comp"]["f1:key-base-only"] = dict(cfg=["cap-1.5", 5, 1.0, 0.28, 0.5], note="key on W only",
                                           cos_with_PP_K3_median=float(np.median(np.einsum(
                                               "ij,ij->i", b1, key_scale(fpp, W3)[0]) / np.maximum(np.linalg.norm(
                                                   b1, axis=1) * np.linalg.norm(key_scale(fpp, W3)[0], axis=1), 1e-300))))
    rec("f1:key-base-only|K3", b1, None, "key scales W only, residual unkeyed")
    rec("f1:key-base-only|K3+sh", b1, SHARP_Q, "key scales W only + sharpening")
    # controls: the PP key permuted over targets (must lose), and a V1 sanity (PP|V1 vs PP|K3 = 415's claim)
    rng = np.random.default_rng(418)
    rec("CTRL:PP|K3shuf", key_scale(fpp, W3[rng.permutation(len(W3))])[0], None, "CONTROL: key permuted")
    # retro rows (no key)
    for name, cfg in RETRO.items():
        final, parts = comp(cfg)
        res["comp"][name] = dict(cfg=list(cfg), fractions=parts["fractions"], s=parts["s"])
        if name == "R:319p2(p0.5)":
            assert np.array_equal(final, fpp), "retro 319p2 must equal PP"
        rec(name, final, None, "historical board step, no key")
    res["peak_rss_gib"] = E062.peak_rss_gib()
    tag = "comp418x" if extra else "comp418"
    if extra:   # known answer: the reference rows must reproduce the main grid's exactly
        main = json.loads((WORK / f"comp418-{held}.json").read_text())["arms"]
        for nm in ("PP|none", "PP|K3", "PP|K3+sh", "PP|V1", "PP|V1+sh", "CTRL:PP|K3shuf"):
            d = max(abs(res["arms"][nm][k] - main[nm][k]) for k in ("cos_agg", "pds", "pds_lam0.5"))
            log(f"KA-main {nm}: |d| {d:.1e}")
            assert d == 0.0, f"KA-main FAILED {nm}"
    (WORK / f"{tag}-{held}.json").write_text(json.dumps(res, indent=1))
    np.savez_compressed(WORK / f"{'arrays418x' if extra else 'arrays418'}-{held}.npz", names=np.asarray(list(arrays)),
                        d=np.asarray((truth * truth).sum(1)), W3=W3,
                        **{f"V1__{k}": v for k, v in V1.items()},
                        **{f"{n}__{k}": v for n, a in arrays.items() for k, v in a.items()})
    log(f"done; peak RSS {res['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    run(sys.argv[1], extra="--extra" in sys.argv[2:])
