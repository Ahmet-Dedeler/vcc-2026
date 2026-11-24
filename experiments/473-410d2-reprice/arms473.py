"""473b — 228's per-gene gain m on 446b's bulk composition, PLOLO, one fold per process (420's machinery, unchanged).

Base = `stack|K3+sh1.75` (lib420: stack H1 x2 + CD4 x1.5, the 410 key clip (0.3, 3), sharpening q 1.75 row-renorm on
each challenge frame's gate), scored as DESIGNED and as REALISED by 414's ballast emulation at KMAX 16 (the carriers keep
PP unkeyed, as in 446b). m = the fold's 228 gain (`gain-<fold>.npz`, gain_reg, t 0.25: m = 1 + t (g/median - 1)); on
kolf2 that m IS the builder's pinned artifact (414 price410d). Arms (one variable each vs the base):
  A1      m clipped [1/2, 2] x the keyed bulk BEFORE sharpening, one constant restoring the energy on the gated
          non-panel genes summed over the three frames (the builder's `...rg` rule)
  A2      unclipped m, same placement      A3  m clipped [1/3, 3]
  A1post  clipped m AFTER sharpening (per frame), each row's gated norm restored
  C1      A1 with m permuted over genes (CONTROL: must lose)

KNOWN ANSWERS FIRST (abort on failure, before any arm number):
  KA-414g  my gain code on 228's composition: `Kinf+gf@0.25` == arms414-<fold>.json (cos_agg, pds), 1e-12
  KA-PP    compose_arm(defaults) == sens415.compose_pp bit-for-bit; K3 weights == 418's saved W3
  KA-420   `stack|K3+sh1.75` design and realised (KMAX 16) == comp420-<fold>.json, |d| < 1e-12

    uv run python experiments/473-410d2-reprice/arms473.py <fold>     # k562_gwps | orion_hct116 | orion_hek293t | kolf2
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/420-bulk-stack"))
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import lib420  # noqa: E402
from lib420 import (COMPS, E062, FOLDS, PROD_OPS, TS, Donors, Operator, betas_for, build_fold, cd4_for,  # noqa: E402
                    cd4_frozen_and_partial, compose_arm, compose_pp, de_profiles, deflate, design_and_realised, frames,
                    key_scale, pre_normrel, profile_from_effect, realise_k, score, shared_energy_key, sharpen,
                    source_residual, union_rows, weights)
import lodo415  # noqa: E402
from arms228 import multiplier  # noqa: E402
from comp420 import cmp, score_rows  # noqa: E402

W420 = ROOT / "experiments/420-bulk-stack/work"
W414 = ROOT / "experiments/414-bulk-reprice/work"
W418 = ROOT / "experiments/418-bulk-composition/work"
WORK = HERE / "work"
Q = 1.75
TOL = 1e-12


def frob_preserve(E, M):
    return M * (np.linalg.norm(E) / max(np.linalg.norm(M), 1e-300))


def run(held: str) -> None:
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
    res = dict(fold=held, n_targets=len(targets), sources=sources, ka={}, arms={})
    arrays = {}

    # ---------------------------------------------------------------- the fold's 228 gain
    z = np.load(ROOT / f"experiments/228-theta-gain-plolo/work/gain-{held}.npz", allow_pickle=True)
    live = np.asarray(z["live"], dtype=bool)
    greg = np.asarray(z["gain_reg"], dtype=np.float64)
    m = multiplier(greg, live, 0.25)
    res["m_q"] = np.quantile(m, [0, .01, .1, .5, .9, .99, 1]).tolist()
    res["m_frac_clipped_2"] = float(np.mean((m < 0.5) | (m > 2.0)))

    # ---------------------------------------------------------------- KA-414g (228 composition, Kinf, gain on the final row)
    cd4_old = cd4_for(targets)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in lodo415.UNION])
    base228, _, _, _ = lodo415.compose(op, pre_cached, cd4_old, h1, wcols, pp_pos)
    key410, carried410, _ = shared_energy_key(targets, sources, op, wcols, cols=None)
    Kinf = key_scale(base228, weights(key410, carried410, 0.3, 1e9)[0])[0]
    s_g, _ = score(frob_preserve(Kinf, Kinf * m[None, :]), fold, truth, betas)
    r414 = json.loads((W414 / f"arms414-{held}.json").read_text())["arms"]["Kinf+gf@0.25"]
    d = max(abs(s_g["cos_agg"] - r414["cos_agg"]), abs(s_g["pds"] - r414["pds"]))
    log(f"KA-414g Kinf+gf@0.25 |d| {d:.1e} (cos_agg {s_g['cos_agg']:.6f}, pds {s_g['pds']:.6f})")
    assert d < TOL, f"KA-414g FAILED {held}: {d}"
    res["ka"]["ka414g"] = d
    del base228, Kinf, cd4_old

    # ---------------------------------------------------------------- PP + stack composition (420's path) + KA-PP
    dn = Donors(targets, sources)
    frozen, partial = cd4_frozen_and_partial(targets)
    cd4, _ = deflate(frozen, 5)
    pre228 = dn.pre(sources, False)
    pre228[:, pp_pos] = 0.0
    nocov = sum(dn.A[s].astype(int) for s in sources) == 0
    pre228[nocov] = pre_cached[nocov]
    fb = pre228.copy()
    pre = pre_normrel(dn, sources, "cap", -1.5, fallback_rows=fb)
    f415, _, _, radd415, _ = compose_pp(op, pre, cd4, partial, h1, wcols, pp_pos)
    whit = op.apply(pre * lib420.ALPHA, wcols, lib420.GAMMA, lib420.LAM)
    fPP, pPP = compose_arm(whit, cd4, partial, h1, wcols, pp_pos)
    assert np.array_equal(fPP, f415) and np.array_equal(pPP["radd"], radd415), "KA-PP compose"
    mc, gh = COMPS["stack"]
    fST, pST = compose_arm(whit, cd4, partial, h1, wcols, pp_pos, m_cd4=mc, g_h1=gh)
    W3 = weights(key410, carried410, 0.3, 3.0)[0]
    z418 = np.load(W418 / f"arrays418-{held}.npz")
    assert np.array_equal(W3, z418["W3"]), "KA K3 weights"
    log("KA-PP and K3 weights PASS")
    del whit, f415, pre, pre228, fb, dn
    P_DE = de_profiles(fPP)
    bulk = key_scale(fST, W3)[0]

    gates = [g & nonpanel for (_, g) in frames().values()]

    def rg_const(M):
        before = sum(float((bulk[:, g] ** 2).sum()) for g in gates)
        after = sum(float(((bulk * M[None, :])[:, g] ** 2).sum()) for g in gates)
        return float(np.sqrt(before / after))

    def rec(name, B=None, post=None, note=""):
        if post is None:
            dr = design_and_realised(B, P_DE, q=Q, kmaxes=(16,))
        else:                       # A1post: sharpen per frame, multiply, restore each row's gated norm, then render
            dr = {}
            for ctx, (b, gate) in frames().items():
                F = sharpen(bulk, Q, gate, mode="row")
                n0 = np.linalg.norm(F, axis=1)
                F = F * post[None, :]
                F = F * (n0 / np.maximum(np.linalg.norm(F, axis=1), 1e-300))[:, None]
                des, rea = np.empty_like(F), np.empty_like(F)
                K = np.empty(F.shape[0], int)
                for i in range(F.shape[0]):
                    p_bk = profile_from_effect(F[i], b, gate, "minnorm")
                    des[i] = np.log1p(TS * p_bk) - b
                    p_re, _, K[i], _, _ = realise_k(P_DE[ctx][i], p_bk, 16)
                    rea[i] = np.log1p(TS * p_re) - b
                dr[ctx] = dict(design=des, realised={16: rea}, stats={16: dict(K=K)})
        sd, ad = score_rows({c: dr[c]["design"] for c in "ABC"}, scored, truth, betas, nonpanel)
        sr, ar = score_rows({c: dr[c]["realised"][16] for c in "ABC"}, scored, truth, betas, nonpanel)
        K = np.concatenate([dr[c]["stats"][16]["K"] for c in "ABC"])
        sr["frac_capped"] = float(np.mean(K >= 16))
        res["arms"][name] = dict(note=note, design=sd, realised={"16": sr})
        arrays[name] = ad
        arrays[f"{name}@R16"] = ar
        log(f"{name:34s} design cos {sd['cos_agg']:.6f} pds@.5 {sd['pds_lam0.5']:.6f} | R16 cos {sr['cos_agg']:.6f} "
            f"pds {sr['pds']:.6f} pds@.5 {sr['pds_lam0.5']:.6f} A {np.mean([sr['AB'][c]['A'] for c in 'ABC']):.4f} "
            f"capped {100 * sr['frac_capped']:.0f}%")
        return res["arms"][name]

    # ---------------------------------------------------------------- KA-420 (the base, before any arm)
    e = rec("stack|K3+sh1.75", bulk, note="446b's bulk composition; KA-420")
    c420 = json.loads((W420 / f"comp420-{held}.json").read_text())["arms"]["stack|K3+sh1.75"]
    dD = cmp(e["design"], c420["design"])
    dR = cmp(e["realised"]["16"], c420["realised"]["16"])
    log(f"KA-420 stack|K3+sh1.75 design |d| {dD:.1e}, realised(16) |d| {dR:.1e}")
    # comp420 stored some designed rows via lib418.score and others via its row-scorer (its KA-REAL: equal to 1e-12)
    assert dD < TOL and dR < TOL, f"KA-420 FAILED {held}: {dD} {dR}"
    res["ka"]["ka420"] = dict(design=dD, realised=dR)
    log("ALL KNOWN ANSWERS PASS")

    # ---------------------------------------------------------------- arms
    rng = np.random.default_rng(473)
    mc2 = np.clip(m, 0.5, 2.0)
    for name, M, note in (("A1 clip2 pre rg", mc2, "m clip [1/2,2], pre-sharpen, render-gated constant"),
                          ("A2 unclipped pre rg", m, "m unclipped, pre-sharpen, render-gated"),
                          ("A3 clip3 pre rg", np.clip(m, 1 / 3, 3.0), "m clip [1/3,3], pre-sharpen, render-gated"),
                          ("C1 clip2 permuted", mc2[rng.permutation(mc2.size)], "CONTROL: A1's m permuted over genes")):
        c = rg_const(M)
        res.setdefault("constants", {})[name] = c
        rec(name, bulk * M[None, :] * c, note=note)
    rec("A1post clip2 post-sharpen", post=mc2, note="m clip [1/2,2] after sharpening, row gated norm restored")

    res["peak_rss_gib"] = E062.peak_rss_gib()
    (WORK / f"arms473-{held}.json").write_text(json.dumps(res, indent=1, default=float))
    np.savez_compressed(WORK / f"arrays473-{held}.npz", names=np.asarray(list(arrays)), d=np.asarray((truth * truth).sum(1)),
                        **{f"{n}__{k}": v for n, a in arrays.items() for k, v in a.items()})
    log(f"done; peak RSS {res['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    run(sys.argv[1])
