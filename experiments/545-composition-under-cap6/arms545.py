"""545 (PLOLO) — the stack's composition re-scanned under the cap-6 LODO key (495): CD4 multiplier m_cd4 (x CD4_F 0.34;
production 1.5, i.e. VCC_CD4_TARGET_F 0.51) x H1 union g_h1 (production 0.56 = 2 x 0.28, VCC_H1_UNION_G), each keyed by
the fold's LODO key at (0.5, 6) (and the production mix at (0.5, 3) for the known answer), q 1.75, realised 16, against
PP's unkeyed carriers (the DE channel is untouched). KNOWN ANSWERS: the production mix at (0.5, 3) == comp420's
`stack|V1_5_3+sh1.75`, and at (0.5, 6) == 540b's `V1_0.5_6`, |d| < 1e-12.

(540b's docstring follows.)
540b (PLOLO) — the cap extension of 540 (caps 4 ... uncapped; same code, CLIPS changed).

540 (PLOLO) — the LODO key's clip (floor, cap) on 446bj2's stack at q 1.75, one fold per process (lib420, 420's path).

The fold's LODO-V1 key vector kv1 is computed as comp420 does (lodo_bases + lodo_key on PP; KA: weights(kv1, 0.3, 3) ==
418's saved `V1__PP`), then each clip w = weights(kv1, carried>=2, lo, hi) keys the fold's unkeyed stack composition (fST),
which is sharpened at q 1.75 on each frame's gate (row mode), rendered (minnorm) and realised at KMAX 16 against PP's
unkeyed carriers. KNOWN ANSWER (abort unless |d| < 1e-12): clip (0.5, 3) == comp420's `stack|V1_5_3+sh1.75` (designed and
realised 16).

    uv run python experiments/545-composition-under-cap6/arms545.py <fold>        (~5.5 GB peak, ~20 s)
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
from lib420 import (COMPS, E062, FOLDS, PROD_OPS, TS, Donors, Operator, betas_for, build_fold,  # noqa: E402
                    cd4_frozen_and_partial, compose_arm, compose_pp, de_profiles, deflate, frames, key_scale,
                    pre_normrel, profile_from_effect, realise_k, sharpen, source_residual, union_rows, weights)
import lodo415  # noqa: E402
from comp420 import cmp, score_rows  # noqa: E402

W420 = ROOT / "experiments/420-bulk-stack/work"
W418 = ROOT / "experiments/418-bulk-composition/work"
WORK = HERE / "work"
TOL = 1e-12
CLIPS = [(0.5, 3.0), (0.5, 4.0), (0.5, 6.0), (0.5, 8.0), (0.5, 10.0), (0.5, 15.0), (0.5, 1e9), (0.3, 6.0), (0.4, 6.0),
         (0.3, 1e9)]


def tag(lo, hi):
    return f"V1_{lo:g}_{hi:g}"


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
    res = dict(fold=held, n_targets=len(targets), ka={}, arms={}, keys={})
    arrays = {}
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in lodo415.UNION])
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
    fST, _ = compose_arm(whit, cd4, partial, h1, wcols, pp_pos, m_cd4=mc, g_h1=gh)
    import lib418 as _l418
    d_rows = {s_: _l418.donor_rows(s_, targets)[0] for s_ in sources}
    A_ = {s_: dn.A[s_] for s_ in sources}
    carried_all = sum(A_[s_].astype(int) for s_ in sources) >= 2
    bases = _l418.lodo_bases(op, dn, sources, "cap", -1.5, pre, wcols)
    kv1 = np.nan_to_num(_l418.lodo_key(fPP, pPP, bases, d_rows, A_, carried_all, sources, pp_pos), nan=0.0)
    assert np.array_equal(weights(kv1, carried_all, 0.3, 3.0)[0], np.load(W418 / f"arrays418-{held}.npz")["V1__PP"]), "KA V1"
    del bases, d_rows, f415, pre, pre228, fb, dn
    P_DE = de_profiles(fPP)
    log("KA-PP and KA V1 (== 418's V1__PP) PASS")

    def rows_for(bulk):
        dr = {}
        for ctx, (b, gate) in frames().items():
            F = sharpen(bulk, 1.75, gate, mode="row")
            des, rea = np.empty_like(F), np.empty_like(F)
            K = np.empty(F.shape[0], int)
            for i in range(F.shape[0]):
                p_bk = profile_from_effect(F[i], b, gate, "minnorm")
                des[i] = np.log1p(TS * p_bk) - b
                p_re, _, K[i], _, _ = realise_k(P_DE[ctx][i], p_bk, 16)
                rea[i] = np.log1p(TS * p_re) - b
            dr[ctx] = dict(design=des, realised=rea, K=K)
        return dr

    c420 = json.loads((W420 / f"comp420-{held}.json").read_text())["arms"]
    b540 = json.loads((ROOT / f"experiments/540-lodo-clip-scan/work/arms540b-{held}.json").read_text())["arms"]
    W3 = weights(kv1, carried_all, 0.5, 3.0)[0]
    W6 = weights(kv1, carried_all, 0.5, 6.0)[0]
    grid = [(1.5, 0.56, 3.0)] + [(mc, gh, 6.0) for mc in (1.0, 1.25, 1.5, 1.75, 2.0, 2.5) for gh in (0.42, 0.56, 0.70, 0.84)]
    for mc, gh, hi in grid:
        fC, partsC = compose_arm(whit, cd4, partial, h1, wcols, pp_pos, m_cd4=mc, g_h1=gh)
        dr = rows_for(key_scale(fC, W3 if hi == 3.0 else W6)[0])
        sr, ar = score_rows({c: dr[c]["realised"] for c in "ABC"}, scored, truth, betas, nonpanel)
        K = np.concatenate([dr[c]["K"] for c in "ABC"])
        sr["K_mean_emul"] = float(K.mean())
        name = f"cd4x{mc:g}_h1g{gh:g}|cap{hi:g}"
        res["arms"][name] = dict(realised={"16": sr}, fractions=partsC["fractions"])
        arrays[f"{name}@R16"] = ar
        log(f"{name:22s} R16 cos {sr['cos_agg']:.6f} pds {sr['pds']:.6f} K {K.mean():.2f} fr {({k: round(v, 3) for k, v in partsC['fractions'].items()})}")
        if (mc, gh) == (1.5, 0.56):
            ref = c420["stack|V1_5_3+sh1.75"]["realised"]["16"] if hi == 3.0 else b540["V1_0.5_6"]["realised"]["16"]
            dd = cmp(sr, ref)
            log(f"KA production mix cap {hi:g}: realised |d| {dd:.1e}")
            assert dd < TOL, f"KA FAILED {held} {name}: {dd}"
            res["ka"][name] = dd
    res["peak_rss_gib"] = E062.peak_rss_gib()
    (WORK / f"arms545-{held}.json").write_text(json.dumps(res, indent=1, default=float))
    np.savez_compressed(WORK / f"arrays545-{held}.npz", names=np.asarray(list(arrays)), d=np.asarray((truth * truth).sum(1)),
                        **{f"{nm}__{k}": v for nm, a in arrays.items() for k, v in a.items()})
    log(f"done; peak RSS {res['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    run(sys.argv[1])
