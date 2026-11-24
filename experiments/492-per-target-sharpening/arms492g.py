"""492g (PLOLO) — per-target sharpening on 446bj2's bulk: the stack keyed by the LODO key at clip (0.5, 3), one fold per process.

arms492.py's path (lib420) with the base moved from `stack|K3+sh1.75` (446b) to `stack|V1_5_3+sh1.75` (446bj2). The fold's
LODO weights are 420's: V1_5_3 = weights(kv1, carried, 0.5, 3) = max(418's saved `V1__PP`, 0.5) (clip(x, .5, 3) =
max(clip(x, .3, 3), .5); non-carried targets stay 1.0), so the LODO bases need not be recomputed. Rows are sharpened with
their own q on each frame's gate (169's `sharpen(mode="row")`, each row's gated norm restored), rendered (minnorm) and
realised by 414's ballast emulation at KMAX 16 against PP's unkeyed carriers.

KNOWN ANSWERS (abort unless |d| < 1e-12): the per-target path at q 1.75 for every row == comp420's `stack|V1_5_3+sh1.75`,
and at q 2 == comp420's `stack|V1_5_3+sh2` (designed and realised 16); KA-PP compose; K3 weights == 418's W3.

Arms (q per tercile low, mid, high; the tercile key is the fold's V1_5_3 w, the floor guard sends w <= 0.5 to LOW):
  P2g (2.5, 2.0, 1.75) | P2 without the guard | S_P2 (strength tercile: the unkeyed stack row's non-panel norm) |
  L30g (3.0, 1.75, 1.75) | P2g_K3key (the K3 tercile, guard 0.3, on the LODO bulk) | CTRL = P2g's q permuted (seed 492).

    uv run python experiments/492-per-target-sharpening/arms492g.py <fold>        (~5.4 GB peak, ~15 s)
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
                    pre_normrel, profile_from_effect, realise_k, shared_energy_key, sharpen, source_residual,
                    union_rows, weights)
import lodo415  # noqa: E402
from comp420 import cmp, score_rows  # noqa: E402

W420 = ROOT / "experiments/420-bulk-stack/work"
W418 = ROOT / "experiments/418-bulk-composition/work"
WORK = HERE / "work"
TOL = 1e-12
P2, L30 = (2.5, 2.0, 1.75), (3.0, 1.75, 1.75)
BASE = "stack|V1_5_3+sh1.75"


def tercile(v, spec, floor=None):
    """492's knob: np.quantile cut points + np.digitize; with `floor`, targets at the key's floor go to LOW."""
    t = np.digitize(v, np.quantile(v, [1 / 3, 2 / 3]))
    if floor is not None:
        t = np.where(v <= floor + 1e-12, 0, t)
    return np.asarray([spec[k] for k in t], float), t


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
    res = dict(fold=held, n_targets=len(targets), ka={}, arms={})
    arrays = {}
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in lodo415.UNION])
    key410, carried410, _ = shared_energy_key(targets, sources, op, wcols, cols=None)
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
    z418 = np.load(W418 / f"arrays418-{held}.npz")
    W3 = weights(key410, carried410, 0.3, 3.0)[0]
    assert np.array_equal(W3, z418["W3"]), "KA K3 weights"
    WV5 = np.maximum(np.asarray(z418["V1__PP"], np.float64), 0.5)        # 420's V1_5_3 (see the docstring)
    del whit, f415, pre, pre228, fb, dn
    P_DE = de_profiles(fPP)
    bulk = key_scale(fST, WV5)[0]
    strength = np.linalg.norm(fST[:, nonpanel], axis=1)                    # 492d's strength (unkeyed stack row norm)
    log("KA-PP and K3 weights PASS; LODO V1_5_3 weights from 418's V1__PP")

    def rows_for(qv):
        dr = {}
        for ctx, (b, gate) in frames().items():
            F = np.empty_like(bulk)
            for q in np.unique(qv):
                sel = qv == q
                F[sel] = sharpen(bulk[sel], float(q), gate, mode="row")
            des, rea = np.empty_like(F), np.empty_like(F)
            K = np.empty(F.shape[0], int)
            for i in range(F.shape[0]):
                p_bk = profile_from_effect(F[i], b, gate, "minnorm")
                des[i] = np.log1p(TS * p_bk) - b
                p_re, _, K[i], _, _ = realise_k(P_DE[ctx][i], p_bk, 16)
                rea[i] = np.log1p(TS * p_re) - b
            dr[ctx] = dict(design=des, realised=rea, K=K)
        return dr

    def rec(name, qv, note=""):
        dr = rows_for(qv)
        sd, ad = score_rows({c: dr[c]["design"] for c in "ABC"}, scored, truth, betas, nonpanel)
        sr, ar = score_rows({c: dr[c]["realised"] for c in "ABC"}, scored, truth, betas, nonpanel)
        K = np.concatenate([dr[c]["K"] for c in "ABC"])
        sr["frac_capped"] = float(np.mean(K >= 16))
        sr["K_mean_emul"] = float(K.mean())
        res["arms"][name] = dict(note=note, design=sd, realised={"16": sr})
        arrays[name] = ad
        arrays[f"{name}@R16"] = ar
        log(f"{name:22s} R16 cos {sr['cos_agg']:.6f} pds {sr['pds']:.6f} pds@.5 {sr['pds_lam0.5']:.6f} capped {100 * sr['frac_capped']:.0f}%")
        return res["arms"][name]

    c420 = json.loads((W420 / f"comp420-{held}.json").read_text())["arms"]
    n = len(targets)
    for name, q in ((BASE, 1.75), ("stack|V1_5_3+sh2", 2.0)):
        e = rec(name, np.full(n, q), note="known answer vs comp420")
        dD, dR = cmp(e["design"], c420[name]["design"]), cmp(e["realised"]["16"], c420[name]["realised"]["16"])
        log(f"KA {name}: design |d| {dD:.1e}, realised |d| {dR:.1e}")
        assert dD < TOL and dR < TOL, f"KA FAILED {held} {name}: {dD} {dR}"
        res["ka"][name] = dict(design=dD, realised=dR)
    log("ALL KNOWN ANSWERS PASS")
    qp2g, tg = tercile(WV5, P2, 0.5)
    qp2, tn = tercile(WV5, P2)
    qs, ts = tercile(strength, P2)
    rec("P2g", qp2g, note="q by LODO V1_5_3 tercile (2.5, 2.0, 1.75), floor guard 0.5")
    rec("P2 (no guard)", qp2, note="the same without the floor guard")
    rec("S_P2", qs, note="q by strength tercile (2.5, 2.0, 1.75)")
    rec("L30g", tercile(WV5, L30, 0.5)[0], note="q by LODO tercile (3.0, 1.75, 1.75), floor guard 0.5")
    rec("P2g_K3key", tercile(W3, P2, 0.3)[0], note="q by K3 tercile, guard 0.3, on the LODO bulk")
    rec("CTRL P2g permuted", qp2g[np.random.default_rng(492).permutation(n)], note="CONTROL")
    res["frac_floor_V1_5_3"] = float(np.mean(WV5 <= 0.5 + 1e-12))
    res["tercile_counts"] = dict(P2g=np.bincount(tg, minlength=3).tolist(), P2=np.bincount(tn, minlength=3).tolist(),
                                 S=np.bincount(ts, minlength=3).tolist())
    res["guard_changes_assignment"] = bool(not np.array_equal(qp2g, qp2))
    res["spearman_w_strength"] = float(np.corrcoef(np.argsort(np.argsort(WV5)), np.argsort(np.argsort(strength)))[0, 1])
    res["peak_rss_gib"] = E062.peak_rss_gib()
    log(f"floor share {res['frac_floor_V1_5_3']:.3f}; terciles {res['tercile_counts']}; guard changes it: "
        f"{res['guard_changes_assignment']}")
    (WORK / f"arms492g-{held}.json").write_text(json.dumps(res, indent=1, default=float))
    np.savez_compressed(WORK / f"arrays492g-{held}.npz", names=np.asarray(list(arrays)), d=np.asarray((truth * truth).sum(1)),
                        **{f"{nm}__{k}": v for nm, a in arrays.items() for k, v in a.items()})
    log(f"done; peak RSS {res['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    run(sys.argv[1])
