"""492 — per-target sharpening exponent by the 410 key tercile, PLOLO, one fold per process (473's path, lib420).

Base = `stack|K3+sh1.75`. Arms set q per target from the fold's K3 weight tercile (see PREREGISTRATION.md). Rows are
sharpened with their own q (169's `sharpen(mode="row")` on each frame's gate, each row's gated norm restored), rendered
(minnorm) and realised by 414's ballast emulation at KMAX 16 against PP's unkeyed carriers.
KNOWN ANSWERS: the per-target path at q 1.75 for every row == comp420's `stack|K3+sh1.75` (design, realised 16), and at
q 2 == comp420's `stack|K3+sh2`, |d| < 1e-12; then 473's KA-PP / K3-weight checks.

    uv run python experiments/492-per-target-sharpening/arms492.py <fold>
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
ARMS = {"P1": (2.0, 1.75, 1.75), "P2": (2.5, 2.0, 1.75), "P3": (3.0, 2.25, 1.75)}
# 492b follow-up grid (written after 492's preregistered read; reported as exploratory)
ARMS_B = {"L25": (2.5, 1.75, 1.75), "L30": (3.0, 1.75, 1.75), "L40": (4.0, 1.75, 1.75), "L30M20": (3.0, 2.0, 1.75)}


def qvec_of(w, spec):
    edges = np.quantile(w, [1 / 3, 2 / 3])
    t = np.digitize(w, edges)            # 0 low, 1 mid, 2 high
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
    W3 = weights(key410, carried410, 0.3, 3.0)[0]
    assert np.array_equal(W3, np.load(W418 / f"arrays418-{held}.npz")["W3"]), "KA K3 weights"
    if MODE_E:
        import lib418 as _l418
        d_rows = {s_: _l418.donor_rows(s_, targets)[0] for s_ in sources}
        A_ = {s_: dn.A[s_] for s_ in sources}
        carried_all = sum(A_[s_].astype(int) for s_ in sources) >= 2
        bases = _l418.lodo_bases(op, dn, sources, "cap", -1.5, pre, wcols)
        kv1 = _l418.lodo_key(fPP, pPP, bases, d_rows, A_, carried_all, sources, pp_pos)
        WV1 = weights(np.nan_to_num(kv1, nan=0.0), carried_all, 0.3, 3.0)[0]
        assert np.array_equal(WV1, np.load(W418 / f"arrays418-{held}.npz")["V1__PP"]), "KA V1 weights"
        del bases, d_rows
    del whit, f415, pre, pre228, fb, dn
    P_DE = de_profiles(fPP)
    bulk = key_scale(fST, W3)[0]
    log("KA-PP and K3 weights PASS")

    def rows_for(qv, B=None):
        B = bulk if B is None else B
        dr = {}
        for ctx, (b, gate) in frames().items():
            F = np.empty_like(bulk)
            for q in np.unique(qv):
                sel = qv == q
                F[sel] = sharpen(B[sel], float(q), gate, mode="row")
            des, rea = np.empty_like(F), np.empty_like(F)
            K = np.empty(F.shape[0], int)
            for i in range(F.shape[0]):
                p_bk = profile_from_effect(F[i], b, gate, "minnorm")
                des[i] = np.log1p(TS * p_bk) - b
                p_re, _, K[i], _, _ = realise_k(P_DE[ctx][i], p_bk, 16)
                rea[i] = np.log1p(TS * p_re) - b
            dr[ctx] = dict(design=des, realised=rea, K=K)
        return dr

    def rec(name, qv, note="", B=None):
        dr = rows_for(qv, B)
        sd, ad = score_rows({c: dr[c]["design"] for c in "ABC"}, scored, truth, betas, nonpanel)
        sr, ar = score_rows({c: dr[c]["realised"] for c in "ABC"}, scored, truth, betas, nonpanel)
        K = np.concatenate([dr[c]["K"] for c in "ABC"])
        sr["frac_capped"] = float(np.mean(K >= 16))
        res["arms"][name] = dict(note=note, design=sd, realised={"16": sr})
        arrays[name] = ad
        arrays[f"{name}@R16"] = ar
        log(f"{name:24s} R16 cos {sr['cos_agg']:.6f} pds {sr['pds']:.6f} pds@.5 {sr['pds_lam0.5']:.6f} capped {100 * sr['frac_capped']:.0f}%")
        return res["arms"][name]

    c420 = json.loads((W420 / f"comp420-{held}.json").read_text())["arms"]
    n = len(targets)
    if MODE_E:
        bulkV = key_scale(fST, WV1)[0]
        e = rec("stack|V1+sh1.75", np.full(n, 1.75), B=bulkV)
        dD, dR = cmp(e["design"], c420["stack|V1+sh1.75"]["design"]), cmp(e["realised"]["16"], c420["stack|V1+sh1.75"]["realised"]["16"])
        log(f"KA stack|V1+sh1.75: design |d| {dD:.1e}, realised |d| {dR:.1e}")
        assert dD < TOL and dR < TOL, f"KA V1 FAILED {held}: {dD} {dR}"
        res["ka"]["stack|V1+sh1.75"] = dict(design=dD, realised=dR)
        def tg(v, spec, floor=None):
            t = np.digitize(v, np.quantile(v, [1 / 3, 2 / 3]))
            if floor is not None:
                t = np.where(v <= floor + 1e-12, 0, t)
            return np.asarray([spec[k] for k in t], float)
        strength = np.linalg.norm(fST[:, nonpanel], axis=1)
        rec("V1 P2g", tg(WV1, (2.5, 2.0, 1.75), floor=0.3), B=bulkV)
        rec("V1 S_P2", tg(strength, (2.5, 2.0, 1.75)), B=bulkV)
        rec("V1 CTRL P2g permuted", tg(WV1, (2.5, 2.0, 1.75), floor=0.3)[np.random.default_rng(4925).permutation(n)], B=bulkV)
        res["frac_floor_V1"] = float(np.mean(WV1 <= 0.3 + 1e-12))
        res["peak_rss_gib"] = E062.peak_rss_gib()
        (WORK / f"arms492e-{held}.json").write_text(json.dumps(res, indent=1, default=float))
        np.savez_compressed(WORK / f"arrays492e-{held}.npz", names=np.asarray(list(arrays)), d=np.asarray((truth * truth).sum(1)),
                            **{f"{nm}__{k}": v for nm, a in arrays.items() for k, v in a.items()})
        log(f"done (E); peak RSS {res['peak_rss_gib']:.2f} GiB")
        return
    for name, q in (("stack|K3+sh1.75", 1.75), ("stack|K3+sh2", 2.0)):
        e = rec(name, np.full(n, q), note="known answer vs comp420")
        dD, dR = cmp(e["design"], c420[name]["design"]), cmp(e["realised"]["16"], c420[name]["realised"]["16"])
        log(f"KA {name}: design |d| {dD:.1e}, realised |d| {dR:.1e}")
        assert dD < TOL and dR < TOL, f"KA FAILED {held} {name}: {dD} {dR}"
        res["ka"][name] = dict(design=dD, realised=dR)
    log("ALL KNOWN ANSWERS PASS")
    if MODE_D:
        # 492d: the same q values assigned by other per-target properties (exploratory)
        strength = np.linalg.norm(fST[:, nonpanel], axis=1)          # unkeyed stack row norm = response strength
        DSPEC = {"S_L30": ("strength", (3.0, 1.75, 1.75)), "Sinv_L30": ("-strength", (3.0, 1.75, 1.75)),
                 "S_P2": ("strength", (2.5, 2.0, 1.75)), "W_L30": ("w", (3.0, 1.75, 1.75))}
        if MODE_F:          # 492f: the strong tercile's dose under the strength assignment (exploratory)
            DSPEC = {"S_P2h15": ("strength", (2.5, 2.0, 1.5)), "S_P2h125": ("strength", (2.5, 2.0, 1.25)),
                     "S_P3h15": ("strength", (3.0, 2.25, 1.5)), "S_P2": ("strength", (2.5, 2.0, 1.75))}
        for name, (prop, spec) in DSPEC.items():
            v = {"strength": strength, "-strength": -strength, "w": W3}[prop]
            qv, _ = qvec_of(v, spec)
            rec(name, qv, note=f"q by {prop} tercile {spec}")
        res["spearman_w_strength"] = float(np.corrcoef(np.argsort(np.argsort(W3)), np.argsort(np.argsort(strength)))[0, 1])
    if MODE_C:
        for name, (lo, spec) in ARMS_C.items():
            Wlo = weights(key410, carried410, lo, 3.0)[0]
            qv, _ = qvec_of(W3, spec)                     # the tercile stays the K3 one (same ranking)
            rec(name, qv, note=f"key floor {lo}, q by K3 tercile {spec}", B=key_scale(fST, Wlo)[0])
    for name, spec in ({} if (MODE_C or MODE_D) else (ARMS_B if MODE_B else ARMS)).items():
        qv, t = qvec_of(W3, spec)
        rec(name, qv, note=f"q by K3 tercile {spec}")
    if MODE_C or MODE_D:
        pass
    elif not MODE_B:
        qv2, _ = qvec_of(W3, ARMS["P2"])
        rec("CTRL P2 permuted", qv2[np.random.default_rng(492).permutation(n)], note="CONTROL")
    else:
        qv3, _ = qvec_of(W3, ARMS_B["L30"])
        rec("CTRL L30 permuted", qv3[np.random.default_rng(4922).permutation(n)], note="CONTROL")
    res["tercile_counts"] = np.bincount(qvec_of(W3, ARMS["P2"])[1], minlength=3).tolist()
    res["peak_rss_gib"] = E062.peak_rss_gib()
    (WORK / (f"arms492f-{held}.json" if MODE_F else f"arms492d-{held}.json" if MODE_D else f"arms492c-{held}.json" if MODE_C else f"arms492b-{held}.json" if MODE_B else f"arms492-{held}.json")).write_text(json.dumps(res, indent=1, default=float))
    np.savez_compressed(WORK / (f"arrays492f-{held}.npz" if MODE_F else f"arrays492d-{held}.npz" if MODE_D else f"arrays492c-{held}.npz" if MODE_C else f"arrays492b-{held}.npz" if MODE_B else f"arrays492-{held}.npz"), names=np.asarray(list(arrays)), d=np.asarray((truth * truth).sum(1)),
                        **{f"{nm}__{k}": v for nm, a in arrays.items() for k, v in a.items()})
    log(f"done; peak RSS {res['peak_rss_gib']:.2f} GiB")


MODE_B = "--b" in sys.argv
MODE_C = "--c" in sys.argv
MODE_D = "--d" in sys.argv
MODE_E = "--e" in sys.argv
MODE_F = "--f" in sys.argv
if MODE_F:
    MODE_D = True
# 492c exploratory: the key floor with the low-tercile sharpening (tercile from the K3 weights throughout)
ARMS_C = {"F5": (0.5, (1.75, 1.75, 1.75)), "F4L30": (0.4, (3.0, 1.75, 1.75)), "F5L30": (0.5, (3.0, 1.75, 1.75)),
          "F5P2": (0.5, (2.5, 2.0, 1.75)), "F2L30": (0.2, (3.0, 1.75, 1.75))}

if __name__ == "__main__":
    run(sys.argv[1])
