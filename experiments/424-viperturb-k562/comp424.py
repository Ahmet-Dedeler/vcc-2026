"""424 — VIPerturb-seq K562 in the donor base, priced on PLOLO with 420's instruments, one fold per process.

Every arm is 417a's bulk (PP composition, 410 key clip (0.3, 3) HELD FIXED at the fold's value, sharpening q 1.25,
ballast realisation KMAX 16) with only the donor base changed (lib424.arm_view). The carriers (DE design) are the
arm's own unkeyed, unsharpened PP rows, as a builder that changes the base would make them.

KNOWN ANSWERS FIRST (abort on failure, before any arm number):
  KA-420   my base `PP|K3+sh1.25` (design and realised K16) and `PP|K3` == comp420-<fold>.json, |d| = 0 on cos_agg, pds,
           pds lam 0.5 / 1 and the per-frame render A / B; the fold's K3 weights == 418's saved W3.
  KA-VIEW  arm_view's machinery with VIPerturb switched off reproduces the base `pre` bit-for-bit.
CONTROL: `CTRL:shuf` ((i) with VIPerturb rows permuted across targets) must lose.
SECONDARY: `i:extra` with the 410 key REFIT on sources + VIPerturb ("what a key refit would change").

    uv run python comp424.py <fold>        # k562_gwps | orion_hct116 | orion_hek293t | kolf2
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
import time  # noqa: E402
from types import SimpleNamespace  # noqa: E402

import numpy as np  # noqa: E402

from lib424 import (ALPHA, ARMS, E062, TAG, VIP_PATH, FOLDS, GAMMA, LAM, PROD_OPS, Q, UNION, VIP, W418, W420, WORK, Donors,  # noqa: E402
                    Operator, arm_stats, arm_view, betas_for, build_fold, cd4_frozen_and_partial, cmp, compose_arm,
                    de_profiles, deflate, design_and_realised, key_scale, pre_normrel, score_rows, shared_energy_key,
                    source_residual, union_rows, vip_full_rows, weights)

TOL = 0.0


def run(held: str) -> None:
    t0 = time.time()
    log = lambda *a: print(f"[{held}]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    fold = build_fold(held)
    scored = np.asarray(fold["scored"])
    pp_pos = fold["panel_positions"]
    wcols = np.asarray(fold["wcols"])
    truth_full = np.asarray(fold["truth"], dtype=np.float64)
    truth = truth_full[:, scored]
    pre_cached = np.asarray(fold["pre"], dtype=np.float64)
    targets = list(fold["targets"])
    sources = list(FOLDS[held]["sources"])
    op = Operator(PROD_OPS / "whitener_union_pruned.npz" if FOLDS[held]["op_kind"] == "production" else FOLDS[held]["op"])
    nonpanel = np.ones(truth_full.shape[1], bool)
    nonpanel[pp_pos] = False
    betas = betas_for(held, nonpanel)
    vip_ok = held not in ("k562_gwps",)         # VIPerturb is K562: excluded when K562 is the held truth
    res = dict(fold=held, n_targets=len(targets), sources=sources, vip_used=vip_ok, ka={}, arms={}, stats={})
    arrays = {}

    dn = Donors(targets, sources)
    frozen, partial = cd4_frozen_and_partial(targets)
    cd4, _ = deflate(frozen, 5)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    pre228 = dn.pre(sources, False)
    pre228[:, pp_pos] = 0.0
    nocov = sum(dn.A[s].astype(int) for s in sources) == 0
    pre228[nocov] = pre_cached[nocov]
    fb = pre228
    pre_base = pre_normrel(dn, sources, "cap", -1.5, fallback_rows=fb)
    key410, carried410, _ = shared_energy_key(targets, sources, op, wcols, cols=None)
    W3 = weights(key410, carried410, 0.3, 3.0)[0]
    z418 = np.load(W418 / f"arrays418-{held}.npz")
    assert np.array_equal(W3, z418["W3"]), "KA: K3 weights != 418's W3"
    ref420 = json.loads((W420 / f"work/comp420-{held}.json").read_text())["arms"]

    def rec(name, pre, w=W3, note=""):
        whit = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
        final, parts = compose_arm(whit, cd4, partial, h1, wcols, pp_pos)
        P_DE = de_profiles(final)
        bulk = key_scale(final, w)[0]
        dr = design_and_realised(bulk, P_DE, q=Q, kmaxes=(16,))
        sd, ad = score_rows({c: dr[c]["design"] for c in "ABC"}, scored, truth, betas, nonpanel)
        sr, ar = score_rows({c: dr[c]["realised"][16] for c in "ABC"}, scored, truth, betas, nonpanel)
        K = np.concatenate([dr[c]["stats"][16]["K"] for c in "ABC"])
        sr["frac_capped"], sr["K_mean"] = float(np.mean(K >= 16)), float(np.mean(K))
        res["arms"][name] = dict(design=sd, realised={"16": sr}, fractions=parts["fractions"], s=parts["s"], note=note)
        arrays[name] = ad
        arrays[f"{name}@R16"] = ar
        log(f"{name:28s} cos {sd['cos_agg']:.6f} pds@.5 {sd['pds_lam0.5']:.6f} | R16 cos {sr['cos_agg']:.6f} "
            f"pds@.5 {sr['pds_lam0.5']:.6f} pds {sr['pds']:.6f} capped {100 * sr['frac_capped']:.0f}% | "
            f"f_cd4 {parts['fractions']['cd4']:.4f} s {parts['s']:.4f}")
        return final

    # ---------------------------------------------------------------- KA-420: the base reproduces 420 bit-for-bit
    rec("PP|K3+sh1.25", pre_base)
    d_real = cmp(res["arms"]["PP|K3+sh1.25"]["realised"]["16"], ref420["PP|K3+sh1.25"]["realised"]["16"])
    d_des = cmp(res["arms"]["PP|K3+sh1.25"]["design"], ref420["PP|K3+sh1.25"]["design"])
    log(f"KA-420 PP|K3+sh1.25 vs comp420: realised K16 |d| {d_real:.1e}, design |d| {d_des:.1e}")
    assert d_real == TOL and d_des == TOL, "KA-420 FAILED"
    res["ka"]["ka420"] = dict(realised=d_real, design=d_des)
    # ---------------------------------------------------------------- KA-VIEW: the view machinery with VIP off == base
    view0 = SimpleNamespace(F=dict(dn.F), A=dict(dn.A), med=dict(dn.med), pmcorr=dn.pmcorr)
    assert np.array_equal(pre_normrel(view0, sources, "cap", -1.5, fallback_rows=fb), pre_base), "KA-VIEW"
    log("KA-VIEW PASS; ALL KNOWN ANSWERS PASS")

    effv, availv, vaudit = vip_full_rows()
    res["vip_audit"] = vaudit
    if vip_ok:
        rng = np.random.default_rng(424)
        for arm in ARMS:
            view, subset = arm_view(dn, targets, sources, arm, rng=rng)
            res["stats"][arm] = arm_stats(dn, targets, sources, arm)
            pre = pre_normrel(view, subset, "cap", -1.5, fallback_rows=fb)
            changed = np.flatnonzero(np.any(pre != pre_base, axis=1))
            res["stats"][arm]["rows_changed"] = int(changed.size)
            rec(f"{arm}|K3+sh1.25", pre, note=json.dumps(res["stats"][arm]))
        # secondary: (i) with the 410 key refit on sources + VIPerturb (same fold operator)
        keyv, carv, _ = shared_energy_key(targets, sources + [VIP], op, wcols, cols=None)
        Wv = weights(keyv, carv, 0.3, 3.0)[0]
        res["stats"]["key_refit"] = dict(carried_before=int(carried410.sum()), carried_after=int(carv.sum()),
                                         spearman_w=float(np.corrcoef(np.argsort(np.argsort(W3)), np.argsort(np.argsort(Wv)))[0, 1]))
        view, subset = arm_view(dn, targets, sources, "i:extra")
        rec("i:extra|K3refit+sh1.25", pre_normrel(view, subset, "cap", -1.5, fallback_rows=fb), w=Wv,
            note="(i) + 410 key refit on sources + VIPerturb")
        rec("PP|K3refit+sh1.25", pre_base, w=Wv, note="base + the refit key only (isolates the key change)")
    else:
        log("VIPerturb excluded on this fold (held line is K562): every arm == base by construction")
    res["peak_rss_gib"] = E062.peak_rss_gib()
    (WORK / f"comp424{TAG}-{held}.json").write_text(json.dumps(res, indent=1, default=float))
    np.savez_compressed(WORK / f"arrays424{TAG}-{held}.npz", names=np.asarray(list(arrays)), d=np.asarray((truth * truth).sum(1)),
                        **{f"{n}__{k}": v for n, a in arrays.items() for k, v in a.items()})
    log(f"done; {len(res['arms'])} arms; peak RSS {res['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    run(sys.argv[1])
