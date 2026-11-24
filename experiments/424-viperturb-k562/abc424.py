"""424 on the REAL A/B/C panel: render-space energies of every VIPerturb arm (designed and ballast-realised, KMAX 16),
the input 420's r417a mse route needs (an arm's energy = 417a's measured energy x realised(arm) / realised(417a design)).

Production donors K562-GWPS + Orion HCT116 + HEK293T, production operator union_pruned (e2f65ecf), 410's K3 artifact
(d7f76a22) held fixed, q 1.25, carriers = the arm's own unkeyed unsharpened PP rows.

KNOWN ANSWERS (abort on failure): KA-PP (the PP constants reproduce the 410b build meta, compose_arm == compose_pp);
KA-ABC420 my base `PP|K3+sh1.25` design and realised-K16 energies == abc420.json (|d| = 0).

    uv run python abc424.py      # -> work/abc424.json
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

from lib424 import (ALPHA, ARMS, E062, TAG, VIP_PATH, GAMMA, LAM, MAIN, Q, ROOT, UNION, VIP, W420, WORK, Donors, Operator,  # noqa: E402
                    arm_stats, arm_view, cd4_frozen_and_partial, compose_arm, compose_pp, de_profiles, deflate,
                    design_and_realised, energy, key_scale, pre_normrel, shared_energy_key, source_residual, union_rows,
                    weights)
import lib420  # noqa: E402
import sens415  # noqa: E402


def main():
    t0 = time.time()
    log = lambda *a: print("[abc424]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    targets = list(sens415.PANEL)
    pos = {g: i for i, g in enumerate(sens415.OFFICIAL)}
    pp_pos = np.asarray([pos[g] for g in targets], dtype=np.int64)
    nonpanel = np.ones(len(sens415.OFFICIAL), bool)
    nonpanel[pp_pos] = False
    wcols = np.asarray(np.load(WORK / "plolo/fold-kolf2.npz")["wcols"])
    op = Operator(ROOT / "data/donors/whitener_union_pruned.npz", name="union_pruned")
    dn = Donors(targets, MAIN)
    frozen, partial = cd4_frozen_and_partial(targets)
    cd4, daud = deflate(frozen, 5)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    pre = pre_normrel(dn, MAIN, "cap", -1.5)
    f415, _, c415, _, _ = compose_pp(op, pre, cd4, partial, h1, wcols, pp_pos)
    M = sens415.META
    prod = dict(before_operator_median_norm=M["whitening"]["before_operator_median_norm"],
                operator_preserve_scale=M["whitening"]["operator_preserve_scale"],
                before_all_median_norm=M["whitening"]["before_all_median_norm"],
                after_all_median_norm=M["whitening"]["after_all_median_norm"], frac_cd4=M["cd4"]["fraction"],
                frac_h1=M["h1_union"]["residual_frame"]["applied_fraction"],
                frac_partial=M["cd4_partial"]["residual_frame"]["applied_fraction"],
                s=M["cd4"]["post_add_preserve"]["post_add_global_preserve_scale"])
    worst = max(abs(c415[k] / v - 1) for k, v in prod.items())
    assert worst < 1e-5 and daud["basis_u_k_sha256"] == M["cd4"]["deflation"]["basis_u_k_sha256"], "KA-PP"
    W3 = lib420.abc_weights("K3", targets)
    ref = json.loads((W420 / "work/abc420.json").read_text())["designs"]["PP|K3+sh1.25"]
    res = dict(ka=dict(pp_worst_rel_constant=worst), designs={}, stats={})

    def run(name, pre_arm, w=W3):
        whit = op.apply(pre_arm * ALPHA, wcols, GAMMA, LAM)
        final, parts = compose_arm(whit, cd4, partial, h1, wcols, pp_pos)
        P_DE = de_profiles(final)
        bulk = key_scale(final, w)[0]
        dr = design_and_realised(bulk, P_DE, q=Q, kmaxes=(16,), zero_cols=pp_pos)
        e = dict(design={c: energy(dr[c]["design"], nonpanel) for c in "ABC"},
                 realised={"16": {c: energy(dr[c]["realised"][16], nonpanel) for c in "ABC"}},
                 fractions=parts["fractions"], s=parts["s"],
                 frac_capped=float(np.mean(np.concatenate([dr[c]["stats"][16]["K"] for c in "ABC"]) >= 16)))
        res["designs"][name] = e
        log(f"{name:26s} A des {[round(e['design'][c]['A'], 4) for c in 'ABC']} real16 "
            f"{[round(e['realised']['16'][c]['A'], 4) for c in 'ABC']} vs 417a-design "
            f"{[round(e['realised']['16'][c]['A'] / ref['realised']['16'][c]['A'], 4) for c in 'ABC']} "
            f"capped {100 * e['frac_capped']:.1f}% f_cd4 {parts['fractions']['cd4']:.4f}")
        return e

    e0 = run("PP|K3+sh1.25", pre)
    d = max(max(abs(e0["design"][c][x] / ref["design"][c][x] - 1), abs(e0["realised"]["16"][c][x] / ref["realised"]["16"][c][x] - 1))
            for c in "ABC" for x in ("A", "B"))
    log(f"KA-ABC420 base vs abc420 PP|K3+sh1.25: max rel |d| {d:.1e}")
    assert d == 0.0, "KA-ABC420 FAILED"
    res["ka"]["abc420"] = d
    rng = np.random.default_rng(424)
    for arm in ARMS:
        view, subset = arm_view(dn, targets, list(MAIN), arm, rng=rng)
        res["stats"][arm] = arm_stats(dn, targets, list(MAIN), arm)
        run(f"{arm}|K3+sh1.25", pre_normrel(view, subset, "cap", -1.5))
    # secondary: the 410 key refit on the three donors + VIPerturb (KA: the refit machinery reproduces the artifact key)
    art = np.load(ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.3_hi3.npz")
    k0, c0, _ = shared_energy_key(targets, list(MAIN), op, wcols, cols=None)
    ka_key = float(np.nanmax(np.abs(np.where(c0, k0, 0) - np.where(c0, art["key"], 0))))
    w0 = weights(k0, c0, 0.3, 3.0)[0]
    log(f"KA-KEY 410 key recomputed on MAIN vs artifact: max |d key| {ka_key:.1e}; w equal {np.array_equal(w0, W3)}")
    res["ka"]["key_recompute"] = dict(max_abs=ka_key, w_equal=bool(np.array_equal(w0, W3)))
    if np.array_equal(w0, W3):
        kv, cv, _ = shared_energy_key(targets, list(MAIN) + [VIP], op, wcols, cols=None)
        Wv = weights(kv, cv, 0.3, 3.0)[0]
        view, subset = arm_view(dn, targets, list(MAIN), "i:extra")
        run("i:extra|K3refit+sh1.25", pre_normrel(view, subset, "cap", -1.5), w=Wv)
        run("PP|K3refit+sh1.25", pre, w=Wv)
        res["stats"]["key_refit"] = dict(carried_before=int(c0.sum()), carried_after=int(cv.sum()))
    res["peak_rss_gib"] = E062.peak_rss_gib()
    (WORK / f"abc424{TAG}.json").write_text(json.dumps(res, indent=1, default=float))
    log(f"done; peak RSS {res['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main()
