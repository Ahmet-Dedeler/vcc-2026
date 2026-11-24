"""418 on the REAL A/B/C panel: every composition arm built on the shipped composition (PP, 300 targets, 3 production
donors, production operator union_pruned) with the SHIPPED 410 key (pert_scale_410_lo0.3_hi3.npz, sha d7f76a22...)
applied as apply_pert_scale does, rendered in each challenge frame (plolo.deliver, minnorm), with and without 417a's
sharpening. No truth exists here: this reads what the cliff route needs (render-space A_c / B_c ratios vs 410b's and
417a's bulk) and what each knob does to the effect-space row (direction vs scale).

KNOWN ANSWERS (abort on failure): KA-PP (the PP constants reproduce the 410b build meta to < 1e-5 relative, deflation
energy and U_k sha exact, as sens415.abc); the key npz sha256 is 410b's pinned d7f76a22...; the unmodified composition
equals sens415.compose_pp bit-for-bit.

    uv run python abc418.py
"""
from __future__ import annotations

import hashlib
import json
import time

import numpy as np

from lib418 import (ALPHA, E062, GAMMA, LAM, ROOT, SHARP_Q, UNION, WORK, Donors, Operator, cd4_frozen_and_partial,
                    compose_arm, compose_pp, deflate, deliver, frames, key_on_base_only, key_scale, pre_normrel,
                    sharpen, source_residual, union_rows)
from comp418 import ARMS, EXTRA, NR, RETRO
import sens415

KEY = ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.3_hi3.npz"
KEY_SHA = "d7f76a22f62ad6991d7220ce8a65218db38b2d71dcfe3c3781ec6a8f70409e0e"
MAIN = ("k562_gwps", "orion_hct116", "orion_hek293t")


def render_AB(bulk, nonpanel, q=None):
    out = {}
    for ctx, (b, gate) in frames().items():
        F = bulk if q is None else sharpen(bulk, q, gate, mode="row")
        dl = deliver(F, b, gate, "minnorm")
        H = dl[:, nonpanel]
        Hc = H - H.mean(0, keepdims=True)
        out[ctx] = dict(A=float((H * H).sum(1).mean()), B=float((Hc * Hc).sum(1).mean()))
    return out


def shape(final, ref):
    """What a knob does to the effect-space rows vs PP: fraction of rows whose DIRECTION moved, the median cosine over
    those rows, the energy-weighted cosine, and the Frobenius ratio (a pure scale knob moves only the last)."""
    n1, n0 = np.linalg.norm(final, axis=1), np.linalg.norm(ref, axis=1)
    cos = np.einsum("ij,ij->i", final, ref) / np.maximum(n1 * n0, 1e-300)
    moved = (1 - cos) > 1e-9
    return dict(frac_rows_direction_moved=float(moved.mean()),
                median_cos_moved_rows=float(np.median(cos[moved])) if moved.any() else 1.0,
                energy_weighted_cos=float((np.einsum("ij,ij->i", final, ref)).sum() / np.sqrt((n1 ** 2).sum() * (n0 ** 2).sum())),
                frob_ratio=float(np.linalg.norm(final) / np.linalg.norm(ref)),
                row_norm_ratio_q=np.quantile(n1 / np.maximum(n0, 1e-300), [0.1, 0.5, 0.9]).tolist())


def main():
    t0 = time.time()
    log = lambda *a: print("[abc]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    targets = list(sens415.PANEL)
    pos = {g: i for i, g in enumerate(sens415.OFFICIAL)}
    pp_pos = np.asarray([pos[g] for g in targets], dtype=np.int64)
    nonpanel = np.ones(len(sens415.OFFICIAL), bool)
    nonpanel[pp_pos] = False
    wcols = np.asarray(np.load(WORK / "plolo/fold-kolf2.npz")["wcols"])
    op = Operator(ROOT / "data/donors/whitener_union_pruned.npz", name="union_pruned")
    dn = Donors(targets, MAIN)
    frozen, partial = cd4_frozen_and_partial(targets)
    cd4 = {0: frozen}
    cd4[5], daud = deflate(frozen, 5)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    # ---- KA
    sha = hashlib.sha256(KEY.read_bytes()).hexdigest()
    assert sha == KEY_SHA, f"key sha {sha} != 410b's pinned {KEY_SHA}"
    z = np.load(KEY)
    kt = [str(t) for t in z["targets"]]
    kidx = {t: i for i, t in enumerate(kt)}
    w = np.asarray([float(z["w"][kidx[t]]) for t in targets])
    pre = {"cap-1.5": pre_normrel(dn, MAIN, "cap", -1.5)}
    assert np.array_equal(pre["cap-1.5"], dn.pre(MAIN, True)), "pre_normrel != Donors.pre on A/B/C"
    f415, _, c415, _, _ = compose_pp(op, pre["cap-1.5"], cd4[5], partial, h1, wcols, pp_pos)
    M = sens415.META
    prod = dict(before_operator_median_norm=M["whitening"]["before_operator_median_norm"],
                operator_preserve_scale=M["whitening"]["operator_preserve_scale"],
                before_all_median_norm=M["whitening"]["before_all_median_norm"],
                after_all_median_norm=M["whitening"]["after_all_median_norm"], frac_cd4=M["cd4"]["fraction"],
                frac_h1=M["h1_union"]["residual_frame"]["applied_fraction"],
                frac_partial=M["cd4_partial"]["residual_frame"]["applied_fraction"],
                s=M["cd4"]["post_add_preserve"]["post_add_global_preserve_scale"])
    worst = max(abs(c415[k] / v - 1) for k, v in prod.items())
    log(f"KA-PP worst relative constant error vs 410b meta {worst:.2e}; deflation energy {daud['energy_removed_frac']:.12f} "
        f"(meta {M['cd4']['deflation']['energy_removed_frac']:.12f}), U_k sha {daud['basis_u_k_sha256'][:8]} "
        f"(meta {M['cd4']['deflation']['basis_u_k_sha256'][:8]})")
    assert worst < 1e-5 and daud["basis_u_k_sha256"] == M["cd4"]["deflation"]["basis_u_k_sha256"]
    assert abs(daud["energy_removed_frac"] - M["cd4"]["deflation"]["energy_removed_frac"]) < 1e-9
    whit = {"cap-1.5": op.apply(pre["cap-1.5"] * ALPHA, wcols, GAMMA, LAM)}
    fpp, ppp = compose_arm(whit["cap-1.5"], cd4[5], partial, h1, wcols, pp_pos)
    assert np.array_equal(fpp, f415), "compose_arm(defaults) != compose_pp on A/B/C"
    log("ALL A/B/C KNOWN ANSWERS PASS")
    res = dict(ka=dict(worst_rel_constant=worst, key_sha256=sha), arms={})

    for nr in sorted({v[0] for v in ARMS.values()} | {v[0] for v in RETRO.values()}):
        if nr not in whit:
            mode, exp = NR[nr]
            pre[nr] = pre_normrel(dn, MAIN, mode, exp)
            whit[nr] = op.apply(pre[nr] * ALPHA, wcols, GAMMA, LAM)
    for k in {v[1] for v in ARMS.values()} | {v[1] for v in RETRO.values()}:
        if k not in cd4:
            cd4[k], _ = deflate(frozen, k)
    K3pp = key_scale(fpp, w)[0]
    for name, cfg in {**ARMS, **EXTRA, **RETRO}.items():
        nr, k, m, gh, gp = cfg
        final, parts = (fpp, ppp) if name == "PP" else compose_arm(whit[nr], cd4[k], partial, h1, wcols, pp_pos,
                                                                  m_cd4=m, g_h1=gh, g_part=gp)
        e = dict(cfg=list(cfg), fractions=parts["fractions"], s=parts["s"], shape=shape(final, fpp))
        if name.startswith("R:"):
            e["none"] = render_AB(final, nonpanel)
        else:
            bulk = key_scale(final, w)[0]
            e["shape_keyed"] = shape(bulk, K3pp)
            e["K3"] = render_AB(bulk, nonpanel)
            e["K3+sh"] = render_AB(bulk, nonpanel, SHARP_Q)
            if name == "PP":
                e["none"] = render_AB(final, nonpanel)
        res["arms"][name] = e
        log(f"{name:22s} cd4 {parts['fractions']['cd4']:.4f} h1 {parts['fractions']['h1_union']:.4f} part "
            f"{parts['fractions']['cd4_partial']:.4f} s {parts['s']:.4f} | rows moved {e['shape']['frac_rows_direction_moved']:.2f} "
            f"median cos {e['shape']['median_cos_moved_rows']:.5f} ew-cos {e['shape']['energy_weighted_cos']:.5f} "
            f"frob {e['shape']['frob_ratio']:.4f}" + (f" | A(K3) {[round(e['K3'][c]['A'], 4) for c in 'ABC']}" if "K3" in e else ""))
    b1 = key_on_base_only(whit["cap-1.5"], ppp, w, pp_pos, fpp)
    res["arms"]["f1:key-base-only"] = dict(shape_keyed=shape(b1, K3pp), K3=render_AB(b1, nonpanel),
                                           **{"K3+sh": render_AB(b1, nonpanel, SHARP_Q)})
    log("f1 done")
    res["peak_rss_gib"] = E062.peak_rss_gib()
    (WORK / "abc418.json").write_text(json.dumps(res, indent=1))
    log(f"done; peak RSS {res['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main()
