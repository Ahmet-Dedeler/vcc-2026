"""420 on the REAL A/B/C panel: render-space energies of every bulk design (composition x key x sharpening q), both as
DESIGNED (418's route: plolo.deliver of the design) and as REALISED by the ballast (414's emulation, KMAX 16/24/32), with
the ballast block statistics (fraction of blocks at the cap, clip mass, cos(realised, design)).

KNOWN ANSWERS (abort on failure, before any design number is written):
  KA-PP     the PP constants reproduce the 410b build meta (< 1e-5 rel), deflation energy and U_k sha exact, and
            compose_arm(defaults) == sens415.compose_pp bit-for-bit (abc418's KA).
  KA-KEY    every key artifact's sha is the pinned one; the new V1 (0.5,3) artifact is 415's key re-clipped (asserted).
  KA-418    my design energies of PP|none, PP|K3, PP|K3+sh, b:h1x2|K3+sh, a:cd4x1.5|K3+sh == abc418.json (1e-12 rel).
  KA-414b   realise_k(kmax 16) == ballast414.realise bit-for-bit on every block of PP|K3+sh; PP|K3 and PP|K3+sh block
            stats == ballast418.json.
VALIDATION of the energy routes against the four MEASURED artifacts (render-space A_c from the shipped h5ad files):
  395 (PP|none, welded, no ballast), 410b (PP|K3 ballast), 410c (PP|Kinf ballast), 417a (PP|K3+sh ballast).

    uv run python abc420.py            # -> work/abc420.json
    uv run python abc420.py --level    # per-context scale maps for re-levelling (after the main run) -> work/abc420_level.json
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np

from lib420 import (LB, ALPHA, COMP418, COMPS, E062, GAMMA, KEYS, KMAXES, LAM, QS, ROOT, SHARP_Q, UNION, W414, W418, WORK,
                    Donors, Operator, abc_weights, cd4_frozen_and_partial, compose_arm, compose_pp, de_profiles, deflate,
                    design_and_realised, energy, frames, key_scale, make_v1_53, pre_normrel, realise414, realise_k, rowcos,
                    sharpen, source_residual, union_rows)
import sens415
from plolo import TS, profile_from_effect

MAIN = ("k562_gwps", "orion_hct116", "orion_hek293t")
MEASURED = {  # render-space A_c/B_c of the shipped artifacts (414's decomp414 / cliffcheck414, read-only)
    "395": ("decomp", "395"), "410b": ("decomp", "410b"),
    "410c": ("cliff", "the-kettle-applied-for-tenure.h5ad"), "417a": ("cliff", "the-pigeon-passed-peer-review.h5ad"),
}


def measured():
    dec = json.loads((W414 / "decomp414.json").read_text())
    cc = json.loads((W414 / "cliffcheck414.json").read_text())
    out = {}
    for k, (src, key) in MEASURED.items():
        out[k] = ({c: dict(A=dec[key][c]["A"], B=dec[key][c]["B"]) for c in "ABC"} if src == "decomp"
                  else {c: dict(A=cc[key]["per_context"][c]["A"], B=cc[key]["per_context"][c]["B"]) for c in "ABC"})
    return out


def setup():
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
    assert np.array_equal(pre, dn.pre(MAIN, True)), "pre_normrel != Donors.pre on A/B/C"
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
    assert abs(daud["energy_removed_frac"] - M["cd4"]["deflation"]["energy_removed_frac"]) < 1e-9, "KA-PP deflation"
    whit = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    fpp, _ = compose_arm(whit, cd4, partial, h1, wcols, pp_pos)
    assert np.array_equal(fpp, f415), "KA-PP: compose_arm(defaults) != compose_pp on A/B/C"
    finals = {"PP": fpp}
    for name, (m, g) in COMPS.items():
        if name != "PP":
            finals[name], parts = compose_arm(whit, cd4, partial, h1, wcols, pp_pos, m_cd4=m, g_h1=g)
            assert not np.array_equal(finals[name], fpp)
    del op
    return dict(targets=targets, pp_pos=pp_pos, nonpanel=nonpanel, finals=finals, ka_pp_worst=worst)


def stats_summary(st, cosr):
    K = st["K"]
    return dict(K_median=float(np.median(K)), K_mean=float(np.mean(K)), frac_capped=float(np.mean(K >= st["kmax"])),
                clip_mass_q=np.quantile(st["clip_mass"], [0.5, 0.9, 1.0]).tolist(),
                clip_genes_max=int(np.max(st["clip_genes"])),
                cos_realised_vs_design_q=np.quantile(cosr, [0.01, 0.1, 0.5]).tolist(),
                r_min_q=np.quantile(st["r"], [0.0, 0.1, 0.5]).tolist())


def main():
    t0 = time.time()
    log = lambda *a: print("[abc420]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    S = setup()
    targets, pp_pos, nonpanel, finals = S["targets"], S["pp_pos"], S["nonpanel"], S["finals"]
    v1path, v1sha = make_v1_53()
    W = {k: abc_weights(k, targets) for k in (*KEYS, "Kinf")}
    log(f"KA-PP worst {S['ka_pp_worst']:.2e}; keys loaded (V1_5_3 = {v1path.name} sha {v1sha[:12]}); "
        f"w floor frac: " + ", ".join(f"{k} {np.mean(w <= 0.5 + 1e-12):.3f}" for k, w in W.items()))
    P_DE = de_profiles(finals["PP"])
    abc418 = json.loads((W418 / "abc418.json").read_text())["arms"]
    b418 = json.loads((W418 / "ballast418.json").read_text())

    res = dict(ka={}, designs={}, measured=measured(), v1_53=dict(path=str(v1path.relative_to(ROOT)), sha256=v1sha))

    def run(comp, key, q, kmaxes=KMAXES, tag=None, scale=None):
        bulk = finals[comp] if key is None else key_scale(finals[comp], W[key])[0]
        dr = design_and_realised(bulk, P_DE, q=q, scale=scale, kmaxes=kmaxes, zero_cols=pp_pos)
        e = dict(comp=comp, key=key, q=q, design={}, realised={k: {} for k in kmaxes}, ballast={})
        for ctx in "ABC":
            e["design"][ctx] = energy(dr[ctx]["design"], nonpanel)
            for k in kmaxes:
                e["realised"][k][ctx] = energy(dr[ctx]["realised"][k], nonpanel)
        for k in kmaxes:
            st = {f: np.concatenate([dr[c]["stats"][k][f] for c in "ABC"]) for f in ("K", "r", "clip_mass", "clip_genes")}
            st["kmax"] = k
            cosr = np.concatenate([rowcos(dr[c]["design"][:, nonpanel], dr[c]["realised"][k][:, nonpanel]) for c in "ABC"])
            e["ballast"][k] = stats_summary(st, cosr)
        name = tag or f"{comp}|{key or 'none'}" + ("" if q is None else f"+sh{q:g}")
        res["designs"][name] = e
        b = e["ballast"][kmaxes[0]]
        log(f"{name:28s} A des {[round(e['design'][c]['A'], 4) for c in 'ABC']} real16 "
            f"{[round(e['realised'][kmaxes[0]][c]['A'], 4) for c in 'ABC']} | capped {100 * b['frac_capped']:.1f}% "
            f"cos q1 {b['cos_realised_vs_design_q'][0]:.4f}")
        return e, dr

    # ------------------------------------------------------------------ KA-414b: realise_k(16) == ballast414.realise
    bulk = key_scale(finals["PP"], W["K3"])[0]
    nb = 0
    for ctx, (b, gate) in frames().items():
        F = sharpen(bulk, SHARP_Q, gate, mode="row")
        for i in range(F.shape[0]):
            p_bk = profile_from_effect(F[i], b, gate, "minnorm")
            a1, a2 = realise414(P_DE[ctx][i], p_bk), realise_k(P_DE[ctx][i], p_bk, 16)
            assert np.array_equal(a1[0], a2[0]) and a1[1:] == a2[1:], f"KA-414b realise {ctx} {i}"
            nb += 1
    log(f"KA-414b: realise_k(16) == ballast414.realise on all {nb} blocks of PP|K3+sh")

    # ------------------------------------------------------------------ KA-418 design energies + ballast stats
    worst = 0.0
    for comp, key, q, ref in (("PP", None, None, ("PP", "none")), ("PP", "K3", None, ("PP", "K3")),
                              ("PP", "K3", SHARP_Q, ("PP", "K3+sh")), ("h1x2", "K3", SHARP_Q, ("b:h1x2", "K3+sh")),
                              ("cd4x1.5", "K3", SHARP_Q, ("a:cd4x1.5", "K3+sh"))):
        e, _ = run(comp, key, q)
        r = abc418[ref[0]][ref[1]]
        d = max(abs(e["design"][c][x] / r[c][x] - 1) for c in "ABC" for x in ("A", "B"))
        worst = max(worst, d)
        log(f"KA-418 {comp}|{key}|{q}: design energy vs abc418 max rel |d| {d:.1e}")
        assert d < 1e-12, f"KA-418 FAILED {comp} {key} {q}: {d}"
    for nm, tg in (("PP|K3", "PP|K3"), ("PP|K3+sh1.25", "PP|K3+sh")):
        mine, ref = res["designs"][nm]["ballast"][16], b418[tg]
        d = max(abs(mine["K_median"] - ref["K_median"]), abs(mine["frac_capped"] - ref["frac_K_capped"]),
                *[abs(a - bb) for a, bb in zip(mine["cos_realised_vs_design_q"], ref["cos_realised_vs_design_q"])],
                *[abs(a - bb) for a, bb in zip(mine["clip_mass_q"], ref["clip_mass_q"])])
        log(f"KA-414b ballast stats {nm} vs ballast418 {tg}: max |d| {d:.1e}")
        assert d < 1e-12, f"KA ballast418 FAILED {nm}: {d}"
    res["ka"] = dict(pp_worst_rel_constant=S["ka_pp_worst"], abc418_worst=worst, realise_blocks=nb)
    log("ALL A/B/C KNOWN ANSWERS PASS")

    # ------------------------------------------------------------------ validation: 410c = PP|Kinf through ballast
    run("PP", "Kinf", None)

    # ------------------------------------------------------------------ the grid
    for comp in COMPS:
        for key in KEYS:
            for q in (None, *QS):
                name = f"{comp}|{key}" + ("" if q is None else f"+sh{q:g}")
                if name not in res["designs"]:
                    run(comp, key, q)
        (WORK / "abc420.json").write_text(json.dumps(res, indent=1, default=float))
    res["peak_rss_gib"] = E062.peak_rss_gib()
    (WORK / "abc420.json").write_text(json.dumps(res, indent=1, default=float))
    log(f"done; {len(res['designs'])} designs; peak RSS {res['peak_rss_gib']:.2f} GiB")


def level(names):
    """Per-context scale maps: realised render energy at effect-space scale s (VCC_BULK_CTX_SCALE) for chosen designs."""
    t0 = time.time()
    S = setup()
    targets, pp_pos, nonpanel, finals = S["targets"], S["pp_pos"], S["nonpanel"], S["finals"]
    W = {k: abc_weights(k, targets) for k in KEYS}
    P_DE = de_profiles(finals["PP"])
    grid = (0.94, 0.96, 0.97, 0.98, 0.99, 1.0, 1.01, 1.02, 1.03, 1.04, 1.06, 1.08, 1.10, 1.12, 1.15, 1.18, 1.20)
    out = {}
    for name in names:
        comp, rest = name.split("|")
        key, q = (rest.split("+sh") + [None])[:2]
        q = float(q) if q else None
        bulk = key_scale(finals[comp], W[key])[0]
        rows = {}
        for s in grid:
            dr = design_and_realised(bulk, P_DE, q=q, scale={c: s for c in "ABC"}, kmaxes=(16,), zero_cols=pp_pos)
            rows[s] = {ctx: dict(design=energy(dr[ctx]["design"], nonpanel), realised=energy(dr[ctx]["realised"][16], nonpanel))
                       for ctx in "ABC"}
        out[name] = {str(s): v for s, v in rows.items()}
        print(f"[level] {name}: realised A ratio at s " + ", ".join(
            f"{s:g}: {rows[s]['A']['realised']['A'] / rows[1.0]['A']['realised']['A']:.4f}" for s in grid) +
            f" ({time.time() - t0:.0f}s)", flush=True)
    prev = json.loads((WORK / "abc420_level.json").read_text()) if (WORK / "abc420_level.json").exists() else {}
    prev.update(out)
    (WORK / "abc420_level.json").write_text(json.dumps(prev, indent=1))


MEAN_LIB = 7_919_901.0 / 384      # 410b's median carrier total over 384 cells (ballast414.S_C) -> mean carrier library


def leak(names):
    """DE-channel leak of the ballast, in units of 410's full WELD (board: -0.0171 overall, nmae +0.054 raw).

    A block's per-cell arithmetic mean CPM (what nmae's log2FC and the DE means read) is
        mean = p_de + L_K (p_re - p_de),   L_K = (K + t_c / LB) / 400,   t_c = MEAN_LIB (400 - K)
    (carriers ~ p_de per cell; K ballast cells of LB counts carry S p_re - t_c p_de). A welded arm has L = 1 and p_re = the
    keyed design. leak = mean_i L_i |log2 mean_i - log2 p_de_i|_1 over carrier-live genes, divided by the same for the
    410 weld (PP|K3 unsharpened, L = 1). Also reports the ballast share of the rank sample, mean K / 400."""
    S = setup()
    targets, pp_pos, nonpanel, finals = S["targets"], S["pp_pos"], S["nonpanel"], S["finals"]
    W = {k: abc_weights(k, targets) for k in KEYS}
    P_DE = de_profiles(finals["PP"])
    out = {}

    def dev(comp, key, q, kmax, weld=False, scale=None):
        bulk = key_scale(finals[comp], W[key])[0]
        tot, ks, n = 0.0, [], 0
        for ctx, (b, gate) in frames().items():
            F = bulk if q is None else sharpen(bulk, q, gate, mode="row")
            if scale is not None:
                F = F * float(scale[ctx])
            for i in range(F.shape[0]):
                p_de = P_DE[ctx][i]
                p_bk = profile_from_effect(F[i], b, gate, "minnorm")
                live = p_de * 1e6 > 5.0                      # the DE gene filter: > 5 CPM in the carriers' design
                if weld:
                    L, p_re, K = 1.0, p_bk, 0
                else:
                    p_re, r, K, cm, cg = realise_k(p_de, p_bk, kmax)
                    L = (K + MEAN_LIB * (400 - K) / LB) / 400.0
                mean = p_de + L * (p_re - p_de)
                tot += float(np.abs(np.log2(np.maximum(mean[live], 1e-12)) - np.log2(p_de[live])).sum())
                ks.append(K)
                n += 1
        return tot / n, float(np.mean(ks))

    ref, _ = dev("PP", "K3", None, 16, weld=True)
    out["weld410_L1"] = ref
    prev = json.loads((WORK / "abc420_leak.json").read_text()) if (WORK / "abc420_leak.json").exists() else {}
    for name in names:
        nm, _, sc = name.partition(":")
        scale = {kv.split("=")[0]: float(kv.split("=")[1]) for kv in sc.split(",")} if sc else None
        comp, rest = nm.split("|")
        key, q = (rest.split("+sh") + [None])[:2]
        q = float(q) if q else None
        out[name] = {}
        for k in KMAXES:
            d, km = dev(comp, key, q, k, scale=scale)
            out[name][str(k)] = dict(leak_vs_410weld=d / ref, mean_K=km, ballast_rank_share=km / 400)
        print(f"[leak] {name:22s} " + "  ".join(f"K{k}: leak {out[name][str(k)]['leak_vs_410weld']:.4f} meanK "
                                                  f"{out[name][str(k)]['mean_K']:.2f}" for k in KMAXES), flush=True)
    assert "weld410_L1" not in prev or abs(prev["weld410_L1"] - ref) < 1e-12, "leak reference changed"
    prev.update(out)
    (WORK / "abc420_leak.json").write_text(json.dumps(prev, indent=1))


if __name__ == "__main__":
    if "--leak" in sys.argv:
        leak([a for a in sys.argv[1:] if a != "--leak"])
    elif "--level" in sys.argv:
        level([a for a in sys.argv[1:] if a != "--level"])
    else:
        main()
