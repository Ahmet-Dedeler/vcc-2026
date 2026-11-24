"""418: can the 410b ballast construction REALISE the top composition arms on the REAL A/B/C panel?

414's `ballast414.realise` (the builder's `ballast_count`: K = min(16, ceil(s_c (1/r_min - 1)/995000)), genes whose
ballast-inflated total falls below the carriers' own totals are clipped), with the carriers' design = the shipped
composition UNKEYED (PP, 395's design; the carriers keep it in a bulk-only arm) and the ballast design = each arm keyed
by the shipped 410 weights (and sharpened for the 417a variants). Reported per arm: K median, the fraction of blocks at
K = 16 (the cap), clipped mass, and the median render cosine between the realised and the designed block.
KNOWN ANSWER (loose; the emulation omits the carriers' own sign/emission design): PP|K3 should sit near the 410b build's
ballast_K median 2 / max 16 and clip genes max 8.

    uv run python ballast418.py
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np

from lib418 import (ALPHA, E062, GAMMA, LAM, ROOT, SHARP_Q, UNION, WORK, Donors, Operator, cd4_frozen_and_partial,
                    compose_arm, deflate, frames, key_scale, pre_normrel, sharpen, source_residual, union_rows)
import sens415
sys.path.insert(0, str(ROOT / "experiments/414-bulk-reprice"))
from ballast414 import realise  # noqa: E402
from plolo import TS, profile_from_effect  # noqa: E402

KEY = ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.3_hi3.npz"
MAIN = ("k562_gwps", "orion_hct116", "orion_hek293t")
ARMS = {"PP": (1.0, 0.28), "a:cd4x1.25": (1.25, 0.28), "a:cd4x1.5": (1.5, 0.28), "a:cd4x1.75": (1.75, 0.28),
        "b:h1x1.5": (1.0, 0.42), "b:h1x2": (1.0, 0.56), "b:h1x2.5": (1.0, 0.70)}


def main():
    t0 = time.time()
    targets = list(sens415.PANEL)
    pos = {g: i for i, g in enumerate(sens415.OFFICIAL)}
    pp_pos = np.asarray([pos[g] for g in targets], dtype=np.int64)
    wcols = np.asarray(np.load(WORK / "plolo/fold-kolf2.npz")["wcols"])
    op = Operator(ROOT / "data/donors/whitener_union_pruned.npz", name="union_pruned")
    dn = Donors(targets, MAIN)
    frozen, partial = cd4_frozen_and_partial(targets)
    cd4, _ = deflate(frozen, 5)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    z = np.load(KEY)
    kidx = {str(t): i for i, t in enumerate(z["targets"])}
    w = np.asarray([float(z["w"][kidx[t]]) for t in targets])
    whit = op.apply(pre_normrel(dn, MAIN, "cap", -1.5) * ALPHA, wcols, GAMMA, LAM)
    del op
    fpp, _ = compose_arm(whit, cd4, partial, h1, wcols, pp_pos)
    abc = json.loads((WORK / "abc418.json").read_text())["arms"]
    P_DE = {ctx: [profile_from_effect(fpp[i], b, gate, "minnorm") for i in range(fpp.shape[0])]
            for ctx, (b, gate) in frames().items()}
    out = {}
    for name, (m, gh) in ARMS.items():
        final, _ = compose_arm(whit, cd4, partial, h1, wcols, pp_pos, m_cd4=m, g_h1=gh)
        bulk = key_scale(final, w)[0]
        for q in (None, SHARP_Q):
            tag = f"{name}|K3" + ("+sh" if q else "")
            d = dict(r_min=[], K=[], clip_mass=[], clip_genes=[], cos_real=[])
            A_des = {}
            for ctx, (b, gate) in frames().items():
                F = bulk if q is None else sharpen(bulk, q, gate, mode="row")
                H = []
                for i in range(F.shape[0]):
                    p_bk = profile_from_effect(F[i], b, gate, "minnorm")
                    p_re, r, k, cm, cg = realise(P_DE[ctx][i], p_bk)
                    hd, hr = np.log1p(TS * p_bk) - b, np.log1p(TS * p_re) - b
                    hd[pp_pos] = 0.0
                    hr[pp_pos] = 0.0
                    H.append(hd)
                    d["r_min"].append(r); d["K"].append(k); d["clip_mass"].append(cm); d["clip_genes"].append(cg)
                    d["cos_real"].append(float(hd @ hr / np.sqrt((hd @ hd) * (hr @ hr))))
                A_des[ctx] = float(np.mean([h @ h for h in H]))
            ref = abc[name]["K3+sh" if q else "K3"]
            ka = max(abs(A_des[c] / ref[c]["A"] - 1) for c in "ABC")    # the design render == abc418's (same path)
            assert ka < 1e-9, f"KA design render != abc418 {tag}: {ka}"
            out[tag] = dict(K_median=float(np.median(d["K"])), K_max=int(np.max(d["K"])),
                            frac_K_capped=float(np.mean(np.asarray(d["K"]) >= 16)),
                            r_min_q=np.quantile(d["r_min"], [0.0, 0.1, 0.5]).tolist(),
                            clip_mass_q=np.quantile(d["clip_mass"], [0.5, 0.9, 1.0]).tolist(),
                            clip_genes_max=int(np.max(d["clip_genes"])),
                            cos_realised_vs_design_q=np.quantile(d["cos_real"], [0.01, 0.1, 0.5]).tolist())
            e = out[tag]
            print(f"[ballast] {tag:18s} K median {e['K_median']:.0f} max {e['K_max']} capped {100 * e['frac_K_capped']:.1f}% "
                  f"| r_min q0/10/50 {np.round(e['r_min_q'], 3).tolist()} | clip mass q50/90/max "
                  f"{e['clip_mass_q'][0]:.1e}/{e['clip_mass_q'][1]:.1e}/{e['clip_mass_q'][2]:.1e} | cos(realised, design) "
                  f"q1/10/50 {np.round(e['cos_realised_vs_design_q'], 5).tolist()} ({time.time() - t0:.0f}s)", flush=True)
    out["peak_rss_gib"] = E062.peak_rss_gib()
    (WORK / "ballast418.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
