"""414: can the 410b ballast construction REALISE each bulk design? (builder `ballast_count`, K <= 16 cells of 995k)

Per block the builder sets K = min(16, ceil(s_c (1/r_min - 1) / 995000)), r_min = min_g profile_bulk / profile_de over live
genes, s_c = the carriers' total (~7.92M, 410b meta median); the block's totals are alloc(S, profile_bulk) with
S = s_c + K * 995000, and genes where that is below the carriers' own totals are CLIPPED at the carriers' totals.
So a design whose r_min < 1/(1 + 16*995000/s_c) = 0.332 cannot be realised exactly. Here, on PLOLO, profile_de = the
unscaled base (the carriers' design, as in 410b) and profile_bulk = each arm; we report r_min, K, the clipped mass, and
re-score pds / cos_agg on the REALISED (clipped) bulk against the design. Known answer: a design equal to the base has
r_min = 1, K = 1, no clipping, and scores bit-identically to arms414's base.

    uv run python ballast414.py
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from arms414 import (ALPHA, CD4_F, E062, FOLDS, GAMMA, LAM, PROD_OPS, UNION, UNION_G, Operator,  # noqa: E402
                     build_fold, cd4_for, endpoints, frames, key_scale, merged, sharpen,
                     shared_energy_key, solve_fraction, solve_fraction_carrying, source_residual, union_rows, weights)
from plolo import TS, profile_from_effect  # noqa: E402

ROOT = HERE.parents[1]
WORK = HERE / "work"
S_C, LB, KMAX = 7_919_901.0, 995_000.0, 16


def realise(p_de, p_bulk):
    live = p_de > 1e-9
    r = float(np.min(p_bulk[live] / p_de[live]))
    K = int(min(KMAX, max(1, int(np.ceil(S_C * (1.0 / max(r, 1e-6) - 1.0) / LB)))))
    S = S_C + K * LB
    want = S * p_bulk
    floor = S_C * p_de
    clipped = want < floor
    tot = np.where(clipped, floor, want)
    return tot / tot.sum(), r, K, float((floor - want)[clipped].sum() / S), int(clipped.sum())


def run_fold(held):
    t0 = time.time()
    fold = build_fold(held)
    fold["scored"] = np.asarray(fold["scored"])
    wcols = np.asarray(fold["wcols"])
    truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
    pre = np.asarray(fold["pre"], dtype=np.float64)
    targets = list(fold["targets"])
    op = Operator(PROD_OPS / "whitener_union_pruned.npz" if FOLDS[held]["op_kind"] == "production"
                  else Path(FOLDS[held]["op"]))
    cd4 = cd4_for(targets)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    whitened = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    f_cd4, _ = solve_fraction(whitened, cd4, CD4_F)
    f_u, _ = solve_fraction_carrying(whitened, h1, UNION_G)
    base, _ = E062.add_post_whiten_residual(whitened, merged(cd4, f_cd4, h1, f_u), 1.0, wcols)
    base[:, fold["panel_positions"]] = 0.0
    key, carried, _ = shared_energy_key(targets, fold["sources"], op, wcols, cols=None)
    K3 = key_scale(base, weights(key, carried, 0.3, 3.0)[0])[0]
    Ki = key_scale(base, weights(key, carried, 0.3, 1e9)[0])[0]
    K5 = key_scale(base, weights(key, carried, 0.5, 1e9)[0])[0]
    arms = {"base": (base, None), "K3": (K3, None), "Kinf": (Ki, None), "K0.5_inf": (K5, None),
            "Kinf+sh@1.15": (Ki, 1.15), "Kinf+sh@1.25": (Ki, 1.25), "Kinf+sh@1.5": (Ki, 1.5),
            "K0.5_inf+sh@1.25": (K5, 1.25), "K3+sh@1.15": (K3, 1.15)}
    ref = json.loads((WORK / f"arms414-{held}.json").read_text())["arms"]
    ref.update({k: v for k, v in json.loads((WORK / f"extra414-{held}.json").read_text())["arms"].items() if k not in ref})
    out = dict(fold=held, arms={})
    P_DE = {ctx: [profile_from_effect(base[i], b, gate, "minnorm") for i in range(base.shape[0])]
            for ctx, (b, gate) in frames().items()}
    for name, (E, q) in arms.items():
        per_des, per_real, diag = [], [], dict(r_min=[], K=[], clip_mass=[], clip_genes=[])
        for ctx, (b, gate) in frames().items():
            F = E if q is None else sharpen(E, q, gate, mode="row")
            d_des, d_real = np.empty_like(F), np.empty_like(F)
            for i in range(F.shape[0]):
                p_de = P_DE[ctx][i]
                p_bk = profile_from_effect(F[i], b, gate, "minnorm")
                p_re, r, k, cm, cg = realise(p_de, p_bk)
                d_des[i] = np.log1p(TS * p_bk) - b
                d_real[i] = np.log1p(TS * p_re) - b
                diag["r_min"].append(r); diag["K"].append(k); diag["clip_mass"].append(cm); diag["clip_genes"].append(cg)
            per_des.append(endpoints(d_des[:, fold["scored"]], truth))
            per_real.append(endpoints(d_real[:, fold["scored"]], truth))
        s = dict(pds_design=float(np.mean([e["pds"] for e in per_des])),
                 pds_realised=float(np.mean([e["pds"] for e in per_real])),
                 cos_agg_design=float(np.mean([np.sqrt(max(0, 1 - e["sse"])) for e in per_des])),
                 cos_agg_realised=float(np.mean([np.sqrt(max(0, 1 - e["sse"])) for e in per_real])),
                 r_min_q=np.quantile(diag["r_min"], [0.0, 0.1, 0.5, 0.9]).tolist(),
                 frac_K_capped=float(np.mean(np.asarray(diag["K"]) >= KMAX)),
                 K_median=float(np.median(diag["K"])),
                 clip_mass_q=np.quantile(diag["clip_mass"], [0.5, 0.9, 1.0]).tolist(),
                 clip_genes_median=float(np.median(diag["clip_genes"])))
        d = abs(s["pds_design"] - ref[name]["pds"])
        assert d < 1e-12, f"KA design == arms414 FAILED {held} {name} {d}"
        out["arms"][name] = s
        print(f"[{held}] {name:18s} pds design {s['pds_design']:.5f} realised {s['pds_realised']:.5f}  cos design "
              f"{s['cos_agg_design']:.6f} realised {s['cos_agg_realised']:.6f}  r_min q0/10/50 "
              f"{np.round(s['r_min_q'][:3], 3).tolist()}  K capped {100 * s['frac_K_capped']:.0f}%  clip mass q50/90 "
              f"{s['clip_mass_q'][0]:.2e}/{s['clip_mass_q'][1]:.2e} ({time.time() - t0:.0f}s)", flush=True)
    (WORK / f"ballast414-{held}.json").write_text(json.dumps(out, indent=1))
    del op


def main():
    for held in sys.argv[1:] or ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2"):
        if not (WORK / f"ballast414-{held}.json").exists():
            run_fold(held)


if __name__ == "__main__":
    main()
