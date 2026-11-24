"""500 follow-up: production-scale read-back, ballast dust-K, post-boost jackknife.

    uv run python experiments/500-cloud-bulk-review/repro_followup.py
"""
from __future__ import annotations

import importlib.util
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]
for k in list(os.environ):
    if k.startswith("VCC_"):
        os.environ.pop(k)
os.environ["VCC_PERT_SCALE"] = "lodo415_lo0.5_hi3"
os.environ["VCC_PERT_SCALE_BALLAST"] = "1"

name = "make_submission_pm_cd4"
spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / "make_submission_pm_cd4.py")
M = importlib.util.module_from_spec(spec)
sys.modules[name] = M
spec.loader.exec_module(M)

from cell_eval2.moments import jackknife_correction  # noqa: E402

TS = 50_000.0
G = 4096
N = 400
N_PERT = 12
N_PANEL = 64
Q = 1.75
SEED = 500
# final scales after lodo415_lo0.5_hi3's energy constant (repro_bulk.py): floor, median, cap,
# and cap times the ctx-scale bound 1.4.
SCALES = (0.305, 0.609, 1.000, 1.272, 1.827, 2.558)


def cosine(a, b):
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    return float(a @ b / (na * nb)) if na and nb else float("nan")


def effect(profile, b, mask):
    return (np.log1p(TS * profile) - b)[mask]


def main():
    rng = np.random.default_rng(SEED)
    raw = rng.lognormal(-1.2, 1.8, size=G)
    p = raw / raw.sum()
    b = np.log1p(TS * p)
    panel = np.zeros(G, dtype=bool)
    panel[:N_PANEL] = True
    libs = np.clip(rng.lognormal(np.log(20109), 0.35, size=N), 3275, 52420).astype(np.int64)
    out = {"scales": SCALES}

    # ----- read-back at the production final-scale grid
    rb = {str(s): [] for s in SCALES}
    for i in range(N_PERT):
        h = rng.normal(0.0, 0.12, size=G) + rng.standard_t(3, size=G) * 0.04
        h[panel] *= 0.25
        h[i % N_PANEL] = -1.6
        gate = (1e6 * np.expm1(b) / np.expm1(b).sum() > M.MIN_CPM) & (b > 0)
        for s in SCALES:
            gated = np.where(gate, h * s, 0.0)
            sharp = M.bulk_sharpen(gated, Q, gate, panel)
            prof = np.clip(np.expm1(b + sharp), 0.0, None) / TS
            prof = prof / prof.sum()
            got = np.log1p(TS * prof) - b
            on = gate & ~panel
            fixed = M.mass_consistent_effect(sharp, b, gate, mode="minnorm")
            prof_f = np.clip(np.expm1(b + fixed), 0.0, None) / TS
            prof_f = prof_f / prof_f.sum()
            got_f = np.log1p(TS * prof_f) - b
            rb[str(s)].append({
                "cos": cosine(got[on], sharp[on]),
                "cos_fix": cosine(got_f[on], fixed[on]),
                "clipped": int(np.sum(on & (b + sharp <= 0))),
                "l2_ratio": float(np.linalg.norm(got[on]) / max(np.linalg.norm(sharp[on]), 1e-300)),
                "mass_ratio": float(np.clip(np.expm1(b + sharp), 0, None).sum() / TS),
            })
    rb_sum = {}
    for s, rows in rb.items():
        rb_sum[s] = {
            "cos_min": min(r["cos"] for r in rows),
            "cos_p10": float(np.quantile([r["cos"] for r in rows], 0.1)),
            "cos_median": float(np.median([r["cos"] for r in rows])),
            "cos_fix_min": min(r["cos_fix"] for r in rows),
            "cos_fix_median": float(np.median([r["cos_fix"] for r in rows])),
            "clipped_median": float(np.median([r["clipped"] for r in rows])),
            "l2_ratio_median": float(np.median([r["l2_ratio"] for r in rows])),
            "mass_ratio_median": float(np.median([r["mass_ratio"] for r in rows])),
        }
    out["readback_by_final_scale"] = rb_sum
    print("readback", json.dumps(rb_sum), flush=True)

    # ----- ballast: dust gene that sets K, deterministic cosine, Chat after boost
    rows = []
    for i in range(6):
        h = rng.normal(0.0, 0.12, size=G) + rng.standard_t(3, size=G) * 0.04
        h[panel] *= 0.25
        h[i] = -1.6
        de = M.profile_from_effect(h, b)
        # production final scales: floor and cap
        for s, tag in ((0.305, "floor"), (0.609, "median"), (1.827, "cap")):
            bulk = M.profile_from_effect(h * s, b, sharpen_q=Q, protect=panel)
            live = de > 1e-9
            ratio = np.where(live, bulk / np.maximum(de, 1e-300), np.inf)
            j = int(np.argmin(ratio))
            k, rmin = M.ballast_count(de, bulk, libs, N)
            # thresholds that ignore dust
            def k_at(thr):
                m = de > thr
                if not m.any():
                    return None
                rr = float(np.min(bulk[m] / de[m]))
                sc = float(libs[: N - M.BALLAST_KMAX].sum())
                need = sc * (1.0 / max(rr, 1e-6) - 1.0)
                return int(min(M.BALLAST_KMAX, max(1, math.ceil(need / M.BALLAST_LB)))), rr

            k_cpm5, r_cpm5 = k_at(5e-6)  # CPM 5 at a 1e6 library == fraction 5e-6
            draws = np.vstack([rng.multinomial(int(libs[c]), de) for c in range(N)])
            carriers = draws[: N - k]
            full, diag = M.ballast_rows(carriers, de, bulk, k, seed=7 + i)
            tot = full.sum(0).astype(np.float64)
            got = tot / tot.sum()
            mask = ~panel
            # deterministic allocation cosine: replay ballast_rows' target before the multinomial,
            # using the delivered carrier totals (the noise the ratio mode is designed to keep)
            ce = float(jackknife_correction(sparse.csr_matrix(full.astype(np.float64)),
                                            np.zeros(full.shape[0], np.intp), 1, TS)[0])
            # 400 equal-depth cells from the INTENDED bulk profile, control-like libraries
            ref = np.vstack([rng.multinomial(int(libs[c]), bulk) for c in range(N)])
            ce_ref = float(jackknife_correction(sparse.csr_matrix(ref.astype(np.float64)),
                                                np.zeros(N, np.intp), 1, TS)[0])
            rec = {
                "tag": tag, "scale": s, "K": k, "rmin": rmin,
                "argmin_de": float(de[j]), "argmin_bulk": float(bulk[j]),
                "argmin_cpm": float(1e6 * de[j]),
                "K_if_cpm5": k_cpm5, "rmin_if_cpm5": r_cpm5,
                "clip_genes": diag["ballast_clip_genes"], "clip_mass": diag["ballast_clip_mass"],
                "trim": diag["ballast_trim"], "jk_ratio_builder": diag["ballast_jk_ratio"],
                "cos_vs_bulk": cosine(effect(got, b, mask), effect(bulk, b, mask)),
                "cos_vs_de": cosine(effect(got, b, mask), effect(de, b, mask)),
                "Chat": ce, "C400_bulk": ce_ref, "Chat_over_C400": ce / ce_ref,
            }
            rows.append(rec)
            print(f"ballast {tag} i={i} K={k} (cpm5 would be {k_cpm5}) rmin={rmin:.4g} "
                  f"argmin_cpm={rec['argmin_cpm']:.4g} clip={diag['ballast_clip_mass']} "
                  f"cos={rec['cos_vs_bulk']:.4f} Chat/C400={rec['Chat_over_C400']:.3f}", flush=True)
    out["ballast_scales"] = rows

    # same carriers, K in {1, 4, 16}: how much does an inflated K rotate the bulk?
    h = rng.normal(0.0, 0.12, size=G)
    h[panel] *= 0.25
    h[0] = -1.6
    de = M.profile_from_effect(h, b)
    bulk = M.profile_from_effect(h * 1.827, b, sharpen_q=Q, protect=panel)
    draws = np.vstack([rng.multinomial(int(libs[c]), de) for c in range(N - 16)])
    sweep = []
    mask = ~panel
    for k in (1, 2, 4, 8, 16):
        full, diag = M.ballast_rows(draws, de, bulk, k, seed=1)
        got = full.sum(0).astype(np.float64)
        got /= got.sum()
        sweep.append({
            "K": k, "clip_mass": diag["ballast_clip_mass"], "clip_genes": diag["ballast_clip_genes"],
            "trim": diag["ballast_trim"],
            "cos_vs_bulk": cosine(effect(got, b, mask), effect(bulk, b, mask)),
            "cos_vs_de": cosine(effect(got, b, mask), effect(de, b, mask)),
            "ballast_fraction": (k * M.BALLAST_LB) / float(full.sum()),
        })
    out["k_sweep_at_cap_fixed_carriers"] = sweep
    print("k_sweep", json.dumps(sweep), flush=True)

    dest = Path(__file__).resolve().parent / "work" / "repro_followup.json"
    dest.write_text(json.dumps(out, indent=1))
    print(f"wrote {dest}", flush=True)


if __name__ == "__main__":
    main()
