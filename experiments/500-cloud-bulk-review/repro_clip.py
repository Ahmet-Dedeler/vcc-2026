"""500: is the floor-scale cosine loss sampling noise, and does a clip-aware L2 restore fix read-back?

    uv run python experiments/500-cloud-bulk-review/repro_clip.py
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np

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

TS = 50_000.0
G = 4096
N = 400
N_PANEL = 64
Q = 1.75


def cosine(a, b):
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    return float(a @ b / (na * nb)) if na and nb else float("nan")


def effect(profile, b, mask):
    return (np.log1p(TS * profile) - b)[mask]


def clip_aware_sharpen(gated, q, gate, protect, b, n0):
    """Sharpen, then keep every non-panel coordinate inside the renderer's open clip and put the
    L2 budget only on coordinates that survive it. Rescaling never pushes a coordinate back through
    the floor: the scale is capped by the most-negative realizable coordinate."""
    h = M.bulk_sharpen(gated, q, gate, protect)
    on = gate & ~protect
    h = np.where(on, h, gated)
    for _ in range(8):
        floor = np.where(on, -b, -np.inf)
        h = np.maximum(h, floor)
        room = np.ones(h.shape, dtype=np.float64)
        neg = on & (h < -1e-12)
        # scaling h by alpha keeps h >= -b on negative coords iff alpha <= (-b)/|h| = b/(-h) wait:
        # h is negative, floor is -b (b>0). h >= -b. alpha>1 makes h more negative. The max alpha
        # on a negative coordinate is (-b) / h  if h<0? h is negative so (-b)/h is negative. Use
        # alpha <= (-b) / h when h < 0, i.e. alpha <= b / (-h).
        with np.errstate(divide="ignore", invalid="ignore"):
            room = np.where(neg, (-floor) / np.maximum(-h, 1e-300), np.inf)
        cap = float(np.min(room[on])) if on.any() else 1.0
        cur = float(np.linalg.norm(h[on]))
        if cur <= 0:
            break
        alpha = min(n0 / cur, cap if np.isfinite(cap) else n0 / cur)
        if abs(alpha - 1.0) < 1e-9:
            break
        h = np.where(on, h * alpha, h)
    return h


def main():
    rng = np.random.default_rng(500)
    raw = rng.lognormal(-1.2, 1.8, size=G)
    p = raw / raw.sum()
    b = np.log1p(TS * p)
    panel = np.zeros(G, dtype=bool)
    panel[:N_PANEL] = True
    libs = np.clip(rng.lognormal(np.log(20109), 0.35, size=N), 3275, 52420).astype(np.int64)
    mask = ~panel
    noise_rows = []
    clip_rows = []
    for i in range(8):
        h = rng.normal(0.0, 0.12, size=G) + rng.standard_t(3, size=G) * 0.04
        h[panel] *= 0.25
        h[i % N_PANEL] = -1.6
        de = M.profile_from_effect(h, b)
        # deterministic carriers: largest remainder, so the only error is integer rounding
        for s, tag in ((0.305, "floor"), (0.609, "median"), (1.827, "cap")):
            bulk = M.profile_from_effect(h * s, b, sharpen_q=Q, protect=panel)
            k, rmin = M.ballast_count(de, bulk, libs, N)
            det = np.zeros((N - k, G), dtype=np.int64)
            for c in range(N - k):
                raw_c = de * int(libs[c])
                det[c] = np.floor(raw_c)
                rem = int(libs[c] - det[c].sum())
                if rem > 0:
                    det[c, np.argsort(-(raw_c - det[c]), kind="stable")[:rem]] += 1
            full, diag = M.ballast_rows(det, de, bulk, k, seed=3)
            got = full.sum(0).astype(np.float64)
            got /= got.sum()
            eb, ea, ed = effect(bulk, b, mask), effect(got, b, mask), effect(de, b, mask)
            noise_rows.append({
                "tag": tag, "K": int(k), "rmin": float(rmin),
                "cos_vs_bulk": cosine(ea, eb),
                "cos_vs_de": cosine(ea, ed),
                "effect_l2": float(np.linalg.norm(eb)),
                "resid_l2": float(np.linalg.norm(ea - eb)),
                "clip_mass": diag["ballast_clip_mass"],
            })

        # clip-aware read-back at the cap and at cap*1.4
        gate = (1e6 * np.expm1(b) / max(float(np.expm1(b).sum()), 1e-12) > M.MIN_CPM) & (b > 0)
        for s in (1.827, 2.558, 3.0):
            gated = np.where(gate, h * s, 0.0)
            on = gate & ~panel
            n0 = float(np.linalg.norm(gated[on]))
            sharp = M.bulk_sharpen(gated, Q, gate, panel)
            aware = clip_aware_sharpen(gated, Q, gate, panel, b, n0)

            def pack(hh):
                prof = np.clip(np.expm1(b + hh), 0.0, None) / TS
                prof = prof / prof.sum()
                got = np.log1p(TS * prof) - b
                return {
                    "cos": cosine(got[on], hh[on]),
                    "cos_vs_presharp": cosine(got[on], gated[on]),
                    "l2_ratio_vs_presharp": float(np.linalg.norm(got[on]) / max(n0, 1e-300)),
                    "clipped": int(np.sum(on & (b + hh <= 1e-9))),
                }

            clip_rows.append({"scale": s, "stock": pack(sharp), "clip_aware": pack(aware)})

    def summ(tag):
        xs = [r for r in noise_rows if r["tag"] == tag]
        return {
            "cos_bulk_min": min(r["cos_vs_bulk"] for r in xs),
            "cos_bulk_median": float(np.median([r["cos_vs_bulk"] for r in xs])),
            "resid_over_effect_median": float(np.median([r["resid_l2"] / r["effect_l2"] for r in xs])),
            "effect_l2_median": float(np.median([r["effect_l2"] for r in xs])),
            "clip_mass_max": max(r["clip_mass"] for r in xs),
        }

    by_scale = {}
    for s in (1.827, 2.558, 3.0):
        xs = [r for r in clip_rows if r["scale"] == s]
        by_scale[str(s)] = {
            "stock_cos_min": min(r["stock"]["cos"] for r in xs),
            "stock_cos_median": float(np.median([r["stock"]["cos"] for r in xs])),
            "stock_cos_vs_presharp_median": float(np.median([r["stock"]["cos_vs_presharp"] for r in xs])),
            "aware_cos_min": min(r["clip_aware"]["cos"] for r in xs),
            "aware_cos_median": float(np.median([r["clip_aware"]["cos"] for r in xs])),
            "aware_cos_vs_presharp_median": float(np.median([r["clip_aware"]["cos_vs_presharp"] for r in xs])),
            "stock_l2_median": float(np.median([r["stock"]["l2_ratio_vs_presharp"] for r in xs])),
            "aware_l2_median": float(np.median([r["clip_aware"]["l2_ratio_vs_presharp"] for r in xs])),
        }
    out = {"deterministic_ballast": {t: summ(t) for t in ("floor", "median", "cap")}, "clip_aware": by_scale}
    print(json.dumps(out, indent=1), flush=True)
    dest = Path(__file__).resolve().parent / "work" / "repro_clip.json"
    dest.write_text(json.dumps({"summary": out, "noise_rows": noise_rows, "clip_rows": clip_rows}, indent=1))
    print(f"wrote {dest}", flush=True)


if __name__ == "__main__":
    main()
