#!/usr/bin/env python3
"""625c — can a module flip be READ through final-entry noise? Signal size per module in the delivery-faithful fold model.

Ahmet approved decoding D/E/F's own module programme through final entries (637 §4, 2026-10-02). Feasibility: per-module
flip effects x_M on nmae (the primary decode member) against the noise each entry carries. Model: 480/612b's
delivery-faithful `small_base+mde` shape, P = BETA·cons + s·add, with the October hook = key-gated identity hybrid
(weak-key half: sign(S), S = z(P_C4 - .5) + 0.5 z(I), Flex refs; strong half: key sign). Flipping module M = P_C4 -> 1 - P_C4
on M's genes (625's fixed gene sets, mapped by symbol), recompute S and s. x_M = mean-target nmae(P_M) - nmae(P_base), raw.
Designs: the preregistered 16 modules, and 8 CLEAN4 octiles (module // 2). Noise per decoded module: an entry carrying
σ_e (overall units) gives σ_x = 2σ_e·sqrt(n)/n for an n-row Hadamard (y = S - 2Δ; x = Hᵀy/n): 16 rows -> σ_e/2, 8 -> σ_e/√2.
Fused (637: a bulk switch per entry) σ_e ≈ 0.0002 overall; pure prior-only entries σ_e = 0 (shared draws).
Overall-equivalent = Δraw nmae × 1.648/6 (192's INV) × board multiplier m (1, 1.5). Run with cwd 480's dir. < 2 GB.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "experiments/480-nmae-instrument")); sys.path.insert(0, str(ROOT / "experiments/537-arm-pipeline"))
sys.path.insert(0, str(ROOT / "experiments/622-identity-track"))
import fold480 as F  # noqa: E402
from common480 import load_donor_rows, nmae_rows, load_h1  # noqa: E402
import foldsign537 as FS  # noqa: E402
import hyb622 as Hy  # noqa: E402

C4ALL = ["k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2"]
C4B = {"h1": C4ALL, "orion_hct116": C4ALL, "orion_hek293t": C4ALL, "k562_gwps": ["nadig_jurkat", "nadig_hepg2"]}
MODS = np.load(HERE / "work/modules625.npz", allow_pickle=True)
MOD = dict(zip(MODS["genes"].astype(str), MODS["module"].astype(int)))
h = load_h1(); G = np.isfinite(h["T"])
L0 = np.nan_to_num(np.asarray(np.load(ROOT / "experiments/426-mde-magnitude/work/h1/stats-base_s0_bps.npz")["lfc"], float))
c3 = np.nan_to_num(np.nanmean(np.stack([load_donor_rows(d, h["targets"], h["genes"], pc_cpm=1.0)[0] for d in F.DON]), 0))
BETA = float(np.median(np.abs(L0[G])) / np.median(np.abs(c3[G])))
prof = Hy.profiles()
INV = 1.648 / 6
SIG_E = 0.0002
out = {"BETA": BETA}
for name in C4B:
    S = F.load_surface(name)
    T, cnt, genes, targets = S["T"], S["cnt"], [str(g) for g in S["genes"]], S["targets"]
    cons0 = np.nan_to_num(np.nanmean(np.stack([load_donor_rows(d, targets, np.array(genes), pc_cpm=1.0)[0] for d in S["donors"]]), 0))
    P4 = FS.clean4_P(np.array(genes), C4B[name])
    sc = np.sign(cons0); add = F.mde_add(cnt)[None, :]
    weak = np.zeros(T.shape, bool)
    for i in range(T.shape[0]):
        a = np.abs(cons0[i]); nz = a > 0
        weak[i] = (~nz) | (a < (np.quantile(a[nz], 0.5) if nz.any() else 0.0))
    base = BETA * cons0
    refs = [k for k in ("ctxA", "ctxB", "ctxC", "h1") if k != name]
    zi = Hy.identity(prof, name, refs, genes)
    sd = np.nanstd(P4 - 0.5)
    mod = np.array([MOD.get(g, -1) for g in genes])

    def pred(P4x):
        zc = np.where(np.isfinite(P4x), (P4x - 0.5) / sd, np.nan)
        Sv = np.where(np.isfinite(zc), zc + 0.5 * np.nan_to_num(zi), np.nan_to_num(zi))
        sI = np.where(np.isfinite(Sv) & (Sv != 0), np.where(Sv > 0, -1.0, 1.0), 0.0)
        return base + np.where(weak & (sI[None, :] != 0), sI[None, :], sc) * add
    n0 = float(np.nanmean(nmae_rows(pred(P4), T)))
    res = {"nmae_base": n0, "covered_genes": int((mod >= 0).sum())}
    for design, K, mfun in (("16mod", 16, lambda m: m), ("8oct", 8, lambda m: np.where(m >= 0, m // 2, -1))):
        mm = mfun(mod)
        x = np.zeros(K)
        for k in range(K):
            P4x = np.where(mm == k, 1.0 - P4, P4)
            x[k] = float(np.nanmean(nmae_rows(pred(P4x), T))) - n0
        sig_x = 2 * SIG_E * np.sqrt(K) / K
        ov = np.abs(x) * INV
        res[design] = dict(x_raw=x.tolist(), median_abs_overall_m1=float(np.median(ov)), max_abs_overall_m1=float(ov.max()),
                           sigma_x_fused_overall=float(sig_x),
                           n_detectable_fused_m1=int((ov > 2 * sig_x).sum()), n_detectable_fused_m15=int((1.5 * ov > 2 * sig_x).sum()),
                           n_helpful_flips=int((x < 0).sum()), sum_helpful_overall_m1=float(-x[x < 0].sum() * INV))
        r = res[design]
        print(f"{name:14s} {design:5s} | |x_M| overall (m=1): median {r['median_abs_overall_m1']:.6f}, max {r['max_abs_overall_m1']:.6f} | "
              f"fused σ_x {sig_x:.6f}: detectable {r['n_detectable_fused_m1']}/{K} (m 1.5: {r['n_detectable_fused_m15']}/{K}) | "
              f"helpful flips {r['n_helpful_flips']}, their total +{r['sum_helpful_overall_m1']:.5f} overall", flush=True)
    out[name] = res
(HERE / "work/snr625.json").write_text(json.dumps(out, indent=1))
