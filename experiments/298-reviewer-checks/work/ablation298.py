"""298 item 2 — the FIT of hallmark identity alone, against random matched gene sets.

297 row 1c reads 290's control C3 as an ablation: it reports the DROP when the nine hallmark
projections are replaced by random matched sets inside the FULL 18-column model (+0.3174 ->
+0.2943, i.e. random sets recover 93 %).  That is not the number Wang et al. 2026's Hallmark
contrast disagrees with.  This refits three restricted models against the SAME target
(`remainder_half2`), the SAME component folds and the SAME component block bootstrap:

  (i)   HALLMARK-ONLY   the nine `pe_<HALLMARK>` projections + `pe_prolif_abs`      (10 columns)
  (ii)  RANDOM-MATCHED  the same 10 columns built from random gene sets of matched size, 20 draws
                        (ctrl290.py's construction, same seed, same draw sequence)
  (iii) SPECTRAL-ONLY   gen_pc5, gen_pc20, gen_cos, gen_abscos, breadth, logmag     (6 columns)

and reports (i) - (ii) with its interval.

Code COPIED, not imported: `lib290.py` is loaded by PATH (as 290's own scripts do) so the ridge,
the fold construction and the bootstrap are literally 290's objects; the bootstrap driver below is
`lib290.block_boot` (lib284.py:167-236) with the per-replicate rho kept so a paired difference
against a SET of 20 draws can be formed.

    PYTHONDONTWRITEBYTECODE=1 uv run python ablation298.py

Writes work/ablation298.json, work/ablation298.log.
"""
from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("~/Code/virtual-cell-challenge-2026")
E290 = ROOT / "experiments/290-remainder-response-features"
WORK = ROOT / "experiments/298-reviewer-checks/work"
sys.path.insert(0, str(E290))
import lib290 as L  # noqa: E402

T0 = time.time()
LOG: list[str] = []
N_DRAW = 20

HALL = [f"{c}" for c in L.PE_COLS] + ["pe_prolif_abs"]          # (i)  10 columns
SPEC = ["gen_pc5", "gen_pc20", "gen_cos", "gen_abscos", "breadth", "logmag"]   # (iii) 6 columns


def P(msg: str = "") -> None:
    line = f"[{time.time() - T0:7.1f}s] {msg}" if msg else ""
    print(line, flush=True)
    LOG.append(line)


# ---- lib290.block_boot (lib284.py:167-236), per-replicate rho kept ------------------------------
def block_boot_rhos(y, preds: dict, groups, rng, n: int = L.N_BOOT) -> tuple[dict, np.ndarray, list]:
    groups = np.asarray(groups)
    uniq = np.array(sorted(set(groups.tolist())))
    idx_of = {g: np.flatnonzero(groups == g) for g in uniq}
    keys = list(preds)
    R = np.full((n, len(keys)), np.nan)
    for b in range(n):
        pick = rng.integers(0, len(uniq), size=len(uniq))
        sel = np.concatenate([idx_of[uniq[p]] for p in pick])
        yy = y[sel]
        for j, k in enumerate(keys):
            R[b, j] = L.spearman(yy, preds[k][sel])
    point = {k: float(L.spearman(y, preds[k])) for k in keys}
    out = {k: dict(point=point[k], lo=float(np.nanpercentile(R[:, j], 2.5)),
                   hi=float(np.nanpercentile(R[:, j], 97.5)), sd=float(np.nanstd(R[:, j])))
           for j, k in enumerate(keys)}
    return out, R, keys


def iv(v: np.ndarray, point: float) -> dict:
    return dict(point=float(point), lo=float(np.nanpercentile(v, 2.5)),
                hi=float(np.nanpercentile(v, 97.5)), sd=float(np.nanstd(v)),
                p_le_0=float(np.mean(v <= 0)))


def main() -> None:
    h = L.check_prereg()
    P(f"290 prereg {h} MATCH (lib290 loaded by path; ridge/folds/bootstrap are 290's objects)")

    # ---------------- the unit, target and folds, rebuilt exactly as ctrl290.py:56-67 ------------
    held = L.load_held()
    X21 = held[L.FEATS].to_numpy(float)
    comp = np.load(L.SPLITS284, allow_pickle=True)["comp"]
    fold = L.component_folds(comp, 5)
    ref = np.load(L.REM284, allow_pickle=True)
    ok = np.asarray(ref["ok"])
    y = np.asarray(ref["r2"])
    u = np.flatnonzero(ok)
    yu, foldu, compu = y[u], fold[u], comp[u]
    feat = pd.read_csv(L.WORK / "feat290.csv")
    assert list(feat["gene"].astype(str)) == list(held["gene"].astype(str))

    pred = pd.read_csv(L.WORK / "pred290.csv")
    assert list(pred["gene"].astype(str)) == list(held["gene"].astype(str)[u])
    d_y = float(np.nanmax(np.abs(pred["remainder_half2"].to_numpy(float) - yu)))
    d_f = int(np.abs(pred["fold"].to_numpy() - foldu).max())
    d_c = int(np.abs(pred["component"].to_numpy() - compu).max())
    P(f"unit {len(u)} genes, {len(set(compu.tolist()))} components, "
      f"fold sizes {np.bincount(foldu, minlength=5).tolist()}")
    P(f"  pred290.csv agreement: target {d_y:.1e} (CSV round-trip), fold {d_f}, component {d_c}")
    assert d_y < 1e-12 and d_f == 0 and d_c == 0

    def Xof(cols: list[str]) -> np.ndarray:
        return feat.loc[u, [c if c == "on_axis" else f"{c}_h1" for c in cols]].to_numpy(float)

    # ================= KA: 290's primary cell, refitted =========================================
    Xfull = Xof(L.FEATS290)
    p_full, _ = L.oof_ridge(Xfull, yu, foldu)
    rho_full = L.spearman(yu, p_full)
    P(f"\nKA  290 primary cell (full {Xfull.shape[1]}-column model, half-1 features -> "
      f"remainder_half2): {rho_full:+.6f}  (290 reported +0.317385)  "
      f"delta {abs(rho_full - 0.31738512466156):.2e}")
    assert abs(rho_full - 0.31738512466156) < 1e-9, "KA FAILED — abort (AGENTS.md:1375-1382)"
    p_ann, _ = L.oof_ridge(X21[u], yu, foldu)
    P(f"KA  annotation-21 ridge: {L.spearman(yu, p_ann):+.6f}  (290 reported -0.056642)")

    # ================= (i) HALLMARK-ONLY and (iii) SPECTRAL-ONLY ================================
    Xh = Xof(HALL)
    p_hall, i_hall = L.oof_ridge(Xh, yu, foldu)
    Xs = Xof(SPEC)
    p_spec, i_spec = L.oof_ridge(Xs, yu, foldu)
    P(f"\n(i)   HALLMARK-ONLY  ({Xh.shape[1]} cols)  rho {L.spearman(yu, p_hall):+.4f}")
    P(f"(iii) SPECTRAL-ONLY  ({Xs.shape[1]} cols)  rho {L.spearman(yu, p_spec):+.4f}")

    # ================= (ii) RANDOM MATCHED SETS — ctrl290.py:36-99 ==============================
    per_src, genes_held, keys, axis = L.scan_sources()
    lines = [per_src[k]["line"] for k in keys]
    nG, nA = len(genes_held), len(axis)
    h1 = [i for i, ln in enumerate(lines) if ln in L.HALF1_LINES]
    _, _, _, psizes = L.programme_directions(axis)
    sizes = [psizes[p]["n_on_axis"] for p in L.PROGRAMMES]
    prolif_size = len({g for p in L.PROLIF for g in L.hallmark_sets()[p]} & set(axis))
    P(f"\nhalf-1 sources {len(h1)}; axis {nA}; matched sizes {sizes}; proliferation {prolif_size}")
    assert sizes == [104, 102, 96, 44, 43, 49, 101, 79, 38] and prolif_size == 228, \
        "matched sizes differ from ctrl290 — abort"

    idxs, Us = {}, {}
    for i in h1:
        k = keys[i]
        idx, U, _, _, _ = L.load_unit(per_src[k]["rel"], axis, genes_held)
        idxs[k], Us[k] = idx, U
    P(f"loaded {len(h1)} half-1 banks; peak RSS "
      f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30:.2f} GiB")

    rng = np.random.default_rng(L.SEED + 91)          # ctrl290.py:73 — the SAME draw sequence
    draws, p_rand = [], {}
    for d in range(N_DRAW):
        Wr = np.zeros((len(L.PROGRAMMES), nA))
        for j, s in enumerate(sizes):
            Wr[j, rng.choice(nA, size=s, replace=False)] = 1.0
            Wr[j] /= np.linalg.norm(Wr[j])
        wp = np.zeros(nA)
        wp[rng.choice(nA, size=prolif_size, replace=False)] = 1.0
        wp /= np.linalg.norm(wp)

        PE = np.full((len(L.PROGRAMMES), len(h1), nG), np.nan, np.float32)
        PP = np.full((len(h1), nG), np.nan, np.float32)
        for q, i in enumerate(h1):
            k = keys[i]
            proj = Us[k] @ Wr.T.astype(np.float32)
            for j in range(len(L.PROGRAMMES)):
                PE[j, q, idxs[k]] = proj[:, j]
            PP[q, idxs[k]] = Us[k] @ wp.astype(np.float32)
        with np.errstate(invalid="ignore"):
            pe = np.stack([np.nanmean(PE[j], axis=0) for j in range(len(L.PROGRAMMES))], axis=1)
            pab = np.nanmean(np.abs(PP), axis=0)
        Xd = np.column_stack([pe[u], pab[u][:, None]])          # (ii) the SAME 10 columns as (i)
        pr, _ = L.oof_ridge(Xd, yu, foldu)
        p_rand[f"rand{d:02d}"] = pr
        draws.append(dict(draw=d, rho=float(L.spearman(yu, pr))))
        P(f"  random draw {d:2d}: hallmark-only rho {draws[-1]['rho']:+.4f}")
    del Us, idxs

    r_draw = np.array([x["rho"] for x in draws])
    P(f"\n(ii)  RANDOM-MATCHED hallmark-only: {r_draw.mean():+.4f} +- {r_draw.std():.4f} "
      f"[{r_draw.min():+.4f}, {r_draw.max():+.4f}]; "
      f"{int((r_draw >= L.spearman(yu, p_hall)).sum())}/{N_DRAW} reach (i)")

    # ================= the component block bootstrap =============================================
    preds = dict(hallmark_only=p_hall, spectral_only=p_spec, full=p_full, annotation=p_ann,
                 **p_rand)
    ci, R, keys_b = block_boot_rhos(yu, preds, compu, np.random.default_rng(L.SEED + 11))
    kidx = {k: j for j, k in enumerate(keys_b)}
    rand_cols = [kidx[f"rand{d:02d}"] for d in range(N_DRAW)]

    # (i) - (ii): paired inside each bootstrap replicate, (ii) = mean over the 20 draws
    d_mean = R[:, kidx["hallmark_only"]] - np.nanmean(R[:, rand_cols], axis=1)
    pt_mean = ci["hallmark_only"]["point"] - float(r_draw.mean())
    delta_mean = iv(d_mean, pt_mean)
    # the conservative reading: (i) minus the BEST of the 20 draws, same replicate
    d_best = R[:, kidx["hallmark_only"]] - np.nanmax(R[:, rand_cols], axis=1)
    delta_best = iv(d_best, ci["hallmark_only"]["point"] - float(r_draw.max()))
    # spectral - hallmark, and full - spectral (what carries the signal)
    d_sh = R[:, kidx["spectral_only"]] - R[:, kidx["hallmark_only"]]
    delta_spec_hall = iv(d_sh, ci["spectral_only"]["point"] - ci["hallmark_only"]["point"])
    d_fs = R[:, kidx["full"]] - R[:, kidx["spectral_only"]]
    delta_full_spec = iv(d_fs, ci["full"]["point"] - ci["spectral_only"]["point"])

    P("\ncomponent block bootstrap (2000, 437 components):")
    for k in ("hallmark_only", "spectral_only", "full", "annotation"):
        v = ci[k]
        P(f"  {k:16s} {v['point']:+.4f} [{v['lo']:+.4f}, {v['hi']:+.4f}]")
    P(f"  {'random-matched':16s} {r_draw.mean():+.4f} (mean of 20 draws; per-draw sd "
      f"{r_draw.std():.4f})")
    P(f"\n  (i) - (ii) mean of draws  {delta_mean['point']:+.4f} "
      f"[{delta_mean['lo']:+.4f}, {delta_mean['hi']:+.4f}]  P(<=0) {delta_mean['p_le_0']:.3f}")
    P(f"  (i) - (ii) BEST draw      {delta_best['point']:+.4f} "
      f"[{delta_best['lo']:+.4f}, {delta_best['hi']:+.4f}]  P(<=0) {delta_best['p_le_0']:.3f}")
    P(f"  (iii) - (i) spectral-hall {delta_spec_hall['point']:+.4f} "
      f"[{delta_spec_hall['lo']:+.4f}, {delta_spec_hall['hi']:+.4f}]")
    P(f"  full - (iii)              {delta_full_spec['point']:+.4f} "
      f"[{delta_full_spec['lo']:+.4f}, {delta_full_spec['hi']:+.4f}]")

    out = dict(
        ka=dict(full_rho=float(rho_full), full_ref=0.31738512466156,
                full_delta=abs(rho_full - 0.31738512466156),
                annotation_rho=float(L.spearman(yu, p_ann)), annotation_ref=-0.056642284791368284),
        unit=dict(n=len(u), n_components=int(len(set(compu.tolist()))),
                  fold_sizes=np.bincount(foldu, minlength=5).tolist()),
        columns=dict(hallmark_only=HALL, spectral_only=SPEC, full=L.FEATS290),
        matched_sizes=dict(zip(L.PROGRAMMES, sizes)), prolif_size=prolif_size,
        ridge_alpha=dict(hallmark={str(k): v for k, v in i_hall.items()},
                         spectral={str(k): v for k, v in i_spec.items()}),
        fits={k: ci[k] for k in ("hallmark_only", "spectral_only", "full", "annotation")},
        random_matched=dict(n_draws=N_DRAW, mean=float(r_draw.mean()), sd=float(r_draw.std()),
                            lo=float(r_draw.min()), hi=float(r_draw.max()), draws=draws,
                            n_reaching_hallmark=int((r_draw >= ci["hallmark_only"]["point"]).sum()),
                            per_draw_ci={f"rand{d:02d}": ci[f"rand{d:02d}"] for d in range(N_DRAW)}),
        delta_hallmark_minus_random_mean=delta_mean,
        delta_hallmark_minus_random_best=delta_best,
        delta_spectral_minus_hallmark=delta_spec_hall,
        delta_full_minus_spectral=delta_full_spec,
        seconds=time.time() - T0,
        peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30)
    (WORK / "ablation298.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "ablation298.log").write_text("\n".join(LOG))
    P(f"\npeak RSS {out['peak_rss_gib']:.2f} GiB; wrote work/ablation298.json")


if __name__ == "__main__":
    main()
