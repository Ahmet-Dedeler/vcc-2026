"""316 -- how shape-only is the parity arm on the SCORED (non-panel) columns at PRODUCTION's panel share?

The norm match is on the FULL row (production's `MB:1464-1467`, and the mandate's formula), so the
parity arm's full-axis delivered norm is 1.000000 by construction.  308's harness rows are panel-ZEROED,
so for 308 full-axis and non-panel are the same number; in production they are not, and the panel is
33 % of 313's fixture axis against 300/18,533 = 1.62 % of production's.  This sweeps 313's own fixture
construction over axis width, running the REAL builder functions each time, and reports the delivered
non-panel norm ratio of the parity arm and of the shipped arm against the incumbent.

  PYTHONDONTWRITEBYTECODE=1 uv run python experiments/316-parity-arm/scale316.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E313 = ROOT / "experiments/313-normrel-parity"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(E313))

import fixture313 as F  # noqa: E402
import parity313 as P313  # noqa: E402

EXP = -1.5


def widen(n_genes: int) -> dict:
    """`fixture313.build()` at a different axis width: the same seed, the same 20 targets and 20 panel
    columns, the same coverage FRACTIONS, the same three special rows, the same PerturbMatch window."""
    rng = np.random.default_rng(F.SEED)
    genes = [f"G{i:03d}" for i in range(n_genes)]
    tg = F.targets()
    T = F.N_TARGETS
    shared = rng.normal(0.0, 0.35, size=(T, n_genes))
    scale = n_genes / F.N_GENES
    cover_bounds = {n: (int(round(lo * scale)), int(round(hi * scale)))
                    for n, (lo, hi) in F.GENE_COVER.items()}
    out: dict[str, np.ndarray] = {}
    for name in F.DONORS:
        lo, hi = cover_bounds[name]
        cover = np.arange(lo, hi)
        keep = [t for t in tg if t not in F.PERT_MISSING[name]]
        delta = np.zeros((len(keep), cover.size))
        for r, t in enumerate(keep):
            ti = tg.index(t)
            row = 0.55 * shared[ti, cover] + rng.normal(0.0, 0.30, size=cover.size)
            if lo <= ti < hi:
                row[ti - lo] = -2.4 + rng.normal(0.0, 0.4)
            delta[r] = row
        if name == F.BIG_ROW[0]:
            delta[keep.index(F.BIG_ROW[1])] *= 8.0
        if name == F.SMALL_ROW[0]:
            delta[keep.index(F.SMALL_ROW[1])] *= 0.05
        if name == F.NEAR_ZERO[0]:
            r = keep.index(F.NEAR_ZERO[1])
            other = np.delete(np.arange(len(keep)), r)
            delta[r] = delta[other].mean(0)
        out[f"delta__{name}"] = delta
        out[f"genes__{name}"] = np.asarray([genes[i] for i in cover])
        out[f"perts__{name}"] = np.asarray(keep)
        out[f"b_ctrl__{name}"] = np.abs(rng.normal(2.0, 0.5, size=cover.size))
    pm_lo, pm_hi = (int(round(x * scale)) for x in F.PM_GENES)
    pm_cover = np.arange(pm_lo, pm_hi)
    pm_keep = [t for t in tg if t not in F.PM_MISSING]
    pm_delta = np.zeros((len(pm_keep), pm_cover.size))
    for r, t in enumerate(pm_keep):
        ti = tg.index(t)
        pm_delta[r] = 0.45 * shared[ti, pm_cover] + rng.normal(0.0, 0.45, size=pm_cover.size)
    out["pm_delta"] = pm_delta
    out["pm_genes"] = np.asarray([genes[i] for i in pm_cover])
    out["pm_perts"] = np.asarray(pm_keep)
    out["official"] = np.asarray(genes)
    out["targets"] = np.asarray(tg)
    out["donors"] = np.asarray(F.DONORS)
    out["panel_positions"] = np.arange(T, dtype=np.int64)
    return out


def one(n_genes: int, MB) -> dict:
    fx = widen(n_genes)
    targets = [str(t) for t in fx["targets"]]
    official = [str(g) for g in fx["official"]]
    panel = np.asarray(fx["panel_positions"], dtype=np.int64)
    model = P313.build_model(fx)
    pm = P313.build_pm(fx, official)
    ALPHA, BETA = float(MB.ALPHA), float(MB.BETA)
    uni = MB.donor_weights(model)
    op = MB.production_operator_columns(model, [official[i] for i in panel])
    npn = np.asarray([c for c in range(n_genes) if c not in set(panel.tolist())], dtype=np.int64)

    def run(exponent, pmmode):
        MB.NORMREL_EXP = exponent
        MB.NORMREL_PM = pmmode
        handled: set[str] = set()
        pre = MB.build_v9_pre_whiten_effects(
            model, targets, alpha=ALPHA, beta=BETA, weights=uni,
            pm=pm, operator_columns=op, pm_handled_out=handled)
        out, _ = MB.apply_perturbmatch_k562_replacement(
            pre, model, pm, targets, alpha=ALPHA, weights=uni, already_blended=handled)
        return out

    inc = run(None, None)
    shipped = run(EXP, None)
    parity = run(EXP, "blend")

    def ratios(A, cols):
        na = np.linalg.norm(A[:, cols], axis=1)
        nb = np.linalg.norm(inc[:, cols], axis=1)
        return na / np.maximum(nb, 1e-300)

    full_all = np.arange(n_genes)
    e_pan = np.linalg.norm(inc[:, panel], axis=1) ** 2
    e_all = np.linalg.norm(inc, axis=1) ** 2
    r_par_np, r_shp_np = ratios(parity, npn), ratios(shipped, npn)
    r_par_f, r_shp_f = ratios(parity, full_all), ratios(shipped, full_all)
    return dict(
        n_genes=int(n_genes), panel_share_of_columns=float(panel.size / n_genes),
        median_panel_energy_share=float(np.median(e_pan / np.maximum(e_all, 1e-300))),
        parity_full=[float(np.median(r_par_f)), float(r_par_f.min()), float(r_par_f.max())],
        parity_nonpanel=[float(np.median(r_par_np)), float(r_par_np.min()), float(r_par_np.max())],
        shipped_full=[float(np.median(r_shp_f)), float(r_shp_f.min()), float(r_shp_f.max())],
        shipped_nonpanel=[float(np.median(r_shp_np)), float(r_shp_np.min()), float(r_shp_np.max())],
        parity_matrix_nonpanel=float(np.linalg.norm(parity[:, npn]) / np.linalg.norm(inc[:, npn])),
        shipped_matrix_nonpanel=float(np.linalg.norm(shipped[:, npn]) / np.linalg.norm(inc[:, npn])),
        parity_matrix_full=float(np.linalg.norm(parity) / np.linalg.norm(inc)),
        shipped_matrix_full=float(np.linalg.norm(shipped) / np.linalg.norm(inc)),
    )


def main() -> None:
    import make_submission_pm_cd4 as MB  # noqa: E402
    widths = (60, 120, 300, 600, 1235)
    out = [one(g, MB) for g in widths]
    hdr = (f"{'genes':>6s} {'panel%c':>7s} {'panelE%':>8s} | "
           f"{'parity non-panel [med, min, max]':>36s} | {'shipped non-panel [med, min, max]':>36s} | "
           f"{'par|M|':>7s} {'shp|M|':>7s}")
    lines = [hdr]
    for r in out:
        lines.append(
            f"{r['n_genes']:6d} {100 * r['panel_share_of_columns']:7.2f} "
            f"{100 * r['median_panel_energy_share']:8.2f} | "
            f"{r['parity_nonpanel'][0]:11.6f} {r['parity_nonpanel'][1]:11.6f} {r['parity_nonpanel'][2]:11.6f} | "
            f"{r['shipped_nonpanel'][0]:11.6f} {r['shipped_nonpanel'][1]:11.6f} {r['shipped_nonpanel'][2]:11.6f} | "
            f"{r['parity_matrix_nonpanel']:7.4f} {r['shipped_matrix_nonpanel']:7.4f}")
    lines.append("")
    lines.append("full-axis delivered norm (the quantity the norm match controls):")
    for r in out:
        lines.append(f"  {r['n_genes']:6d}  parity [{r['parity_full'][1]:.7f}, {r['parity_full'][2]:.7f}]  "
                     f"shipped [{r['shipped_full'][1]:.6f}, {r['shipped_full'][2]:.6f}]  "
                     f"whole matrix parity {r['parity_matrix_full']:.7f} shipped {r['shipped_matrix_full']:.6f}")
    text = "\n".join(lines)
    print(text)
    (WORK / "scale316.json").write_text(json.dumps(out, indent=2))
    (WORK / "scale316.log").write_text(text + "\n")


if __name__ == "__main__":
    main()
