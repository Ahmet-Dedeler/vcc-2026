"""204 Phase 0 — reproduce 194's achieved-repression table by an INDEPENDENT route.

194 measured a target's achieved repression from the cell level:

    r_T = sum_c obs_c(T) / sum_c exp_c(T),      exp_c(T) = p_ctrl(T) * depth_c

summed over the target's own perturbed cells, where `obs_c(T)` is that cell's count of
its own target gene and `depth_c` its library size.  That ratio is *identical* to a
closed form on the shipped pseudobulk table, because both numerator and denominator are
pooled sums:

    sum_c obs_c / (p_ctrl * sum_c depth_c)  =  f_T(T) / f_ctrl(T)

and the builder stores exactly those two fractions:

    delta[T, g] = log1p(TS * f_T(g)) - log1p(TS * f_ctrl(g))
    b_ctrl[g]   = log1p(TS * f_ctrl(g))
 => r_T        = expm1(delta[T, T] + b_ctrl[T]) / expm1(b_ctrl[T])

So the shipped tables alone determine every donor's achieved repression, with no cell
streaming at all.  This script computes that, and checks it against 194's stored
per-cell `_kd_target.npy` on 194's own restricted row sets.  Two entirely different code
paths (one streaming 10M cells, one three lines of algebra on a 40 MB matrix) reading the
same number is the known-answer row for the whole experiment.

    uv run python repress204.py
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)

TS = 5e4
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
E194 = ROOT / "experiments/194-donor-reestimation/work"

# 194 RESULTS.md section 1, the numbers this must reproduce
TARGET_194 = {
    "k562_gwps": dict(median_residual=0.155, frac_below_0p2=0.476),
    "orion_hct116": dict(median_residual=0.300, frac_below_0p2=0.214),
    "orion_hek293t": dict(median_residual=0.510, frac_below_0p2=0.044),
}


def sha256_file(p: Path, cap: int | None = None) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while True:
            b = f.read(1 << 22)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def load_donor(name: str) -> dict:
    d = ROOT / "data/donors" / name
    genes = [g.strip() for g in (d / "genes.txt").read_text().splitlines() if g.strip()]
    perts = [p.strip() for p in (d / "perts.txt").read_text().splitlines() if p.strip()]
    return dict(dir=d, genes=genes, perts=perts,
                gpos={g: i for i, g in enumerate(genes)},
                ppos={p: i for i, p in enumerate(perts)})


def residual_closed_form(name: str) -> tuple[list[str], np.ndarray]:
    """r_T for every row of the shipped table; NaN where the target gene is off-axis."""
    b = load_donor(name)
    b_ctrl = np.load(b["dir"] / "b_ctrl.npy").astype(np.float64)
    delta = np.load(b["dir"] / "delta.npy", mmap_mode="r")
    ctrl_lin = np.expm1(b_ctrl)
    r = np.full(len(b["perts"]), np.nan)
    for i, t in enumerate(b["perts"]):
        g = b["gpos"].get(t)
        if g is None or ctrl_lin[g] <= 0:
            continue
        r[i] = float(np.expm1(float(delta[i, g]) + b_ctrl[g])) / float(ctrl_lin[g])
    return b["perts"], r


def main() -> None:
    out: dict = {"inputs": {}, "operators": {}, "donors": {}}

    # ---- rule 11: a path and a sha for every operator on the table ---------------------
    for p in sorted((ROOT / "data/donors").glob("whitener*.npz")):
        with np.load(p) as z:
            shp = {k: list(np.shape(z[k])) for k in z.files}
        out["operators"][str(p.relative_to(ROOT))] = dict(
            sha256=sha256_file(p), shapes=shp)
        print(f"[op] {p.relative_to(ROOT)}  sha {out['operators'][str(p.relative_to(ROOT))]['sha256'][:8]}  {shp}",
              flush=True)

    for name in DONORS:
        d = ROOT / "data/donors" / name
        out["inputs"][name] = {
            f: sha256_file(d / f)[:16] for f in ("delta.npy", "b_ctrl.npy")}

    # ---- the closed-form repression table ---------------------------------------------
    for name in DONORS:
        perts, r = residual_closed_form(name)
        np.save(WORK / f"resid-{name}.npy", r)
        (WORK / f"resid-{name}.txt").write_text("\n".join(perts))
        fin = np.isfinite(r)
        rec = dict(n_rows=len(perts), n_with_readout=int(fin.sum()),
                   median_residual=float(np.nanmedian(r)),
                   p10=float(np.nanpercentile(r, 10)),
                   p90=float(np.nanpercentile(r, 90)),
                   frac_below_0p2=float(np.nanmean(r[fin] < 0.2)),
                   mean_residual=float(np.nanmean(r)))

        # ---- known-answer: 194's restricted row set, its own per-cell measurement ------
        sd = E194 / f"sums-{name}-restricted"
        if sd.exists():
            meta = json.loads((sd / "meta.json").read_text())
            names194 = list(meta["names"])
            kd = np.load(sd / "_kd_target.npy")[0][:-1]      # drop the control group
            pos = {p: i for i, p in enumerate(perts)}
            idx = np.array([pos[t] for t in names194])
            mine = r[idx]
            both = np.isfinite(mine) & np.isfinite(kd)
            rec["known_answer_194"] = dict(
                n_rows=len(names194), n_both=int(both.sum()),
                median_194=float(np.nanmedian(kd)),
                median_mine=float(np.nanmedian(mine)),
                frac_below_0p2_194=float(np.nanmean(kd[np.isfinite(kd)] < 0.2)),
                frac_below_0p2_mine=float(np.nanmean(mine[np.isfinite(mine)] < 0.2)),
                max_abs_diff=float(np.max(np.abs(mine[both] - kd[both]))),
                median_abs_diff=float(np.median(np.abs(mine[both] - kd[both]))),
                corr=float(np.corrcoef(mine[both], kd[both])[0, 1]))
        rec["target_194_RESULTS"] = TARGET_194[name]
        out["donors"][name] = rec
        print(f"[{name}] median r {rec['median_residual']:.4f} "
              f"(194 restricted: {rec.get('known_answer_194', {}).get('median_mine', float('nan')):.4f} "
              f"vs 194's own {rec.get('known_answer_194', {}).get('median_194', float('nan')):.4f}) "
              f"frac<0.2 {rec['frac_below_0p2']:.4f}", flush=True)

    (WORK / "repress204.json").write_text(json.dumps(out, indent=2))
    print("wrote work/repress204.json")


if __name__ == "__main__":
    main()
