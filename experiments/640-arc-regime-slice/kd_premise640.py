"""640 premise check: do the October-like draws carry weaker truth knockdown than Arc's A/B/C panel, per held-out fold?

Truth knockdown per target = 204's residual (perturbed / control expression of the target in the held-out bank's own
pseudobulk; 0 = full repression). Arc's published rule (2026-10-01 blog): every panel target has >= 80 % median
knockdown in every line, i.e. residual <= 0.2.
    uv run python experiments/640-arc-regime-slice/kd_premise640.py
"""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
W438 = ROOT / "experiments/438-october-stack-repricing/work"
W204 = ROOT / "experiments/204-repression-weighted-mixing/work"
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")


def resid(bank):
    names = (W204 / f"resid-{bank}.txt").read_text().split()
    r = np.load(W204 / f"resid-{bank}.npy")
    assert len(names) == r.size, bank
    return dict(zip(names, r))


def main():
    R = {b: resid(b) for b in CANCER}
    out = {}
    for tag in ("abc", "real", "efp", "efr"):
        for f in CANCER:
            p = W438 / tag / f"fold-{f}-held.npz"
            if not p.exists():
                continue
            tg = [str(t) for t in np.load(p, allow_pickle=True)["targets"]]
            v = np.array([R[f].get(t, np.nan) for t in tg])
            m = np.isfinite(v)
            row = dict(n=len(tg), n_kd=int(m.sum()), median=float(np.median(v[m])),
                       frac_le_0_2=float(np.mean(v[m] <= 0.2)), frac_le_0_3=float(np.mean(v[m] <= 0.3)),
                       frac_gt_0_5=float(np.mean(v[m] > 0.5)), q25=float(np.quantile(v[m], .25)),
                       q75=float(np.quantile(v[m], .75)))
            out[f"{tag}|{f}"] = row
            print(f"{tag:5s} {f:14s} n {row['n']:3d} kd {row['n_kd']:3d}  median resid {row['median']:.3f}  "
                  f"IQR [{row['q25']:.3f}, {row['q75']:.3f}]  <=0.2 {row['frac_le_0_2']:.2f}  <=0.3 {row['frac_le_0_3']:.2f}  "
                  f">0.5 {row['frac_gt_0_5']:.2f}")
    (Path(__file__).parent / "work").mkdir(exist_ok=True)
    (Path(__file__).parent / "work/kd_premise640.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
