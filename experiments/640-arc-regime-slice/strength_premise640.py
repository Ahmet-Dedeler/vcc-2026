"""640b premise: is Arc's A/B/C panel stronger (truth response energy) than 438's uniform random October-like draws?

Arc built the panel by binning the H1 pilot's DE counts (negligible / subtle / strong) and sampling each bin for diverse
response modes (FINDINGS H.1); 199/238's draws sample the eligible public pool uniformly. Strength here = the held-out
truth's per-target response energy d (438's fold arrays, the cos_agg denominator), compared within each fold.
    uv run python experiments/640-arc-regime-slice/strength_premise640.py
"""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
W438 = ROOT / "experiments/438-october-stack-repricing/work"
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")


def d_of(tag, f):
    Z = np.load(W438 / tag / f"arrays438-{f}-held.npz")
    tg = [str(t) for t in np.load(W438 / tag / f"fold-{f}-held.npz", allow_pickle=True)["targets"]]
    return tg, np.asarray(Z["d"], float)


def main():
    out = {}
    for f in CANCER:
        tg0, d0 = d_of("abc", f)
        pooled = np.log10(d0)
        qs = np.quantile(pooled, [1 / 3, 2 / 3])
        print(f"== {f}: abc terciles of log10 d at {qs.round(3)}")
        for tag in ("abc", "real", "efp", "efr"):
            tg, d = d_of(tag, f)
            ld = np.log10(d)
            frac = [float(np.mean(ld <= qs[0])), float(np.mean((ld > qs[0]) & (ld <= qs[1]))), float(np.mean(ld > qs[1]))]
            top = np.sort(d)[::-1]
            share10 = float(top[: max(1, len(d) // 10)].sum() / d.sum())
            out[f"{tag}|{f}"] = dict(n=len(d), median_log10_d=float(np.median(ld)), tercile_fracs=frac, top10pct_energy_share=share10)
            print(f"  {tag:4s} n {len(d):3d} median log10 d {np.median(ld):+.3f}  abc-tercile fractions weak/mid/strong "
                  f"{frac[0]:.2f}/{frac[1]:.2f}/{frac[2]:.2f}  top-10% energy share {share10:.2f}")
    (Path(__file__).parent / "work/strength_premise640.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
