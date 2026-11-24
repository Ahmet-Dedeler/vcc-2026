"""Source-only diagnostic: which per-target source property predicts reach damage?

Exploratory analysis of ALREADY-UNBLINDED experiment 149 (HEK293T) and 150
(HCT116) frames. Both used the identical fixed lambda-1 dual-moment emitter, so
per-target reach deltas are comparable across the two contexts with no dose
confound. Any hypothesis this produces MUST be preregistered before it is tested
on a new surface.

Reads only:
  - work/source-build-v1/effects.npz          (source-derived, no held expression)
  - work/scoring-v1/{incumbent,candidate}/run/results.csv   (already unblinded)
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CTX = {
    "HEK293T (149)": ROOT / "experiments/149-dual-moment-hek293t",
    "HCT116 (150)": ROOT / "experiments/150-dual-moment-score-level-hct116",
}
REACH = "de_wilcoxon_direction_reach_raw"


def per_target(path: Path, metric: str) -> pd.Series:
    frame = pd.read_csv(path)
    rows = frame[frame["metric"] == metric]
    return rows.set_index("perturbation")["value"].astype(float)


def features(effects_path: Path) -> pd.DataFrame:
    with np.load(effects_path, allow_pickle=False) as z:
        targets = z["targets"].astype(str)
        inc = z["incumbent"].astype(np.float64)
        dirn = z["direction"].astype(np.float64)
    relax = dirn - inc
    inc_n = np.linalg.norm(inc, axis=1)
    rel_n = np.linalg.norm(relax, axis=1)
    dir_n = np.linalg.norm(dirn, axis=1)

    def safe(a, b):
        return np.divide(a, b, out=np.zeros_like(a), where=b > 0)

    cos_inc_dir = safe((inc * dirn).sum(1), inc_n * dir_n)
    cos_rel_inc = safe((relax * inc).sum(1), rel_n * inc_n)
    both = (inc != 0) & (relax != 0)
    agree = ((np.sign(relax) == np.sign(inc)) & both).sum(1)
    frac_sign_agree = safe(agree.astype(np.float64), both.sum(1).astype(np.float64))
    # concentration of the incumbent effect: share of squared mass in top 1% genes
    k = max(1, inc.shape[1] // 100)
    top = np.sort(inc ** 2, axis=1)[:, -k:].sum(1)
    conc = safe(top, (inc ** 2).sum(1))
    return pd.DataFrame(
        {
            "incumbent_norm": inc_n,
            "relax_norm": rel_n,
            "relax_rel_mag": safe(rel_n, inc_n),
            "cos_incumbent_direction": cos_inc_dir,
            "cos_relax_incumbent": cos_rel_inc,
            "frac_relax_sign_agrees": frac_sign_agree,
            "incumbent_top1pct_mass": conc,
        },
        index=pd.Index(targets, name="perturbation"),
    )


def main() -> None:
    tables = {}
    for name, base in CTX.items():
        scoring = base / "work/scoring-v1"
        inc = per_target(scoring / "incumbent/run/results.csv", REACH)
        cand = per_target(scoring / "candidate/run/results.csv", REACH)
        feats = features(base / "work/source-build-v1/effects.npz")
        common = feats.index.intersection(inc.index).intersection(cand.index)
        table = feats.loc[common].copy()
        table["reach_incumbent"] = inc.loc[common]
        table["reach_candidate"] = cand.loc[common]
        table["reach_delta"] = table["reach_candidate"] - table["reach_incumbent"]
        tables[name] = table

        print(f"\n{'=' * 70}\n{name}  —  {len(table)} targets\n{'=' * 70}")
        print(f"  mean reach: incumbent {table['reach_incumbent'].mean():+.5f}"
              f"  candidate {table['reach_candidate'].mean():+.5f}"
              f"  delta {table['reach_delta'].mean():+.5f}")
        print(f"  targets improved: {(table['reach_delta'] > 0).sum()}/{len(table)}")
        print(f"\n  Spearman corr( feature , reach_delta ):")
        for col in feats.columns:
            rho = table[col].corr(table["reach_delta"], method="spearman")
            pear = table[col].corr(table["reach_delta"], method="pearson")
            print(f"    {col:<28} rho={rho:+.3f}  r={pear:+.3f}")

    # cross-context comparison on shared targets
    a, b = tables["HEK293T (149)"], tables["HCT116 (150)"]
    shared = a.index.intersection(b.index)
    print(f"\n{'=' * 70}\nCROSS-CONTEXT — {len(shared)} shared targets\n{'=' * 70}")
    print(f"  reach_delta HEK  mean {a.loc[shared,'reach_delta'].mean():+.5f}")
    print(f"  reach_delta HCT  mean {b.loc[shared,'reach_delta'].mean():+.5f}")
    rho = a.loc[shared, "reach_delta"].corr(b.loc[shared, "reach_delta"], method="spearman")
    print(f"  Spearman corr(reach_delta_HEK, reach_delta_HCT) = {rho:+.3f}")
    print("    -> if near zero, per-target reach response does NOT transfer at all")
    for col in ["relax_rel_mag", "cos_incumbent_direction", "frac_relax_sign_agrees"]:
        r = a.loc[shared, col].corr(b.loc[shared, col], method="spearman")
        print(f"  corr({col}) across contexts = {r:+.3f}")

    out = ROOT / "analysis/reach_signflip_diagnostic.csv"
    pd.concat(tables, names=["context"]).to_csv(out)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
