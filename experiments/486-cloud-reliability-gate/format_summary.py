"""Print the evaluation tables from work/results/summary.json. No refitting."""

from __future__ import annotations

import json
from pathlib import Path


def main() -> None:
    report = json.loads((Path(__file__).resolve().parent / "work" / "results" / "summary.json").read_text())
    print(
        f"lines learned>equal (ew CI>0): {report['n_lines_learned_beats_equal_ew']}/{report['n_lines']}"
    )
    print(
        f"lines learned>fixed (ew CI>0): {report['n_lines_learned_beats_fixed_ew']}/{report['n_lines']}"
    )
    print(
        f"lines fixed>equal (ew CI>0): {report['n_lines_fixed_beats_equal_ew']}/{report['n_lines']}"
    )
    for fold in report["folds"]:
        print(f"\n== {fold['held_out']} n={fold['n_targets']} {fold['seconds']:.0f}s ==")
        print(f"  coef {fold['fit']['coef']} b={fold['fit']['intercept']:.3f} nfev={fit_nfev(fold)}")
        for method, block in fold["metrics"]["methods"].items():
            ew = block["ew"]
            uw = block["uw"]
            rk = block["rank"]
            print(
                f"  {method:12s} ew={ew['point']:.4f} [{ew['ci95'][0]:.4f},{ew['ci95'][1]:.4f}] "
                f"uw={uw['point']:.4f} rank={rk['point']:.4f} top1={block['top1']['point']:.4f}"
            )
        print("  gains:")
        for name, metrics in fold["metrics"]["gains"].items():
            g = metrics["ew"]
            print(
                f"    {name} ew={g['point']:+.4f} [{g['ci95'][0]:+.4f},{g['ci95'][1]:+.4f}] "
                f"frac+={g['frac_positive']:.2f}"
            )
        drops = sorted(fold["ablation"], key=lambda r: -r["drop_vs_learned"])
        print("  ablation drops:", ", ".join(f"{r['feature']} {r['drop_vs_learned']:+.4f}" for r in drops))


def fit_nfev(fold: dict) -> int:
    return int(fold["fit"]["nfev"])


if __name__ == "__main__":
    main()
