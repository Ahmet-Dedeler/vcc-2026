"""810 -- the decomposition tables (markdown) from analyse810's per-fold summaries: means over folds per surface set and tier,
with the number of folds where a quantity is positive. Local, light.
    uv run python experiments/810-strong-direction/table810.py > experiments/810-strong-direction/work/table810.md
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import analyse810 as A  # noqa: E402

SETS = {"BULK-W (18 L folds)": lambda s: s.startswith("L"),
        "BULK-ABC cancer (3)": lambda s: s.startswith("ABC") and "kolf2" not in s,
        "ABC/kolf2": lambda s: s == "ABC__kolf2"}
TIERS = {0: "weak", 1: "mid", 2: "strong"}


def m(rows, k, fmt="{:+.4f}", pos=False):
    v = np.asarray([r[k] for r in rows if r.get(k) is not None], float)
    if v.size == 0:
        return "—"
    s = fmt.format(v.mean())
    if pos:
        s += f" ({int((v > 0).sum())}/{v.size})"
    return s


def main():
    res, PT = A.main()
    out = []
    for name, sel in SETS.items():
        stems = [s for s in res if sel(s)]
        if not stems:
            continue
        out.append(f"\n### {name}\n")
        out.append("| tier | X_real | share of a / d | cos_real med | cos_fin med | contrib base / cd4 / h1 / part | energy base / cd4 / h1 / part "
                   "| H1 carried | med cos base / cd4 / h1 / pre / U / src0 / src1 | donor agreement src01 (raw01) | cos(fin,U) |")
        out.append("|---|---|---|---|---|---|---|---|---|---|---|")
        for t in (0, 1, 2):
            rows = [res[s][t] for s in stems if t in res[s]]
            out.append(f"| {TIERS[t]} | {m(rows, 'X_real')} | {m(rows, 'share_a', '{:.2f}')} / {m(rows, 'share_d', '{:.2f}')} "
                       f"| {m(rows, 'med_cos_real')} | {m(rows, 'med_cos_fin')} "
                       f"| {m(rows, 'contrib_base')} / {m(rows, 'contrib_cd4')} / {m(rows, 'contrib_h1')} / {m(rows, 'contrib_part')} "
                       f"| {m(rows, 'energy_base', '{:.3f}')} / {m(rows, 'energy_cd4', '{:.3f}')} / {m(rows, 'energy_h1', '{:.3f}')} / {m(rows, 'energy_part', '{:.3f}')} "
                       f"| {m(rows, 'carry_h1', '{:.2f}')} "
                       f"| {m(rows, 'med_cos_base')} / {m(rows, 'med_cos_cd4')} / {m(rows, 'med_cos_h1')} / {m(rows, 'med_cos_pre')} / {m(rows, 'med_cos_U')} / {m(rows, 'med_cos_src0')} / {m(rows, 'med_cos_src1')} "
                       f"| {m(rows, 'med_cos_src01')} ({m(rows, 'med_cos_raw01')}) | {m(rows, 'med_cos_fin_U', '{:.3f}')} |")
        out.append("\nDirection-only ablations at the realised norms (ΔX_tier, composition level; (k/N) = folds > 0):\n")
        out.append("| tier | Xc_fin | − base | − cd4 | − h1 | − part | U only | pre only | src0 only | src1 only | blend U 0.25 | 0.5 | 1 | base→pre 0.5 | base→pre 1 | LS oracle parts | LS oracle all |")
        out.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for t in (0, 1, 2):
            rows = [res[s][t] for s in stems if t in res[s]]
            ks = ["dXc_no_base", "dXc_no_cd4", "dXc_no_h1", "dXc_no_part", "dXc_only_U", "dXc_only_pre", "dXc_only_src0",
                  "dXc_only_src1", "dXc_blendU_b0.25", "dXc_blendU_b0.5", "dXc_blendU_b1", "dXc_swapbase_pre_b0.5", "dXc_swapbase_pre_b1", "dXc_oracle_parts", "dXc_oracle_all"]
            out.append(f"| {TIERS[t]} | {m(rows, 'Xc_fin')} | " + " | ".join(m(rows, k, pos=True) for k in ks) + " |")
        out.append("\nKnockdown (raw log-ratio of the target's own gene; median over the tier's targets, mean over folds):\n")
        out.append("| tier | kd sources | kd held line | pds (literal, incumbent) |")
        out.append("|---|---|---|---|")
        for t in (0, 1, 2):
            rows = [res[s][t] for s in stems if t in res[s]]
            out.append(f"| {TIERS[t]} | {m(rows, 'med_kd_src', '{:+.3f}')} | {m(rows, 'med_kd_held', '{:+.3f}')} | {m(rows, 'pds', '{:.4f}')} |")
    print("\n".join(out))


if __name__ == "__main__":
    main()
