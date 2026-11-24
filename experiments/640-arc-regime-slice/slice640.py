"""640: re-price 438's October decisions on knockdown-matched target subsets.

Arc's 2026-10-01 blog: every panel target has >= 80 % median on-target knockdown in every line (median 89 %), and the
same sgRNA sequences are used in all six lines (2025: Replogle's K562 day-8 library protospacers). 438's October-like
draws are scored on public truths (K562 GWPS, Orion HCT116 / HEK293T held out) whose per-target knockdown varies widely.
Question: do 438's verdicts (above all CD4 x1.5 losing on every draw) hold on the targets that look like Arc's regime?

Subsets (each fold separately; resid = 204's perturbed/control expression of the target in that bank's pseudobulk):
  truth_strong   held-out truth's own resid <= that fold's median over the panel (the stronger-knockdown half)
  truth_weak     the complement (resid > fold median)
  truth_le03     truth resid <= 0.3 (K562 / HCT116 folds only have enough)
  k562_le02      K562 GWPS resid <= 0.2 (Arc-guide proxy: same library as Arc's 2025 constructs), all folds
  k562_gt02      K562 GWPS resid > 0.2 (incl. no K562 measurement)
  no_tail        truth resid <= 0.5 (drop only the weak tail the draws carry ~2x more of)
Prices: 438's levelled route, k = 1, cancer-fold mean (robust438.mean_price), energies fixed at the full panel.
    uv run python experiments/640-arc-regime-slice/slice640.py
"""
from __future__ import annotations

import sys
sys.dont_write_bytecode = True
from pathlib import Path
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/438-october-stack-repricing"))
from price438 import CANCER, REF, WORK, Panel  # noqa: E402
from robust438 import STEPS  # noqa: E402

W204 = ROOT / "experiments/204-repression-weighted-mixing/work"
OUT = Path(__file__).resolve().parent / "work"


def resid(bank):
    names = (W204 / f"resid-{bank}.txt").read_text().split()
    r = np.load(W204 / f"resid-{bank}.npy")
    return dict(zip(names, r))


R = {b: resid(b) for b in CANCER}


def targets(tag, f):
    return [str(t) for t in np.load(WORK / tag / f"fold-{f}-held.npz", allow_pickle=True)["targets"]]


def keep_mask(tag, f, rule):
    tg = targets(tag, f)
    tv = np.array([R[f].get(t, np.nan) for t in tg])
    kv = np.array([R["k562_gwps"].get(t, np.nan) for t in tg])
    med = np.nanmedian(tv)
    if rule == "all":
        m = np.ones(len(tg), bool)
    elif rule == "truth_strong":
        m = np.isfinite(tv) & (tv <= med)
    elif rule == "truth_weak":
        m = np.isfinite(tv) & (tv > med)
    elif rule == "truth_le03":
        m = np.isfinite(tv) & (tv <= 0.3)
    elif rule == "no_tail":
        m = ~(np.isfinite(tv) & (tv > 0.5))
    elif rule == "k562_le02":
        m = np.isfinite(kv) & (kv <= 0.2)
    elif rule == "k562_gt02":
        m = ~(np.isfinite(kv) & (kv <= 0.2))
    else:
        raise ValueError(rule)
    return m


def panel(tag, rule):
    kw = dict(full_tag="efr_s1" if tag == "efr" else None, energy="fold" if tag == "efr" else "full")
    P = Panel(tag, **kw)
    n = {}
    for f in P.folds:
        m = keep_mask(tag, f, rule)
        P.keep[f] = None if m.all() else m
        n[f] = int(m.sum())
    P._st = {}
    return P, n


def mean_price(P, arm, ref):
    p = [P.price(f, arm, 1.0, ref, True)["ov"] for f in P.folds]
    return float(np.mean(p)), int(sum(v > 0 for v in p)), p


def main():
    rules = ("all", "truth_strong", "truth_weak", "no_tail", "truth_le03", "k562_le02", "k562_gt02")
    out = {}
    lines = ["| step | panel | " + " | ".join(rules) + " |", "|---|---|" + "---|" * len(rules)]
    ntab = {}
    for lab, arm, ref in STEPS:
        for tag in ("abc", "real", "efp", "efr"):
            cells = []
            for rule in rules:
                try:
                    P, n = panel(tag, rule)
                except FileNotFoundError:
                    cells.append("n/a")
                    continue
                if len(P.folds) < 3 or not all(f"{arm}:L" in P.R[f] for f in P.folds):
                    cells.append("n/a")
                    continue
                if min(n.values()) < 15:
                    cells.append(f"(n<15: {min(n.values())})")
                    continue
                mp, npos, per = mean_price(P, arm, ref)
                out.setdefault(lab, {}).setdefault(tag, {})[rule] = dict(mean=mp, n_pos=npos, per_fold=dict(zip(P.folds, per)), n=n)
                ntab[(tag, rule)] = n
                cells.append(f"{mp:+.5f} {npos}/3")
            lines.append(f"| {lab} | {tag} | " + " | ".join(cells) + " |")
            print(lines[-1], flush=True)
    lines.append("")
    lines.append("Targets kept per fold (k562 / hct116 / hek293t):")
    for (tag, rule), n in sorted(ntab.items()):
        lines.append(f"- {tag} {rule}: " + " / ".join(str(n[f]) for f in CANCER))
    OUT.mkdir(exist_ok=True)
    (OUT / "slice640.json").write_text(json.dumps(out, indent=1))
    (OUT / "slice640.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
