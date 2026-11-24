#!/usr/bin/env python3
"""712 (PREREGISTRATION.md): donor-strength tail ratios per panel and fold line, as a proxy for the LODO key's fold cap hi.

Strength of a target in line h = L2 norm of h's bank delta row over the bank's genes minus the panel's own genes.
Tail ratio = P 96.647(strength) / median(strength) over the panel's covered targets (hi's order statistic); P90 / median too.
Validity: Spearman(proxy, 522's hi546) over the 9 real / efp / efr folds >= 0.6, and A/B/C's proxy below all 9 on the same
line. L7 / L8 are sealed and never read.
usage: uv run python experiments/712-arc-panel-tails/tail712.py
"""
from __future__ import annotations

import csv
import glob
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/438-october-stack-repricing"))
import numpy as np  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import lib438 as L  # noqa: E402

P_PROD = 96.64697522769474
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t")
R238 = ROOT / "experiments/238-essential-free-panel/work"


def read_targets(path: Path) -> list[str]:
    with open(path) as f:
        return [r["target_gene"] for r in csv.DictReader(f) if r["target_gene"] != "non-targeting"]


def panels() -> dict[str, list[str]]:
    P = {"abc": read_targets(ROOT / "data/challenge/pert_counts.csv")}
    h1 = []
    for f in sorted(glob.glob(str(ROOT / "data/external/vcc2025/*/pert_counts_*.csv"))):
        h1 += read_targets(Path(f))
    P["h1_2025"] = sorted(set(h1))
    P["real"] = read_targets(ROOT / "experiments/199-def-runbook/work/def-bundle-r/pert_counts.csv")
    P["efp"] = read_targets(R238 / "ef-bundle/pert_counts.csv")
    P["efr"] = read_targets(R238 / "ef-bundle-r/pert_counts.csv")
    for k in range(1, 7):                                   # L7 / L8 are sealed: never read
        P[f"L{k}"] = read_targets(ROOT / f"experiments/801-lab-folds/work/L{k}/bundle/pert_counts.csv")
    return P


def main() -> None:
    P = panels()
    banks = {}
    for h in FOLDS:
        b = L.bank(h)
        banks[h] = (b, np.load(b.delta_path, mmap_mode="r"), {g: i for i, g in enumerate(b.genes)},
                    {q: i for i, q in enumerate(b.perts)})
    out = {"panels": {k: len(v) for k, v in P.items()}, "tails": {}}
    print("panel     n   | " + " | ".join(f"{h:>13s} cov  P96.6/med  P90/med" for h in FOLDS), flush=True)
    for name, tg in P.items():
        row = {}
        for h in FOLDS:
            b, delta, gpos, ppos = banks[h]
            cols = np.array(sorted(i for g, i in gpos.items() if g not in set(tg)))
            st = np.array([float(np.linalg.norm(np.asarray(delta[ppos[t]], np.float64)[cols])) for t in tg if t in ppos])
            med = float(np.median(st))
            row[h] = dict(covered=int(len(st)), tail=float(np.percentile(st, P_PROD) / med),
                          p90=float(np.percentile(st, 90) / med), median=med)
        out["tails"][name] = row
        print(f"{name:8s} {len(tg):4d} | " + " | ".join(f"{row[h]['covered']:17d}  {row[h]['tail']:9.2f}  {row[h]['p90']:7.2f}"
                                                         for h in FOLDS), flush=True)
    # validity: 522's hi546 on real / efp / efr
    hi, px = [], []
    for p in ("real", "efp", "efr"):
        for h in FOLDS:
            hi.append(json.loads((ROOT / f"experiments/522-october-bulk-key/work/dump522-{p}-{h}.json").read_text())["hi546"])
            px.append(out["tails"][p][h]["tail"])
    rho = float(spearmanr(px, hi).statistic)
    abc_below = all(out["tails"]["abc"][h]["tail"] < min(out["tails"][p][h]["tail"] for p in ("real", "efp", "efr"))
                    for h in FOLDS)
    valid = rho >= 0.6 and abc_below
    out["validity"] = dict(spearman_proxy_vs_hi546=rho, abc_below_all_9=abc_below, valid=valid,
                           hi546=hi, proxy=px)
    print(f"validity: Spearman(proxy, hi546) over 9 folds = {rho:+.3f} (>= 0.6?), A/B/C below all of real/efp/efr on every "
          f"line: {abc_below} -> {'VALID' if valid else 'NO READ'}", flush=True)
    if valid:
        draws = [p for p in P if p not in ("abc", "h1_2025")]
        light = sum(out["tails"]["h1_2025"][h]["tail"] <= out["tails"]["abc"][h]["tail"]
                    and out["tails"]["h1_2025"][h]["tail"] < min(out["tails"][p][h]["tail"] for p in draws) for h in FOLDS)
        inside = sum(min(out["tails"][p][h]["tail"] for p in draws) <= out["tails"]["h1_2025"][h]["tail"]
                     <= max(out["tails"][p][h]["tail"] for p in draws) for h in FOLDS)
        verdict = "ARC-LIGHT" if light >= 2 else ("A/B/C-UNUSUAL" if inside >= 2 else "MIXED")
        out["reading"] = dict(lines_light=light, lines_inside_draws=inside, verdict=verdict)
        print(f"reading: H1 2025 light on {light}/3 lines, inside the draws' range on {inside}/3 -> {verdict}", flush=True)
    (HERE / "work" / "tail712.json").write_text(json.dumps(out, indent=1))
    print("wrote work/tail712.json")


if __name__ == "__main__":
    main()
