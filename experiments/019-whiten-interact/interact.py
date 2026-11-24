"""019 — does whitening change the rest of the recipe?

Experiment 017 measured whitening on top of the incumbent configuration: `gene_norm`
treatment, donors at 4:1:1. Both of those were chosen *without* whitening in the picture, and
both plausibly interact with it:

- **`gene_norm`** divides each gene by its mean |effect| across perturbations. That is the
  *diagonal* of the same covariance whitening inverts properly. If whitening subsumes it,
  `gene_norm` should become neutral or harmful once whitening is on -- and carrying a
  redundant transform costs a floor parameter and a percentile heuristic for nothing.

- **Donor weights.** Experiment 010 found equal weighting *worse than a single donor*,
  because averaging smooths predictions toward each other and destroys the discrimination
  `pds` measures. Whitening attacks that same smoothing directly. So the penalty for
  averaging should shrink, and the optimal mix should move toward equal weight.

Both are cheap to settle and both change what gets built.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from vcc2026.donor_transfer import Donor, _norm_scale   # noqa: E402
from vcc2026.whiten import Whitener                      # noqa: E402

DONORS = {
    "k562_gwps": "data/donors/k562_gwps",
    "orion_hct116": "data/donors/orion_hct116",
    "orion_hek293t": "data/donors/orion_hek293t",
}


def pds_cosine(P, T):
    Pn = P / (np.linalg.norm(P, axis=1, keepdims=True) + 1e-12)
    Tn = T / (np.linalg.norm(T, axis=1, keepdims=True) + 1e-12)
    D = 1.0 - Pn @ Tn.T
    n = D.shape[1]
    sd = np.diag(D)[:, None]
    return float((1.0 - ((D < sd).sum(1) + 0.5 * ((D == sd).sum(1) - 1)) / (n - 1)).mean())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--treats", nargs="+", default=["gene_norm", "raw"])
    ap.add_argument("--weights", nargs="+",
                    default=["1:0", "8:1", "4:1", "2:1", "1:1", "1:2"])
    ap.add_argument("--gammas", type=float, nargs="+", default=[0.0, 0.55])
    ap.add_argument("--basis-perts", type=int, default=6000)
    args = ap.parse_args()

    raw_bank = {n: Donor.load(p) for n, p in DONORS.items()}
    vcc = list(pd.read_csv("data/challenge/pert_counts.csv")["target_gene"])

    results = {}
    for treat in args.treats:
        banks = {}
        for name, d in raw_bank.items():
            _, s = d.decompose()
            if treat == "gene_norm":
                s = (s / _norm_scale(s)).astype(np.float32)
            gl = [str(g) for g in d.genes]
            banks[name] = dict(s=s, raw=d.decompose()[1],
                               gpos={g: i for i, g in enumerate(gl)}, gset=set(gl),
                               idx={p: i for i, p in enumerate(d.perts)})
        common = [g for g in vcc if all(g in b["idx"] for b in banks.values())]
        rng = np.random.default_rng(0)

        folds = {}
        for held in banks:
            others = [n for n in banks if n != held]
            shared = sorted(set.intersection(*[banks[n]["gset"] for n in [held] + others]))
            feat = [g for g in shared if g not in set(common)]
            col = {n: np.fromiter((banks[n]["gpos"][g] for g in feat), np.int64, len(feat))
                   for n in [held] + others}
            T = np.vstack([banks[held]["raw"][banks[held]["idx"][p]][col[held]]
                           for p in common]).astype(np.float64)
            sc = {n: float(np.median(np.linalg.norm(banks[n]["s"][:, col[n]], axis=1))) or 1.0
                  for n in others}
            ref = float(np.median(list(sc.values())))
            parts, Bl = {}, []
            for n in others:
                parts[n] = (ref / sc[n]) * np.vstack(
                    [banks[n]["s"][banks[n]["idx"][p]][col[n]] for p in common]
                ).astype(np.float64)
                M = banks[n]["s"][:, col[n]]
                if M.shape[0] > args.basis_perts:
                    M = M[np.sort(rng.choice(M.shape[0], args.basis_perts, False))]
                Bl.append((ref / sc[n]) * M.astype(np.float64))
            print(f"[{treat}] fitting whitener for fold {held} "
                  f"({sum(b.shape[0] for b in Bl)} x {len(feat)})", flush=True)
            folds[held] = dict(T=T, parts=parts, others=others,
                               W=Whitener.fit(np.vstack(Bl)))

        for wspec in args.weights:
            ws = [float(x) for x in wspec.split(":")]
            for g in args.gammas:
                per = {}
                for h, f in folds.items():
                    w = dict(zip(f["others"], ws))
                    tot = sum(w.values())
                    P = sum((w[n] / tot) * f["parts"][n] for n in f["others"])
                    per[h] = pds_cosine(f["W"].apply(P, g), f["T"])
                macro = float(np.mean(list(per.values())))
                key = f"{treat:9} w={wspec:5} gamma={g:g}"
                results[key] = {**per, "macro": macro}
                fs = "  ".join(f"{k.replace('orion_','')[:9]}={v:.4f}" for k, v in per.items())
                print(f"  {key:32} macro PDS {macro:.4f}   [{fs}]", flush=True)

    best = max(results, key=lambda k: results[k]["macro"])
    print(f"\nbest: {best} -> {results[best]['macro']:.4f}")
    Path("experiments/019-whiten-interact/results.json").write_text(
        json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
