"""022 — re-test the graveyard under whitening.

Equal donor weighting was *the worst* configuration in experiment 010 and is *the best* in
experiment 019. Nothing about the data changed; whitening changed. That is worth taking
seriously as a general warning, because several other rejected ideas were rejected with the
**same mechanism** attached to them:

| rejected idea | the reason given | is it the hubness story? |
|---|---|---|
| equal donor weighting | "smooths predictions toward the shared response, killing discrimination" | yes — and it flipped |
| empirical-Bayes shrinkage | "+10.6% cosine, -0.017 on the board. Wrong objective." | yes — cosine up, discrimination down |
| per-gene inverse-variance pooling | "over-weights whichever donor expresses a gene highly" | yes — abundance is the shared component |
| low-rank truncation | "the transferable signal is high-dimensional" | no — truncation keeps the shared part |
| basal-similarity donor weighting | "concentrates weight on the wrong donor" | no — an accuracy failure, not a geometry one |

The first three all describe *the dominant shared programme drowning the perturbation-specific
signal*, which is precisely what whitening removes. If the mechanism is right, they should
move toward neutral or positive; the last two should not move at all. That asymmetry is the
actual test — a result where everything improves would mean whitening is just adding slack.

Also swept: how many donor perturbations the covariance is estimated from. 6,000 per donor was
an arbitrary cap.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from vcc2026.donor_transfer import Donor, _norm_scale, eb_shrink   # noqa: E402
from vcc2026.whiten import Whitener                                 # noqa: E402

DONORS = {"k562_gwps": "data/donors/k562_gwps",
          "orion_hct116": "data/donors/orion_hct116",
          "orion_hek293t": "data/donors/orion_hek293t"}
CELLS = {n: json.loads((Path(p) / "meta.json").read_text())["cells_per_pert"]["median"]
         for n, p in DONORS.items()}


def pds(pred, true):
    P = pred / (np.linalg.norm(pred, axis=1, keepdims=True) + 1e-12)
    T = true / (np.linalg.norm(true, axis=1, keepdims=True) + 1e-12)
    D = 1.0 - P @ T.T
    n = D.shape[1]
    sd = np.diag(D)[:, None]
    return float((1.0 - ((D < sd).sum(1) + 0.5 * ((D == sd).sum(1) - 1)) / (n - 1)).mean())


raw_donor = {n: Donor.load(p) for n, p in DONORS.items()}
vcc = list(pd.read_csv("data/challenge/pert_counts.csv")["target_gene"])

# Three variants of the per-donor effect matrix, all target-blind.
VARIANTS = ["gene_norm", "eb", "gene_norm+eb"]
banks = {}
for name, d in raw_donor.items():
    _, s = d.decompose()
    gl = [str(g) for g in d.genes]
    v = {}
    v["gene_norm"] = (s / _norm_scale(s)).astype(np.float32)
    v["eb"] = eb_shrink(s, d.b_ctrl, CELLS[name])
    v["gene_norm+eb"] = (eb_shrink(s, d.b_ctrl, CELLS[name]) / _norm_scale(s)).astype(np.float32)
    banks[name] = dict(v=v, raw=s, b_ctrl=d.b_ctrl,
                       gpos={g: i for i, g in enumerate(gl)}, gset=set(gl),
                       idx={p: i for i, p in enumerate(d.perts)})

common = [g for g in vcc if all(g in b["idx"] for b in banks.values())]
print(f"panel {len(common)}\n")
rng = np.random.default_rng(0)
GAMMAS = [0.0, 0.65]
results = {}

for held in banks:
    others = [n for n in banks if n != held]
    shared_g = sorted(set.intersection(*[banks[n]["gset"] for n in [held] + others]))
    feat = [g for g in shared_g if g not in set(common)]
    col = {n: np.fromiter((banks[n]["gpos"][g] for g in feat), np.int64, len(feat))
           for n in [held] + others}
    T = np.vstack([banks[held]["raw"][banks[held]["idx"][p]][col[held]]
                   for p in common]).astype(np.float64)

    # Poisson precision per (donor, gene) from the donor's own control profile -- the
    # inverse-variance weights that lost by 0.020 in experiment 015.
    prec = {}
    for n in others:
        e = np.expm1(banks[n]["b_ctrl"][col[n]])
        pg = e / max(e.sum(), 1e-9)
        prec[n] = np.maximum(CELLS[n] * 20000.0 * pg, 1e-3)

    for variant in VARIANTS:
        sc = {n: float(np.median(np.linalg.norm(banks[n]["v"][variant][:, col[n]], axis=1)))
              or 1.0 for n in others}
        ref = float(np.median(list(sc.values())))
        parts = {n: (ref / sc[n]) * np.vstack(
            [banks[n]["v"][variant][banks[n]["idx"][p]][col[n]] for p in common]
        ).astype(np.float64) for n in others}

        combos = {
            "equal": 0.5 * (parts[others[0]] + parts[others[1]]),
            "primary4:1": (4 * parts[others[0]] + parts[others[1]]) / 5,
        }
        # inverse-variance pooling, per gene
        wsum = prec[others[0]] + prec[others[1]]
        combos["invvar"] = ((prec[others[0]] * parts[others[0]]
                             + prec[others[1]] * parts[others[1]]) / wsum)
        # magnitude-discarding: sign only
        combos["signonly"] = np.sign(combos["equal"])

        if variant == "gene_norm":
            Bl = []
            for n in others:
                M = banks[n]["v"][variant][:, col[n]]
                if M.shape[0] > 6000:
                    M = M[np.sort(rng.choice(M.shape[0], 6000, False))]
                Bl.append((ref / sc[n]) * M.astype(np.float64))
            W = Whitener.fit(np.vstack(Bl))
            print(f"fold {held:14} whitener fitted", flush=True)
            # low-rank truncation of the equal-weight prediction, before whitening
            U, S, Vt = np.linalg.svd(combos["equal"], full_matrices=False)
            for r in (50, 200):
                combos[f"rank{r}"] = (U[:, :r] * S[:r]) @ Vt[:r]
        for cname, P in combos.items():
            for g in GAMMAS:
                results.setdefault(f"{variant:12} {cname:10} g={g:g}", {})[held] = \
                    pds(W.apply(P, g), T)

print(f"\n{'config':34} " + "  ".join(f"{k.replace('orion_','')[:9]:>10}" for k in banks)
      + f" {'macro':>8}")
out = {}
for k, v in sorted(results.items()):
    macro = float(np.mean(list(v.values())))
    out[k] = {**v, "macro": macro}
    print(f"{k:34} " + "  ".join(f"{v[n]:10.4f}" for n in banks) + f" {macro:8.4f}")
best = max(out, key=lambda k: out[k]["macro"])
print(f"\nbest: {best} -> {out[best]['macro']:.4f}")
Path("experiments/022-graveyard/results.json").write_text(json.dumps(out, indent=2))
