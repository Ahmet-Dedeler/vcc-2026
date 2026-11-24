"""Experiment 010 — leave-one-donor-out validation, scored as PDS.

Why the previous validation stack was wrong, in two independent ways:

**Wrong target biology.** Everything was tuned on K562 -> H1 hESC. H1 is embryonic; the 2026
contexts are six lines "from different tissues", almost certainly somatic. A +20% transfer
cosine on H1 produced a *negative* pds move on A/B/C.

**Wrong objective.** Mean cosine and 300-way PDS ranking are not the same thing. Averaging
donors smooths every prediction toward the shared response, which raises mean alignment while
*reducing* discrimination between perturbations — and discrimination is exactly what PDS
measures. A change can improve cosine and lower pds at the same time. That is consistent with
what v3 did.

This harness fixes both. Each donor line is held out in turn and predicted from the others,
which is somatic-to-somatic transfer across different tissues — structurally the same task as
A/B/C — and it is scored with PDS itself, implemented to the vcc2026 spec:

  - profiles b = log1p(TS * P / sum P), effects delta = b_p - b_ctrl, per line
  - cosine distance, `exclusion_scope: panel` — ALL panel target genes removed from the
    feature space, not just the perturbation's own
  - PDS_p = 1 - k_p/(n-1) with midrank tie handling
  - each held-out line scored separately, then macro-averaged

The held-out line contributes nothing to weights, normalisation or shrinkage.

    uv run python experiments/010-lodo/lodo.py --panel 300
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vcc2026.donor_transfer import Donor, eb_shrink, _norm_scale  # noqa: E402

DONORS = {
    "k562_gwps": "data/donors/k562_gwps",
    "orion_hct116": "data/donors/orion_hct116",
    "orion_hek293t": "data/donors/orion_hek293t",
}


def pds_cosine(pred: np.ndarray, true: np.ndarray) -> np.ndarray:
    """PDS per the vcc2026 spec: 1 - normalised rank of the true match, midrank ties.

    `pred[i]` and `true[i]` are the predicted and measured effect for the same perturbation,
    already restricted to the scored feature space (panel target genes removed).
    """
    P = pred / (np.linalg.norm(pred, axis=1, keepdims=True) + 1e-12)
    T = true / (np.linalg.norm(true, axis=1, keepdims=True) + 1e-12)
    D = 1.0 - P @ T.T                      # (n_pred, n_true) cosine distances
    n = D.shape[1]
    self_d = np.diag(D)[:, None]
    lt = (D < self_d).sum(1)
    eq = (D == self_d).sum(1)
    k = lt + 0.5 * (eq - 1)                # midrank: a tied block shares its span's average
    return 1.0 - k / (n - 1)


def load(treat: str, cells: dict) -> dict:
    out = {}
    for name, path in DONORS.items():
        if not Path(path).exists():
            continue
        d = Donor.load(path)
        m, s = d.decompose()
        if treat == "eb":
            s = eb_shrink(s, d.b_ctrl, cells[name])
        elif treat == "gene_norm":
            s = (s / _norm_scale(s)).astype(np.float32)
        gl = [str(g) for g in d.genes]
        out[name] = dict(s=s, genes=np.asarray(gl, dtype=object),
                         gpos={g: i for i, g in enumerate(gl)},   # O(1) gene lookup
                         gset=set(gl),
                         idx={p: i for i, p in enumerate(d.perts)}, b_ctrl=d.b_ctrl)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", type=int, default=300)
    ap.add_argument("--treat", nargs="+", default=["raw", "gene_norm", "eb"])
    ap.add_argument("--weights", nargs="+", default=["1:1", "2:1", "1:0", "0:1"],
                    help="relative weights of the two remaining donors, in DONORS order")
    ap.add_argument("--basal", action="store_true",
                    help="also test weighting donors by basal-profile similarity to the "
                         "held-out line -- the only context signal the Challenge provides")
    ap.add_argument("--basal-temp", type=float, nargs="+", default=[5.0, 10.0, 20.0])
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    cells = {n: json.loads((Path(p) / "meta.json").read_text())["cells_per_pert"]["median"]
             for n, p in DONORS.items() if Path(p).exists()}
    print(f"donors: {list(cells)}  cells/pert {cells}")

    # Panel: the real VCC targets, restricted to those every line carries so the held-out
    # ranking problem has the same size for each fold.
    vcc = [g for g in pd.read_csv("data/challenge/pert_counts.csv")["target_gene"]]
    banks = load("raw", cells)
    common = [g for g in vcc if all(g in b["idx"] for b in banks.values())]
    rng = np.random.default_rng(args.seed)
    if len(common) > args.panel:
        common = list(np.asarray(common)[np.sort(rng.choice(len(common), args.panel, False))])
    print(f"panel: {len(common)} perturbations present in all {len(banks)} lines\n")

    results: dict[str, dict[str, float]] = {}
    for treat in args.treat:
        banks = load(treat, cells)
        specs = list(args.weights)
        if args.basal:
            specs += [f"basal{t:g}" for t in args.basal_temp]
        for wspec in specs:
            per_fold = {}
            for held in banks:
                others = [n for n in banks if n != held]
                if wspec.startswith("basal"):
                    # Weight each predictor by how similar its basal control profile is to
                    # the held-out line's. This is the only context signal the Challenge
                    # gives us, so if it works here it is usable on A/B/C.
                    temp = float(wspec[5:])
                    sh = sorted(set.intersection(*[banks[n]["gset"] for n in [held] + others]))
                    cx = {n: np.fromiter((banks[n]["gpos"][g] for g in sh), np.int64, len(sh))
                          for n in [held] + others}
                    tb = banks[held]["b_ctrl"][cx[held]]
                    sims = {}
                    for n in others:
                        ob = banks[n]["b_ctrl"][cx[n]]
                        sims[n] = float(ob @ tb / (np.linalg.norm(ob) * np.linalg.norm(tb)))
                    mx = max(sims.values())
                    w = {n: float(np.exp(temp * (sims[n] - mx))) for n in others}
                else:
                    ws = [float(x) for x in wspec.split(":")]
                    if len(ws) != len(others):
                        continue
                    w = dict(zip(others, ws))
                if sum(w.values()) == 0:
                    continue

                # Shared feature space: genes the held-out line and every predictor carry.
                shared = sorted(set.intersection(*[banks[n]["gset"] for n in [held] + others]))
                # exclusion_scope: panel -- drop EVERY panel target gene from the space
                panel_set = set(common)
                feat = [g for g in shared if g not in panel_set]
                # O(1) dict lookups; the previous list(...).index(g) was O(n_genes) per gene
                # and rebuilt the list each time -- ~270M ops per fold.
                col = {n: np.fromiter((banks[n]["gpos"][g] for g in feat), np.int64, len(feat))
                       for n in [held] + others}

                T = np.vstack([banks[held]["s"][banks[held]["idx"][p]][col[held]]
                               for p in common])
                # Equalise donor magnitude so weights mean what they say, using ONLY the
                # predictor lines -- the held-out line must not touch normalisation.
                scale = {n: float(np.median(np.linalg.norm(banks[n]["s"][:, col[n]], axis=1)))
                         or 1.0 for n in others}
                ref = float(np.median([scale[n] for n in others]))
                tot = sum(w.values())
                P = np.zeros_like(T)
                for n in others:
                    P += (w[n] / tot) * (ref / scale[n]) * np.vstack(
                        [banks[n]["s"][banks[n]["idx"][p]][col[n]] for p in common])
                per_fold[held] = float(pds_cosine(P, T).mean())

            if per_fold:
                macro = float(np.mean(list(per_fold.values())))
                key = f"{treat} w={wspec}"
                results[key] = {**per_fold, "macro": macro}
                folds = "  ".join(f"{k.replace('orion_','')[:9]}={v:.4f}"
                                  for k, v in per_fold.items())
                print(f"  {key:22} macro PDS {macro:.4f}   [{folds}]", flush=True)

    best = max(results, key=lambda k: results[k]["macro"])
    print(f"\nbest: {best} -> macro PDS {results[best]['macro']:.4f}")
    print("(PDS 0.5 = no discriminating information; competition raw pds for v1 was 0.7383)")
    Path("experiments/010-lodo").mkdir(parents=True, exist_ok=True)
    Path("experiments/010-lodo/results.json").write_text(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
