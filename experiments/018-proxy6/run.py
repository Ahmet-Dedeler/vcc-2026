"""018 driver — score whitening (and any other effect transform) on all six metrics.

Experiment 017 found spectral flattening worth +0.047 macro PDS, still climbing at the edge
of the grid. PDS is one metric of six, and whitening changes the *relative magnitudes* of
genes inside each prediction, which is exactly what the four DE metrics read. So the
question this answers is not "does whitening help pds" -- 017 settled that -- but "what does
it cost everywhere else".

    uv run python experiments/018-proxy6/run.py --gammas 0 0.5 1 1.5 2 3
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
sys.path.insert(0, str(ROOT / "experiments" / "018-proxy6"))
from vcc2026.donor_transfer import Donor, _norm_scale        # noqa: E402
from vcc2026.whiten import Whitener                          # noqa: E402
from proxy6 import pds_cosine, topk_mask, nmae_curve, de_metrics, scale  # noqa: E402

DONORS = {
    "k562_gwps": "data/donors/k562_gwps",
    "orion_hct116": "data/donors/orion_hct116",
    "orion_hek293t": "data/donors/orion_hek293t",
}


def build_folds(weights: str, treat: str = "gene_norm", basis_perts: int = 6000) -> dict:
    banks = {}
    for name, path in DONORS.items():
        d = Donor.load(path)
        _, s = d.decompose()
        if treat == "gene_norm":
            s = (s / _norm_scale(s)).astype(np.float32)
        gl = [str(g) for g in d.genes]
        banks[name] = dict(s=s, gpos={g: i for i, g in enumerate(gl)}, gset=set(gl),
                           idx={p: i for i, p in enumerate(d.perts)})
    # Truth is always the raw held-out delta -- the scorer does not normalise the reference.
    raw = {n: Donor.load(p).decompose()[1] for n, p in DONORS.items()}

    vcc = list(pd.read_csv("data/challenge/pert_counts.csv")["target_gene"])
    common = [g for g in vcc if all(g in b["idx"] for b in banks.values())]
    folds = {}
    for held in banks:
        others = [n for n in banks if n != held]
        shared = sorted(set.intersection(*[banks[n]["gset"] for n in [held] + others]))
        feat = [g for g in shared if g not in set(common)]
        col = {n: np.fromiter((banks[n]["gpos"][g] for g in feat), np.int64, len(feat))
               for n in [held] + others}
        T = np.vstack([raw[held][banks[held]["idx"][p]][col[held]] for p in common])
        ws = dict(zip(others, [float(x) for x in weights.split(":")]))
        sc = {n: float(np.median(np.linalg.norm(banks[n]["s"][:, col[n]], axis=1))) or 1.0
              for n in others}
        ref = float(np.median(list(sc.values())))
        parts = {}
        for n in others:
            parts[n] = (ref / sc[n]) * np.vstack(
                [banks[n]["s"][banks[n]["idx"][p]][col[n]] for p in common]).astype(np.float64)
        P = sum((ws[n] / sum(ws.values())) * parts[n] for n in others)
        # Donor-basis whitener, exactly as experiments 017/019 validated it: fitted on the
        # predictor lines' full perturbation catalogue, blind to the held-out line.
        rng = np.random.default_rng(0)
        Bl = []
        for n in others:
            M = banks[n]["s"][:, col[n]]
            if M.shape[0] > basis_perts:
                M = M[np.sort(rng.choice(M.shape[0], basis_perts, False))]
            Bl.append((ref / sc[n]) * M.astype(np.float64))
        folds[held] = dict(P=P, T=T.astype(np.float64), parts=parts, others=others,
                           W=Whitener.fit(np.vstack(Bl)))
    return folds, common





def geo_conf(parts: dict, others: list) -> np.ndarray:
    """Experiment 013's winning per-gene confidence: sqrt(|A|*|B|) zeroed on sign disagreement."""
    A, B = parts[others[0]], parts[others[1]]
    return np.sqrt(np.abs(A) * np.abs(B)) * (np.sign(A) == np.sign(B))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--weights", default="1:1")
    ap.add_argument("--gammas", type=float, nargs="+", default=[0.0, 0.5, 1.0, 1.5, 2.0, 3.0])
    ap.add_argument("--k-real", type=int, nargs="+", default=[500])
    ap.add_argument("--k-pred-mult", type=float, nargs="+", default=[1.0])
    ap.add_argument("--alphas", type=float, nargs="+",
                    default=[0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0])
    ap.add_argument("--call-conf", default="absP", choices=["absP", "geo"],
                    help="signal that decides WHICH genes we call DE. 'absP' is what a "
                         "minimum-variance emitter produces implicitly; 'geo' is exp 013's "
                         "sqrt(|A||B|) x sign-agreement, deliverable via per-gene emission "
                         "variance without touching the mean profile.")
    args = ap.parse_args()

    folds, common = build_folds(args.weights)
    print(f"panel {len(common)}  weights {args.weights}  "
          f"P{list(folds.values())[0]['P'].shape}\n")
    alphas = np.array(args.alphas)
    out = {}

    hdr = (f"{'config':26} {'pds':>7} {'nmae*':>7} {'a*':>5} {'fid':>7} {'jac':>7} "
           f"{'reach':>7} {'rch+':>7} | {'s_pds':>6} {'s_nm':>6} {'s_fid':>6} "
           f"{'s_jac':>6} {'s_rch':>6} {'avg5':>7}")
    print(hdr)
    for K in args.k_real:
        for mult in args.k_pred_mult:
            for g in args.gammas:
                acc = {k: [] for k in ("pds", "nmae", "astar", "fid", "jac", "reach", "rchp")}
                for h, f in folds.items():
                    P = f["W"].apply(f["P"], g)
                    T = f["T"]
                    R = topk_mask(T, K)
                    acc["pds"].append(pds_cosine(P, T))
                    # Normalise P to unit median row-norm so alpha means the same thing
                    # across gammas -- whitening changes the overall scale by orders of
                    # magnitude and an un-normalised alpha grid would not be comparable.
                    Pn = P / (np.median(np.linalg.norm(P, axis=1)) + 1e-12)
                    Tn = T / (np.median(np.linalg.norm(T, axis=1)) + 1e-12)
                    cur = nmae_curve(Pn, Tn, R, alphas)
                    j = int(np.argmin(cur))
                    acc["nmae"].append(cur[j]); acc["astar"].append(alphas[j])
                    gc = geo_conf(f["parts"], f["others"])
                    call = P if args.call_conf == "absP" else np.sign(P) * gc
                    m = de_metrics(call, T, R, max(1, int(K * mult)))
                    acc["fid"].append(m["fid"]); acc["jac"].append(m["jac"])
                    acc["reach"].append(m["reach"])
                    m2 = de_metrics(call, T, R, max(1, int(K * mult)), conf=gc)
                    acc["rchp"].append(m2["reach"])
                v = {k: float(np.mean(x)) for k, x in acc.items()}
                s = dict(pds=scale("pds", v["pds"]), nmae=scale("nmae", v["nmae"]),
                         fid=scale("fid", v["fid"]), jac=scale("jac", v["jac"]),
                         reach=scale("reach", v["rchp"]))
                avg5 = float(np.mean(list(s.values())))
                key = f"K={K} kp={mult:g} g={g:g} {args.call_conf}"
                out[key] = {**v, **{f"s_{k}": x for k, x in s.items()}, "avg5": avg5}
                print(f"{key:26.26} {v['pds']:7.4f} {v['nmae']:7.4f} {v['astar']:5.2f} "
                      f"{v['fid']:7.4f} {v['jac']:7.4f} {v['reach']:7.4f} {v['rchp']:7.4f} | "
                      f"{s['pds']:6.3f} {s['nmae']:6.3f} {s['fid']:6.3f} {s['jac']:6.3f} "
                      f"{s['reach']:6.3f} {avg5:7.4f}", flush=True)

    Path("experiments/018-proxy6/results.json").write_text(json.dumps(out, indent=2))
    best = max(out, key=lambda k: out[k]["avg5"])
    print(f"\nbest avg5: {best} -> {out[best]['avg5']:.4f}")


if __name__ == "__main__":
    main()
