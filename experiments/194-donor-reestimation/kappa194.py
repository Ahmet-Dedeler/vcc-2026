"""194 — recompute 159's donor-table route to kappa^2 under each estimator.

159 §3.3 route 2:

    tau = kappa^2 = cos_same_target(K562-gwps <-> Orion)  /  cos_same_target(K562-gwps
                                                              <-> K562-essential)

The denominator is the SAME cell line measured in a DIFFERENT library, so it carries the
two screens' noise attenuation but no line difference; the numerator adds the line
difference on top.  159 measured 0.0897 / 0.2627 = 0.341 on the shipped tables, which is
the known-answer row here.

The pair statistic is 159's own `donorpairs.pair`, imported rather than reimplemented, so
the protocol (3,000-target cap, gene-axis intersection, shared-target coordinates dropped,
residual space, permutation null) is identical by construction.  Only the tables differ.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/159-transfer-ceiling"))

import donorpairs as DP  # noqa: E402

PAIRS = (("k562_gwps", "k562_essential"),          # same line, different library
         ("k562_gwps", "orion_hct116"),
         ("k562_gwps", "orion_hek293t"),
         ("orion_hct116", "orion_hek293t"))


def load_table(donor: str, variant: str | None) -> dict | None:
    p = (ROOT / "data/donors" / donor) if variant is None else \
        (WORK / "tables" / donor / variant.replace("@", "_at_"))
    if not (p / "delta.npy").exists():
        return None
    delta = np.load(p / "delta.npy")
    perts = [l.strip() for l in (p / "perts.txt").read_text().splitlines() if l.strip()]
    genes = [l.strip() for l in (p / "genes.txt").read_text().splitlines() if l.strip()]
    assert delta.shape == (len(perts), len(genes)), (donor, variant, delta.shape)
    return dict(delta=np.nan_to_num(delta.astype(np.float32)), perts=perts, genes=genes)


def run(variant: str | None) -> dict:
    rng = np.random.default_rng(DP.SEED)
    tabs = {d: load_table(d, variant) for d in
            {x for p in PAIRS for x in p}}
    out: dict[str, dict] = {}
    for a, b in PAIRS:
        if tabs[a] is None or tabs[b] is None:
            continue
        r = DP.pair(tabs[a], tabs[b], rng)
        if r is None:
            continue
        out[f"{a}|{b}"] = dict(n_perts=r["n_perts"], n_genes=r["n_genes"],
                               resid_same=r["resid"]["same"],
                               resid_diff=r["resid"]["diff"],
                               resid_same_permuted=r["resid"]["same_permuted"],
                               full_same=r["full"]["same"],
                               full_diff=r["full"]["diff"],
                               within_A=r["resid"]["within_A"],
                               within_B=r["resid"]["within_B"])
    den = out.get("k562_gwps|k562_essential", {}).get("resid_same")
    for key in ("k562_gwps|orion_hct116", "k562_gwps|orion_hek293t"):
        if den and key in out:
            out[key]["tau_route2"] = out[key]["resid_same"] / den
    return dict(variant=variant or "SHIPPED", pairs=out,
                tau_mean=float(np.mean([v["tau_route2"] for v in out.values()
                                        if "tau_route2" in v]))
                if den else float("nan"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", nargs="*", default=["v0"])
    ap.add_argument("--shipped", action="store_true")
    ap.add_argument("--out", default=str(WORK / "kappa.json"))
    args = ap.parse_args()
    rows = []
    if args.shipped:
        rows.append(run(None))
    for v in args.variants:
        rows.append(run(v))
    for r in rows:
        print(f"\n=== {r['variant']}  tau_mean {r['tau_mean']:.4f} "
              f"(kappa {np.sqrt(max(r['tau_mean'], 0)):.4f})")
        for k, v in r["pairs"].items():
            print(f"  {k:<32} n={v['n_perts']:>5} g={v['n_genes']:>6} "
                  f"same {v['resid_same']:+.4f} diff {v['resid_diff']:+.4f} "
                  f"perm {v['resid_same_permuted']:+.4f}"
                  + (f"  tau {v['tau_route2']:.4f}" if "tau_route2" in v else ""))
    Path(args.out).write_text(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
