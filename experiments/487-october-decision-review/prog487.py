"""487 decision 6: what the programme/gate evidence licenses on a panel with more design-weak targets (481 F3/F4).

On H1 (467's production-emitter runs, 3 seeds) block-prior arms share their emission draws, so a target's DE members
depend only on the prior it gets (468: exact split). Any gate rule is therefore an exact per-target mix of two runs:
  base     = 435's CLEAN4 run,           ungated = 462's Arc-programme run (h1_arcfull).
KA: 467's gated33 run must equal, target by target, ungated on strength >= 2 and base on strength <= 1.
Per strength class (0 / 1 / >= 2; strength = 467's n_gt_0.2), the mean per-target delta (programme - CLEAN4) of each DE
member, then the expected per-context overall-equivalent (467's spans) on a panel with class shares pi:
  sum_c pi_c * mean_c(delta_member) / span / 6,
for three rules: ungated; le1 (strength <= 1 keeps CLEAN4); eq0 (only strength 0 keeps CLEAN4 = what q33 does once more
than 33 % of targets have strength <= 1: efr, and the Arc-style stand-in in 100 % of 481's draws).
Panels: A/B/C (453 dump, builder), the Arc-style stand-in (481 knnmap: the 2025 H1 panel's non-essential targets), efr.
Targets are the 67 non-essential h1-test targets (the 2026 panel is essential-free) and, for reference, all 85.
Class-stratified target bootstrap (4,000).
usage: uv run python experiments/487-october-decision-review/prog487.py -> work/prog487.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
F = "per-target-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse.csv"
M = {"fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
     "nmae": "de_wilcoxon_lfc_nmae", "jac": "de_wilcoxon_sig_jaccard"}
SPANS = {"fid": 0.30, "reach": 0.90, "nmae": -0.60, "jac": 0.37}
ESS = set((ROOT / "data/donors/k562_essential/perts.txt").read_text().split())
PANELS = {  # class shares (strength 0, 1, >= 2)
    "A/B/C (453 builder)": (0.070, 0.207, 0.723),
    "Arc-style stand-in (481 knn, 2025 H1 non-essential)": (0.135, 0.261, 0.604),
    "efr (465 builder)": (0.227, 0.313, 0.460),
    "H1 non-essential (467)": None,   # its own shares, for the KA of the reweighting
}


def load(d):
    df = pd.concat([pd.read_csv(d / F.format(s=s)).assign(seed=s) for s in (0, 1, 2)])
    return df.pivot_table(index=["seed", "perturbation"], columns="metric", values="value")


def main():
    base = load(ROOT / "experiments/435-zs-on-coarse/work/h1")
    ung = load(ROOT / "experiments/462-broad16-lab-signature/work/h1_arcfull")
    g33 = load(ROOT / "experiments/467-strength-gated-programme/work/h1_gated33")
    st = json.loads((ROOT / "experiments/467-strength-gated-programme/work/strength467.json").read_text())
    order = st["targets_in_order"]
    strength = pd.Series({t: st["per_target"][str(i)]["n_gt_0.2"] for i, t in enumerate(order)})
    targets = sorted(base.loc[0].index)
    s = strength.reindex(targets).values
    cls = np.where(s == 0, 0, np.where(s == 1, 1, 2))
    ess = np.array([t in ESS for t in targets])
    # KA: gated33 == ungated on strength >= 2, == base on strength <= 1 (per target, every member, every seed)
    ka_bad = 0
    for seed in (0, 1, 2):
        B, U, G = base.loc[seed].reindex(targets), ung.loc[seed].reindex(targets), g33.loc[seed].reindex(targets)
        for col in M.values():
            exp = np.where(cls >= 2, U[col].values, B[col].values)
            ka_bad += int(np.sum(~np.isclose(G[col].values, exp, equal_nan=True, rtol=0, atol=1e-12)))
    print(f"KA exact split (gated33 = ungated on >=2, base on <=1): mismatches {ka_bad}")
    out = {"ka_mismatches": ka_bad, "class_counts_all": np.bincount(cls, minlength=3).tolist(),
           "class_counts_nonessential": np.bincount(cls[~ess], minlength=3).tolist()}
    # per-target deltas averaged over seeds
    D = sum((ung.loc[sd].reindex(targets) - base.loc[sd].reindex(targets)) for sd in (0, 1, 2)) / 3
    rng = np.random.default_rng(4876)

    def _m(v):
        v = v[np.isfinite(v)]
        return float(v.mean()) if v.size else 0.0

    def member_means(idx_by_class):
        return {c: {k: _m(D[col].values[ix]) for k, col in M.items()} for c, ix in idx_by_class.items()}

    def value(means, pi, rule):
        on = {"ungated": (1, 1, 1), "le1": (0, 0, 1), "eq0": (0, 1, 1)}[rule]
        return sum(pi[c] * on[c] * means[c][k] / SPANS[k] for c in range(3) for k in M) / 6

    res = {}
    for sub, mask in (("non-essential", ~ess), ("all", np.ones(len(targets), bool))):
        idx = {c: np.flatnonzero(mask & (cls == c)) for c in range(3)}
        means = member_means(idx)
        own = tuple(len(idx[c]) / sum(len(v) for v in idx.values()) for c in range(3))
        panels = {(k if v is not None else f"H1 {sub} (own shares)"): (v if v is not None else own) for k, v in PANELS.items()}
        boots = []
        for _ in range(4000):
            bi = {c: ix[rng.integers(0, len(ix), len(ix))] if len(ix) else ix for c, ix in idx.items()}
            boots.append(member_means(bi))
        r = {"class_n": [len(idx[c]) for c in range(3)],
             "class_means_overall_equiv": {c: round(sum(means[c][k] / SPANS[k] for k in M) / 6, 5) for c in range(3)},
             "class_means_members": {c: {k: round(v, 4) for k, v in means[c].items()} for c in range(3)}}
        for pname, pi in panels.items():
            row = {}
            for rule in ("ungated", "le1", "eq0"):
                v = value(means, pi, rule)
                bs = np.array([value(b, pi, rule) for b in boots])
                row[rule] = {"point": round(v, 5), "ci95": [round(float(np.quantile(bs, .025)), 5), round(float(np.quantile(bs, .975)), 5)],
                             "p_le_0": round(float((bs <= 0).mean()), 3)}
            d_le1_eq0 = np.array([value(b, pi, "le1") - value(b, pi, "eq0") for b in boots])
            row["le1 - eq0"] = {"point": round(value(means, pi, "le1") - value(means, pi, "eq0"), 5),
                                "ci95": [round(float(np.quantile(d_le1_eq0, .025)), 5), round(float(np.quantile(d_le1_eq0, .975)), 5)]}
            r[pname] = row
        res[sub] = r
        print(f"\n== {sub} targets: class n (0/1/>=2) {r['class_n']}; per-class overall-equiv of programme vs CLEAN4: "
              f"{r['class_means_overall_equiv']}")
        for c in range(3):
            print(f"   class {c}: " + " ".join(f"{k} {v:+.4f}" for k, v in r["class_means_members"][c].items()))
        for pname in panels:
            row = r[pname]
            print(f"  {pname:52s} " + "  ".join(f"{rule} {row[rule]['point']:+.4f} [{row[rule]['ci95'][0]:+.4f},{row[rule]['ci95'][1]:+.4f}]"
                                                 for rule in ("ungated", "le1", "eq0"))
                  + f"  le1-eq0 {row['le1 - eq0']['point']:+.4f} [{row['le1 - eq0']['ci95'][0]:+.4f},{row['le1 - eq0']['ci95'][1]:+.4f}]")
    out["results"] = res
    (HERE / "work/prog487.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
