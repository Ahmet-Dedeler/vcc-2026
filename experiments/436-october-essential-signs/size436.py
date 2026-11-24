"""436 (c): October sizing. How many targets of each panel does each essential screen cover, the union, and how many
of the challenge's genes carry a donor sign on a covered target. Panels: Arc's real A/B/C panel; 199's pessimistic
(`def-bundle`, any donor) and realistic (`def-bundle-r`, K562 n Orion n CD4) synthetic D/E/F draws (both carry ~30
essential genes, FINDINGS 8.1g); 238's essential-free draws (`ef-bundle` = efp, `ef-bundle-r` = efr), which match
Arc's panel on both observed properties (disjoint from A/B/C, 0 on the K562-essential list).
Also the chance expectation of A/B/C's union coverage under a random draw from 199's pools. Writes work/size436.json."""
from __future__ import annotations

import csv
import json
from math import comb
from pathlib import Path

import numpy as np

import signs436 as S6

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PANELS = {
    "abc (Arc, real)": ROOT / "data/challenge",
    "199 pess (def-bundle)": ROOT / "experiments/199-def-runbook/work/def-bundle",
    "199 real (def-bundle-r)": ROOT / "experiments/199-def-runbook/work/def-bundle-r",
    "238 efp (ef-bundle)": ROOT / "experiments/238-essential-free-panel/work/ef-bundle",
    "238 efr (ef-bundle-r)": ROOT / "experiments/238-essential-free-panel/work/ef-bundle-r",
}


def read_col(p, col):
    with open(p) as f:
        return [row[col] for row in csv.DictReader(f)]


genes = read_col(ROOT / "data/challenge/gene_names.csv", "gene_name")
ess_list = set(S6._lines(S6.ESS_DIRS["k562_essential"] / "perts.txt"))
screens = {nm: set(S6._lines(d / "perts.txt")) for nm, d in S6.ESS_DIRS.items()}
union = set().union(*screens.values())
screen_genes = {nm: set(S6._lines(d / "genes.txt")) for nm, d in S6.ESS_DIRS.items()}
panel_genes = {nm: set(S6._lines(d / "genes.txt")) for nm, d in S6.PANEL_DIRS.items()}
gset = set(genes)
out = {"union_size": len(union), "union_on_k562ess_list": len(union & ess_list),
       "screen_sizes": {nm: len(v) for nm, v in screens.items()},
       "screen_genes_on_challenge_axis": {nm: len(v & gset) for nm, v in screen_genes.items()},
       "panel_donor_genes_on_challenge_axis": {nm: len(v & gset) for nm, v in panel_genes.items()},
       "panels": {}}
print(f"union of the four essential screens: {len(union)} targets ({len(union & ess_list)} on the K562-essential list)")
print("genes each screen measures on the challenge axis:", out["screen_genes_on_challenge_axis"],
      "| panel donors:", out["panel_donor_genes_on_challenge_axis"])
for nm, d in PANELS.items():
    tg = read_col(d / "pert_counts.csv", "target_gene")
    assert len(tg) == 300, (nm, len(tg))
    row = {nm2: int(sum(t in v for t in tg)) for nm2, v in screens.items()}
    cov_t = [t for t in tg if t in union]
    row["union"] = len(cov_t)
    row["union_on_k562ess_list"] = int(sum(t in ess_list for t in cov_t))
    row["union_not_on_list"] = len(cov_t) - row["union_on_k562ess_list"]
    row["panel_on_k562ess_list"] = int(sum(t in ess_list for t in tg))
    ntab = [sum(t in v for v in screens.values()) for t in cov_t]
    row["screens_per_covered_target"] = {int(k): int(v) for k, v in zip(*np.unique(ntab, return_counts=True))} if ntab else {}
    # genes carrying a donor sign on a covered target: the challenge-axis genes any covering table measures
    # (all mode: the essential screens covering it + any production donor covering it)
    if cov_t:
        s_all, info_all = S6.donor_signs(np.array(cov_t), np.array(genes), mode="all")
        s_ess, _ = S6.donor_signs(np.array(cov_t), np.array(genes), mode="ess")
        row["signed_genes_per_covered_target_all"] = {"median": float(np.median((s_all != 0).sum(1))),
                                                      "min": int((s_all != 0).sum(1).min()), "max": int((s_all != 0).sum(1).max())}
        row["signed_genes_per_covered_target_ess"] = {"median": float(np.median((s_ess != 0).sum(1)))}
        row["covered_not_on_list"] = [t for t in cov_t if t not in ess_list]
    out["panels"][nm] = row
    print(f"{nm:25s} " + " ".join(f"{k} {row[k]:3d}" for k in (*screens, "union")) +
          f" | union on K562-ess list {row['union_on_k562ess_list']}, off it {row['union_not_on_list']}; panel on list "
          f"{row['panel_on_k562ess_list']}; screens/target {row['screens_per_covered_target']}; signed genes/target (all) "
          f"{row.get('signed_genes_per_covered_target_all', {}).get('median', 0):.0f}")

# chance expectation of A/B/C's union coverage (199's pools, from 238's bundle record) and the rule-of-three bound
b238 = json.loads((ROOT / "experiments/238-essential-free-panel/work/bundle238.json").read_text())
pools = {}
for key in ("efp", "efr"):
    rec = b238[key]
    pools[key] = {k: rec.get(k) for k in ("pool", "eligible_with_essential", "eligible_essential_free")}
out["pools_238"] = pools
# A random 300 from a pool of N genes that contains u union genes has E = 300 u / N. The union's overlap with each pool is
# not recorded, so bound it with the union's size (every union gene is a K562-screened, expressed gene).
for N in sorted({v["eligible_with_essential"] for v in pools.values() if v["eligible_with_essential"]}):
    e = 300 * len(union) / N
    out.setdefault("abc_chance", {})[str(N)] = {"expected_union_hits_upper": e, "p_zero_poisson": float(np.exp(-e))}
    print(f"A/B/C: a random 300 of {N} pool genes would carry <= {e:.1f} union targets (P(0) = {np.exp(-e):.1e}); observed 0")
# and on an essential-free pool: only the union's OFF-list genes can land, as they do on 238's efp/efr (9 each)
off = len(union - ess_list)
for key, v in pools.items():
    N = v["eligible_essential_free"]
    if N:
        e = 300 * off / N
        out.setdefault("abc_chance_essfree", {})[key] = {"pool": N, "union_off_list": off, "expected_upper": e,
                                                         "p_zero_poisson": float(np.exp(-e))}
        print(f"A/B/C vs an essential-free pool ({key}, {N} genes): <= {e:.1f} off-list union targets expected "
              f"(P(0) = {np.exp(-e):.1e}); observed 0")
# if Arc drew 600 targets and split them at random 300/300, P(A/B/C gets none of k union genes) = C(600-k,300)/C(600,300)
out["split_600"] = {str(k): comb(600 - k, 300) / comb(600, 300) for k in (1, 3, 5, 9, 15, 30)}
print("random 600 -> 300/300 split, P(A/B/C carries 0 of k union genes):",
      {k: f"{v:.1e}" for k, v in out["split_600"].items()})
# rule of three: 0 of 300 on A/B/C -> 95 % upper bound on the per-target rate 3/300 -> <= 3 covered D/E/F targets
out["rule_of_three_upper_targets"] = 3
(HERE / "work/size436.json").write_text(json.dumps(out, indent=1))
