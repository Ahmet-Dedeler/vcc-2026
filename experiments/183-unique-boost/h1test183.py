"""183 — the same unique-component boost, priced on h1-test's delivered space.

Second surface, exactly 173's pipeline (which is 168d's): a15-era effect matrix `H` recovered
from the 095 artifact, gated at `min_cpm 10`, mass-fixed (`minnorm`), rendered, read back, and
scored against h1-test's own pseudobulk truth.  Donors are `k562_gwps` + `orion_hek293t` only —
`data/donors/h1_2025` covers 85/85 of these targets and is pure oracle leakage; it never enters.

Known-answer row: alpha = 0 must reproduce 168d/173's base `pds` 0.86134.

Caveat carried forward, not buried: h1-test's panel is 85 rows against the challenge's 300, so
its centroid is a noisier estimate of the shared direction AND its `pds` decoy pool is 3.4x
smaller.  Absolute values do not transfer; only the sign and rough shape of the alpha response.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "experiments/168-sharpening"))
sys.path.insert(0, os.path.join(REPO, "experiments/156-delivered-geometry"))
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
sys.path.insert(0, HERE)
from h1_tiebreak import group_profiles, H1, CTRL, TS                   # noqa: E402
from massfix import mass_fix                                           # noqa: E402
from screen import endpoints, render                                   # noqa: E402
from uboost import norm_check, shared_energy, unique_boost             # noqa: E402

ALPHAS = (0.00, 0.25, 0.50, 0.75, 1.00)
BASE_PDS_168D = 0.86134            # 173's reproduced known answer

# ---- truth, exactly as 168d/173 built it ------------------------------------------------
tg, genes_t, acc_t = group_profiles(H1)
ci = int(np.flatnonzero(tg == CTRL)[0])
prof_t = np.log1p(TS * acc_t / acc_t.sum(1, keepdims=True))
b_ctrl = prof_t[ci]
keep_t = np.asarray([i for i in range(len(tg)) if i != ci])
T_all = prof_t[keep_t] - b_ctrl[None, :]
targets = tg[keep_t]

# ---- design effect H --------------------------------------------------------------------
E = np.load(os.path.join(REPO, "experiments/169-sharpening-six-member/work/h1test-effect.npz"),
            allow_pickle=True)
assert np.array_equal(np.asarray([str(g) for g in E["genes"]]),
                      np.asarray([str(g) for g in genes_t])), "gene axes disagree"
tpos = {str(t): i for i, t in enumerate(E["targets"])}
have = np.asarray([i for i, t in enumerate(targets) if str(t) in tpos])
T = T_all[have]
targets = targets[have]
H = E["H"][np.asarray([tpos[str(t)] for t in targets])].astype(np.float64)
gate = E["gate"].astype(bool)
gpos = {str(g): i for i, g in enumerate(genes_t)}
drop = {gpos[str(t)] for t in targets if str(t) in gpos}
mask = np.asarray([i for i in range(len(genes_t)) if i not in drop])

print(f"h1-test: {len(targets)} targets, {int(gate.sum())} gated genes, "
      f"{len(mask)} scored genes")

Hg = np.where(gate[None, :], H, 0.0)
D0 = render(mass_fix(Hg, b_ctrl, gate), b_ctrl)      # the delivered incumbent

worst_norm = 0.0
rows = []
for alpha in ALPHAS:
    Dpost = unique_boost(D0, alpha)
    worst_norm = max(worst_norm, norm_check(D0, Dpost))
    Hp = unique_boost(Hg, alpha)
    worst_norm = max(worst_norm, norm_check(Hg, Hp))
    variants = dict(
        post=Dpost,
        re=render(mass_fix(np.where(gate[None, :], Dpost, 0.0), b_ctrl, gate), b_ctrl),
        pre=render(mass_fix(np.where(gate[None, :], Hp, 0.0), b_ctrl, gate), b_ctrl),
    )
    for place, M in variants.items():
        P, Q = M[:, mask], T[:, mask]
        e = endpoints(P, Q)
        g = shared_energy(P)
        rows.append(dict(alpha=alpha, place=place, **e, **g))
        print(f"  a={alpha:.2f} {place:<5} pds {e['pds']:.6f}  cos {e['cos']:.6f}  "
              f"pair cos {g['mean_pairwise_cos']:+.6f}  centroid E "
              f"{g['centroid_energy_frac']:.5f}  med |h| {g['median_row_norm']:.6f}")

print(f"\nworst |‖h'‖/‖h‖ − 1| = {worst_norm:.3e}  (want < 1e-12)")
assert worst_norm < 1e-12, "the transform is not norm-preserving"

R = {(r["alpha"], r["place"]): r for r in rows}
base = R[(0.00, "pre")]["pds"]
print(f"\nKNOWN-ANSWER: alpha=0 pds {base:.5f} vs 168d/173's {BASE_PDS_168D:.5f}  "
      f"|d| {abs(base - BASE_PDS_168D):.2e}")
assert abs(base - BASE_PDS_168D) < 1e-3, "alpha=0 does not reproduce the incumbent"

print(f"\n{'=' * 96}\nh1-test, delivered space, deltas vs alpha=0")
print(f"  {'alpha':>6} {'pds post':>10} {'pds re':>10} {'pds pre':>10} "
      f"{'d pds pre':>10} {'d cos %':>9} {'pair cos':>10} {'centroid E':>11}")
for alpha in ALPHAS:
    p = {pl: R[(alpha, pl)] for pl in ("post", "re", "pre")}
    dc = 100.0 * (p["pre"]["cos"] / R[(0.00, "pre")]["cos"] - 1.0)
    print(f"  {alpha:>6.2f} {p['post']['pds']:10.6f} {p['re']['pds']:10.6f} "
          f"{p['pre']['pds']:10.6f} {p['pre']['pds'] - base:+10.6f} {dc:+9.2f} "
          f"{p['pre']['mean_pairwise_cos']:+10.6f} "
          f"{p['pre']['centroid_energy_frac']:11.5f}")

np.savez_compressed(os.path.join(HERE, "work", "h1test183.npz"),
                    alphas=np.asarray(ALPHAS),
                    pds=np.asarray([[R[(a, pl)]["pds"] for pl in ("post", "re", "pre")]
                                    for a in ALPHAS]),
                    cos=np.asarray([[R[(a, pl)]["cos"] for pl in ("post", "re", "pre")]
                                    for a in ALPHAS]),
                    targets=targets)
