#!/usr/bin/env python3
"""362 step 1 — build modified designs + design-stage weighted sign accuracy.

Reads (read-only): stage345.npz, effects.npz::cap15, jurkat_common.npz.
Writes into THIS directory only: work/designs362.npz, work/design_report.json.
No emission, no scorer. Fast (<1 min), low memory.
"""
from __future__ import annotations

import json
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)

STAGE = ROOT / "experiments/345-common-program-structure/work/stage345.npz"
EFF = ROOT / "experiments/323-current-family-validation/work/effects.npz"
JUR = ROOT / "experiments/350-same-line-priors-data/work/jurkat_common.npz"

st = np.load(STAGE, allow_pickle=True)
Cmean = np.asarray(st["Cmean"], dtype=np.float64)
Dfull = np.asarray(st["Dfull"], dtype=np.float64)
gate = np.asarray(st["gate"], dtype=bool)
elig = np.asarray(st["elig"], dtype=np.int64)
R = np.asarray(st["R"], dtype=np.float64)
genes_uni = np.asarray(st["genes"]).astype(str)
targets = np.asarray(st["targets"]).astype(str)
T, G = Dfull.shape
assert Cmean.shape == (T, G) and R.shape == (T, G) and len(elig) == 58

e = np.load(EFF, allow_pickle=True)
cap15 = np.asarray(e["cap15"], dtype=np.float64)
genes_eff = np.asarray(e["genes"]).astype(str)
targets_eff = np.asarray(e["targets"]).astype(str)
assert (targets_eff == targets).all(), "target order drift effects vs stage345"
assert cap15.shape == (len(targets), len(genes_eff))

j = np.load(JUR, allow_pickle=True)
genes_j = np.asarray(j["genes"]).astype(str)
vec_j = np.asarray(j["vector"], dtype=np.float64)
assert vec_j.shape == genes_j.shape

# ---- universe -> effects mapping by name ----
pos_eff = {g: i for i, g in enumerate(genes_eff)}
col = np.array([pos_eff.get(g, -1) for g in genes_uni])
miss = col < 0
missing_genes = sorted(genes_uni[miss].tolist())
print(f"universe {G}, mapped {int((~miss).sum())}, missing {int(miss.sum())}: {missing_genes}")

# sanity: Dfull must equal cap15 mapped (345's construction)
check = cap15[:, np.where(miss, 0, col)].copy()
check[:, miss] = 0.0
dmax = float(np.max(np.abs(check - Dfull)))
print(f"max|mapped cap15 - Dfull| = {dmax:.3e}")
assert dmax == 0.0, "Dfull is not the name-mapped cap15"

# ---- Jurkat vector onto the universe by name ----
pos_j = {g: i for i, g in enumerate(genes_j)}
Juni = np.array([vec_j[pos_j[g]] if g in pos_j else 0.0 for g in genes_uni])
print(f"Jurkat universe coverage: {int((Juni != 0).sum())}/{G} "
      f"(unmapped genes get 0)")


def pooled_scale(D, C):
    # 345's exact rule (run345b.py::pooled_scale): median over ALL pooled eligible
    # gate pairs, zeros included (Cmean is 0 where no other target has gate signal;
    # transfer vectors are 0 where uncovered). No zero filtering.
    ad = np.concatenate([np.abs(D[ti, gate[ti]]) for ti in elig])
    ac = np.concatenate([np.abs(C[ti, gate[ti]]) if C.ndim == 2
                         else np.abs(C[gate[ti]]) for ti in elig])
    return float(np.median(ad) / np.median(ac)), float(np.median(ad))


s_oracle, medD = pooled_scale(Dfull, Cmean)
medC = float(np.median(np.concatenate([np.abs(Cmean[ti, gate[ti]]) for ti in elig])))
print(f"s_oracle = {s_oracle:.4f} (345 recorded 2.7549), medD={medD:.4f} medC={medC:.4f}")

J2 = np.repeat(Juni[None, :], T, axis=0)
s_generic, _ = pooled_scale(Dfull, J2)
medJ = float(np.median(np.concatenate([np.abs(Juni[gate[ti]]) for ti in elig])))
print(f"s_generic = {s_generic:.4f}, medJ Pool={medJ:.4f}")

# ---- arms on the universe ----
rng = np.random.default_rng(362)
pi = rng.permutation(G)
Cshuf = Cmean[:, pi]

arms_uni = {
    "cap15": Dfull.copy(),
    "oracle_a05": Dfull + 0.5 * s_oracle * Cmean,
    "oracle_a1": Dfull + 1.0 * s_oracle * Cmean,
    "generic_a1": Dfull + 1.0 * s_generic * J2,
    "shuffled_a1": Dfull + 1.0 * s_oracle * Cshuf,
}
n0 = np.linalg.norm(Dfull, axis=1)
n1 = np.linalg.norm(arms_uni["oracle_a1"], axis=1)
print(f"normpreserve factors n0/n1: min {float((n0/n1).min()):.4f} "
      f"median {float(np.median(n0/n1)):.4f} max {float((n0/n1).max()):.4f}")
arms_uni["oracle_a1_normpreserve"] = arms_uni["oracle_a1"] * (n0 / n1)[:, None]


def wacc(sign_mat):
    out = []
    for ti in elig:
        m_ = gate[ti]
        l = R[ti, m_]
        w = np.abs(l)
        ok = (sign_mat[ti, m_] == np.sign(l))
        out.append(float((w * ok).sum() / w.sum()))
    return float(np.mean(out)), [float(v) for v in out]


report = {"s_oracle": s_oracle, "s_generic": s_generic, "medD": medD, "medC": medC,
          "medJ_pool": medJ, "shuffle_seed": 362, "missing_universe_genes": missing_genes,
          "dmax_mapped_check": dmax,
          "expect_345": {"a0": 0.6002, "a05": 0.6527, "a1": 0.6564}}
design_w = {}
for name, U in arms_uni.items():
    w, per = wacc(np.sign(U))
    design_w[name] = {"wacc": w, "per_target": per}
    print(f"design wacc {name}: {w:.4f}")
report["design_wacc"] = design_w

# known-answer check vs 345's grid
assert abs(design_w["cap15"]["wacc"] - 0.6002) < 0.002, design_w["cap15"]["wacc"]
assert abs(design_w["oracle_a05"]["wacc"] - 0.6527) < 0.002, design_w["oracle_a05"]["wacc"]
assert abs(design_w["oracle_a1"]["wacc"] - 0.6564) < 0.002, design_w["oracle_a1"]["wacc"]
print("345 grid reproduced within 0.002 on all three anchors")

# ---- lift to the full effects axis (universe genes only; rest bit-identical) ----
order = ["cap15", "oracle_a05", "oracle_a1", "oracle_a1_normpreserve",
         "generic_a1", "shuffled_a1"]
full = {}
for name in order:
    M = cap15.copy()
    M[:, col[~miss]] = arms_uni[name][:, ~miss]
    full[name] = M

np.savez_compressed(WORK / "designs362.npz", genes=genes_eff, targets=targets, **full)
(WORK / "design_report.json").write_text(json.dumps(report, indent=1) + "\n")
print("wrote work/designs362.npz + work/design_report.json")
