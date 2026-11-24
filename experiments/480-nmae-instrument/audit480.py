"""480 resume-2 audit: recompute RESULTS.md's H1 headline numbers with fresh code (does NOT import common480).

Independent re-implementation of the gate, the nmae formula (metrics BRIEF §7), the donor conversion (inverse of
build_donor_table.py's log1p(5e4 * P / sum P)) and the constructions quoted in RESULTS.md §1, §2, §4.6 and §5. Every
row carries the value RESULTS.md quotes; the script exits non-zero if any row misses its tolerance.

    uv run python audit480.py      # -> work/audit480.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TS = 5e4

# ---- truth: the scorer's gate on h1-test (323's reference table) -------------------------------------------------
ref = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=True)
keep = ref["keep"].astype(bool)
GENES = ref["genes"].astype(str)[keep]
TARGETS = ref["targets"].astype(str)
q, lfc = ref["q"].astype(float), ref["lfc"].astype(float)
gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & (GENES[None, :] != TARGETS[:, None])
scored = gate.sum(1) >= 10
gate &= scored[:, None]


def per_target(pred: np.ndarray) -> np.ndarray:
    """NMAE_p for every scored target (NaN elsewhere); non-finite predictions read as 0."""
    P = np.where(np.isfinite(pred), pred, 0.0)
    num = (np.abs(P - np.where(gate, lfc, 0.0)) * gate).sum(1)
    den = (np.abs(np.where(gate, lfc, 0.0))).sum(1)
    out = np.full(len(TARGETS), np.nan)
    out[scored] = num[scored] / den[scored]
    return out


def score(pred: np.ndarray) -> float:
    return float(np.nanmean(per_target(pred)))


def donor_log2(name: str, pc: float) -> np.ndarray:
    """Same-target donor rows on (TARGETS x GENES): log2((P_pert + pc) / (P_ctrl + pc)), P in pseudobulk CPM.
    NaN where the donor has no row for the target or no column for the gene."""
    base = ROOT / "data/donors" / name
    perts = [s.strip() for s in open(base / "perts.txt")]
    genes = [s.strip() for s in open(base / "genes.txt")]
    prow = {p: i for i, p in enumerate(perts)}
    gcol: dict[str, int] = {}
    for i, g in enumerate(genes):
        gcol.setdefault(g, i)
    col = np.array([gcol.get(g, -1) for g in GENES])
    have = col >= 0
    delta = np.load(base / "delta.npy", mmap_mode="r")
    b = np.load(base / "b_ctrl.npy").astype(float)[col[have]]
    ctrl_cpm = np.expm1(b) / TS * 1e6
    out = np.full((len(TARGETS), len(GENES)), np.nan)
    for k, t in enumerate(TARGETS):
        if t in prow:
            d = np.asarray(delta[prow[t]], float)[col[have]]
            pert_cpm = np.maximum(np.expm1(b + d) / TS * 1e6, 0.0)
            with np.errstate(divide="ignore", invalid="ignore"):
                v = np.log2((pert_cpm + pc) / (ctrl_cpm + pc))
            out[k, have] = np.where(np.isfinite(v), v, np.nan)
    return out


def boot(d: np.ndarray, n: int = 4000, seed: int = 0):
    d = d[np.isfinite(d)]
    bs = d[np.random.default_rng(seed).integers(0, d.size, (n, d.size))].mean(1)
    return [float(d.mean()), float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))]


rows = []


def check(label: str, got: float, want: float, tol: float):
    ok = abs(got - want) <= tol
    rows.append(dict(row=label, quoted=want, recomputed=got, tol=tol, ok=bool(ok)))
    print(f"{'ok ' if ok else 'BAD'} {label:62s} quoted {want:+.5f}  recomputed {got:+.5f}")


print("scored targets", int(scored.sum()), "gate cells", int(gate.sum()))
check("scored targets (58)", float(scored.sum()), 58, 0)
check("zero prediction", score(np.zeros(gate.shape)), 1.0, 1e-12)
check("oracle (truth itself)", score(np.where(gate, lfc, 0.0)), 0.0, 1e-12)

L_base = np.load(ROOT / "experiments/426-mde-magnitude/work/h1/stats-base_s0_bps.npz", allow_pickle=True)
L_crs = np.load(ROOT / "experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz",
                allow_pickle=True)
assert (L_base["targets"].astype(str) == TARGETS).all() and (L_crs["targets"].astype(str) == TARGETS).all()
base = np.nan_to_num(L_base["lfc"].astype(float))
ours = np.nan_to_num(L_crs["lfc"].astype(float))
assert base.shape == gate.shape == ours.shape
check("KA 426 base_s0_bps", score(base), 0.95237, 5e-6)
check("shipped coarse MDE (435)", score(ours), 0.94270, 5e-6)

# §5 oracles on the shipped construction
T = np.where(gate, lfc, 0.0)
sg = np.sign(ours)
check("§5 oracle magnitudes at our signs", score(sg * np.abs(T)), 0.8442, 5e-5)
# RESULTS.md first quoted 0.6891 (a rounding slip; the value is 0.68916): corrected to 0.6892 on 2026-09-26
check("§5 oracle signs at our magnitudes", score(np.sign(T) * np.abs(ours)), 0.6892, 5e-5)
check("§5 abstain where our sign is wrong", score(np.where(sg == np.sign(T), ours, 0.0)), 0.8132, 5e-5)

# §2 donor constructions (scale 1)
D1 = {d: donor_log2(d, 1.0) for d in ["k562_gwps", "orion_hct116", "orion_hek293t"]}
check("§2 (b) K562 scale 1, pc1", score(D1["k562_gwps"]), 1.0052, 5e-5)
check("§2 (b) K562 scale 1, naive20 (delta/ln2)", score(donor_log2("k562_gwps", 20.0)), 0.9705, 5e-5)
check("§2 (c) Orion HCT116 scale 1, pc1", score(D1["orion_hct116"]), 1.0487, 5e-5)
check("§2 (c) Orion HEK293T scale 1, pc1", score(D1["orion_hek293t"]), 1.0171, 5e-5)
with np.errstate(all="ignore"):
    cons = np.nan_to_num(np.nanmean(np.stack(list(D1.values())), 0))
    orion = np.nanmean(np.stack([D1["orion_hct116"], D1["orion_hek293t"]]), 0)
    fam = np.nan_to_num(np.nanmean(np.stack([D1["k562_gwps"], orion]), 0))
check("§2 (c) consensus of 3, scale 1, pc1", score(cons), 0.9632, 5e-5)
check("§2 (e) family-balanced, scale 1, pc1", score(fam), 0.9611, 5e-5)
check("§2 (d) our sign x K562 magnitude, pc1", score(sg * np.abs(np.nan_to_num(D1["k562_gwps"]))), 1.0208, 5e-5)
d_b = boot(per_target(D1["k562_gwps"]) - per_target(ours))
check("§2 Δ (b) K562 pc1 − (a) coarse, all 58", d_b[0], 0.0625, 5e-5)

# §4.6 rider at k 0.2 (mixed) on h1-test
d_r = per_target(ours + 0.2 * cons) - per_target(ours)
check("§4.6 rider ours + 0.2·cons (pc1), Δ", float(np.nanmean(d_r)), -0.0030, 5e-5)
d_r1 = per_target(ours + 0.1 * cons) - per_target(ours)
check("§4.6 rider ours + 0.1·cons (pc1), Δ", float(np.nanmean(d_r1)), -0.0024, 5e-5)

# §5.2 gate-wide |l|-weighted sign accuracy of our shipped construction
w = np.abs(T)
acc = float((w * (sg == np.sign(T))).sum() / w.sum())
check("§5.2 pooled |l|-weighted gate sign accuracy", acc, 0.568, 5e-4)

json.dump(dict(rows=rows, boot_b_minus_a=d_b, rider02_boot=boot(d_r)), open(HERE / "work/audit480.json", "w"),
          indent=1)
bad = [r for r in rows if not r["ok"]]
print(f"{len(rows) - len(bad)}/{len(rows)} rows reproduce")
sys.exit(1 if bad else 0)
