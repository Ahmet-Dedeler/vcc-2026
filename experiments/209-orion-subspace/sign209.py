"""209 — the ORDERING channel: consensus sign accuracy on the h1-test significance-selected
pool, so the subspace arm is directly comparable to 175's rank-1 arm (0.5697 -> 0.5770).

Same instrument as 175's addendum, unchanged: `Surface("164/de164.a15.npz")`, pool =
`real_sig & adj`, donor rows from `DonorHandle.rows` mean-centred per gene exactly as
production's `Donor.decompose`, no scaling, no whitening.  The only thing that changes is
the projection applied to the two Orion donors' rows.

The subspace basis is estimated on the 18,533-gene OFFICIAL axis (`subspace209.py`); the
surface has its own gene axis, so the basis is restricted to the shared symbols and
re-orthonormalised.  The retained norm per direction is reported -- if it were small the
comparison would be meaningless.

TWO consensuses are reported, and the difference between them is the mandate's point:
  `all3`     k562 + hct116 + hek293t, 175's own object and the one the known-answer row is on
  `h1surf`   k562 + hek293t ONLY -- the donors `h1-test` actually carries (AGENTS.md §6.4
             trap 1).  **That surface holds out HCT116, one of the two lines whose shared
             excess is being removed, so it cannot see the mechanism** -- it can only see a
             2-D projection removed from a single donor.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/153-de-bench"))
sys.path.insert(0, str(ROOT / "experiments/172-sign-confidence-key"))
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
from plolo import OFFICIAL  # noqa: E402
from debench import Surface  # noqa: E402
from donors172 import DonorHandle  # noqa: E402

SEED = 20260903
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
ORION = ("orion_hct116", "orion_hek293t")


def basis_on_axis(name: str, r: int | None, genes: np.ndarray) -> tuple[np.ndarray, dict]:
    """Restrict an official-axis basis to a surface gene axis and re-orthonormalise."""
    if name == "blob2":
        z = np.load(WORK / "blob-axis.npz")
        V = np.asarray(z["basis"], dtype=np.float64)
        cols = np.arange(len(OFFICIAL))
    else:
        z = np.load(WORK / f"subspace-{name}.npz")
        V = np.asarray(z["basis"], dtype=np.float64)
        cols = np.asarray(z["cols"])
        if r is not None:
            V = V[:r]
    full = np.zeros((V.shape[0], len(OFFICIAL)))
    full[:, cols] = V
    pos = {g: i for i, g in enumerate(OFFICIAL)}
    take = np.asarray([pos.get(str(g), -1) for g in genes])
    keep = take >= 0
    sub = np.zeros((full.shape[0], len(genes)))
    sub[:, keep] = full[:, take[keep]]
    retained = np.linalg.norm(sub, axis=1)
    Q, _ = np.linalg.qr(sub.T)
    return Q.T, dict(retained_norm=retained.tolist(), n_shared_genes=int(keep.sum()))


def random_basis(dim: int, source: str, seed: int, genes: np.ndarray) -> np.ndarray:
    z = np.load(WORK / f"subspace-{source}.npz")
    cols = np.asarray(z["cols"])
    keepset = {str(OFFICIAL[c]) for c in cols}
    mask = np.asarray([str(g) in keepset for g in genes])
    rng = np.random.default_rng(seed)
    G = rng.standard_normal((int(mask.sum()), dim))
    Q, _ = np.linalg.qr(G)
    out = np.zeros((dim, len(genes)))
    out[:, mask] = Q.T
    return out


def main() -> None:
    S = Surface(str(ROOT / "experiments/164-callset-artefact/de164.a15.npz"))
    pool = S.real_sig & S.adj
    sT = np.sign(S.real_lfc)
    genes = np.asarray([str(g) for g in S.genes])
    targets = [str(p) for p in S.perts]

    D, own_pc = {}, {}
    rng = np.random.default_rng(0)                      # 175's seed, unchanged
    for n in DONORS:
        h = DonorHandle(n)
        Sd, carried = h.rows(targets, genes)
        D[n] = (Sd.astype(np.float64), carried)
        tcols, dcols = h.project(genes)
        sub = np.sort(rng.choice(h.n_perts, min(4000, h.n_perts), replace=False))
        M = np.zeros((sub.size, len(genes)))
        for r, i in enumerate(sub):
            M[r, tcols] = (np.asarray(h.delta[i], dtype=np.float64) - h.m)[dcols]
        _, _, Vt = np.linalg.svd(M, full_matrices=False)
        own_pc[n] = Vt[:1]
        del M
    allc = np.all([c for _, c in D.values()], axis=0)
    m3 = pool & allc[:, None] & np.all([np.sign(Sd) != 0 for Sd, _ in D.values()], axis=0)
    m2 = pool & (D["k562_gwps"][1] & D["orion_hek293t"][1])[:, None] \
        & (np.sign(D["k562_gwps"][0]) != 0) & (np.sign(D["orion_hek293t"][0]) != 0)

    def acc(key, m):
        return float(np.mean(key[m] == sT[m]))

    bases: dict[str, tuple] = {}
    for tag, r in (("C", 2), ("C", 1), ("C", 5), ("WIDE", 2), ("CblobDeflated", 3)):
        bases[f"{tag}_r{r}"] = basis_on_axis(tag, r, genes)
    bases["blob2"] = basis_on_axis("blob2", None, genes)
    for s in (1, 2, 3):
        bases[f"perm_s{s}"] = (random_basis(2, "C", SEED + s, genes),
                               dict(retained_norm=[1.0, 1.0]))

    arms = [("identity", None, 0.0, ())]
    for lam in (0.25, 0.5, 0.75, 1.0):
        arms.append((f"subC_r2_l{lam:g}", "C_r2", lam, ORION))
    for key in ("C_r1", "C_r5", "WIDE_r2", "CblobDeflated_r3", "blob2"):
        arms.append((f"{key}_l1", key, 1.0, ORION))
    for s in (1, 2, 3):
        arms.append((f"permC_s{s}_l1", f"perm_s{s}", 1.0, ORION))

    out = {"n_pool_all3": int(m3.sum()), "n_pool_h1surf": int(m2.sum()), "arms": {}}
    print(f"pool: all-three {int(m3.sum())}, h1-surface donors (k562+hek293t) "
          f"{int(m2.sum())}", flush=True)
    for key in bases:
        print(f"  basis {key:18s} retained norms "
              + " ".join(f"{v:.3f}" for v in bases[key][1]["retained_norm"][:4]), flush=True)

    # 175's own rank-1 rows, as the known-answer comparator
    for k in (0, 1):
        defl = {n: (Sd - (Sd @ own_pc[n][:k].T) @ own_pc[n][:k] if k else Sd)
                for n, (Sd, _) in D.items()}
        cons3 = sum(defl.values()) / 3
        cons2 = (defl["k562_gwps"] + defl["orion_hek293t"]) / 2
        out["arms"][f"rank1_own_k{k}"] = dict(
            all3=acc(np.sign(cons3), m3), h1surf=acc(np.sign(cons2), m2),
            note="175's per-donor own-top-PC deflation, all three donors")

    for name, bkey, lam, donors in arms:
        defl = {}
        for n, (Sd, _) in D.items():
            if lam > 0 and n in donors and bkey is not None:
                V = bases[bkey][0]
                defl[n] = Sd - lam * ((Sd @ V.T) @ V)
            else:
                defl[n] = Sd
        cons3 = sum(defl.values()) / 3
        cons2 = (defl["k562_gwps"] + defl["orion_hek293t"]) / 2
        out["arms"][name] = dict(all3=acc(np.sign(cons3), m3),
                                 h1surf=acc(np.sign(cons2), m2),
                                 per_donor={n: acc(np.sign(defl[n]), m3) for n in D})

    print(f"\n{'arm':22s} {'all3':>8s} {'h1surf':>8s} {'K562':>8s} {'HCT116':>8s} "
          f"{'HEK293T':>8s}")
    base3 = out["arms"]["identity"]["all3"]
    base2 = out["arms"]["identity"]["h1surf"]
    for name in out["arms"]:
        r = out["arms"][name]
        pd = r.get("per_donor", {})
        print(f"{name:22s} {r['all3']:8.4f} {r['h1surf']:8.4f} "
              + " ".join(f"{pd.get(n, float('nan')):8.4f}" for n in DONORS)
              + f"   d_all3 {r['all3'] - base3:+.4f}  d_h1surf {r['h1surf'] - base2:+.4f}")
    out["baseline"] = dict(all3=base3, h1surf=base2)
    out["basis_audit"] = {k: v[1] for k, v in bases.items()}
    (WORK / "sign209.json").write_text(json.dumps(out, indent=2))
    print("\nwrote work/sign209.json")


if __name__ == "__main__":
    main()
