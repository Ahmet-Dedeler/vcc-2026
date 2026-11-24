"""159 -- reconcile the harness and the donor table, which disagree 100x on the DECOY level.

Same two cell lines (Orion HCT116, Orion HEK293T), same space, same estimator:

  thinned harness   (132 targets, 150 cells, 18,533-gene challenge axis)
      same +0.1881   diff +0.1790   within-A +0.2476   within-B +0.4175
  donor tables      (3,000 targets, all cells, 35,584-gene Orion axis)
      same +0.0306   diff +0.0016   within-A +0.0025   within-B +0.0067

One of those is an artefact.  The candidates are (i) the gene axis, (ii) the target set,
(iii) the harness's SUBSAMPLED control pool, (iv) the harness's 150-cell depth.  This walks
from one to the other one variable at a time, so the answer is not a guess.

`Donor.decompose()`-style centring is NOT applied anywhere here; the donor tables are read as
`delta.npy` off disk.
"""
from __future__ import annotations

import json
import os

import numpy as np

from reliability import load, prof, unit

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")


def offdiag_mean(S):
    n = S.shape[0]
    return float(S[~np.eye(n, dtype=bool)].mean())


def report(tag, A, B, extra=""):
    a, b = unit(A), unit(B)
    S = a @ b.T
    same = float(np.diag(S).mean())
    diff = offdiag_mean(S)
    print(f"  {tag:<52} n={A.shape[0]:4d} g={A.shape[1]:6d}  "
          f"same {same:+.4f}  diff {diff:+.4f}  excess {same - diff:+.5f}  "
          f"ratio {same / diff if abs(diff) > 1e-9 else float('nan'):7.3f}  "
          f"wA {offdiag_mean(a @ a.T):+.4f}  wB {offdiag_mean(b @ b.T):+.4f} {extra}")
    return {"tag": tag, "n": int(A.shape[0]), "g": int(A.shape[1]), "same": same,
            "diff": diff, "excess": same - diff,
            "ratio": same / diff if abs(diff) > 1e-9 else None,
            "within_A": offdiag_mean(a @ a.T), "within_B": offdiag_mean(b @ b.T)}


def donor(rel):
    p = os.path.join(REPO, rel)
    return (np.load(os.path.join(p, "delta.npy")),
            [l.strip() for l in open(os.path.join(p, "perts.txt")) if l.strip()],
            [l.strip() for l in open(os.path.join(p, "genes.txt")) if l.strip()])


def main() -> None:
    rows = []
    hA, hB = load("hct116"), load("hek293t")
    shared = sorted(set(hA["targets"]) & set(hB["targets"]))
    chal = list(hA["genes"])
    print(f"harness pair: {len(shared)} shared targets, {len(chal)} challenge genes")

    dA, pA, gA = donor("data/donors/orion_hct116")
    dB, pB, gB = donor("data/donors/orion_hek293t")
    iA = {p: i for i, p in enumerate(pA)}
    iB = {p: i for i, p in enumerate(pB)}
    jA = {g: i for i, g in enumerate(gA)}
    jB = {g: i for i, g in enumerate(gB)}

    print("\nA. DONOR TABLES, walking the gene axis and the target set to the harness's")
    common_t = [t for t in shared if t in iA and t in iB]
    axes = {
        "orion axis (35,584 shared genes)": [g for g in gA if g in jB],
        "challenge axis (18,533)": [g for g in chal if g in jA and g in jB],
    }
    tsets = {"3000 random targets": None, "the harness's 120 targets": common_t}
    rng = np.random.default_rng(20260902)
    allc = sorted(set(pA) & set(pB))
    rand3000 = [allc[i] for i in sorted(rng.choice(len(allc), 3000, replace=False))]
    tsets["3000 random targets"] = rand3000

    for tlabel, ts in tsets.items():
        for alabel, gs in axes.items():
            drop = set(ts)
            gg = [g for g in gs if g not in drop]
            X = dA[np.asarray([iA[t] for t in ts])][:, np.asarray([jA[g] for g in gg])]
            Y = dB[np.asarray([iB[t] for t in ts])][:, np.asarray([jB[g] for g in gg])]
            rows.append(report(f"donor | {tlabel} | {alabel}",
                               X.astype(np.float64), Y.astype(np.float64)))

    print("\nB. HARNESS, at its own depth and with its own control pool")
    kt, kc = hA["kt"], hA["kc"]
    ia = np.asarray([int(np.flatnonzero(hA["targets"] == s)[0]) for s in shared])
    ib = np.asarray([int(np.flatnonzero(hB["targets"] == s)[0]) for s in shared])
    # scored genes: challenge axis, minus the panel target coordinates, minus anything the
    # donor tables do not carry (so the control-swap arm sees the identical gene set)
    mask = np.asarray([i for i in range(len(chal))
                       if chal[i] not in set(shared) and chal[i] in jA and chal[i] in jB])

    cA = prof(hA["csum"].sum(0))[mask]
    cB = prof(hB["csum"].sum(0))[mask]
    EA = prof(hA["tsum"][ia].sum(1))[:, mask] - cA[None, :]
    EB = prof(hB["tsum"][ib].sum(1))[:, mask] - cB[None, :]
    rows.append(report("harness | 120 targets | challenge axis | own control", EA, EB))

    # (iii) swap in the DONOR's control profile, keeping the harness's 150-cell target cells.
    # If the shared blob is a property of the harness's control POOL, this removes it.
    bA = np.load(os.path.join(REPO, "data/donors/orion_hct116/b_ctrl.npy"))
    bB = np.load(os.path.join(REPO, "data/donors/orion_hek293t/b_ctrl.npy"))
    selA = np.asarray([jA[chal[i]] for i in mask])
    selB = np.asarray([jB[chal[i]] for i in mask])
    EA2 = prof(hA["tsum"][ia].sum(1))[:, mask] - bA[selA][None, :]
    EB2 = prof(hB["tsum"][ib].sum(1))[:, mask] - bB[selB][None, :]
    rows.append(report("harness | 120 targets | challenge axis | DONOR control", EA2, EB2,
                       "<- control pool swapped"))

    # (iv) is the blob just the harness control being SHALLOWER/DEEPER?  Compare the harness
    # control profile with the donor control profile on the same axis.
    for lbl, hc, dc in (("HCT116", cA, bA[selA]), ("HEK293T", cB, bB[selB])):
        d = hc - dc
        print(f"    {lbl}: ||b_ctrl(harness) - b_ctrl(donor)|| = {np.linalg.norm(d):.4f}  "
              f"(||b_ctrl|| {np.linalg.norm(hc):.1f})  "
              f"cos to harness mean effect "
              f"{float(d @ (EA if lbl == 'HCT116' else EB).mean(0) / max(np.linalg.norm(d) * np.linalg.norm((EA if lbl == 'HCT116' else EB).mean(0)), 1e-9)):+.4f}")

    # (v) same 120 targets, harness cells, but the DONOR's own effect rows -- isolates
    # "which cells" from "which control".
    X = dA[np.asarray([iA[t] for t in shared])][:, selA].astype(np.float64)
    Y = dB[np.asarray([iB[t] for t in shared])][:, selB].astype(np.float64)
    rows.append(report("donor | 120 targets | challenge axis | donor control", X, Y))

    print("\nC. the mean effect vector -- the object the shared blob lives in")
    for lbl, E, D in (("HCT116", EA, X), ("HEK293T", EB, Y)):
        me, md = E.mean(0), D.mean(0)
        print(f"    {lbl}: ||mean effect||  harness {np.linalg.norm(me):7.3f}   "
              f"donor {np.linalg.norm(md):7.3f}   "
              f"cos(harness mean, donor mean) "
              f"{float(me @ md / max(np.linalg.norm(me) * np.linalg.norm(md), 1e-9)):+.4f}   "
              f"median ||row|| harness {np.median(np.linalg.norm(E, axis=1)):7.3f} "
              f"donor {np.median(np.linalg.norm(D, axis=1)):7.3f}")
    print(f"    cos(HCT116 harness mean, HEK293T harness mean) = "
          f"{float(EA.mean(0) @ EB.mean(0) / (np.linalg.norm(EA.mean(0)) * np.linalg.norm(EB.mean(0)))):+.4f}")

    print("\nD. after removing each line's own mean effect (the blob), harness vs donor")
    rows.append(report("harness | 120 targets | mean-removed",
                       EA - EA.mean(0), EB - EB.mean(0)))
    rows.append(report("donor   | 120 targets | mean-removed",
                       X - X.mean(0), Y - Y.mean(0)))

    with open(os.path.join(WORK, "reconcile.json"), "w") as fh:
        json.dump(rows, fh, indent=1)
    print(f"\nwrote {os.path.join(WORK, 'reconcile.json')}")


if __name__ == "__main__":
    main()
