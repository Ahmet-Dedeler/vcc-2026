"""201 gate 0 (scoring half) — read 185's `rc_hist_a` row off the tables `ka201.py` wrote.

Every number here must reproduce 185 section 3b/3c/3d/3e for the SAME arm, computed with
185's own instruments (`debench.Surface`, `six185.walk`, `six185.members`).  A miss means
the re-key is not 185's and the preregistration does not apply to it.

Four arms carry known-answer rows: `a15` and `a15_b` (the two incumbent twins whose |z|
multisets the rc arms consume), `rc_hybrid` (the arm 191 built and the board scored, so its
row proves 201's `call_set_first` switch did not move it), and `rc_hist_a` with its twin.

Also reported, because this is the h1-test preview of the challenge-side falsifier 3: the
CALL-SET overlap between `rc_hist_a` and `rc_hybrid`, which is the variable the slot buys.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
sys.path.insert(0, f"{ROOT}/experiments/185-rank-controlled-emission")
from debench import Surface  # noqa: E402
import six185 as SIX  # noqa: E402

P = os.path.join(WORK, "de201.%s.npz")
ARMS = [a for a in ("a15", "a15_b", "rc_hybrid", "rc_hist_a", "rc_hist_a_seq",
                    "rc_hist_a_b") if os.path.exists(P % a)]

# 185 section 1 / 3b / 3c / 3d, verbatim.  Tolerances are 191's.
EXPECT = {
    ("a15", "reach"): (0.18283, 2e-5),
    ("a15", "n_pred"): (6034, 0.5),
    ("a15", "jac"): (0.06824, 5e-5),
    ("a15", "nmae"): (0.96659, 5e-5),
    ("a15", "fid"): (0.50701, 5e-5),
    ("a15", "P1"): (0.7317, 5e-4),
    ("a15", "H_call"): (-0.391, 2e-3),
    ("a15", "E_call"): (+0.486, 2e-3),
    ("a15_b", "reach"): (0.17537, 2e-5),
    ("a15_b", "n_pred"): (6031, 0.5),
    ("a15_b", "jac"): (0.06834, 5e-5),
    ("a15_b", "nmae"): (0.96641, 5e-5),
    ("a15_b", "P1"): (0.6707, 5e-4),
    ("rc_hybrid", "reach"): (0.21681, 2e-5),
    ("rc_hybrid", "n_pred"): (6043, 0.5),
    ("rc_hybrid", "jac"): (0.06827, 5e-5),
    ("rc_hybrid", "nmae"): (0.96670, 5e-5),
    ("rc_hybrid", "fid"): (0.50693, 5e-5),
    ("rc_hybrid", "P1"): (0.8780, 5e-4),
    ("rc_hybrid", "H_call"): (+0.984, 2e-3),
    ("rc_hybrid", "E_call"): (+0.007, 5e-3),
    # ---- the arm itself (the mandate's known-answer row) --------------------------------
    ("rc_hist_a", "reach"): (0.23278, 2e-5),
    ("rc_hist_a", "n_pred"): (6040, 0.5),
    ("rc_hist_a", "jac"): (0.05932, 5e-5),
    ("rc_hist_a", "fid"): (0.50808, 5e-5),
    ("rc_hist_a", "nmae"): (0.96670, 5e-5),
    ("rc_hist_a", "recall"): (0.616, 5e-4),
    ("rc_hist_a", "P1"): (0.8780, 5e-4),
    ("rc_hist_a", "P3"): (0.7851, 5e-4),
    ("rc_hist_a", "P10"): (0.6914, 5e-4),
    ("rc_hist_a", "H_univ"): (+0.999, 2e-3),
    ("rc_hist_a", "H_call"): (+0.996, 2e-3),
    ("rc_hist_a", "E_call"): (+0.076, 5e-3),
    ("rc_hist_a_b", "reach"): (0.23272, 2e-5),
    ("rc_hist_a_b", "n_pred"): (6032, 0.5),
    ("rc_hist_a_b", "jac"): (0.05930, 5e-5),
    ("rc_hist_a_b", "nmae"): (0.96670, 5e-5),
    ("rc_hist_a_b", "P1"): (0.8780, 5e-4),
}


def call_sets(path):
    """Our BH call set per perturbation, self gene excluded, exactly `debench`'s."""
    table = np.load(path, allow_pickle=True)
    genes = np.asarray(table["genes"]).astype(str)
    perts = np.asarray(table["perts"]).astype(str)
    padj = table["pred_padj"]
    sig = padj < 0.05
    for i, target in enumerate(perts):
        own = np.flatnonzero(genes == target)
        if own.size:
            sig[i, own[0]] = False
    return sig


def main() -> None:
    meta = np.load(os.path.join(
        ROOT, "experiments/185-rank-controlled-emission/work/ctrl_meta.npz"
    ), allow_pickle=True)
    expr = np.log1p(meta["cpm"][meta["keep"]])
    z15 = np.load(P % "a15", allow_pickle=True)
    absH = z15["absH"].astype(np.float64)

    rows: list[tuple[str, str, float, float, float, bool]] = []
    out: dict[str, dict] = {}
    print(f"{'arm':<16}{'n_pred':>7}{'reach':>10}{'fid':>9}{'jac':>9}{'nmae':>9}"
          f"{'recall':>8}{'P(1)':>8}{'P(3)':>8}{'P(10)':>8}{'H_univ':>9}{'H_call':>9}"
          f"{'E_call':>9}")
    for arm in ARMS:
        surface = Surface(P % arm)
        members = SIX.members(surface)
        _reach, pk, _ks = SIX.walk(surface)
        table = np.load(P % arm, allow_pickle=True)
        spear = SIX.spearmans(table, surface, absH, expr)
        record = dict(members, P1=pk[1], P3=pk[3], P10=pk[10],
                      H_call=spear["H_call"], E_call=spear["E_call"],
                      H_univ=spear["H_univ"], E_univ=spear["E_univ"])
        out[arm] = record
        print(f"{arm:<16}{record['n_pred']:>7}{record['reach']:>10.5f}{record['fid']:>9.5f}"
              f"{record['jac']:>9.5f}{record['nmae']:>9.5f}{record['recall']:>8.4f}"
              f"{pk[1]:>8.4f}{pk[3]:>8.4f}{pk[10]:>8.4f}"
              f"{spear['H_univ']:>+9.3f}{spear['H_call']:>+9.3f}{spear['E_call']:>+9.3f}")
        for key, (want, tol) in EXPECT.items():
            if key[0] != arm:
                continue
            got = float(record[key[1]])
            rows.append((arm, key[1], got, want, tol, abs(got - want) <= tol))

    print("\nKNOWN-ANSWER ROWS against 185 (the mandate's gate 0)")
    bad = []
    for arm, name, got, want, tol, ok in rows:
        print(f"  {arm:<16}{name:<8}{got:>11.5f}   185 {want:<10.5f} tol {tol:<8g} "
              f"{'PASS' if ok else 'FAIL'}")
        if not ok:
            bad.append(f"{arm}/{name}")

    if "rc_hist_a" in out and "rc_hist_a_seq" in out:
        print("\nSEEDING: 185's two-stream bookkeeping vs production's single stream")
        for name in ("reach", "jac", "nmae", "fid"):
            a, b = out["rc_hist_a"][name], out["rc_hist_a_seq"][name]
            print(f"  {name:<8}{a:>11.5f}{b:>11.5f}   |diff| {abs(a - b):.5f}"
                  f"   (185 section 3e twin spread on this family: 0.00007)")

    overlap = {}
    if "rc_hist_a" in out and "rc_hybrid" in out:
        print("\nCALL-SET OVERLAP on h1-test (the challenge-side falsifier 3's preview)")
        left, right = call_sets(P % "rc_hist_a"), call_sets(P % "rc_hybrid")
        inter = (left & right).sum(1)
        union = (left | right).sum(1)
        jacc = inter / np.maximum(union, 1)
        share = inter / np.maximum(left.sum(1), 1)
        a15sets = call_sets(P % "a15")
        share15 = (left & a15sets).sum(1) / np.maximum(left.sum(1), 1)
        hy15 = (right & a15sets).sum(1) / np.maximum(right.sum(1), 1)
        overlap = dict(
            jaccard_median=float(np.median(jacc)), jaccard_mean=float(jacc.mean()),
            share_of_rc_hist_a_in_rc_hybrid=float(np.median(share)),
            share_of_rc_hist_a_in_a15=float(np.median(share15)),
            share_of_rc_hybrid_in_a15=float(np.median(hy15)),
            n_pred_rc_hist_a=float(np.median(left.sum(1))),
            n_pred_rc_hybrid=float(np.median(right.sum(1))),
        )
        for key, value in overlap.items():
            print(f"  {key:<34}{value:>12.4f}")

    with open(os.path.join(WORK, "score201.json"), "w") as handle:
        json.dump(dict(members=out, known=[list(r) for r in rows], overlap=overlap),
                  handle, default=float)
    if bad:
        sys.exit(f"ABORT (rule 10 known-answer): {bad}")
    print("\nall known-answer rows PASS")


if __name__ == "__main__":
    main()
