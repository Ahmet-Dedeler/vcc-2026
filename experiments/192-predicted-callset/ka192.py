"""192 Phase 0 — known-answer rows.  ABORT if any fails.

Two families, both mandatory:

A. 188's POOL RECONSTRUCTION.  188 §Phase-3 / `pool188.py` reports that the h1-test
   reference's own BH<0.05 significant set is "reproduced at 0.7219" by `|lfc|*sqrt(cpm)`,
   0.6406 by `|h|`, 0.3796 by `|lfc|`, at a median pool size of 89.

   THE STATISTIC, stated exactly (188's own definition, read off `pool188.py:99-106`):

       for each perturbation p with n_real_p = |{g : padj_p,g < 0.05}| > 0,
           T_p(key) = the top n_real_p genes of row p by `key`
           o_p      = |T_p ∩ SIG_p| / n_real_p
       overlap(key) = mean_p o_p          (unweighted over the 85 perturbations)

   Because |T_p| = |SIG_p| = n_real_p exactly, o_p is simultaneously the recall and the
   precision of the key at that row's own true size.  It is NOT a rank correlation and NOT
   an AUC.  The keys use the REFERENCE'S OWN `lfc` and the control CPM -- so 0.7219 is an
   ORACLE number, the ceiling of the family, not a buildable score.  That distinction is the
   whole subject of this experiment and it is recorded here at the top.

B. h1-test production row: a15 reach 0.18283, median n_pred 6035 (raw BH) / 6034 (debench
   convention), fid 0.50701, jac 0.06824 (185 §1, §3b).  Donors `k562_gwps` + `orion_hek293t`
   only; `data/donors/h1_2025` never enters (it covers 85/85 of h1-test's targets).
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
E169 = os.path.join(REPO, "experiments/169-sharpening-six-member")
E185 = os.path.join(REPO, "experiments/185-rank-controlled-emission")
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
from debench import Surface  # noqa: E402

TS = 50_000.0
ALPHA = 0.05


def lfc_to_h(L, cpm):
    """The scorer's profile space.  h = log1p(TS*p_pert) - log1p(TS*p_ctrl), exact."""
    base = TS * 1e-6 * np.asarray(cpm, dtype=np.float64)
    return np.log1p(base[None, :] * np.exp2(np.nan_to_num(L))) - np.log1p(base)[None, :]


def topn_mask(key, n_per_row):
    """Top-n_p per row by `key`; ties broken by numpy's partition order (188's convention)."""
    out = np.zeros(key.shape, bool)
    for i in range(key.shape[0]):
        n = int(n_per_row[i])
        if n <= 0:
            continue
        n = min(n, key.shape[1])
        out[i, np.argpartition(-key[i], n - 1)[:n]] = True
    return out


def overlap(key, sig, n):
    """188's statistic: mean over rows of |top-n_p ∩ SIG_p| / n_p."""
    t = topn_mask(key, n)
    return float(np.mean([(t[i] & sig[i]).sum() / max(n[i], 1)
                          for i in range(len(n)) if n[i] > 0]))


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    fails, out = [], {}

    # ------------------------------------------------------------------ A
    print("=" * 100)
    print("A.  188's POOL RECONSTRUCTION (ORACLE keys: the reference's own lfc x the control CPM)")
    E = np.load(os.path.join(E169, "work/h1test-effect.npz"), allow_pickle=True)
    pos = {str(g): i for i, g in enumerate(E["genes"])}
    s = Surface(os.path.join(E169, "de169.q1.00_mf.npz"))
    cpm = np.array([E["cpm"][pos[g]] for g in s.genes], dtype=np.float64)
    R = np.nan_to_num(s.real_lfc)
    Rh = lfc_to_h(R, cpm)
    sig = s.real_sig & np.isfinite(s.real_lfc)
    n = sig.sum(1)

    keys = dict(abs_lfc=np.abs(R), abs_h=np.abs(Rh), wald=np.abs(R) * np.sqrt(cpm)[None, :],
                cpm_only=np.broadcast_to(cpm[None, :], R.shape).copy())
    want = dict(abs_lfc=0.3796, abs_h=0.6406, wald=0.7219, cpm_only=0.1512)
    print(f"    {len(s.perts)} perturbations x {len(s.genes)} gated genes; n_real per row: "
          f"median {np.median(n):.0f}, mean {n.mean():.0f}, min {n.min()}, max {n.max()}")
    for lab, val, exp, tol in [("median n_real", float(np.median(n)), 89.0, 0.5),
                               ("mean n_real", float(n.mean()), 765.0, 0.5)]:
        ok = abs(val - exp) <= tol
        print(f"    {lab:<22}{val:>10.4f}   expect {exp:<10.4f} {'PASS' if ok else 'FAIL'}")
        if not ok:
            fails.append(lab)
    out["n_real"] = dict(median=float(np.median(n)), mean=float(n.mean()),
                         min=int(n.min()), max=int(n.max()),
                         vector=[int(x) for x in n])
    for k, key in keys.items():
        o = overlap(key, sig, n)
        out[f"oracle_overlap_{k}"] = o
        ok = abs(o - want[k]) <= 5e-4
        print(f"    overlap({k:<9}){o:>10.4f}   expect {want[k]:<10.4f} {'PASS' if ok else 'FAIL'}")
        if not ok:
            fails.append(f"overlap {k}")

    # the random baseline for THIS statistic, which 188 never reported
    rng = np.random.default_rng(0)
    rnd = float(np.mean([overlap(rng.random(R.shape), sig, n) for _ in range(5)]))
    expct = float(np.mean([n[i] / R.shape[1] for i in range(len(n)) if n[i] > 0]))
    out["oracle_overlap_random"] = rnd
    out["oracle_overlap_expected_random"] = expct
    print(f"    overlap(random)  {rnd:>10.4f}   analytic mean n_p/G = {expct:.4f}  "
          f"-> wald is {want['wald']/expct:.1f}x random on this statistic")

    # ------------------------------------------------------------------ B
    print(f"\n{'=' * 100}")
    print("B.  h1-test production row (a15), 185 §1/§3b")
    s15 = Surface(os.path.join(E185, "work/de185.a15.npz"))
    nr, npd = s15.counts()
    rows = [("a15 reach", s15.reach(), 0.18283, 2e-5),
            ("a15 median n_pred (debench)", float(np.median(npd)), 6034.0, 0.5),
            ("a15 fid", s15.fid(), 0.50701, 2e-5),
            ("a15 jac", s15.jac(), 0.06824, 2e-5),
            ("a15 nmae", s15.nmae(), 0.96659, 2e-5),
            ("h1 n_real median (debench)", float(np.median(nr)), 89.0, 0.5)]
    z15 = np.load(os.path.join(E185, "work/de185.a15.npz"), allow_pickle=True)
    absH = z15["absH"].astype(np.float64)
    rows.append(("|H| key reach", s15.reach(key=absH), 0.23533, 2e-5))
    for lab, got, exp, tol in rows:
        ok = abs(got - exp) <= tol
        print(f"    {lab:<32}{got:>12.5f}   expect {exp:<12.5f} {'PASS' if ok else 'FAIL'}")
        out[lab.replace(" ", "_")] = float(got)
        if not ok:
            fails.append(lab)

    # raw-BH n_pred (6035) -- debench drops the row's own target gene, hence 6034
    raw_np = ((z15["pred_padj"] < ALPHA) & np.isfinite(z15["pred_padj"])).sum(1)
    ok = abs(float(np.median(raw_np)) - 6035.0) <= 0.5
    print(f"    {'a15 median n_pred (raw BH)':<32}{np.median(raw_np):>12.5f}   "
          f"expect {6035.0:<12.5f} {'PASS' if ok else 'FAIL'}")
    if not ok:
        fails.append("raw n_pred")
    out["a15_median_n_pred_rawBH"] = float(np.median(raw_np))

    # ------------------------------------------------------------------ leakage audit hook
    print(f"\n{'=' * 100}")
    print("C.  Donor / surface target overlap (AGENTS.md h1-test trap 1)")
    for d in ("k562_gwps", "orion_hek293t", "h1_2025", "orion_hct116", "rpe1",
              "k562_essential", "nadig_jurkat"):
        p = os.path.join(REPO, "data/donors", d, "perts.txt")
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            dp = {ln.strip() for ln in fh if ln.strip()}
        ov = len(dp & set(s15.perts.tolist()))
        # NB: high coverage is expected and permitted for the two PRODUCTION donors -- a15
        # predicts h1-test's targets FROM k562_gwps + orion_hek293t, which is the design.
        # The leak AGENTS.md names is `h1_2025`: H1's OWN perturbation responses.
        print(f"    {d:<16} covers {ov:>3}/85 of h1-test's targets"
              + ("   <-- SAME-LINE ORACLE, never used on this surface" if d == "h1_2025"
                 else "   (production donor)" if d in ("k562_gwps", "orion_hek293t") else ""))
        out[f"donor_overlap_{d}"] = ov

    with open(os.path.join(WORK, "ka192.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"\nwrote {WORK}/ka192.json")
    if fails:
        sys.exit(f"ABORT (rule 10 / mandate Phase 0): {fails}")
    print("\nALL PHASE-0 KNOWN-ANSWER ROWS PASS.")


if __name__ == "__main__":
    main()
