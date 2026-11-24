"""173 — read the four DE members (sign-repair arms) off the h1-test tables, and convert them to board units.

`reach`, `fid`, `jac` and `nmae` for every sharpening arm, on the ONE surface validated for
`reach` (157 §1: h1-test reproduces both board raw `reach` points to within 0.014; repaired
HCT116 misses a15's by 2.6x).

RULE 9 IS IN FORCE AND IS HANDLED EXPLICITLY
--------------------------------------------
`fid`, `jac` and `reach`'s normalisation all carry `n_real`, which is nothing like the
challenge's on any harness.  So:

  * `reach` -- read at native BH (`f = 1`, the board-calibrated setting) AND under 157's
    PROPORTIONAL power matching, `pool_p = the reference's top-round(f*n_real_p)`, which
    raises power while preserving the per-perturbation heterogeneity the metric averages
    over.  `debench.Surface(n_real_target=N)` is NOT used: 157 §3a measured that flattening
    every perturbation to a common `N` deletes the small-`n_real` tercile that carries 98%
    of our reach and collapses the statistic 23x for a 2.7x pool change.
  * `fid` -- reported as PRECISION `mean_p(k/n_pred)`, which carries no `n_real` at all.
    Raw h1-test `fid` is printed but flagged; 164 §3.3 measured it moving the OPPOSITE way
    from precision purely through the yield cap.
  * `jac` -- raw h1-test `jac` is unreadable (its union carries `n_real`); the overlap RATIO
    is what transfers and it is rebuilt onto the board's own `n_real`/`n_pred`.

BOARD CALIBRATION FOR `reach`: the a8 <-> a15 pair is the known-answer row.  h1-test moves
0.2029 -> 0.1828 (-0.0201) where the board moved 0.1890 -> 0.1811 (-0.0079), so the h1-test
delta is multiplied by that ratio.  164 used the same construction and called it 0.428.

No cell-eval2, no build, no submission.
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
from debench import ANCHOR, Surface  # noqa: E402

ALPHA = 0.05
FACTORS = (1.0, 2.0, 3.5)
ORDER = ["KA_a15_raw", "base_mf", "base_mf_b", "flip_mf", "flip_mf_b", "zero_mf", "permflip_mf"]
# 157 §1, both measured at f = 1 on this surface, against the board's own pair.
BOARD_REACH = dict(a8=0.1890, a15=0.1811)
H1_REACH = dict(a8=0.2029, a15=0.1828)
# nmae has no debench anchor; AGENTS.md §4 gives b 1.0009-1.0017, r 0.369-0.431.
NMAE_B, NMAE_R = 1.0013, 0.400


class PowerSurface(Surface):
    """157's proportional power matching: `pool_p = reference top-round(f * n_real_p)`."""

    def __init__(self, path, factor=1.0):
        super().__init__(path)
        z = np.load(path, allow_pickle=True)
        rp = np.where(np.isfinite(z["real_p"]), z["real_p"], np.inf)
        base = (z["real_padj"] < ALPHA) & np.isfinite(z["real_padj"])
        sig = np.zeros(rp.shape, bool)
        for i in range(rp.shape[0]):
            k = int(round(factor * int(base[i].sum())))
            k = max(0, min(k, int(np.isfinite(z["real_p"][i]).sum())))
            if k:
                sig[i, np.argsort(rp[i], kind="stable")[:k]] = True
        self.real_sig = sig & np.isfinite(z["real_p"])
        for p, j in enumerate(self.self_idx):
            if j >= 0:
                self.real_sig[p, j] = False


def precision(s):
    """`mean_p(k / n_pred)` on the adjudicable call set -- the `n_real`-free half of `fid`."""
    out = []
    for p in range(len(s.perts)):
        sel = np.flatnonzero(s.pred_sig[p] & s.adj[p])
        if sel.size == 0:
            continue
        out.append(float((np.sign(s.pred_lfc[p, sel]) ==
                          np.sign(s.real_lfc[p, sel])).sum()) / sel.size)
    return float(np.mean(out))


def overlap_ratio(s):
    """`mean_p(|pred & real| / n_real)` -- recall on the reference's set, no union, so the
    only `n_real` dependence is the one the board's own `n_real` then supplies."""
    out = []
    for p in range(len(s.perts)):
        nr = int(s.real_sig[p].sum())
        if nr == 0:
            continue
        out.append(int((s.real_sig[p] & s.pred_sig[p]).sum()) / nr)
    return float(np.mean(out))


def main() -> None:
    paths = {os.path.basename(p)[len("de173."):-4]: p
             for p in sorted(glob.glob(os.path.join(HERE, "de173.*.npz")))}
    arms = [a for a in ORDER if a in paths] + [a for a in paths if a not in ORDER]
    if not arms:
        sys.exit("no de173.*.npz -- run realize169.py first")

    rec = {}
    print(f"{'arm':<14}{'n_pred':>8}{'reach f1':>10}{'reach f2':>10}{'reach f3.5':>11}"
          f"{'precision':>11}{'fid raw':>9}{'jac raw':>9}{'recall':>8}{'nmae':>9}"
          f"{'pds':>9}{'sign':>8}")
    for a in arms:
        s = Surface(paths[a])
        nr, npd = s.counts()
        row = dict(n_pred=float(np.median(npd)), n_real=float(np.median(nr)),
                   precision=precision(s), fid=s.fid(), jac=s.jac(), nmae=s.nmae(),
                   pds=s.pds(), sign=s.sign_accuracy(), recall=overlap_ratio(s))
        for f in FACTORS:
            row[f"reach_f{f}"] = PowerSurface(paths[a], f).reach()
        rec[a] = row
        print(f"{a:<14}{row['n_pred']:>8.0f}{row['reach_f1.0']:>10.5f}"
              f"{row['reach_f2.0']:>10.5f}{row['reach_f3.5']:>11.5f}"
              f"{row['precision']:>11.5f}{row['fid']:>9.5f}{row['jac']:>9.5f}"
              f"{row['recall']:>8.5f}{row['nmae']:>9.5f}{row['pds']:>9.5f}"
              f"{row['sign']:>8.5f}", flush=True)

    # ---------------------------------------------------------------- gates
    print(f"\n--- known-answer gates ---")
    ok = True
    if "KA_a15_raw" in rec:
        g = rec["KA_a15_raw"]
        for label, got, want, tol in (("median n_pred", g["n_pred"], 6035, 0.5),
                                      ("reach f=1", g["reach_f1.0"], 0.18283, 5e-4)):
            good = abs(got - want) <= tol
            ok &= good
            print(f"  164/157 a15 {label:<16} got {got:.5f}  expect {want:.5f}  "
                  f"{'PASS' if good else 'FAIL'}")
    if "base_mf" in rec and "base_mf_b" in rec and "KA_a15_raw" in rec:
        seed = abs(rec["base_mf"]["reach_f1.0"] - rec["base_mf_b"]["reach_f1.0"])
        rt = abs(rec["base_mf"]["reach_f1.0"] - rec["KA_a15_raw"]["reach_f1.0"])
        print(f"  recovery round trip: |q1.00 - KA_a15_raw| reach = {rt:.5f} against a "
              f"re-seed spread of {seed:.5f}  "
              f"{'PASS' if rt <= max(3 * seed, 0.01) else 'FAIL'}")
        ok &= rt <= max(3 * seed, 0.01)
        for m in ("precision", "nmae", "pds", "n_pred"):
            print(f"      {m:<10} KA {rec['KA_a15_raw'][m]:.5f}  base_mf {rec['base_mf'][m]:.5f}"
                  f"  seed spread {abs(rec['base_mf'][m] - rec['base_mf_b'][m]):.5f}")
    if not ok:
        print("  *** A GATE FAILED -- everything below is void ***")

    # ---------------------------------------------------------------- board units
    k_reach = ((BOARD_REACH["a8"] - BOARD_REACH["a15"]) /
               (H1_REACH["a8"] - H1_REACH["a15"]))
    print(f"\nboard/h1-test `reach` transfer factor from the a8<->a15 pair: {k_reach:.4f}")
    base = "base_mf"
    b0 = rec[base]
    span_r = ANCHOR["reach"][1] - ANCHOR["reach"][0]
    span_f = ANCHOR["fid"][1] - ANCHOR["fid"][0]
    span_j = ANCHOR["jac"][1] - ANCHOR["jac"][0]
    span_n = NMAE_B - NMAE_R
    print(f"\nDELTAS vs `{base}` (the incumbent a15 + a23 mass fix), in OVERALL board units")
    print(f"{'arm':<14}{'d reach raw':>12}{'x transfer':>11}{'d reach ovr':>12}"
          f"{'d prec':>9}{'d fid ovr':>10}{'d recall':>10}{'d jac ovr':>10}"
          f"{'d nmae':>9}{'d nmae ovr':>11}{'SUM DE':>10}")
    for a in arms:
        r = rec[a]
        dr = r["reach_f1.0"] - b0["reach_f1.0"]
        dr_ovr = dr * k_reach / span_r / 6
        # `fid` on the challenge is precision * min(1, n_pred/n_real); n_pred ~ 4424 vs
        # n_real ~ 240, so the yield cap is saturated and precision IS fid there.
        df_ovr = (r["precision"] - b0["precision"]) / span_f / 6
        # `jac` rebuilt on the board's own counts: J = |A n B| / (n_pred + n_real - |A n B|)
        NR, NP = 240.0, 4424.0
        j0 = b0["recall"] * NR / (NP + NR - b0["recall"] * NR)
        j1 = r["recall"] * NR / (NP + NR - r["recall"] * NR)
        dj_ovr = (j1 - j0) / span_j / 6
        dn = r["nmae"] - b0["nmae"]
        dn_ovr = -dn / span_n / 6
        print(f"{a:<14}{dr:>12.5f}{k_reach:>11.4f}{dr_ovr:>+12.6f}"
              f"{r['precision'] - b0['precision']:>+9.5f}{df_ovr:>+10.6f}"
              f"{r['recall'] - b0['recall']:>+10.5f}{dj_ovr:>+10.6f}"
              f"{dn:>+9.5f}{dn_ovr:>+11.6f}"
              f"{dr_ovr + df_ovr + dj_ovr + dn_ovr:>+10.6f}")

    print(f"\npower-matched `reach` deltas (rule 9: the challenge's n_real ~ 240 sits between "
          f"f=2 and f=3.5 on\nthis surface, 157 §3b) -- RAW, before the transfer factor")
    print(f"{'arm':<14}" + "".join(f"{'f=' + str(f):>12}" for f in FACTORS))
    for a in arms:
        print(f"{a:<14}" + "".join(
            f"{rec[a][f'reach_f{f}'] - b0[f'reach_f{f}']:>+12.5f}" for f in FACTORS))

    with open(os.path.join(WORK, "six173.json"), "w") as fh:
        json.dump(rec, fh, indent=1, default=float)
    print(f"\nwrote {os.path.join(WORK, 'six169.json')}")


if __name__ == "__main__":
    sys.exit(main())
