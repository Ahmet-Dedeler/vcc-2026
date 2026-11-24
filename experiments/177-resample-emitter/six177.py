"""177 — read the DE members off the h1-test tables for the resampling-emitter arms, convert
them to board units, and report the ORDERING diagnostics that are this experiment's mechanism.

`reach`, `fid`, `jac` and `nmae` for every arm, on the ONE surface validated for `reach`
(157 section 1: h1-test reproduces both board raw `reach` points to within 0.014; repaired
HCT116 misses a15's by 2.6x; 165 restricts that validation to EMISSION-LAW changes at a fixed
profile, which is what every arm here is -- `realize177.py` proves the fixity with one sha256
per (arm, target)).

RULE 9 IS IN FORCE AND IS HANDLED EXPLICITLY
--------------------------------------------
`fid`, `jac` and `reach`'s normalisation all carry `n_real`, which is nothing like the
challenge's on any harness.  So:

  * `reach` -- read at native BH (`f = 1`, the board-calibrated setting) AND under 157's
    PROPORTIONAL power matching, `pool_p = the reference's top-round(f*n_real_p)`.
  * `fid` -- reported as PRECISION `mean_p(k/n_pred)`, which carries no `n_real` at all, AND
    the yield cap is priced at BOTH ends of the inferred challenge `n_real` range (240 and
    2,670), because this experiment moves `n_pred` by an order of magnitude and the cap is
    therefore NOT saturated the way it is for every previous emission arm.
  * `jac` -- raw h1-test `jac` is unreadable (its union carries `n_real`); the overlap RATIO
    is what transfers and it is rebuilt onto the board's own `n_real`/`n_pred`.

THE MECHANISM COLUMNS: `rho expr` / `rho |H|` are Spearman(-log10 pred padj, x) INSIDE our own
call set, the exact statistic INBOX-4's falsifier and 171 section 12a reported.  a15 reads
**+0.485 expr / -0.391 |H|**: our confidence order tracks EXPRESSION and ANTI-tracks our own
predicted effect.  `reach |H| ceiling` is the counterfactual INBOX-4 section 2 priced at
+0.036 raw inside the call block (0.21914 against 0.18283).

BOARD CALIBRATION FOR `reach`: the a8 <-> a15 pair is the known-answer row.  h1-test moves
0.2029 -> 0.1828 (-0.0201) where the board moved 0.1890 -> 0.1811 (-0.0079), so the h1-test
delta is multiplied by that ratio (0.393; 164 used the same construction and called it 0.428).

No cell-eval2, no build, no submission.
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np
from scipy.stats import rankdata

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
from debench import ANCHOR, Surface  # noqa: E402

ALPHA = 0.05
FACTORS = (1.0, 2.0, 3.5)
ORDER = ["KA_a15_raw", "base_mf", "base_mf_b", "quasi100", "rs_all", "rs_all_b",
         "rs_sched", "rs_sched_b", "rs_sched_perm", "rs_hkey", "rs_hkey_b",
         "hy_kt2000", "hy_kt2000_b", "hy_kt4000", "hy_kt6000", "hy_kt8000"]
# 157 section 1, both measured at f = 1 on this surface, against the board's own pair.
BOARD_REACH = dict(a8=0.1890, a15=0.1811)
H1_REACH = dict(a8=0.2029, a15=0.1828)
# nmae has no debench anchor; AGENTS.md section 4 gives b 1.0009-1.0017, r 0.369-0.431.
NMAE_B, NMAE_R = 1.0013, 0.400
EFF = os.path.join(REPO, "experiments/169-sharpening-six-member/work/h1test-effect.npz")


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
    """`mean_p(|pred & real| / n_real)` -- recall on the reference's set, no union."""
    out = []
    for p in range(len(s.perts)):
        nr = int(s.real_sig[p].sum())
        if nr == 0:
            continue
        out.append(int((s.real_sig[p] & s.pred_sig[p]).sum()) / nr)
    return float(np.mean(out))


def call_sign(s):
    """Sign accuracy over our own adjudicable call set, pooled -- 169's `sign` column."""
    k = n = 0
    for p in range(len(s.perts)):
        sel = np.flatnonzero(s.pred_sig[p] & s.adj[p])
        k += int((np.sign(s.pred_lfc[p, sel]) == np.sign(s.real_lfc[p, sel])).sum())
        n += sel.size
    return k / max(n, 1)


def rho_pair(s, absH, lcpm):
    """Spearman(-log10 pred padj, .) inside our call set, median over perturbations."""
    re, rh = [], []
    for p in range(len(s.perts)):
        sel = np.flatnonzero(s.pred_sig[p])
        if sel.size < 50:
            continue
        lp = -np.log10(np.maximum(s.pred_padj[p, sel], 1e-300))
        if np.ptp(lp) == 0:
            continue
        re.append(np.corrcoef(rankdata(lp), rankdata(lcpm[sel]))[0, 1])
        rh.append(np.corrcoef(rankdata(lp), rankdata(absH[p, sel]))[0, 1])
    return (float(np.median(re)) if re else np.nan,
            float(np.median(rh)) if rh else np.nan)


def main() -> None:
    paths = {os.path.basename(p)[len("de177."):-4]: p
             for p in sorted(glob.glob(os.path.join(HERE, "de177.*.npz")))}
    arms = [a for a in ORDER if a in paths] + [a for a in paths if a not in ORDER]
    if not arms:
        sys.exit("no de177.*.npz -- run realize177.py first")

    E = np.load(EFF, allow_pickle=True)
    eg = {str(g): i for i, g in enumerate(E["genes"])}
    et = {str(t): i for i, t in enumerate(E["targets"])}

    rec = {}
    print(f"{'arm':<15}{'n_pred':>8}{'reach f1':>10}{'reach f2':>10}{'reach f3.5':>11}"
          f"{'prec':>8}{'fid raw':>9}{'jac raw':>9}{'recall':>8}{'nmae':>9}"
          f"{'sign':>8}{'rho expr':>10}{'rho |H|':>9}")
    for a in arms:
        s = Surface(paths[a])
        gi = np.array([eg[g] for g in s.genes])
        pi = np.array([et[p] for p in s.perts])
        absH = np.abs(E["H"][pi][:, gi])
        lcpm = np.log10(np.maximum(E["cpm"][gi], 1e-9))
        nr, npd = s.counts()
        re, rh = rho_pair(s, absH, lcpm)
        row = dict(n_pred=float(np.median(npd)), n_real=float(np.median(nr)),
                   precision=precision(s), fid=s.fid(), jac=s.jac(), nmae=s.nmae(),
                   pds=s.pds(), sign=call_sign(s), recall=overlap_ratio(s),
                   rho_expr=re, rho_H=rh,
                   reach_H=s.reach(key=absH),
                   reach_sigH=s.reach(key=np.where(s.pred_sig, 1e6 + absH, absH)))
        for f in FACTORS:
            row[f"reach_f{f}"] = PowerSurface(paths[a], f).reach()
        rec[a] = row
        print(f"{a:<15}{row['n_pred']:>8.0f}{row['reach_f1.0']:>10.5f}"
              f"{row['reach_f2.0']:>10.5f}{row['reach_f3.5']:>11.5f}"
              f"{row['precision']:>8.5f}{row['fid']:>9.5f}{row['jac']:>9.5f}"
              f"{row['recall']:>8.5f}{row['nmae']:>9.5f}{row['sign']:>8.5f}"
              f"{re:>+10.3f}{rh:>+9.3f}", flush=True)

    # ---------------------------------------------------------------- gates
    print("\n--- known-answer gates ---")
    ok = True
    if "KA_a15_raw" in rec:
        g = rec["KA_a15_raw"]
        for label, got, want, tol in (("median n_pred", g["n_pred"], 6035, 1.5),
                                      ("reach f=1", g["reach_f1.0"], 0.18283, 5e-4),
                                      ("reach key=|H|", g["reach_H"], 0.23533, 5e-4),
                                      ("reach sig-then-|H|", g["reach_sigH"], 0.21914, 5e-4)):
            good = abs(got - want) <= tol
            ok &= good
            print(f"  164/157/INBOX-4 a15 {label:<20} got {got:.5f}  expect {want:.5f}  "
                  f"{'PASS' if good else 'FAIL'}")
    if "base_mf" in rec:
        g = rec["base_mf"]
        for label, got, want, tol in (("median n_pred", g["n_pred"], 6066, 40),
                                      ("reach f=1", g["reach_f1.0"], 0.18432, 0.01)):
            good = abs(got - want) <= tol
            ok &= good
            print(f"  173 base_mf {label:<27} got {got:.5f}  expect {want:.5f}  "
                  f"{'PASS' if good else 'FAIL'}")
    if "base_mf" in rec and "base_mf_b" in rec and "KA_a15_raw" in rec:
        seed = abs(rec["base_mf"]["reach_f1.0"] - rec["base_mf_b"]["reach_f1.0"])
        rt = abs(rec["base_mf"]["reach_f1.0"] - rec["KA_a15_raw"]["reach_f1.0"])
        good = rt <= max(3 * seed, 0.01)
        ok &= good
        print(f"  recovery round trip: |base_mf - KA_a15_raw| reach = {rt:.5f} against a "
              f"matched-seed spread of {seed:.5f}  {'PASS' if good else 'FAIL'}")
    if not ok:
        print("  *** A GATE FAILED -- everything below is void ***")

    print("\n--- matched-seed twin spreads (the resolution floor) ---")
    for a, bb in (("base_mf", "base_mf_b"), ("rs_all", "rs_all_b"),
                  ("rs_sched", "rs_sched_b"), ("rs_hkey", "rs_hkey_b"),
                  ("hy_kt2000", "hy_kt2000_b")):
        if a in rec and bb in rec:
            print(f"  {a:<12} vs {bb:<14} " + "  ".join(
                f"{m} {abs(rec[a][m]-rec[bb][m]):.5f}"
                for m in ("reach_f1.0", "precision", "nmae", "sign")))

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
    print(f"\nDELTAS vs `{base}` (the incumbent a15 + a23 mass fix), in OVERALL board units.")
    print("`d fid` is given at BOTH ends of the inferred challenge n_real (rule 9): 240 (cap")
    print("saturated for every arm) and 2,670 (cap binds for the low-n_pred arms).")
    print(f"{'arm':<15}{'d reach raw':>12}{'d reach ovr':>12}"
          f"{'d prec':>9}{'d fid240':>10}{'d fid2670':>11}{'d recall':>10}{'d jac ovr':>10}"
          f"{'d nmae':>9}{'d nmae ovr':>11}{'SUM@240':>10}{'SUM@2670':>10}")
    NP_BOARD = 4424.0                       # a15's board-side n_pred (169's convention)
    for a in arms:
        r = rec[a]
        dr = r["reach_f1.0"] - b0["reach_f1.0"]
        dr_ovr = dr * k_reach / span_r / 6
        scale = r["n_pred"] / b0["n_pred"]   # our n_pred moves; carry it to the board's
        df = {}
        for NR in (240.0, 2670.0):
            f1 = r["precision"] * min(1.0, NP_BOARD * scale / NR)
            f0 = b0["precision"] * min(1.0, NP_BOARD / NR)
            df[NR] = (f1 - f0) / span_f / 6
        NR, NP = 240.0, NP_BOARD
        j0 = b0["recall"] * NR / (NP + NR - b0["recall"] * NR)
        j1 = r["recall"] * NR / (NP * scale + NR - r["recall"] * NR)
        dj_ovr = (j1 - j0) / span_j / 6
        dn = r["nmae"] - b0["nmae"]
        dn_ovr = -dn / span_n / 6
        print(f"{a:<15}{dr:>12.5f}{dr_ovr:>+12.6f}"
              f"{r['precision'] - b0['precision']:>+9.5f}{df[240.0]:>+10.6f}"
              f"{df[2670.0]:>+11.6f}{r['recall'] - b0['recall']:>+10.5f}{dj_ovr:>+10.6f}"
              f"{dn:>+9.5f}{dn_ovr:>+11.6f}"
              f"{dr_ovr + df[240.0] + dj_ovr + dn_ovr:>+10.6f}"
              f"{dr_ovr + df[2670.0] + dj_ovr + dn_ovr:>+10.6f}")

    print(f"\nhow much of INBOX-4's inside-block ordering counterfactual each arm REALISES")
    print(f"(ceiling = `reach` with the called genes reordered by |H|, on that arm's own table)")
    print(f"{'arm':<15}{'reach':>10}{'ceil sig|H|':>12}{'gap':>9}{'ceil all |H|':>13}"
          f"{'realised vs base':>18}")
    for a in arms:
        r = rec[a]
        real = r["reach_f1.0"] - b0["reach_f1.0"]
        ceil = b0["reach_sigH"] - b0["reach_f1.0"]
        print(f"{a:<15}{r['reach_f1.0']:>10.5f}{r['reach_sigH']:>12.5f}"
              f"{r['reach_sigH'] - r['reach_f1.0']:>9.5f}{r['reach_H']:>13.5f}"
              f"{(real / ceil if ceil else np.nan):>17.1%}")

    print(f"\npower-matched `reach` deltas (rule 9: the challenge's n_real sits between f=2 and "
          f"f=3.5\non this surface, 157 section 3b) -- RAW, before the transfer factor")
    print(f"{'arm':<15}" + "".join(f"{'f=' + str(f):>12}" for f in FACTORS))
    for a in arms:
        print(f"{a:<15}" + "".join(
            f"{rec[a][f'reach_f{f}'] - b0[f'reach_f{f}']:>+12.5f}" for f in FACTORS))

    with open(os.path.join(WORK, "six177.json"), "w") as fh:
        json.dump(rec, fh, indent=1, default=float)
    print(f"\nwrote {os.path.join(WORK, 'six177.json')}")


if __name__ == "__main__":
    sys.exit(main())
