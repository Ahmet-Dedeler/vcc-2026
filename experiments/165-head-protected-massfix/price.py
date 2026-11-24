"""165 — the six-member price of every mass-correction placement, and the board conversion.

Reads the `work/de.<arm>.s<seed>.npz` tables and `work/profiles.npz` written by `arms.py`.
Deterministic given those files: no emitter, no cell-eval2, no build, no submission.

HOW EACH MEMBER IS PRICED, and which are trusted
------------------------------------------------
`reach`, `fid`, `jac`, `nmae`   `debench.Surface`, on the h1-test tables, at NATIVE power.
                                Rule 9 says the LEVELS of the first three do not transfer.
                                `reach` is converted to the board by CALIBRATING on `mn`,
                                which is a23 and whose board delta is known (-0.007933): the
                                calibration absorbs the harness->board factor AND this
                                artifact's 1.1-1.4x larger mass excess in one number, and it
                                is exactly the ratio the experiment is asking about.
`pds`, `mse`                    profile level, 156's `screen.endpoints` on the DELIVERED
                                effect `log1p(TS*p) - b_ctrl` against the h1-test reference's
                                own pseudobulk effect in the same space, panel targets
                                excluded exactly as `pds_cosine` excludes them. `sse` is the
                                uncorrected `mse` ratio at its own best global amplitude.
                                Both are converted by calibrating on `mn` against a23's
                                +0.001737 / +0.001525.

KNOWN-ANSWER GATES (AGENTS.md 11.1 rule 10's corollary).  ABORTS if either fails:
  a15 raw reach == 0.18283 and a8 raw reach == 0.20291, both from 157 `power_h1test.log` at
  f = 1.0, and themselves within 0.014 of the board's 0.1811 / 0.1890.
"""
from __future__ import annotations

import glob
import os
import re
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
sys.path.insert(0, f"{ROOT}/experiments/156-delivered-geometry")
from debench import Surface                                   # noqa: E402
from screen import endpoints                                  # noqa: E402

TS = 5.0e4
SEEDS = (90_000, 7_000_000)
ARMS = ["a15", "a8", "mn", "mn_half", "cst", "hp_p500", "hp_p2000", "hp_p6000",
        "hp_p9000", "hp_w2000", "eps_keep", "eps_tail", "inveff", "eff", "mn_c10"]
GATE = {"a15": 0.18283, "a8": 0.20291}            # 157 power_h1test.log, f = 1.0
SPAN = {"pds": 0.45208, "fid": 0.29725, "reach": 0.88512, "jac": 0.36689}
# a23's realised board deltas, `A23-BOARD-RESULT.md` (entry IzQVfDpeDWufkPT8Lmd2)
A23 = {"pds": +0.001737, "mse": +0.001525, "nmae": -0.000106,
       "fid": +0.000337, "reach": -0.007933, "jac": +0.000081}


def main() -> None:
    # Every arm is read at EVERY seed it has. `reso.py` shows the two-seed reach estimate is
    # noise-dominated, so the six arms carried to 12 seeds must be priced on all 12.
    seeds = {a: sorted(int(re.search(r"\.s(\d+)\.npz$", f).group(1))
                       for f in glob.glob(f"{WORK}/de.{a}.s*.npz"))
             for a in ARMS if glob.glob(f"{WORK}/de.{a}.s*.npz")}
    have = [a for a in ARMS if a in seeds]
    S = {(a, s): Surface(f"{WORK}/de.{a}.s{s}.npz") for a in have for s in seeds[a]}

    # ---- known-answer gates -----------------------------------------------------------
    print("=" * 104)
    print("(0) KNOWN-ANSWER GATES — 157's h1-test reach at f = 1, which reproduce two board "
          "points to 0.014")
    print("=" * 104)
    bad = []
    for arm, want in GATE.items():
        if arm not in have:
            continue
        got = S[(arm, SEEDS[0])].reach()   # the gate is a SINGLE-seed identity, seed 90000
        ok = abs(got - want) < 5e-4
        print(f"  {arm:<10} reach {got:.5f}  expect {want:.5f}  "
              f"delta {got - want:+.5f}  {'PASS' if ok else 'FAIL'}")
        if not ok:
            bad.append(arm)
    if bad:
        sys.exit(f"ABORT: known-answer reach gate failed for {bad}")

    # ---- the DE members, both seeds ---------------------------------------------------
    print("\n" + "=" * 104)
    print("(1) h1-test DE members, native power, per arm (two emission seeds)")
    print("=" * 104)
    print(f"{'arm':<11}{'seeds':>7}{'n_pred':>8}{'reach':>10}{'d reach':>10}"
          f"{'fid':>9}{'jac':>9}{'nmae':>9}{'pds(lfc)':>10}")
    de = {}
    for a in have:
        sd = seeds[a]
        de[a] = dict(
            reach=float(np.mean([S[(a, s)].reach() for s in sd])),
            fid=float(np.mean([S[(a, s)].fid() for s in sd])),
            jac=float(np.mean([S[(a, s)].jac() for s in sd])),
            nmae=float(np.mean([S[(a, s)].nmae() for s in sd])),
            pds_lfc=float(np.mean([S[(a, s)].pds() for s in sd])),
            n_pred=float(np.median(S[(a, sd[0])].counts()[1])), n=len(sd))
    # a15 is the base for every arm and must be read on the SAME seeds as that arm, because
    # the level of every DE member drifts with the seed set.
    for a in have:
        d = de[a]
        b15 = float(np.mean([S[("a15", s)].reach() for s in seeds[a]]))
        print(f"{a:<11}{d['n']:>7}{d['n_pred']:>8.0f}{d['reach']:>10.5f}"
              f"{d['reach'] - b15:>+10.5f}{d['fid']:>9.5f}{d['jac']:>9.5f}"
              f"{d['nmae']:>9.5f}{d['pds_lfc']:>10.5f}")
    base = de["a15"]

    # ---- profile level: pds, cos, sse --------------------------------------------------
    print("\n" + "=" * 104)
    print("(2) profile level — 156's endpoints on the DELIVERED effect, panel targets "
          "excluded")
    print("=" * 104)
    z = np.load(f"{WORK}/profiles_all.npz", allow_pickle=True)
    b = z["b_ctrl"]
    perts = np.asarray([str(p) for p in z["perts"]])
    genes = np.asarray([str(g) for g in z["genes"]])
    rt = z["real_tot"]
    T = np.log1p(TS * rt / rt.sum(1, keepdims=True)) - b[None, :]
    gpos = {g: i for i, g in enumerate(genes)}
    drop = {gpos[p] for p in perts if p in gpos}
    mask = np.asarray([i for i in range(genes.size) if i not in drop])
    print(f"  {len(perts)} perturbations, {mask.size}/{genes.size} scored genes "
          f"({genes.size - mask.size} panel targets removed)")
    prof = {}
    print(f"\n{'arm':<11}{'pds':>10}{'d pds':>10}{'cos':>10}{'d cos %':>10}"
          f"{'row-row':>10}{'sse*':>11}{'d sse*':>11}")
    ep = {}
    for a in have:
        P = np.log1p(TS * z[f"p_{a}"].astype(np.float64)) - b[None, :]
        prof[a] = P
        ep[a] = endpoints(P[:, mask], T[:, mask])
    e0 = ep["a15"]
    for a in have:
        e = ep[a]
        print(f"{a:<11}{e['pds']:>10.5f}{e['pds'] - e0['pds']:>+10.5f}{e['cos']:>10.5f}"
              f"{100 * (e['cos'] / e0['cos'] - 1):>+9.2f}%{e['rr']:>10.5f}"
              f"{e['sse']:>11.6f}{e['sse'] - e0['sse']:>+11.6f}")

    # ---- how big is each correction? ---------------------------------------------------
    print("\n" + "=" * 104)
    print("(3) the correction itself — size, and where it lands relative to the reach head")
    print("=" * 104)
    a15p = z["p_a15"].astype(np.float64)
    keep = z["keep"]
    hp = np.load(f"{WORK}/head.npz", allow_pickle=True) if \
        os.path.exists(f"{WORK}/head.npz") else None
    print(f"{'arm':<11}{'|dp|/|p| %':>12}{'d log-profile L2':>18}"
          f"{'share of |d| on top-500 by p':>32}")
    for a in have:
        if a in ("a15", "a8"):
            continue
        dp = z[f"p_{a}"].astype(np.float64) - a15p
        rel = np.median(np.linalg.norm(dp, axis=1) / np.linalg.norm(a15p, axis=1))
        dl = np.median(np.linalg.norm(prof[a] - prof["a15"], axis=1))
        share = ""
        if hp is not None:
            H = hp["head"]           # (T, G_full) int index order, ascending pred p
            s = []
            for i in range(dp.shape[0]):
                d = np.abs(prof[a][i] - prof["a15"][i])
                s.append(d[H[i][:500]].sum() / max(d.sum(), 1e-30))
            share = f"{100 * np.mean(s):.1f}%"
        print(f"{a:<11}{100 * rel:>12.4f}{dl:>18.5f}{share:>32}")

    # ---- board conversion --------------------------------------------------------------
    print("\n" + "=" * 104)
    print("(4) BOARD CONVERSION — every member calibrated on `mn`, which IS a23")
    print("=" * 104)
    dr_mn = de["mn"]["reach"] - float(np.mean([S[("a15", s)].reach()
                                               for s in seeds["mn"]]))
    dp_mn = ep["mn"]["pds"] - e0["pds"]
    ds_mn = ep["mn"]["sse"] - e0["sse"]
    dn_mn = de["mn"]["nmae"] - base["nmae"]
    df_mn = de["mn"]["fid"] - base["fid"]
    dj_mn = de["mn"]["jac"] - base["jac"]
    print(f"  calibration (h1-test move for `mn` -> a23's realised board move):")
    print(f"    reach {dr_mn:+.5f} -> {A23['reach']:+.6f}   "
          f"(naive SPAN conversion would give {dr_mn / SPAN['reach']:+.6f}, "
          f"{dr_mn / SPAN['reach'] / A23['reach']:.2f}x)")
    print(f"    pds   {dp_mn:+.5f} -> {A23['pds']:+.6f}")
    print(f"    sse*  {ds_mn:+.6f} -> {A23['mse']:+.6f}")
    for nm, dv, key in (("fid", df_mn, "fid"), ("jac", dj_mn, "jac"),
                        ("nmae", dn_mn, "nmae")):
        print(f"    {nm:<5} {dv:+.5f} -> {A23[key]:+.6f}")

    def cal(x, x_mn, board):
        if abs(x_mn) < 1e-12:
            return 0.0
        return x / x_mn * board

    print(f"\n{'arm':<11}{'pds':>11}{'mse':>11}{'nmae':>11}{'fid':>11}{'reach':>11}"
          f"{'jac':>11}{'SUM/6':>12}")
    rows = {}
    for a in have:
        if a in ("a15", "a8"):
            continue
        m = dict(
            pds=cal(ep[a]["pds"] - e0["pds"], dp_mn, A23["pds"]),
            mse=cal(ep[a]["sse"] - e0["sse"], ds_mn, A23["mse"]),
            nmae=cal(de[a]["nmae"] - base["nmae"], dn_mn, A23["nmae"]),
            fid=cal(de[a]["fid"] - base["fid"], df_mn, A23["fid"]),
            reach=cal(de[a]["reach"] - float(np.mean([S[("a15", s)].reach()
                                                      for s in seeds[a]])),
                      dr_mn, A23["reach"]),
            jac=cal(de[a]["jac"] - base["jac"], dj_mn, A23["jac"]),
        )
        rows[a] = m
        tot = sum(m.values()) / 6.0
        print(f"{a:<11}{m['pds']:>+11.6f}{m['mse']:>+11.6f}{m['nmae']:>+11.6f}"
              f"{m['fid']:>+11.6f}{m['reach']:>+11.6f}{m['jac']:>+11.6f}{tot:>+12.6f}")
    print(f"\n  `mn` reproduces a23's realised overall {sum(A23.values()) / 6:+.6f} "
          f"by construction — it is the calibration point, not a prediction.")
    np.savez_compressed(f"{WORK}/priced.npz",
                        arms=np.asarray(have),
                        **{f"de_{k}": np.asarray([de[a][k] for a in have])
                           for k in ("reach", "fid", "jac", "nmae", "pds_lfc", "n_pred")},
                        **{f"ep_{k}": np.asarray([ep[a][k] for a in have])
                           for k in ("pds", "cos", "rr", "sse")})
    print(f"\nwrote {WORK}/priced.npz")


if __name__ == "__main__":
    main()
