"""188 — the four DE members for this experiment's own h1-test arms, in board units.

Conventions are 169 §4's, verbatim, so its five-seed `q1.00_mf` ladder is the matched control:
`reach` at native BH with 157's board transfer factor from the a8<->a15 pair, `fid` read as
PRECISION (no `n_real`), `jac` rebuilt on the board's own counts from recall, `nmae` from the
emitted tables.  Rule 9 is in force throughout.

ARMS PRICED HERE
  `g10_q1.15_s*` / `g20_q1.15_s*`  — 169's h-space arm with the transform restricted to the
      CPM>10 / CPM>20 sub-gate, i.e. to production's own gate width.  Five and three matched
      seeds.  This is the gate-width falsification test.
  `map_med_s1` / `map_bud_s1`      — the fitted magnitude maps, one seed (their `pds` is
      already NO-GO on PLOLO, so they are priced for completeness under rule 12, not to
      choose between them).
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
E169 = os.path.join(REPO, "experiments/169-sharpening-six-member")
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
sys.path.insert(0, E169)
from debench import ANCHOR, Surface       # noqa: E402
from six169 import overlap_ratio, precision   # noqa: E402

BOARD_REACH = dict(a8=0.1890, a15=0.1811)
H1_REACH = dict(a8=0.2029, a15=0.1828)
NMAE_B, NMAE_R = 1.0013, 0.400
NR, NP = 240.0, 4424.0


def main() -> None:
    paths = {}
    for p in sorted(glob.glob(os.path.join(E169, "de169.*.npz"))):
        paths[os.path.basename(p)[len("de169."):-4]] = p
    for p in sorted(glob.glob(os.path.join(HERE, "de188.*.npz"))):
        paths[os.path.basename(p)[len("de188."):-4]] = p

    LADDER = {
        "q1.15_mf (169, full gate)": (
            ["q1.15_mf", "q1.15_mf_b", "q1.15_mf_s3", "q1.15_mf_s4", "q1.15_mf_s5"],
            ["q1.00_mf", "q1.00_mf_s2", "q1.00_mf_s3", "q1.00_mf_s4", "q1.00_mf_s5"]),
        "g10_q1.15 (CPM>10 sub-gate)": (
            ["g10_q1.15_s1", "g10_q1.15_s2", "g10_q1.15_s3", "g10_q1.15_s4",
             "g10_q1.15_s5"],
            ["q1.00_mf", "q1.00_mf_s2", "q1.00_mf_s3", "q1.00_mf_s4", "q1.00_mf_s5"]),
        "g20_q1.15 (CPM>20 sub-gate)": (
            ["g20_q1.15_s1", "g20_q1.15_s2", "g20_q1.15_s3"],
            ["q1.00_mf", "q1.00_mf_s2", "q1.00_mf_s3"]),
        "q1.15L_mf (169, lfc space)": (
            ["q1.15L_mf", "q1.15L_mf_s2", "q1.15L_mf_s3", "q1.15L_mf_s4", "q1.15L_mf_s5"],
            ["q1.00_mf", "q1.00_mf_s2", "q1.00_mf_s3", "q1.00_mf_s4", "q1.00_mf_s5"]),
        "map_med (fitted median map)": (["map_med_s1"], ["q1.00_mf"]),
        "map_bud (fitted budget map)": (["map_bud_s1"], ["q1.00_mf"]),
    }

    cache = {}

    def read(a):
        if a not in cache:
            s = Surface(paths[a])
            nr, npd = s.counts()
            cache[a] = dict(reach=s.reach(), precision=precision(s), recall=overlap_ratio(s),
                            nmae=s.nmae(), pds=s.pds(), n_pred=float(np.median(npd)),
                            sign=s.sign_accuracy())
        return cache[a]

    k = ((BOARD_REACH["a8"] - BOARD_REACH["a15"]) / (H1_REACH["a8"] - H1_REACH["a15"]))
    sr = ANCHOR["reach"][1] - ANCHOR["reach"][0]
    sf = ANCHOR["fid"][1] - ANCHOR["fid"][0]
    sj = ANCHOR["jac"][1] - ANCHOR["jac"][0]
    sn = NMAE_B - NMAE_R
    print(f"board/h1-test `reach` transfer factor (a8<->a15): {k:.4f};  spans reach {sr:.4f} "
          f"fid {sf:.4f} jac {sj:.4f} nmae {sn:.4f}")
    print(f"\n{'arm':<30}{'n':>3}{'d reach raw':>12}{'sem':>9}{'d reach ovr':>12}"
          f"{'d fid ovr':>11}{'d jac ovr':>11}{'d nmae raw':>12}{'d nmae ovr':>12}"
          f"{'SUM DE':>10}")
    out = {}
    for label, (arm, base) in LADDER.items():
        pairs = [(read(a), read(b)) for a, b in zip(arm, base) if a in paths and b in paths]
        if not pairs:
            print(f"{label:<30}  (missing tables)")
            continue
        dr = np.array([a["reach"] - b["reach"] for a, b in pairs])
        dp = np.array([a["precision"] - b["precision"] for a, b in pairs])
        dc = np.array([a["recall"] - b["recall"] for a, b in pairs])
        dn = np.array([a["nmae"] - b["nmae"] for a, b in pairs])
        n = len(pairs)
        sem = float(dr.std(ddof=1) / np.sqrt(n)) if n > 1 else float("nan")
        j = lambda r: r * NR / (NP + NR - r * NR)     # noqa: E731
        dj = float(np.mean([j(a["recall"]) - j(b["recall"]) for a, b in pairs]))
        r_ovr = float(dr.mean()) * k / sr / 6
        f_ovr = float(dp.mean()) / sf / 6
        j_ovr = dj / sj / 6
        n_ovr = -float(dn.mean()) / sn / 6
        out[label] = dict(n=n, dreach=float(dr.mean()), sem=sem, reach_ovr=r_ovr,
                          fid_ovr=f_ovr, jac_ovr=j_ovr, dnmae=float(dn.mean()),
                          nmae_ovr=n_ovr, sum_de=r_ovr + f_ovr + j_ovr + n_ovr,
                          n_pred=float(np.mean([a["n_pred"] for a, _ in pairs])),
                          pds_lfc=float(np.mean([a["pds"] - b["pds"] for a, b in pairs])))
        print(f"{label:<30}{n:>3}{dr.mean():>+12.5f}{sem:>9.5f}{r_ovr:>+12.6f}"
              f"{f_ovr:>+11.6f}{j_ovr:>+11.6f}{dn.mean():>+12.5f}{n_ovr:>+12.6f}"
              f"{out[label]['sum_de']:>+10.6f}")

    print(f"\nlevels (median n_pred, raw nmae, h1-test DE-table pds — the lfc-space "
          f"statistic Phase 3 shows is the wrong space)")
    for a in sorted(cache):
        c = cache[a]
        print(f"  {a:<16} n_pred {c['n_pred']:>6.0f}  reach {c['reach']:.5f}  "
              f"nmae {c['nmae']:.5f}  prec {c['precision']:.5f}  recall {c['recall']:.5f}  "
              f"pds(lfc) {c['pds']:.5f}  sign {c['sign']:.5f}")
    with open(os.path.join(WORK, "six188.json"), "w") as fh:
        json.dump(dict(arms=out, levels=cache), fh, indent=1, default=float)
    print(f"\nwrote {WORK}/six188.json")


if __name__ == "__main__":
    sys.exit(main())
