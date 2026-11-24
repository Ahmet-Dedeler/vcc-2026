"""234 Phase 1 + Phase 2 — price every key we hold as a SIGN-FLIP arm, exactly.

THE ALGEBRA (why no scoring run is needed for `q_w`)
---------------------------------------------------
A flip arm reverses the sign of a set `S` of (target, gene) coordinates at bit-identical
magnitude.  Its effect on the |t|-weighted sign accuracy is closed-form:

    Dq_w = mean_i [ sum_{j in pool_i & S_i} w_ij (1 - 2 c_ij) / den_i ]

so the arm pays IFF the key selects a set whose weighted sign accuracy is BELOW 0.5.  This
script computes that for every key, every dose and three regions, and converts with 232's
exchange rate (+0.00399 overall per +0.01 `q_w`, k_cos 1.384; +0.00379 at 1.000).

REGIONS
  all         every scored coordinate -- a flip key is NOT a donor-combination rule, so it
              is not confined to the contested set the way 232's vote was.  232's +0.17325
              ceiling is the ceiling for a VOTE; a flip key's ceiling is `1 - q_w` = +0.42.
  contested   >= 2 covering donors AND they disagree in sign (232's disagreement set)
  unanimous   >= 2 covering donors, all agreeing

DOSE.  The selection uses ONLY the key -- never the pool, never the truth.  The pool is the
reference's DE set and is not knowable at prediction time.  Both orientations of every key
are swept (low-key-first and high-key-first) because a "confidence" can point either way,
and the better of the two is reported as FITTED -- deliberately generous, because the
headline result here is a negative one and a negative must be argued at the arm's best.

PHASE 2 (the ceiling).  A perfectly calibrated version of a key flips exactly those bins of
the key where P(correct | key) < 0.5.  In-fold calibration is the ceiling; leave-one-fold-out
calibration is what is actually reachable.  The gap between them separates "the key is weak"
from "the family is weak".

    uv run python screen234.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common234 import (ALL_FOLDS, CANCER, GATE, KA232_CEILING, QW_NEEDED,  # noqa: E402
                       RATE_PER_QW, SEED, WORK, acc_in, mass_share, prep_weights, qw)

POOL = "wald89"
FRACS = (0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20,
         0.30, 0.40, 0.50, 0.70, 1.00)
REGIONS = ("all", "contested", "unanimous")
NBINS = (20, 50, 200)
# keys whose a-priori orientation is "HIGH = we are more likely RIGHT".  The screen sweeps
# both orientations anyway; this is recorded so the fitted-vs-prior distinction is auditable.
PRIOR_HIGH_IS_CONFIDENT = {
    "absH": True, "absfin": True, "abspre": True, "key189": True, "theta": True,
    "expr": True, "cpm": True, "zpred": True, "absH_rank": True, "whitenkept": True,
    "rownorm": True, "ncov": True, "margin_maj": True, "margin_z": True, "coh": True,
    "tdonor": True, "minz": True, "sumz": True,
}


# --------------------------------------------------------------------------- fold loading
class Fold:
    def __init__(self, held: str):
        z = np.load(WORK / f"ex234-{held}.npz", allow_pickle=True)
        self.held = held
        self.truth = np.asarray(z["truth"], dtype=np.float64)
        self.correct = np.asarray(z["correct"], dtype=bool)
        self.sgn_del = np.asarray(z["sgn_del"], dtype=np.int8)
        self.s_mean = np.asarray(z["s_mean"], dtype=np.int8)
        self.disagree = np.asarray(z["disagree"], dtype=bool)
        self.pool = np.asarray(z[f"pool_{POOL}"], dtype=bool)
        self.keys = {k[4:]: np.asarray(z[k], dtype=np.float32)
                     for k in z.files if k.startswith("key_")}
        self.ncov = self.keys["ncov"]
        self.w, self.den = prep_weights(self.truth, self.pool)
        self.live = self.den > 0
        self.q_w = qw(self.sgn_del, self.truth, self.w, self.den)
        self.regions = dict(
            all=np.ones(self.truth.shape, dtype=bool),
            contested=(self.ncov >= 2) & self.disagree,
            unanimous=(self.ncov >= 2) & ~self.disagree)
        # `g` is the per-coordinate contribution to Dq_w if that coordinate is flipped
        self.g = (self.w * (1.0 - 2.0 * self.correct)) / np.where(self.den > 0, self.den,
                                                                 np.inf)[:, None]


def curve(f: Fold, key: np.ndarray, region: str, sign: int) -> np.ndarray:
    """Dq_w as a function of the dose, for `sign * key` ascending (flip lowest first).

    Returns the fold-mean cumulative Dq_w evaluated on `FRACS`.  Rows outside the region
    are pushed to the tail by +inf so a prefix of the sorted row is always inside it.
    """
    A = f.regions[region]
    k = np.where(A, sign * key.astype(np.float64), np.inf)
    order = np.argsort(k, axis=1, kind="stable")
    gs = np.take_along_axis(f.g, order, axis=1)
    cs = np.concatenate([np.zeros((gs.shape[0], 1)), np.cumsum(gs, axis=1)], axis=1)
    m = A.sum(1)
    out = []
    for fr in FRACS:
        idx = np.minimum(np.rint(fr * m).astype(int), m)
        out.append(float(np.mean(cs[np.arange(cs.shape[0]), idx][f.live])))
    return np.asarray(out)


def calibrated_ceiling(f: Fold, key: np.ndarray, nbins: int, ref=None):
    """Phase 2 — the key perfectly calibrated: flip the bins where P(correct|key) < 0.5.

    `ref` is an (edges, flip_decision) pair from OTHER folds (leave-one-fold-out); when it
    is None the calibration is fitted in-fold, which is the CEILING and is optimistic.
    """
    pool = f.pool
    kp = key[pool].astype(np.float64)
    if ref is None:
        edges = np.unique(np.quantile(kp, np.linspace(0, 1, nbins + 1)))
        edges[0], edges[-1] = -np.inf, np.inf
        b = np.digitize(kp, edges[1:-1])
        wq = f.w[pool]
        cq = f.correct[pool].astype(np.float64)
        num = np.bincount(b, weights=wq * cq, minlength=len(edges) - 1)
        den = np.bincount(b, weights=wq, minlength=len(edges) - 1)
        acc = np.where(den > 0, num / np.maximum(den, 1e-300), 1.0)
        flip_bin = acc < 0.5
    else:
        edges, flip_bin = ref
        b = np.digitize(kp, edges[1:-1])
        flip_bin = flip_bin[:len(edges) - 1]
    sel = np.zeros(f.truth.shape, dtype=bool)
    sel[pool] = flip_bin[np.minimum(b, len(flip_bin) - 1)]
    d = float(np.mean((f.g * sel).sum(1)[f.live]))
    return d, (edges, flip_bin), sel


def fit_calibration(folds, key_name: str, nbins: int):
    """Pool the OTHER folds' (key, correctness) into one calibration curve."""
    kp = np.concatenate([f.keys[key_name][f.pool].astype(np.float64) for f in folds])
    wq = np.concatenate([f.w[f.pool] for f in folds])
    cq = np.concatenate([f.correct[f.pool].astype(np.float64) for f in folds])
    edges = np.unique(np.quantile(kp, np.linspace(0, 1, nbins + 1)))
    edges[0], edges[-1] = -np.inf, np.inf
    b = np.digitize(kp, edges[1:-1])
    num = np.bincount(b, weights=wq * cq, minlength=len(edges) - 1)
    den = np.bincount(b, weights=wq, minlength=len(edges) - 1)
    acc = np.where(den > 0, num / np.maximum(den, 1e-300), 1.0)
    return edges, acc < 0.5


# ----------------------------------------------------------------------------- KA-4 / KA-5
def ka_consensus(folds) -> dict:
    """232 §2's contested-mass arithmetic, at the CONSENSUS level, reproduced here."""
    rows = {}
    for f in folds:
        dis = f.disagree & (f.ncov >= 2)
        share = mass_share(dis, f.w, f.den)
        q_all = qw(f.s_mean, f.truth, f.w, f.den)
        wd = np.where(dis, f.w, 0.0)
        dend = wd.sum(1)
        liv = dend > 0
        okm = (f.s_mean != 0) & (f.s_mean == np.sign(f.truth))
        q_dis = float(np.mean((wd * okm).sum(1)[liv] / dend[liv]))
        rows[f.held] = dict(disagree_mass=share, q_mean_all=q_all, q_mean_contested=q_dis,
                            ceiling_vote=share * (1.0 - q_dis))
    cm = {k: float(np.mean([rows[f.held][k] for f in folds if f.held in CANCER]))
          for k in ("disagree_mass", "q_mean_all", "q_mean_contested", "ceiling_vote")}
    ok = (abs(cm["disagree_mass"] - KA232_CEILING["disagree_mass"]) < 5e-4
          and abs(cm["q_mean_contested"] - KA232_CEILING["q_contested"]) < 5e-4
          and abs(cm["ceiling_vote"] - KA232_CEILING["ceiling"]) < 1e-3)
    print(f"[KA-4] 232 §2 consensus arithmetic, 3 cancer folds: disagree mass "
          f"{cm['disagree_mass']:.4f} (232: 0.3627)   sign(mean) on contested "
          f"{cm['q_mean_contested']:.5f} (232: 0.52292)   vote ceiling "
          f"{cm['ceiling_vote']:+.5f} (232: +0.17325)   {'PASS' if ok else 'FAIL'}")
    assert ok, "KA-4 FAILED: 232 §2's contested arithmetic does not reproduce"
    return dict(per_fold=rows, cancer_mean=cm, PASS=True)


# ------------------------------------------------------------------------------- the screen
def main() -> None:
    t0 = time.time()
    folds = {h: Fold(h) for h in ALL_FOLDS}
    C = [folds[h] for h in CANCER]
    audit = {"pool": POOL, "fracs": list(FRACS), "regions": list(REGIONS),
             "rate": RATE_PER_QW, "gate": GATE, "qw_needed": QW_NEEDED}

    print("\nBASE (delivered, wald89):")
    for h in ALL_FOLDS:
        f = folds[h]
        print(f"  {h:<15} q_w {f.q_w:.5f}   pool-oracle headroom 1-q_w {1-f.q_w:+.5f}   "
              f"contested mass {mass_share(f.regions['contested'], f.w, f.den):.4f}   "
              f"acc(contested) {acc_in(f.regions['contested'], f.correct, f.w):.5f}   "
              f"acc(unanimous) {acc_in(f.regions['unanimous'], f.correct, f.w):.5f}")
    audit["base"] = {h: dict(
        q_w=folds[h].q_w, headroom=1 - folds[h].q_w,
        mass=dict((r, mass_share(folds[h].regions[r], folds[h].w, folds[h].den))
                  for r in REGIONS),
        acc=dict((r, acc_in(folds[h].regions[r], folds[h].correct, folds[h].w))
                 for r in REGIONS)) for h in ALL_FOLDS}
    audit["KA_4"] = ka_consensus(list(folds.values()))

    keynames = sorted(C[0].keys)
    audit["keys"] = keynames

    # ---------------- Phase 1: the dose curves
    curves = {}
    for kn in keynames:
        for region in REGIONS:
            for sg in (+1, -1):
                per_fold = {f.held: curve(f, f.keys[kn], region, sg) for f in folds.values()}
                curves[(kn, region, sg)] = per_fold
        print(f"  swept {kn} ({time.time()-t0:.0f}s)", flush=True)

    def cancer_mean(pf):
        return np.mean([pf[h] for h in CANCER], axis=0)

    rows = []
    for kn in keynames:
        for region in REGIONS:
            best = None
            for sg in (+1, -1):
                pf = curves[(kn, region, sg)]
                cm = cancer_mean(pf)
                i = int(np.argmax(cm))
                # leave-one-fold-out dose: pick the frac on the other two cancer folds
                lofo = []
                for h in CANCER:
                    other = np.mean([pf[o] for o in CANCER if o != h], axis=0)
                    lofo.append(pf[h][int(np.argmax(other))])
                cand = dict(key=kn, region=region, orientation=("low_first" if sg > 0
                                                                else "high_first"),
                            prior_orientation=("low_first"
                                               if PRIOR_HIGH_IS_CONFIDENT.get(kn, True)
                                               else "high_first"),
                            best_frac=FRACS[i], dqw_infold=float(cm[i]),
                            dqw_lofo=float(np.mean(lofo)),
                            dqw_per_fold={h: float(pf[h][i]) for h in ALL_FOLDS},
                            curve={h: [float(v) for v in pf[h]] for h in ALL_FOLDS})
                if best is None or cand["dqw_infold"] > best["dqw_infold"]:
                    best = cand
            best["overall_infold"] = {k: best["dqw_infold"] * v
                                      for k, v in RATE_PER_QW.items()}
            best["overall_lofo"] = {k: best["dqw_lofo"] * v for k, v in RATE_PER_QW.items()}
            best["clears_gate"] = bool(best["dqw_lofo"] >= QW_NEEDED["k1.384"])
            rows.append(best)

    # ---------------- Phase 2: the perfectly calibrated ceiling
    ceil_rows = []
    for kn in keynames:
        r = dict(key=kn)
        for nb in NBINS:
            inf = [calibrated_ceiling(f, f.keys[kn], nb)[0] for f in C]
            lo = []
            for h in CANCER:
                ref = fit_calibration([folds[o] for o in CANCER if o != h], kn, nb)
                lo.append(calibrated_ceiling(folds[h], folds[h].keys[kn], nb, ref)[0])
            r[f"ceiling_infold_{nb}"] = float(np.mean(inf))
            r[f"ceiling_lofo_{nb}"] = float(np.mean(lo))
            r[f"ceiling_lofo_{nb}_per_fold"] = {h: float(v) for h, v in zip(CANCER, lo)}
        ceil_rows.append(r)

    # ---------------- the two reference ceilings
    orc = float(np.mean([1.0 - f.q_w for f in C]))
    orc_cont = float(np.mean([
        mass_share(f.regions["contested"], f.w, f.den)
        * (1.0 - acc_in(f.regions["contested"], f.correct, f.w)) for f in C]))
    audit["oracle_full_dqw"] = orc
    audit["oracle_contested_dqw"] = orc_cont

    # ---------------- print
    print(f"\n{'='*118}\nPHASE 1 — every key as a sign-flip arm.  "
          f"GATE: d q_w >= +{QW_NEEDED['k1.384']:.5f} (k1.384) .. "
          f"+{QW_NEEDED['k1.000']:.5f} (k1.000)\n{'='*118}")
    hdr = (f"{'key':<12}{'region':<11}{'orient':<11}{'f*':>7}{'dqw(fit)':>11}"
           f"{'overall':>10}{'dqw(LOFO)':>11}{'overall':>10}{'x short':>9}{'k+':>5}")
    print(hdr)
    print("-" * len(hdr))
    for r in sorted(rows, key=lambda x: -x["dqw_lofo"]):
        sh = QW_NEEDED["k1.384"] / r["dqw_lofo"] if r["dqw_lofo"] > 0 else float("inf")
        npos = sum(1 for h in CANCER if r["dqw_per_fold"][h] > 0)
        print(f"{r['key']:<12}{r['region']:<11}{r['orientation']:<11}"
              f"{r['best_frac']:>7.4f}{r['dqw_infold']:>+11.5f}"
              f"{r['overall_infold']['k1.384']:>+10.5f}{r['dqw_lofo']:>+11.5f}"
              f"{r['overall_lofo']['k1.384']:>+10.5f}"
              f"{(f'{sh:.0f}x' if np.isfinite(sh) else '-'):>9}{npos:>4d}/3")

    print(f"\n{'='*100}\nPHASE 2 — the same keys PERFECTLY CALIBRATED (flip every bin with "
          f"P(correct|key) < 0.5)\n{'='*100}")
    hdr2 = (f"{'key':<12}" + "".join(f"{'infold@'+str(n):>13}" for n in NBINS)
            + "".join(f"{'LOFO@'+str(n):>13}" for n in NBINS))
    print(hdr2)
    print("-" * len(hdr2))
    for r in sorted(ceil_rows, key=lambda x: -x[f"ceiling_lofo_{NBINS[1]}"]):
        print(f"{r['key']:<12}"
              + "".join(f"{r[f'ceiling_infold_{n}']:>+13.5f}" for n in NBINS)
              + "".join(f"{r[f'ceiling_lofo_{n}']:>+13.5f}" for n in NBINS))
    print(f"\n  FULL ORACLE (flip every wrong pool coordinate)      d q_w {orc:+.5f}  "
          f"= {orc*RATE_PER_QW['k1.384']:+.5f} overall")
    print(f"  ORACLE restricted to the CONTESTED set (232's vote ceiling) "
          f"d q_w {orc_cont:+.5f}  = {orc_cont*RATE_PER_QW['k1.384']:+.5f} overall")

    audit["phase1"] = rows
    audit["phase2"] = ceil_rows
    audit["seconds"] = time.time() - t0
    (WORK / "screen234.json").write_text(json.dumps(audit, indent=1))
    print(f"\n-> {WORK/'screen234.json'}  ({audit['seconds']:.0f}s)")


if __name__ == "__main__":
    main()
