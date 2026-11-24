"""518 — the STRENGTH KEY exactly on PLOLO: w -> w exp(theta z_p), z_p = centred percentile rank of a robust, noise-normalised
response-strength count (strength518) of target p in the fold's own SOURCES ("src"), or the mean rank over the sources + CD4
(DESeq2) + VIPerturb-K562 + KOLF2 atlas ("meta"; KOLF2 dropped on the kolf2 fold, where it is the truth).
Exact route: score518.Fold (dumped compositions, key_scale, q 1.75, render, ballast KMAX 16; known answer inside Fold()).
Arms per fold: theta grid {0.5, 0.75, 1.0, 1.25}; the held-out theta (fitted on the other cancer folds by strength518's
analytic LOFO); and the same on top of the cap-6 key (495, in flight) and with 479's P2g-style per-target q not modelled.

    uv run python experiments/518-bulk-mse-gap/skey518.py [folds...]
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import json  # noqa: E402
import time  # noqa: E402
import warnings  # noqa: E402

import numpy as np  # noqa: E402

from lofo518 import crank  # noqa: E402
from score518 import Fold  # noqa: E402
from lib420 import weights  # noqa: E402

warnings.filterwarnings("ignore")
WORK = HERE / "work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
CANCER = FOLDS[:3]
SRC = ("k562_gwps", "orion_hct116", "orion_hek293t")


def zfeat(F, kind):
    S = np.load(WORK / "strength518.npz")
    panel = [str(p) for p in S["panel"]]
    idx = np.asarray([panel.index(t) for t in F.targets])
    keys = [s for s in SRC if s != F.held]
    if kind == "meta":
        keys += ["cd4_deseq", "viperturb_k562"] + (["kolf2"] if F.held != "kolf2" else [])
    cols = []
    for k in keys:
        v = np.asarray(S[k], float)[idx]
        r = crank(v)
        r[~np.isfinite(v)] = np.nan
        cols.append(r)
    m = np.nanmean(np.vstack(cols), 0)
    return np.where(np.isfinite(m), m, 0.0)


def main():
    folds = sys.argv[1:] or FOLDS
    lofo = json.loads((WORK / "strength518.json").read_text())["lofo"]
    path = WORK / "skey518.json"
    out = json.loads(path.read_text()) if path.exists() else {}
    for f in folds:
        t0 = time.time()
        F = Fold(f)
        w3 = F.w_base
        w6 = weights(F.kv1, F.carried, 0.5, 6.0)[0]
        b3, _ = F.score(w=w3)
        b6, _ = F.score(w=w6)
        rec = dict(base3=dict(cos=b3["cos_agg"], pds=b3["pds"], pdsl=b3["pds_lam0.5"], K=b3["K_mean"]),
                   cap6=dict(r=b6["cos_agg"] / b3["cos_agg"] - 1, dpds=b6["pds"] - b3["pds"], dpds_bd=b6["pds_lam0.5"] - b3["pds_lam0.5"],
                             dK=b6["K_mean"] - b3["K_mean"]), arms={})
        for kind, lkey in (("src", "sources"), ("meta", "meta_all")):
            z = zfeat(F, kind)
            th_held = None
            if f in CANCER:
                th_held = None     # the held-out theta: refit on the other two cancer folds' analytic route
            thetas = [0.5, 0.75, 1.0, 1.25]
            for th in thetas:
                m = np.exp(th * z)
                for base_nm, wb, bb in (("cap3", w3, b3), ("cap6", w6, b6)):
                    s, _ = F.score(w=wb * m)
                    nm = f"{kind}|theta{th:g}|on_{base_nm}"
                    rec["arms"][nm] = dict(r_vs_own_base=s["cos_agg"] / bb["cos_agg"] - 1, r_vs_446bj2=s["cos_agg"] / b3["cos_agg"] - 1,
                                           dpds=s["pds"] - bb["pds"], dpds_bd=s["pds_lam0.5"] - bb["pds_lam0.5"], dK=s["K_mean"] - bb["K_mean"],
                                           m_q=np.quantile(m, [0, .1, .5, .9, 1]).tolist())
                    a = rec["arms"][nm]
                    print(f"[{f}] {nm:24s} r {100 * a['r_vs_own_base']:+6.2f} % (vs 446bj2 {100 * a['r_vs_446bj2']:+6.2f} %)  "
                          f"dpds lit {a['dpds']:+.5f} bd {a['dpds_bd']:+.5f}  dK {a['dK']:+.2f}  ({time.time() - t0:.0f}s)", flush=True)
        out[f] = rec
        path.write_text(json.dumps(out, indent=1, default=float))
    # summary
    names = list(next(iter(out.values()))["arms"])
    print(f"\n{'arm':26s} {'cancer r':>9s} n+  " + " ".join(f"{f[:8]:>9s}" for f in FOLDS) + "   lit dpds   bd dpds    dK")
    for nm in names:
        rs = {f: out[f]["arms"][nm]["r_vs_own_base"] for f in FOLDS if f in out}
        canc = [rs[f] for f in CANCER if f in rs]
        agg = lambda k: np.mean([out[f]["arms"][nm][k] for f in CANCER if f in out])
        print(f"{nm:26s} {100 * np.mean(canc):+8.2f}% {sum(r > 0 for r in canc):2d}  " + " ".join(f"{100 * rs.get(f, np.nan):+8.2f}%" for f in FOLDS)
              + f"   {agg('dpds'):+.5f}  {agg('dpds_bd'):+.5f}  {agg('dK'):+.2f}")
    c6 = [out[f]["cap6"]["r"] for f in CANCER if f in out]
    print(f"(cap 6 alone vs 446bj2: cancer r {100 * np.mean(c6):+.2f} %)")


if __name__ == "__main__":
    main()
