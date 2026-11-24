"""518 — exact PLOLO screen of bulk-mse levers on 446bj2's fold bulk (score518.Fold; known answer inside Fold()).

Every arm is ONE variable vs the base V1_0.5_3 (446bj2's bulk: stack composition, LODO key clip (0.5, 3), q 1.75, 10-CPM
render gate, realised at KMAX 16). Reported per fold: r = cos_agg/cos_agg(base) - 1 (at the board's cliff d raw mse =
-0.406 k r), literal d pds, board-emulated d pds (lam 0.5), mean ballast K.

    uv run python experiments/518-bulk-mse-gap/screen518.py <arm-set> [folds...]      arm sets: gate | strength | centre | genecls | all
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

from score518 import Fold, gate_at  # noqa: E402

WORK = HERE / "work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
CANCER = FOLDS[:3]


def strength(F):
    """the keyed stack row norm on non-panel genes (what the key and the energy constant leave each target)."""
    from lib420 import key_scale
    K = key_scale(F.fST, F.w_base)[0]
    return np.linalg.norm(K[:, F.nonpanel], axis=1)


def arms_for(F, which):
    st = strength(F)
    rel = st / np.median(st)
    pre = np.linalg.norm(F.fST[:, F.nonpanel], axis=1)
    prel = pre / np.median(pre)
    A = {}
    if which in ("gate", "all"):
        for m in (2.0, 5.0, 7.0, 15.0, 20.0):
            A[f"gate{m:g}"] = dict(min_cpm=m)
    if which in ("strength", "all"):
        for be in (-0.25, 0.25, 0.5, 1.0):
            A[f"tilt_keyed^{be:+g}"] = dict(w=F.w_base * rel ** be)
        for be in (-0.25, 0.25, 0.5):
            A[f"tilt_prekey^{be:+g}"] = dict(w=F.w_base * prel ** be)
        top = st >= np.quantile(st, 0.9)
        for g in (0.8, 1.25, 1.5):
            A[f"top10x{g:g}"] = dict(w=F.w_base * np.where(top, g, 1.0))
        bot = st <= np.quantile(st, 0.5)
        for g in (0.5, 0.75):
            A[f"bottom50x{g:g}"] = dict(w=F.w_base * np.where(bot, g, 1.0))
    if which in ("centre", "all"):
        A["centre_full"] = dict(centre="effect")
        A["centre_half"] = dict(centre=("shrink", 0.5))
        A["mean_x1.5"] = dict(centre=("shrink", -0.5))
    if which in ("genecls", "all"):
        for gam in (-0.1, 0.1, 0.2, 0.35):
            gm = {}
            for c, (b, g10) in F.fr.items():
                cpm = 1e6 * np.expm1(b) / np.expm1(b).sum()
                med = np.median(cpm[g10])
                gm[c] = np.where(g10, (np.maximum(cpm, 1e-12) / med) ** gam, 1.0)
            A[f"genecpm^{gam:+g}"] = dict(gene_mult=gm)
    if which in ("q", "all"):
        for q in (1.5, 2.0, 2.25):
            A[f"q{q:g}"] = dict(q=q)
    return A


def main():
    which = sys.argv[1]
    folds = sys.argv[2:] or FOLDS
    out = {}
    path = WORK / f"screen518-{which}.json"
    if path.exists():
        out = json.loads(path.read_text())
    for f in folds:
        t0 = time.time()
        F = Fold(f)
        base, _ = F.score()
        rec = dict(base=dict(cos_agg=base["cos_agg"], pds=base["pds"], pds_lam=base["pds_lam0.5"], K=base["K_mean"]), arms={})
        for name, kw in arms_for(F, which).items():
            s, _ = F.score(**kw)
            rec["arms"][name] = dict(r=s["cos_agg"] / base["cos_agg"] - 1, dpds=s["pds"] - base["pds"],
                                     dpds_bd=s["pds_lam0.5"] - base["pds_lam0.5"], K=s["K_mean"], dK=s["K_mean"] - base["K_mean"])
            a = rec["arms"][name]
            print(f"[{f}] {name:22s} r {100 * a['r']:+7.3f} %  dpds lit {a['dpds']:+.5f} bd {a['dpds_bd']:+.5f}  K {a['K']:.2f} "
                  f"({time.time() - t0:.0f}s)", flush=True)
        out[f] = rec
        path.write_text(json.dumps(out, indent=1, default=float))
    names = list(next(iter(out.values()))["arms"])
    print(f"\n{'arm':22s} {'cancer r':>9s} {'n+':>3s} " + " ".join(f"{f[:8]:>9s}" for f in FOLDS) + "   lit dpds (cancer)  bd dpds  dK")
    for nm in names:
        rs = {f: out[f]["arms"][nm]["r"] for f in FOLDS if f in out}
        canc = [rs[f] for f in CANCER if f in rs]
        dp = np.mean([out[f]["arms"][nm]["dpds"] for f in CANCER if f in out])
        db = np.mean([out[f]["arms"][nm]["dpds_bd"] for f in CANCER if f in out])
        dk = np.mean([out[f]["arms"][nm]["dK"] for f in CANCER if f in out])
        print(f"{nm:22s} {100 * np.mean(canc):+8.3f}% {sum(r > 0 for r in canc):3d} " +
              " ".join(f"{100 * rs.get(f, np.nan):+8.3f}%" for f in FOLDS) + f"   {dp:+.5f}  {db:+.5f}  {dk:+.2f}")


if __name__ == "__main__":
    main()
