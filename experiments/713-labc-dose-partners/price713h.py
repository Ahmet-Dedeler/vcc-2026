#!/usr/bin/env python3
"""713 part 2 (science [d80a10], 10-16; the coordinator's fallback slot (2)): H1 x2.5 and x3 at CD4 x1 on LABC, against the
October composition h2c1 under the incumbent key, in score708's exact arithmetic (its A:<design> arm). Registered before any
number (REGISTRATION-H1.json); the dumps are the lab's LABCd rebuild (dump801 --designs h2.5c1,h3c1).

Arms: A:h2.5c1 and A:h3c1 = weights wp (percentile-cap key), composition fin_<design> vs fin (h2c1); q 1.75 row sharpening,
minnorm render, ballast KMAX 16; 522 units, k 1, ret 1, mean over the 3 folds.
Known answers (abort unless they hold): KA0 the incumbent's cos_agg / pds per fold == the lab's stored LABC `pct` arm
(|d| < 1e-12, as part 1); KA1 the array-recomputed price == score708's price for every row (|d| < 1e-12).
CI: paired target bootstrap within folds, 4,000 draws, seed 7132; mean over folds. Best-of-N lower bound = CI low - e(N) sd,
with the lab's N 4 (e 1.029) for the H1/CD4 design family.
Decision rule (registered): HARM on LABC = mean <= -0.0005 with >= 2/3 folds negative. A:h3c1 gets slot (2) only if it is NOT
HARM here AND 710's board read is not HARM. If both points read <= 0 on the mean, slot (3) takes the slot.
usage: uv run python experiments/713-labc-dose-partners/price713h.py        (through gate602; < 3 GB)
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/708-def-panel-repricing"))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

import score708 as S7  # noqa: E402
import fold522 as FD  # noqa: E402
from price713 import agg, price_ag  # noqa: E402  (part 1's array arithmetic, KA'd there)

DUMPS = ROOT / "experiments/801-lab-folds/work/dumps"
META = ROOT / "experiments/801-lab-folds/work/ka-labc/meta"
REG = json.loads((HERE / "REGISTRATION-H1.json").read_text())
FOLDS = S7.FOLDS
NB, SEED, EN = 4000, 7132, 1.029
DESIGNS = ("h2.5c1", "h3c1")


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 23), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    t0 = time.time()
    for f in FOLDS:
        assert sha256(DUMPS / f"dump708-LABCd-{f}.npz") == REG["dumps"][f], f"LABCd/{f}: sha differs from the registration"
    assert sha256(ROOT / REG["scorer"]["path"]) == REG["scorer"]["sha256"], "score708.py differs from the registration"
    per, arrays, dvec = {}, {}, {}
    for f in FOLDS:
        S = S7.Fold(DUMPS / f"dump708-LABCd-{f}.npz")
        assert set(DESIGNS) <= set(S.extra), f"LABCd/{f}: designs {sorted(S.extra)}"
        wp = FD.wclip(S.kv, S.carried, 0.5, FD.pct_hi(S.kv, S.carried, FD.P_PROD))
        inc_s, inc_a, _ = S.score(wp)
        m = json.loads((META / f"LABC__{f}.json").read_text())["arms"]["pct"]
        d0 = max(abs(inc_s["cos_agg"] - m["cos_agg"]), abs(inc_s["pds"] - m["pds"]))
        assert d0 < 1e-12, f"KA0 LABCd/{f}: incumbent vs the lab's pct arm |d| {d0}"
        dvec[f] = (S.truth * S.truth).sum(1)
        allidx = np.arange(S.n)
        per[f], arrays[f] = {}, {"inc": inc_a}
        for d in DESIGNS:
            sa, aa, _ = S.score(wp, fin=S.extra[d])
            p = S7.price(sa, inc_s)
            p_ag = price_ag(agg(aa, dvec[f], allidx), agg(inc_a, dvec[f], allidx))
            assert abs(p - p_ag) < 1e-12, f"KA1 LABCd/{f} A:{d}: {p!r} vs {p_ag!r}"
            per[f][f"A:{d}"] = p
            arrays[f][f"A:{d}"] = aa
        print(f"LABCd/{f}: KA0 |d| {d0:.1e}; " + " | ".join(f"{k} {v:+.5f}" for k, v in per[f].items())
              + f" ({time.time() - t0:.0f}s)", flush=True)
    print("KA0 PASS (incumbent == the lab's LABC pct arm on 3 folds); KA1 PASS (array price == score708 price)", flush=True)
    rng = np.random.default_rng(SEED)
    names = [f"A:{d}" for d in DESIGNS]
    boot = {k: np.empty(NB) for k in names + ["A:h3c1-A:h2.5c1"]}
    for b in range(NB):
        acc = {k: 0.0 for k in boot}
        for f in FOLDS:
            idx = rng.integers(0, len(dvec[f]), len(dvec[f]))
            inc = agg(arrays[f]["inc"], dvec[f], idx)
            pb = {k: price_ag(agg(arrays[f][k], dvec[f], idx), inc) for k in names}
            for k in names:
                acc[k] += pb[k] / len(FOLDS)
            acc["A:h3c1-A:h2.5c1"] += (pb["A:h3c1"] - pb["A:h2.5c1"]) / len(FOLDS)
        for k in boot:
            boot[k][b] = acc[k]
    out = {"registration": "experiments/713-labc-dose-partners/REGISTRATION-H1.json", "units": "522 board units, k 1, ret 1",
           "per_fold": per, "configs": {}, "nb": NB, "seed": SEED, "e_N": EN}
    for k in names + ["A:h3c1-A:h2.5c1"]:
        if k in names:
            v = np.array([per[f][k] for f in FOLDS])
        else:
            v = np.array([per[f]["A:h3c1"] - per[f]["A:h2.5c1"] for f in FOLDS])
        ci = [float(np.percentile(boot[k], 2.5)), float(np.percentile(boot[k], 97.5))]
        sd = float(boot[k].std(ddof=1))
        harm = bool(v.mean() <= -0.0005 and (v < 0).sum() >= 2)
        out["configs"][k] = dict(per_fold=dict(zip(FOLDS, map(float, v))), mean=float(v.mean()), n_pos=int((v > 0).sum()),
                                 n_neg=int((v < 0).sum()), ci95=ci, sd=sd, best_of_n_lb=ci[0] - EN * sd,
                                 harm_on_labc=harm if k in names else None)
        r = out["configs"][k]
        print(f"{k:16s}: mean {r['mean']:+.5f} ({r['n_pos']}/3 > 0) per fold " + " / ".join(f"{x:+.5f}" for x in v)
              + f" CI95 [{ci[0]:+.5f}, {ci[1]:+.5f}] sd {sd:.5f} best-of-N lb {r['best_of_n_lb']:+.5f}"
              + (f" HARM-on-LABC {harm}" if k in names else " (paired step x2.5 -> x3)") + " (522 units)", flush=True)
    a3, a25 = out["configs"]["A:h3c1"], out["configs"]["A:h2.5c1"]
    if a3["mean"] <= 0 and a25["mean"] <= 0:
        verdict = "both points <= 0 on LABC: slot (3) takes the slot"
    elif a3["harm_on_labc"]:
        verdict = "A:h3c1 is HARM on LABC: no slot (2)"
    else:
        verdict = "A:h3c1 not HARM on LABC: slot (2) eligible IF 710's board read is not HARM"
    out["verdict"] = verdict
    print("VERDICT:", verdict, flush=True)
    out["seconds"] = round(time.time() - t0, 1)
    (HERE / "work" / "price713h.json").write_text(json.dumps(out, indent=1))
    print("wrote work/price713h.json")


if __name__ == "__main__":
    main()
