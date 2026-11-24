#!/usr/bin/env python3
"""713 (science [d80a10], 10-16): dose partners of 705's rider (theta) and 707's swap (beta) on LABC, the lab's Arc-drawn
held-out surface (A/B/C's 300 targets on the three cancer held lines with the day-recipe operators; LABC-PROTOCOL.md).
Registered before any number is read (REGISTRATION.json); the lab logs every configuration in its wear ledger.

Arithmetic: score708.py's, exactly (its Fold, zsrc_fold, unit_mean, swap_fin and price are imported, not copied): the
percentile-cap incumbent wp (P 96.647), q 1.75 row sharpening on each frame's gate, minnorm render, ballast KMAX 16, scored
on the fold's truth; price in 522 units, k 1, ret 1: 0.0710 * (X_arm / X_inc - 1) + 0.1245 * (pds_arm - pds_inc).
This is the arithmetic 708 runs on D/E/F's panel on 10-22. The lab's own 705 / 707 rows on LABC use inc = pct+P2g instead,
so they are not this instrument's known answer; its stored `pct` arm is.

Configurations (each against the incumbent; the partners also paired against their family's reference dose):
  s_b0.125 / s_b0.25 (reference) / s_b0.5 : 707's swap on the top z_src tertile at beta 0.125 / 0.25 / 0.5
  r_t0.25 (reference) / r_t0.5             : 522's floor-preserving rider w' = max(0.5, w exp(theta z)) above the floor
Known answers (abort unless they hold):
  KA0 the incumbent's cos_agg and pds on each LABC fold == the lab's stored `pct` arm (ka-labc/meta), |d| < 1e-12;
  KA1 the price recomputed from the per-target arrays (cos_agg per frame = |sum x| / sqrt(sum a * sum d), pds = mean of the
      per-target ranks) == score708's price from its summaries, |d| < 1e-12, for every configuration and fold.
CI: paired target bootstrap within folds (the same resampled targets for arm and incumbent), 4,000 draws, seed 713, the
statistic being the mean over the three folds. Reported, never decided here; the decision rule is in REGISTRATION.json.
usage: uv run python experiments/713-labc-dose-partners/price713.py        (through gate602; < 3 GB)
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
import numpy as np  # noqa: E402

import score708 as S7  # noqa: E402
import fold522 as FD  # noqa: E402  (on sys.path through score708)

DUMPS = ROOT / "experiments/801-lab-folds/work/dumps"
META = ROOT / "experiments/801-lab-folds/work/ka-labc/meta"
REG = json.loads((HERE / "REGISTRATION.json").read_text())
FOLDS = S7.FOLDS
NB, SEED = 4000, 713


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 23), b""):
            h.update(chunk)
    return h.hexdigest()


def agg(arr, d, idx):
    """(X, pds) from score708's per-target arrays on the target subset idx: X = frame mean of |sum x| / sqrt(sum a sum d)."""
    xs = [abs(arr["x"][f][idx].sum()) / np.sqrt(arr["a"][f][idx].sum() * d[idx].sum()) for f in range(arr["x"].shape[0])]
    return float(np.mean(xs)), float(arr["pds"][:, idx].mean())


def price_ag(arm, inc):
    return 0.0710 * (arm[0] / inc[0] - 1.0) + 0.1245 * (arm[1] - inc[1])


def main() -> None:
    t0 = time.time()
    for f in FOLDS:                                          # the registered inputs, byte for byte
        assert sha256(DUMPS / f"dump708-LABC-{f}.npz") == REG["dumps"][f], f"LABC/{f}: dump sha differs from the registration"
    assert sha256(ROOT / REG["scorer"]["path"]) == REG["scorer"]["sha256"], "score708.py differs from the registration"
    with open(ROOT / "data/challenge/pert_counts.csv") as fh:
        panel = [r["target_gene"] for r in csv.DictReader(fh) if r["target_gene"] != "non-targeting"]
    cfg = {"s_b0.125": ("s", 0.125), "s_b0.25": ("s", 0.25), "s_b0.5": ("s", 0.5), "r_t0.25": ("r", 0.25), "r_t0.5": ("r", 0.5)}
    assert set(cfg) == {c["name"] for c in REG["configs"]}, "configs differ from the registration"
    per, arrays, dvec = {}, {}, {}
    for f in FOLDS:
        S = S7.Fold(DUMPS / f"dump708-LABC-{f}.npz")
        wp = FD.wclip(S.kv, S.carried, 0.5, FD.pct_hi(S.kv, S.carried, FD.P_PROD))
        z = S7.zsrc_fold(S.targets, panel, S.sources)
        U = S7.unit_mean(S.targets, S.sources, S.nonpanel)
        have = np.linalg.norm(U, axis=1) > 0
        tier = np.full(S.n, -1)
        r_ = np.argsort(np.argsort(z[have], kind="stable"), kind="stable")
        tier[np.flatnonzero(have)] = np.minimum((r_ * 3) // max(int(have.sum()), 1), 2)
        inc_s, inc_a, _ = S.score(wp)
        m = json.loads((META / f"LABC__{f}.json").read_text())["arms"]["pct"]
        d0 = max(abs(inc_s["cos_agg"] - m["cos_agg"]), abs(inc_s["pds"] - m["pds"]))
        assert d0 < 1e-12, f"KA0 LABC/{f}: incumbent {inc_s['cos_agg']!r}/{inc_s['pds']!r} vs the lab's pct {m['cos_agg']!r}/{m['pds']!r}"
        dvec[f] = (S.truth * S.truth).sum(1)
        allidx = np.arange(S.n)
        per[f], arrays[f] = {}, {"inc": inc_a}
        for name, (kind, v) in cfg.items():
            if kind == "s":
                w, fin = wp, S7.swap_fin(S.fin, U, tier == 2, S.nonpanel, beta=v)
            else:
                w, fin = np.where(wp <= 0.5 + 1e-12, wp, np.maximum(0.5, wp * np.exp(v * z))), None
            sa, aa, _ = S.score(w, fin=fin)
            p = S7.price(sa, inc_s)
            p_ag = price_ag(agg(aa, dvec[f], allidx), agg(inc_a, dvec[f], allidx))
            assert abs(p - p_ag) < 1e-12, f"KA1 LABC/{f} {name}: summary price {p!r} vs array price {p_ag!r}"
            per[f][name] = p
            arrays[f][name] = aa
        print(f"LABC/{f}: KA0 |d| {d0:.1e}; " + " | ".join(f"{k} {v:+.5f}" for k, v in per[f].items())
              + f" ({time.time() - t0:.0f}s)", flush=True)
    print("KA0 PASS (incumbent == the lab's pct arm on 3 folds); KA1 PASS (array price == score708 price, 15 rows)", flush=True)
    # paired target bootstrap: the same resampled targets for every configuration and the incumbent, within each fold
    rng = np.random.default_rng(SEED)
    names = list(cfg)
    pairs = {"s_b0.125-s_b0.25": ("s_b0.125", "s_b0.25"), "s_b0.5-s_b0.25": ("s_b0.5", "s_b0.25"),
             "r_t0.5-r_t0.25": ("r_t0.5", "r_t0.25")}
    boot = {k: np.empty(NB) for k in names + list(pairs)}
    for b in range(NB):
        acc = {k: 0.0 for k in boot}
        for f in FOLDS:
            idx = rng.integers(0, len(dvec[f]), len(dvec[f]))
            inc = agg(arrays[f]["inc"], dvec[f], idx)
            pb = {k: price_ag(agg(arrays[f][k], dvec[f], idx), inc) for k in names}
            for k in names:
                acc[k] += pb[k] / len(FOLDS)
            for k, (a_, r_) in pairs.items():
                acc[k] += (pb[a_] - pb[r_]) / len(FOLDS)
        for k in boot:
            boot[k][b] = acc[k]
    out = {"registration": "experiments/713-labc-dose-partners/REGISTRATION.json", "units": "522 board units, k 1, ret 1",
           "per_fold": per, "configs": {}, "pairs": {}, "nb": NB, "seed": SEED}
    for k in names:
        v = np.array([per[f][k] for f in FOLDS])
        out["configs"][k] = dict(mean=float(v.mean()), n_pos=int((v > 0).sum()), ci95=[float(np.percentile(boot[k], 2.5)),
                                 float(np.percentile(boot[k], 97.5))])
    for k, (a_, r_) in pairs.items():
        v = np.array([per[f][a_] - per[f][r_] for f in FOLDS])
        out["pairs"][k] = dict(mean=float(v.mean()), per_fold={f: float(per[f][a_] - per[f][r_]) for f in FOLDS},
                               n_pos=int((v > 0).sum()), ci95=[float(np.percentile(boot[k], 2.5)),
                                                                float(np.percentile(boot[k], 97.5))],
                               sd=float(boot[k].std(ddof=1)))     # for the lab's best-of-N bound: lb - e(N) * sd
    for k, r in out["configs"].items():
        print(f"{k:10s} vs incumbent: mean {r['mean']:+.5f} ({r['n_pos']}/3 > 0) CI95 [{r['ci95'][0]:+.5f}, {r['ci95'][1]:+.5f}]",
              flush=True)
    for k, r in out["pairs"].items():
        print(f"{k:18s} paired: mean {r['mean']:+.5f} ({r['n_pos']}/3 > 0) CI95 [{r['ci95'][0]:+.5f}, {r['ci95'][1]:+.5f}] "
              f"sd {r['sd']:.5f} (522 units)", flush=True)
    out["seconds"] = round(time.time() - t0, 1)
    np.savez_compressed(HERE / "work" / "boot713.npz", **{k.replace("|", "_"): v for k, v in boot.items()})
    (HERE / "work" / "price713.json").write_text(json.dumps(out, indent=1))
    print("wrote work/price713.json")


if __name__ == "__main__":
    main()
