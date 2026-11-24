"""Known-answer reproduction for 358 via an independent pipeline.

Reproduces 341's pooled means (source mean / cross-line diagonal) for
weighted sign, own cosine and retrieval on the six-line frame, using only
reimplemented primitives (no import of core334/core341/run334).
"""
import hashlib
import json
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[2]
BANKS = ("data/donors/k562_gwps", "data/donors/orion_hct116",
         "data/donors/orion_hek293t", "data/donors/rpe1",
         "data/donors/nadig_jurkat",
         "experiments/166-donor-acquisition/donors/nadig_hepg2")
LINES = ("k562_gwps", "orion_hct116", "orion_hek293t", "rpe1",
         "nadig_jurkat", "nadig_hepg2")


def order(names, salt):
    return sorted(range(len(names)),
                  key=lambda i: hashlib.sha256(
                      (salt + str(names[i])).encode()).hexdigest())


def load(bank, targets, genes):
    p = ROOT / bank
    pi = {t: i for i, t in enumerate(p.joinpath("perts.txt").read_text().splitlines())}
    gi = {g: i for i, g in enumerate(p.joinpath("genes.txt").read_text().splitlines())}
    ti, ix = [pi[t] for t in targets], [gi[g] for g in genes]
    raw = np.load(p / "delta.npy", mmap_mode="r")
    out = np.array(raw[np.ix_(ti, ix)], dtype=np.float64)
    basal = np.array(np.load(p / "b_ctrl.npy")[ix], dtype=np.float64)
    assert np.isfinite(out).all() and np.isfinite(basal).all()
    return out, basal


def rms_direction(a):
    rms = np.sqrt(np.mean(a * a, axis=-1))
    return a / np.maximum(rms[..., None], 1e-12), rms


def unit(a):
    return a / np.maximum(np.linalg.norm(a, axis=-1, keepdims=True), 1e-12)


def endpoints(prediction, truth):
    p, t = unit(prediction), unit(truth)
    s = p @ t.T
    own = np.diag(s).copy()
    beaten = (s < own[:, None]).sum(1)
    ties = (s == own[:, None]).sum(1) - 1
    n = len(truth) - 1
    return {"cosine": own, "retrieval": (beaten + 0.5 * ties) / n}


def sign_endpoints(prediction, truth):
    weight = np.abs(truth)
    correct = np.sign(prediction) == np.sign(truth)
    full = (weight * correct).sum(1) / weight.sum(1)
    top = np.argsort(-weight, axis=1, kind="stable")[:, :64]
    w = np.take_along_axis(weight, top, axis=1)
    c = np.take_along_axis(correct, top, axis=1)
    return {"weighted_sign": full, "top64": (w * c).sum(1) / w.sum(1)}


def main():
    sel = json.loads((ROOT / "experiments/334-text-conditioned-transfer"
                      / "selection.json").read_text())
    targets, genes = sel["targets"], sel["genes"]
    assert len(targets) == 1024 and len(genes) == 512 and len(sel["lines"]) == 6
    assert sel["lines"] == list(LINES), sel["lines"]
    loaded = [load(b, targets, genes) for b in BANKS]
    raw = np.array([a for a, _ in loaded])
    direction, _ = rms_direction(raw)
    del loaded, raw
    perm = order(targets, "334-split|")
    blocks = [np.array(perm[::2]), np.array(perm[1::2])]
    agg = {}
    for outer in range(6):
        donors = [i for i in range(6) if i != outer]
        baseline = direction[donors].mean(0)
        truth = direction[outer]
        for arm in ("source_mean", "crossline_diag"):
            agg.setdefault(arm, {}).setdefault("w", []).append(None)
        out = {a: np.full_like(baseline, np.nan) for a in agg}
        for block in (0, 1):
            train, test = blocks[1 - block], blocks[block]
            num = np.zeros(512)
            den = np.zeros(512)
            for d in donors:
                others = [i for i in donors if i != d]
                x = direction[others][:, train].mean(0)
                y = direction[d][train]
                num += (x * y).sum(0)
                den += (x * x).sum(0)
            prior = float(np.median(den))
            gain = (num + prior) / np.maximum(den + prior, 1e-12)
            out["source_mean"][test] = baseline[test]
            out["crossline_diag"][test] = baseline[test] * gain
        for arm, pred in out.items():
            m = {**endpoints(pred, truth), **sign_endpoints(pred, truth)}
            for k, v in m.items():
                agg[arm].setdefault(k, []).append(v)
    print(f"{'arm':16s} {'w_sign':>10s} {'top64':>10s} {'cos':>10s} {'retr':>10s}")
    means = {}
    for arm in ("source_mean", "crossline_diag"):
        row = {}
        for k in ("weighted_sign", "top64", "cosine", "retrieval"):
            row[k] = float(np.array(agg[arm][k]).mean())
        means[arm] = row
        print(f"{arm:16s} {row['weighted_sign']:10.6f} {row['top64']:10.6f} "
              f"{row['cosine']:10.6f} {row['retrieval']:10.6f}")
    rec = {"source_mean": {"weighted_sign": 0.589398, "cosine": 0.195532,
                           "retrieval": 0.790095, "top64": 0.642517},
           "crossline_diag": {"weighted_sign": 0.589398, "cosine": 0.198735,
                              "retrieval": 0.794614, "top64": 0.642517}}
    ok = True
    for arm in rec:
        for k, v in rec[arm].items():
            d = abs(means[arm][k] - v)
            flag = "MATCH" if d < 5e-6 else "MISS"
            if d >= 5e-6:
                ok = False
            print(f"KA {arm}/{k}: repro {means[arm][k]:.6f} recorded {v:.6f} "
                  f"diff {d:.2e} {flag}")
    print("KA_OVERALL:", "PASS" if ok else "FAIL")


if __name__ == "__main__":
    with threadpool_limits(limits=4):
        main()
