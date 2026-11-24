"""492 — the ballast's nmae cost, calibrated on the board's own q step. The emulated ballast count K (414's realise, KMAX 16,
the same route as abc473/abc492) for 433k's bulk (q 1.5, its factors) and 446b's (q 1.75), whose board step cost nmae
+0.001346 raw (BOARD-RESULT-446b: the deep ballast cells enter the per-cell arithmetic means with K/400 of the weight),
then each 492 arm's K. Cost(arm) ~ +0.001346 * (mean K(arm) - mean K(446b)) / (mean K(446b) - mean K(433k)).
Known answer: the emulated K medians reproduce the builds' metas (433k 4, 446b 6).

    uv run python experiments/492-per-target-sharpening/k492.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo  # noqa: E402

plolo.WORK = HERE / "work" / "plolo"
from plolo import MIN_CPM, profile_from_effect  # noqa: E402

S_C, LB, KMAX = 7_919_901.0, 995_000.0, 16
CTX = ("A", "B", "C")
F446 = {"A": 1.1004, "B": 1.1517, "C": 1.1720}
F433 = {"A": 1.0464, "B": 1.092, "C": 1.114}
SPECS = {"P1": (2.0, 1.75, 1.75), "P2": (2.5, 2.0, 1.75), "P3": (3.0, 2.25, 1.75), "G2": (2.0, 2.0, 2.0)}


def Kof(p_de, p_bulk):
    live = p_de > 1e-9
    r = float(np.min(p_bulk[live] / p_de[live]))
    return int(min(KMAX, max(1, int(np.ceil(S_C * (1.0 / max(r, 1e-6) - 1.0) / LB)))))


def main():
    z = np.load(ROOT / "experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz", allow_pickle=False)
    perts = [str(p) for p in z["perts"]]
    genes = [str(g) for g in z["genes"]]
    bulk = np.asarray(z["bulk"], np.float64)
    carriers = np.asarray(np.load(ROOT / "experiments/417-bulk-sharpen/work/dump_ka417a_g065.npz")["carriers"], np.float64)
    gpos = {g: i for i, g in enumerate(genes)}
    nonpanel = np.ones(len(genes), bool)
    nonpanel[[gpos[t] for t in perts]] = False
    zk = np.load(ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.3_hi3.npz")
    kpos = {str(t): i for i, t in enumerate(zk["targets"])}
    w = np.asarray([float(zk["w"][kpos[t]]) for t in perts])
    ter = np.digitize(w, np.quantile(w, [1 / 3, 2 / 3]))
    QV = {k: np.asarray([v[t] for t in ter], float) for k, v in SPECS.items()}
    QV["CTRL P2 permuted"] = QV["P2"][np.random.default_rng(492).permutation(len(perts))]
    QV["433k"] = np.full(len(perts), 1.5)
    QV["446b"] = np.full(len(perts), 1.75)
    frames = {}
    for c in CTX:
        b = plolo.challenge_bctrl(c)
        cpm = 1e6 * np.expm1(b) / np.expm1(b).sum()
        frames[c] = (b, (cpm > MIN_CPM) & (b > 0))
    P_DE = {c: [profile_from_effect(carriers[i], frames[c][0], frames[c][1], "minnorm") for i in range(len(perts))] for c in CTX}
    out = {}
    for name, qv in QV.items():
        fac = F433 if name == "433k" else F446
        Ks = {}
        for c in CTX:
            b, gate = frames[c]
            G = gate & nonpanel
            F = bulk.copy()
            H = bulk[:, G]
            S = np.sign(H) * np.abs(H) ** qv[:, None]
            S *= (np.linalg.norm(H, axis=1) / np.maximum(np.linalg.norm(S, axis=1), 1e-300))[:, None]
            F[:, G] = S
            F *= fac[c]
            Ks[c] = np.array([Kof(P_DE[c][i], profile_from_effect(F[i], b, gate, "minnorm")) for i in range(len(perts))])
        allK = np.concatenate(list(Ks.values()))
        out[name] = dict(K_mean=float(allK.mean()), K_median=float(np.median(allK)),
                         by_tercile=[float(allK[np.tile(ter, 3) == t].mean()) for t in range(3)])
        print(f"{name:18s} K mean {allK.mean():.3f} median {np.median(allK):.0f}  by key tercile (low/mid/high) "
              f"{[round(x, 2) for x in out[name]['by_tercile']]}", flush=True)
    dk = out["446b"]["K_mean"] - out["433k"]["K_mean"]
    for name in QV:
        out[name]["nmae_cost_raw"] = 0.001346 * (out[name]["K_mean"] - out["446b"]["K_mean"]) / dk
        out[name]["overall_cost"] = -out[name]["nmae_cost_raw"] * 0.2808 / 6 * 6 / 6   # scaled nmae per raw ~ 1/span 0.6; see below
    # overall per raw nmae: 446b's board record: +0.001346 raw nmae = -0.000378 overall -> -0.2808 overall per unit raw
    for name in QV:
        out[name]["overall_cost"] = -0.2808 * out[name]["nmae_cost_raw"]
        print(f"{name:18s} calibrated ballast nmae cost {out[name]['nmae_cost_raw']:+.5f} raw -> overall {out[name]['overall_cost']:+.5f}")
    ka = (out["433k"]["K_median"], out["446b"]["K_median"])
    print(f"KA: emulated K medians 433k {ka[0]:.0f} (meta 4), 446b {ka[1]:.0f} (meta 6)")
    out["ka_medians"] = ka
    (HERE / "work" / "k492.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
