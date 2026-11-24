"""492g — the ballast's nmae cost, re-examined on the board's five bulk-only steps before pricing 479's arms on 446bj2.

k492.py priced the ballast's nmae cost from ONE board step (433k -> 446b: nmae +0.001346 raw per +1.37 emulated mean K).
On the LODO base (446bj2) the per-target arms raise the LOW tercile's ballast far more than on the K3 base (floor 0.5 gives
those rows 1.67x the amplitude), so that single-point model now decides a third of the prize. Five board steps changed
only the bulk (the DE design fixed within each era), which lets the model be tested rather than assumed:

    410b -> 417a  (q 1 -> 1.25)                         nmae raw +0.000094   (before the coarse MDE / zero repair)
    417a -> 417e  (the H1/CD4 stack by LOAD)            nmae raw -0.000147
    433h -> 433k  (q 1.25 -> 1.5, re-levelled)          nmae raw +0.000814   (MDE era: same DE design from 433h on)
    433k -> 446b  (q 1.5 -> 1.75, re-levelled)          nmae raw +0.001346
    446b -> 446bj2 (LODO key, re-levelled)              nmae raw +0.000060

Two cost models per configuration, from 414's ballast emulation (KMAX 16) against the carriers' design profiles:
  K model  : cost ~ emulated mean ballast count K (k492).
  lfc model: the ballast's displacement of the DE channel. Per block the per-cell arithmetic mean CPM is
             m = ((400 - K) p_de + K q_b) / 400, q_b = the ballast cells' composition (realise's tot - S_C p_de, normalised);
             D = mean over blocks of mean over the block's top-100 design genes (|log2 p_de/p_ctrl|, gated, non-panel) of
             |lfc_X - lfc_Y| (lfc with the scorer's epsilon, floored at the zero repair's -5), and the signed change of
             |lfc| and of the distance from the carriers' own design lfc.
Known answer: the emulated K medians vs each build's meta (rc_callset/ballast_K median).

    uv run python experiments/492-per-target-sharpening/ballast492g.py        (~1.5 GB, ~1 min)
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
from plolo import MIN_CPM, profile_from_effect  # noqa: E402

W417 = ROOT / "experiments/417-bulk-sharpen/work"
S_C, LB, KMAX, NCELL, TOPN = 7_919_901.0, 995_000.0, 16, 400, 100
EPS, FLOOR = 1e-15, -5.0     # the scorer's epsilon (1e-9 CPM, in fractions) and the builder's zero-repair floor (log2 -5)


def lfc(m, pc):
    """the scorer's log2 fold change of arithmetic mean CPM, floored where the zero repair floors the carriers."""
    return np.maximum(np.log2((m + EPS) / (pc + EPS)), FLOOR)
CTX = ("A", "B", "C")
ONE = {"A": 1.0, "B": 1.0, "C": 1.0}
F433H = {"A": 1.0088, "B": 1.05, "C": 1.0736}
F433 = {"A": 1.0464, "B": 1.092, "C": 1.114}
F446 = {"A": 1.1004, "B": 1.1517, "C": 1.1720}
FBJ2 = {"A": 1.0962, "B": 1.1461, "C": 1.1681}
# (bulk source, q, factors, meta ballast-K median)
CONFIGS = {
    "410b": ("pp", None, ONE, 2), "417a": ("pp", 1.25, ONE, 2), "417e": ("k3", 1.25, ONE, 2),
    "433h": ("k3", 1.25, F433H, 3), "433k": ("k3", 1.5, F433, 4), "446b": ("k3", 1.75, F446, 6),
    "446bj2": ("lodo", 1.75, FBJ2, 7),
}
STEPS = {("410b", "417a"): 0.936552 - 0.936458, ("417a", "417e"): 0.936405 - 0.936552,
         ("433h", "433k"): 0.895685 - 0.894871, ("433k", "446b"): 0.897031 - 0.895685,
         ("446b", "446bj2"): 0.897091 - 0.897031}


def realise(p_de, p_bulk):
    """ballast414.realise (KMAX 16), returning the block totals too."""
    live = p_de > 1e-9
    r = float(np.min(p_bulk[live] / p_de[live]))
    K = int(min(KMAX, max(1, int(np.ceil(S_C * (1.0 / max(r, 1e-6) - 1.0) / LB)))))
    S = S_C + K * LB
    want = S * p_bulk
    floor = S_C * p_de
    tot = np.where(want < floor, floor, want)
    return tot, K


def tercile(v, spec, floor=None):
    t = np.digitize(v, np.quantile(v, [1 / 3, 2 / 3]))
    if floor is not None:
        t = np.where(v <= floor + 1e-12, 0, t)
    return np.asarray([spec[k] for k in t], float), t


def main():
    t0 = time.time()
    log = lambda *a: print("[ballast492g]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    zl = np.load(W417 / "dump_h1x2_cd4x15_lodo05.npz", allow_pickle=False)
    z3 = np.load(W417 / "dump_h1x2_cd4x15.npz", allow_pickle=False)
    zp = np.load(W417 / "dump_ka417a_g065.npz", allow_pickle=False)
    perts = [str(p) for p in zl["perts"]]
    genes = [str(g) for g in zl["genes"]]
    for z in (z3, zp):
        assert [str(p) for p in z["perts"]] == perts and [str(g) for g in z["genes"]] == genes
    bulks = {"lodo": np.asarray(zl["bulk"], np.float64), "k3": np.asarray(z3["bulk"], np.float64),
             "pp": np.asarray(zp["bulk"], np.float64)}
    carriers = np.asarray(zp["carriers"], np.float64)
    stack_car = np.asarray(zl["carriers"], np.float64)
    gpos = {g: i for i, g in enumerate(genes)}
    nonpanel = np.ones(len(genes), bool)
    nonpanel[[gpos[t] for t in perts]] = False
    n = len(perts)
    frames, P_DE, TOP, LFC_DE = {}, {}, {}, {}
    for c in CTX:
        b = plolo.challenge_bctrl(c)
        pc = np.expm1(b) / np.expm1(b).sum()
        gate = (1e6 * pc > MIN_CPM) & (b > 0)
        frames[c] = (b, gate, pc)
        P_DE[c] = [profile_from_effect(carriers[i], b, gate, "minnorm") for i in range(n)]
        G = np.flatnonzero(gate & nonpanel)
        TOP[c] = np.empty((n, TOPN), int)
        for i in range(n):
            l = np.abs(lfc(P_DE[c][i][G], pc[G]))
            TOP[c][i] = G[np.argsort(-l)[:TOPN]]
        LFC_DE[c] = np.stack([lfc(P_DE[c][i][TOP[c][i]], pc[TOP[c][i]]) for i in range(n)])
    log("frames, carriers' design profiles and top-100 design genes ready")

    def sharpen(B, gate, qv):
        if qv is None:
            return B
        G = gate & nonpanel
        out = B.copy()
        H = B[:, G]
        S = np.sign(H) * np.abs(H) ** np.asarray(qv, float)[:, None]
        out[:, G] = S * (np.linalg.norm(H, axis=1) / np.maximum(np.linalg.norm(S, axis=1), 1e-300))[:, None]
        return out

    def run(B, qv, fac):
        """K per block and the DE channel's per-cell mean lfc on each block's top design genes."""
        K = {}
        L = {}
        for c in CTX:
            b, gate, pc = frames[c]
            F = sharpen(B, gate, qv) * fac[c]
            K[c] = np.empty(n, int)
            L[c] = np.empty((n, TOPN))
            for i in range(n):
                p = profile_from_effect(F[i], b, gate, "minnorm")
                tot, k = realise(P_DE[c][i], p)
                Q = tot - S_C * P_DE[c][i]
                qb = Q / max(Q.sum(), 1e-300)
                m = ((NCELL - k) * P_DE[c][i] + k * qb) / NCELL
                t = TOP[c][i]
                K[c][i] = k
                L[c][i] = lfc(m[t], pc[t])
        return K, L

    out = dict(configs={}, steps={}, arms={})
    R = {}
    for name, (src, q, fac, kmeta) in CONFIGS.items():
        qv = None if q is None else np.full(n, q)
        K, L = run(bulks[src], qv, fac)
        Kall = np.concatenate([K[c] for c in CTX])
        R[name] = (K, L)
        out["configs"][name] = dict(K_mean=float(Kall.mean()), K_median=float(np.median(Kall)), K_meta_median=kmeta,
                                    lam_mean=float(((Kall + S_C / LB) / NCELL).mean()))
        log(f"{name:7s} K mean {Kall.mean():.3f} median {np.median(Kall):.0f} (meta {kmeta})")

    def dist(La, Lb):
        d = np.concatenate([(La[c] - Lb[c]).ravel() for c in CTX])
        sh = np.concatenate([(np.abs(La[c]) - np.abs(Lb[c])).ravel() for c in CTX])
        # toward/away from the carriers' own design lfc (the ballast pulls the per-cell means toward the bulk)
        aw = np.concatenate([(np.abs(La[c] - LFC_DE[c]) - np.abs(Lb[c] - LFC_DE[c])).ravel() for c in CTX])
        return dict(D_abs=float(np.mean(np.abs(d))), D_shrink=float(np.mean(sh)), D_from_design=float(np.mean(aw)))

    for (a, b_), dn in STEPS.items():
        Ka = np.concatenate([R[b_][0][c] for c in CTX]).mean() - np.concatenate([R[a][0][c] for c in CTX]).mean()
        dd = dist(R[b_][1], R[a][1])
        out["steps"][f"{a}->{b_}"] = dict(board_dnmae_raw=dn, dK_mean=float(Ka), **dd)
        log(f"step {a:6s}->{b_:7s}: board dnmae {dn:+.6f} | dK {Ka:+.3f} | D_abs {dd['D_abs']:.5f} D_shrink {dd['D_shrink']:+.5f} "
            f"D_from_design {dd['D_from_design']:+.5f}")

    # arms on the 446bj2 base (the LODO bulk at 446bj2's factors)
    wl = np.asarray(np.load(ROOT / "experiments/420-bulk-stack/work/pert_scale_415pp_lo0.5_hi3.npz")["w"], float)
    strength = np.linalg.norm(stack_car[:, nonpanel], axis=1)
    P2, L30 = (2.5, 2.0, 1.75), (3.0, 1.75, 1.75)
    qp2g, _ = tercile(wl, P2, 0.5)
    arms = {"P2g": qp2g, "S_P2": tercile(strength, P2)[0], "L30g": tercile(wl, L30, 0.5)[0], "G2": np.full(n, 2.0),
            "CTRL P2g permuted": qp2g[np.random.default_rng(492).permutation(n)]}
    Kb, Lb0 = R["446bj2"]
    for name, qv in arms.items():
        K, L = run(bulks["lodo"], qv, FBJ2)
        dK = np.concatenate([K[c] for c in CTX]).mean() - np.concatenate([Kb[c] for c in CTX]).mean()
        dd = dist(L, Lb0)
        out["arms"][name] = dict(dK_mean=float(dK), **dd)
        log(f"arm {name:18s}: dK {dK:+.3f} | D_abs {dd['D_abs']:.5f} D_shrink {dd['D_shrink']:+.5f} D_from_design {dd['D_from_design']:+.5f}")
    # the same arms on 446b's K3 base (492's original pricing base), for the comparison of the two bases' ballast cost
    zk3 = np.load(ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.3_hi3.npz")
    kp = {str(t): i for i, t in enumerate(zk3["targets"])}
    w3 = np.asarray([float(zk3["w"][kp[t]]) for t in perts])
    K, L = run(bulks["lodo"], tercile(w3 := np.asarray([float(np.load(ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.3_hi3.npz")["w"][i]) for i in [{str(t): j for j, t in enumerate(np.load(ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.3_hi3.npz")["targets"])}[t] for t in perts]]), P2, 0.3)[0], FBJ2)
    dK = np.concatenate([K[c] for c in CTX]).mean() - np.concatenate([Kb[c] for c in CTX]).mean()
    out["arms"]["P2g_K3key"] = dict(dK_mean=float(dK), **dist(L, Lb0))
    log(f"arm P2g_K3key         : dK {dK:+.3f} | D_from_design {out['arms']['P2g_K3key']['D_from_design']:+.5f}")
    Kb3, Lb3 = R["446b"]
    for name, qv in {"K3 base: P2g": tercile(w3, P2, 0.3)[0], "K3 base: S_P2": tercile(strength, P2)[0]}.items():
        K, L = run(bulks["k3"], qv, F446)
        dK = np.concatenate([K[c] for c in CTX]).mean() - np.concatenate([Kb3[c] for c in CTX]).mean()
        dd = dist(L, Lb3)
        out["arms"][name] = dict(dK_mean=float(dK), **dd)
        log(f"arm {name:18s}: dK {dK:+.3f} | D_abs {dd['D_abs']:.5f} D_shrink {dd['D_shrink']:+.5f} D_from_design {dd['D_from_design']:+.5f}")
    # the MDE-era fits (433h -> 433k, 433k -> 446b, 446b -> 446bj2 share the DE design): slope through the origin
    era = ("433h->433k", "433k->446b", "446b->446bj2")
    for key in ("dK_mean", "D_from_design", "D_abs"):
        x = np.array([out["steps"][s_][key] for s_ in era])
        y = np.array([out["steps"][s_]["board_dnmae_raw"] for s_ in era])
        sl = float((x * y).sum() / (x * x).sum())
        out.setdefault("fits", {})[key] = dict(slope=sl, resid=[float(v) for v in (y - sl * x)])
        log(f"fit nmae ~ {key}: slope {sl:.5f}; residuals {np.round(y - sl * x, 6).tolist()}")
    for name, a in out["arms"].items():
        a["nmae_raw_K"] = out["fits"]["dK_mean"]["slope"] * a["dK_mean"]
        a["nmae_raw_D"] = out["fits"]["D_from_design"]["slope"] * a["D_from_design"]
        log(f"{name:22s} ballast nmae raw: K model {a['nmae_raw_K']:+.5f}, D model {a['nmae_raw_D']:+.5f} "
            f"-> overall {-0.2808 * a['nmae_raw_K']:+.5f} / {-0.2808 * a['nmae_raw_D']:+.5f}")
    (WORK / "ballast492g.json").write_text(json.dumps(out, indent=1, default=float))
    log("done")


if __name__ == "__main__":
    main()
