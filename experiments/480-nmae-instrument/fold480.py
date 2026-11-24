"""480 step 1 (resume plan): do the panel-donor constructions rank on cancer folds as they do on H1?

Surfaces (truth = the scorer's reference DE rules on a cell harness; relative comparisons only on the folds):
  h1            h1-test (323's reference; donors = K562-GWPS, Orion HCT116, Orion HEK293T)
  orion_hct116  Orion HCT116 thinned harness (gates480.py); donors = K562-GWPS + Orion HEK293T (same-study sibling)
  orion_hek293t Orion HEK293T thinned harness; donors = K562-GWPS + Orion HCT116 (same-study sibling)
  k562_gwps     K562-GWPS panel cells (k562harness480.py -> gates480.py); donors = the two Orion lines (cross-study)
The held-out line's own donor table is never used on its fold.

Constructions (identical code on every surface; sign source = the available donors' consensus sign unless noted):
  zero          all-zero prediction (known answer: nmae exactly 1)
  mde_cons      sign(cons) x shipped MDE law min(0.1/sqrt(cnt), 0.1) for cnt >= 0.3 (information-light, "ours-like")
  cons@1        donor consensus (mean of available donors' log2 ratios, pc 1 CPM) at scale 1 (pack-style)
  <donor>@1     each single donor at scale 1 (pc 1 CPM; naive20 = delta/ln2 also for K562 on H1-like reading)
  cons@k        consensus at its in-sample L1-optimal global k (upper bound for a global scale)
  cons@kLOSO    consensus at the k fitted on the other surfaces (mean of their in-sample k)
  mde+cons@k    mde_cons + k * cons, k in-sample (does donor magnitude add to the MDE law?)
  signorc_mde   oracle sign x MDE law (what the law is worth with perfect signs)
  cons_sign_orc sign(cons) x oracle |l| (the ceiling for this sign source)
On h1 also: a_coarse (ours as shipped), a_base, ours_sign x MDE (no base).
Writes work/fold480.json and work/fold480-rows.npz (per-target rows for every surface and construction).
"""
from __future__ import annotations

import json
import sys

import numpy as np

from common480 import HERE, ROOT, boot_delta, gate_from, load_donor_rows, load_h1, nmae, nmae_rows

DON = ["k562_gwps", "orion_hct116", "orion_hek293t"]
SURF = ["h1", "h1_c150", "h1_thin", "orion_hct116", "orion_hek293t", "k562_gwps"]
LINE = {"h1": "h1", "h1_c150": "h1", "h1_thin": "h1", "orion_hct116": "hct116", "orion_hek293t": "hek293t",
        "k562_gwps": "k562"}


def mde_add(cnt, c=0.1, cap=0.1, mn=0.3):
    return np.where(cnt >= mn, np.minimum(c / np.sqrt(np.maximum(cnt, 1e-3)), cap), 0.0)


def load_surface(name):
    """Returns dict(targets, genes, T, cnt, donors, extra) for one surface."""
    if name == "h1":
        h = load_h1()
        return dict(targets=h["targets"], genes=h["genes"], T=h["T"], cnt=h["cnt"], donors=DON)
    z = np.load(HERE / f"work/gate-{name}.npz", allow_pickle=True)
    targets, genes = z["targets"].astype(str), z["genes"].astype(str)
    T = gate_from(z["q"].astype(float), z["lfc"].astype(float), genes, targets)
    cnt = z["ctrl_cpm"].astype(float) * float(z["lib_med"]) / 1e6
    return dict(targets=targets, genes=genes, T=T, cnt=cnt,
                donors=DON if name.startswith("h1") else [d for d in DON if d != name],
                ncell=z["ncell"], lib_med=float(z["lib_med"]))


def ours_on(targets, genes):
    """Our delivered H1 lfc (426 base, 435 coarse) mapped by name onto another H1 gene axis (0 where absent)."""
    h = load_h1()
    assert (h["targets"] == targets).all()
    col = {g: i for i, g in enumerate(h["genes"])}
    idx = np.array([col.get(g, -1) for g in genes])
    out = []
    for f in ["experiments/426-mde-magnitude/work/h1/stats-base_s0_bps.npz",
              "experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz"]:
        L = np.nan_to_num(np.asarray(np.load(ROOT / f)["lfc"], float))
        M = np.zeros((len(targets), len(genes)))
        M[:, idx >= 0] = L[:, idx[idx >= 0]]
        out.append(M)
    return out


def best_k(P, T, lo=0.0, hi=2.0):
    ks = np.linspace(lo, hi, 81)
    v = [nmae(k * P, T) for k in ks]
    k0 = ks[int(np.argmin(v))]
    ks2 = np.linspace(max(lo, k0 - 0.05), k0 + 0.05, 41)
    v2 = [nmae(k * P, T) for k in ks2]
    return float(ks2[int(np.argmin(v2))])


def sign_acc(S, T):
    """Gate sign accuracy (zero prediction scores 1/2): pooled unweighted and |l|-weighted."""
    G = np.isfinite(T)
    t = T[G]
    s = S[G]
    hit = np.where(s == 0, 0.5, (np.sign(s) == np.sign(t)).astype(float))
    return float(hit.mean()), float((hit * np.abs(t)).sum() / np.abs(t).sum())


K_H1_MDECONS = 0.575   # H1's in-sample k for mde+cons (first run of this script, grid step 0.025)
BETA = None            # set in main() from H1: median |L_base| / median |cons| over H1 gate cells


def main():
    global BETA
    h = load_h1()
    G = np.isfinite(h["T"])
    L0 = np.nan_to_num(np.asarray(np.load(ROOT / "experiments/426-mde-magnitude/work/h1/stats-base_s0_bps.npz")["lfc"],
                                  float))
    c3 = np.nan_to_num(np.nanmean(np.stack([load_donor_rows(d, h["targets"], h["genes"], pc_cpm=1.0)[0]
                                            for d in DON]), 0))
    BETA = float(np.median(np.abs(L0[G])) / np.median(np.abs(c3[G])))
    print(f"BETA (median |our base| / median |cons| on H1 gate cells) = {BETA:.4f}")
    res, rows_all, ks_in = {}, {}, {}
    cache = {}
    for name in SURF:
        S = load_surface(name)
        T, cnt, genes, targets = S["T"], S["cnt"], S["genes"], S["targets"]
        scored = np.isfinite(T).any(1)
        R = {}
        for d in S["donors"]:
            L2, has, okc, _ = load_donor_rows(d, targets, genes, pc_cpm=1.0)
            R[d] = L2
        cons = np.nanmean(np.stack([R[d] for d in S["donors"]]), 0)
        cons0 = np.nan_to_num(cons)
        s = np.sign(cons0)
        add = mde_add(cnt)[None, :]
        P = {"zero": np.zeros_like(T), "mde_cons": s * add, "cons@1": cons0}
        for d in S["donors"]:
            P[f"{d}@1"] = np.nan_to_num(R[d])
        k = best_k(cons0, T)
        ks_in[name] = k
        P["cons@k"] = k * cons0
        # mde + k * cons, k in-sample
        grid = np.linspace(0, 2.5, 101)
        v = [nmae(s * add + kk * cons0, T) for kk in grid]
        k2 = float(grid[int(np.argmin(v))])
        P["mde+cons@k"] = s * add + k2 * cons0
        # a-analog: the H1-calibrated construction (MDE law + consensus at H1's in-sample k), applied unchanged
        P["a_analog"] = s * add + K_H1_MDECONS * cons0
        # small-base MDE calibration: base = BETA * cons (BETA matched to our delivered base's magnitude on H1)
        P["small_base"] = BETA * cons0
        P["small_base+mde"] = BETA * cons0 + s * add
        Tz = np.nan_to_num(T)
        P["signorc_mde"] = np.sign(Tz) * add
        P["cons_sign_orc"] = s * np.abs(Tz)
        if name.startswith("h1"):
            L0, Lc = ours_on(targets, genes)
            P["a_coarse"], P["a_base"] = Lc, L0
            P["ours_sign_mde"] = np.sign(Lc) * add
            P["ours_sign_orc"] = np.sign(Lc) * np.abs(Tz)
        rows = {kk: nmae_rows(v_, T) for kk, v_ in P.items()}
        z0 = float(np.nanmean(rows["zero"]))
        if abs(z0 - 1.0) > 1e-12:
            print(f"KNOWN ANSWER FAILED: zero prediction nmae {z0} on {name}")
            sys.exit(1)
        if name == "h1" and abs(float(np.nanmean(rows["a_base"])) - 0.95237) > 5e-5:
            print("KNOWN ANSWER FAILED: h1 a_base")
            sys.exit(1)
        G = np.isfinite(T)
        info = dict(n_targets=len(targets), scored=int(scored.sum()), gate_cells=int(G.sum()),
                    gate_median=float(np.median(G.sum(1)[scored])), donors=S["donors"], k_cons=k, k_mde_cons=k2,
                    cons_cov=float(np.isfinite(cons)[G].mean()),
                    sign_acc_cons=sign_acc(cons0, T), lib_med=S.get("lib_med"),
                    ncell_med=float(np.median(S["ncell"])) if "ncell" in S else None,
                    nmae={kk: float(np.nanmean(v_)) for kk, v_ in rows.items()})
        if name.startswith("h1"):
            info["sign_acc_ours"] = sign_acc(P["a_coarse"], T)
        res[name] = info
        rows_all[name] = rows
        cache[name] = (cons0, T)
        print(f"\n== {name}: {info['scored']} scored targets, gate median {info['gate_median']:.0f}, donors "
              f"{S['donors']}, cons coverage {info['cons_cov']:.3f}, cons sign acc {info['sign_acc_cons']}")
        for kk, v_ in info["nmae"].items():
            print(f"  {kk:18s} {v_:.4f}")
        print(f"  k_cons {k:.3f}  k(mde+cons) {k2:.3f}")
    # LOSO k: consensus scale fitted on the other surfaces (mean of their in-sample k)
    for name in SURF:
        kl = float(np.mean([ks_in[o] for o in ["h1", "orion_hct116", "orion_hek293t", "k562_gwps"]
                            if LINE[o] != LINE[name]]))
        cons0, T = cache[name]
        rows_all[name]["cons@kLOSO"] = nmae_rows(kl * cons0, T)
        res[name]["k_loso"] = kl
        res[name]["nmae"]["cons@kLOSO"] = float(np.nanmean(rows_all[name]["cons@kLOSO"]))
    # bootstrap deltas vs mde_cons (target bootstrap, 4,000 draws; unit = scored target of that surface)
    for name in SURF:
        b = {}
        for kk, v_ in rows_all[name].items():
            if kk in ("mde_cons", "zero"):
                continue
            b[kk] = boot_delta(rows_all[name]["mde_cons"], v_)
        res[name]["boot_vs_mde_cons"] = b
    for name in SURF:
        r = rows_all[name]
        res[name]["a_analog_vs_cons@1"] = boot_delta(r["cons@1"], r["a_analog"])
        res[name]["smallbase_mde_delta"] = boot_delta(r["small_base"], r["small_base+mde"])
        if name.startswith("h1"):
            res[name]["a_coarse_vs_cons@1"] = boot_delta(r["cons@1"], r["a_coarse"])
            res[name]["real_mde_delta"] = boot_delta(r["a_base"], r["a_coarse"])
    res["_beta"] = BETA
    json.dump(res, open(HERE / "work/fold480.json", "w"), indent=1)
    np.savez_compressed(HERE / "work/fold480-rows.npz",
                        **{f"{n}__{k}": v for n, r in rows_all.items() for k, v in r.items()})
    print("\nΔ vs mde_cons [target bootstrap 95 %]")
    keys = ["cons@1", "cons@k", "cons@kLOSO", "mde+cons@k", "a_analog", "a_coarse"]
    for name in SURF:
        line = []
        for kk in keys:
            if kk in res[name]["boot_vs_mde_cons"]:
                d, lo, hi, _ = res[name]["boot_vs_mde_cons"][kk]
                line.append(f"{kk} {d:+.4f} [{lo:+.4f},{hi:+.4f}]")
        print(f"{name:14s} " + " | ".join(line))
    print("\na_analog - cons@1 (negative = H1-calibrated small-magnitude construction beats pack-style scale 1)")
    for name in SURF:
        d, lo, hi, _ = res[name]["a_analog_vs_cons@1"]
        e, elo, ehi, _ = res[name]["smallbase_mde_delta"]
        print(f"{name:14s} {d:+.4f} [{lo:+.4f},{hi:+.4f}]   | MDE add on small base {e:+.4f} [{elo:+.4f},{ehi:+.4f}]")
    for name in [n for n in SURF if n.startswith("h1")]:
        d, lo, hi, _ = res[name]["a_coarse_vs_cons@1"]
        e, elo, ehi, _ = res[name]["real_mde_delta"]
        print(f"{name} real (a) {d:+.4f} [{lo:+.4f},{hi:+.4f}]   | real MDE add {e:+.4f} [{elo:+.4f},{ehi:+.4f}]")


if __name__ == "__main__":
    main()
