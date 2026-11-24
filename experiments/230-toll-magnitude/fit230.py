"""230 Phase 1c/2 — fit the board `reach` delta against PROFILE-CHANGE MAGNITUDE.

Rule 13 reads the per-gene panel-total digest as a boolean and charges a class constant.  This
file puts a SIZE on that change for every scored arm whose emission law is known, lines the
sizes up against the arm's own realised board `reach` delta, and tests STEP against LINEAR
against SATURATING -- including the out-of-sample retrodiction that a law has to survive.

Magnitudes come from three routes, in decreasing order of directness:

  measured on the artifacts     `cd4f034` vs `rc_callset`, from 214's cached block totals
  measured through the emitter  a8 / a21 / a22 -- emit230.py, real cells, 0 genes differing
  closed form on the renderer   a23's mass fix and 208's Hamilton rounding -- magnitude230.py
  closed form on the effects    `h1union` vs `cd4f034`, from 222's own effect matrices, with
                                `cd4f034`'s cached totals as the known-answer row

    uv run python fit230.py
"""
from __future__ import annotations

import hashlib
import itertools
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
ROOT = HERE.parent.parent
sys.path.insert(0, str(ROOT / "experiments/196-emitted-depth"))
import recon196                                                        # noqa: E402
from magnitude230 import magnitudes                                    # noqa: E402

TS, N, MIN_CPM = 50_000.0, 400, 10.0
CTX = ("A", "B", "C")


def sha16(a) -> str:
    return hashlib.sha256(np.asarray(a, dtype=np.int64).tobytes()).hexdigest()[:16]


def profile_from_effect(effect, b_ctrl, min_cpm=MIN_CPM):
    """`make_submission_pm_cd4.profile_from_effect` at the shipped flags (no scale, no
    gene_scale, no MASS_FIX) -- the three lines that survive when every knob is default."""
    expression = np.expm1(b_ctrl)
    cpm = 1e6 * expression / max(float(expression.sum()), 1e-12)
    gate = (cpm > min_cpm) & (b_ctrl > 0)
    profile = np.clip(np.expm1(b_ctrl + np.where(gate, effect, 0.0)), 0.0, None) / TS
    return profile / float(profile.sum())


# ===========================================================================================
# THE BOARD RECORD.  Every `d_reach` is a scored, deterministic board measurement.
# ===========================================================================================
ARMS = [
    # name          incumbent      d_reach scaled   what changed
    ("a8",          "a2-prune",    +0.021398, "emission law: overdispersion 0 -> 4.0"),
    ("a10",         "a8",          -0.000754, "emission law: spread 0.25->0.15, phi 4->16"),
    ("a21",         "a15",         -0.018482, "emission law: RESAMPLE_X 0.14 (86% resampled)"),
    ("a22",         "a15",         -0.010656, "emission law: resample_mix, the complement"),
    ("a23",         "a15",         -0.007933, "renderer: VCC_MASS_FIX=minnorm"),
    ("rc_hybrid",   "a15",         +0.054940, "emission law: |z| re-assigned by |H|"),
    ("rc_callset",  "rc_hybrid",   +0.022590, "emission law + call set (n_pred 4480->2599)"),
    ("cd4f034",     "rc_callset",  +0.026402, "PROFILE: VCC_CD4_TARGET_F 0.10 -> 0.34"),
    ("h1union",     "cd4f034",     +0.003639, "PROFILE: VCC_H1_UNION_G = 0.28"),
]
# 0.22293924023183273 (h1union, leaderboard 2026-09-04T0524Z) - 0.219300 (cd4f034, 210)
CELL_CHANNEL = {          # emitted per-cell library CV, measured in emit230.py
    "a8": 0.0043, "a10": 0.0187, "a21": 0.3821, "a22": 0.1050, "a23": 0.0187,
    "rc_hybrid": 0.0187, "rc_callset": 0.0187, "cd4f034": 0.0187, "h1union": 0.0187,
}
FRAC_RESAMPLED = {"a21": 0.86, "a22": 0.14}


def measure_pair(Ta_all, Tb_all, live) -> dict:
    rs = [magnitudes(Ta_all[i], Tb_all[i], live) for i in range(Ta_all.shape[0])]
    return {k: float(np.mean([r[k] for r in rs])) for k in
            ("frac_moved", "relL1", "TV", "max_abs_d")} | {
        "digest_same": int(sum(r["digest_same"] for r in rs)), "blocks": len(rs)}


def main() -> None:
    ka, mags = [], {}
    print("=== PHASE 1c — magnitudes measured on the ARTIFACTS ===\n")

    # ---- rc_callset vs a15: the identity 208 gate 14 recorded ---------------------------
    za = np.load(ROOT / "experiments/196-emitted-depth/work/chal-totals-a15.npz",
                 allow_pickle=True)
    a15_blocks = np.asarray(za["blocks"]).astype(str)
    a15_tot = np.asarray(za["totals"]).astype(np.int64)
    zc = np.load(ROOT / "experiments/214-cheap-positives/work/totals-rc_callset.npz",
                 allow_pickle=True)
    z4 = np.load(ROOT / "experiments/214-cheap-positives/work/totals-cd4f034.npz",
                 allow_pickle=True)
    order = {b: i for i, b in enumerate(np.asarray(zc["blocks"]).astype(str))}
    idx = np.array([order[b] for b in a15_blocks])
    rc_tot = np.asarray(zc["totals"]).astype(np.int64)[idx]
    o4 = {b: i for i, b in enumerate(np.asarray(z4["blocks"]).astype(str))}
    cd_tot = np.asarray(z4["totals"]).astype(np.int64)[np.array([o4[b] for b in a15_blocks])]

    same = int(sum(sha16(a15_tot[i]) == sha16(rc_tot[i]) for i in range(900)))
    ka.append(("KA-A  rc_callset totals == a15's (208 gate 14)", f"{same}/900", "900/900",
               same == 900))
    d0 = int(sum(sha16(cd_tot[i]) == sha16(rc_tot[i]) for i in range(900)))
    ka.append(("KA-B  cd4f034 totals != rc_callset's (214/210 f3a)", f"{d0}/900", "0/900",
               d0 == 0))
    ka.append(("KA-C  cd4f034 (A,ABCD1) digest (214 recorded)",
               sha16(cd_tot[int(np.flatnonzero(a15_blocks == "A|ABCD1")[0])]),
               "2cdbc2f0c36ba28b",
               sha16(cd_tot[int(np.flatnonzero(a15_blocks == "A|ABCD1")[0])])
               == "2cdbc2f0c36ba28b"))

    live = np.ones(a15_tot.shape[1], dtype=bool)
    mags["cd4f034"] = measure_pair(rc_tot, cd_tot, live)

    # ---- h1union vs cd4f034, closed form on 222's own effect matrices --------------------
    eff_off = np.load(ROOT / "experiments/222-h1-union-build/work/effects-off.npy")
    eff_pre = np.load(ROOT / "experiments/222-h1-union-build/work/effects-pre.npy")
    eff_on = np.load(ROOT / "experiments/222-h1-union-build/work/effects-on028.npy")
    ka.append(("KA-D  222 effects-off == effects-pre (default path)",
               str(np.array_equal(eff_off, eff_pre)), "True",
               bool(np.array_equal(eff_off, eff_pre))))

    rows_h1, rows_ka = [], []
    for ctx in CTX:
        D = recon196.load_chal(ctx)
        b, r, gate = D["b"], D["r"], D["gate"]
        sel = np.flatnonzero(np.array([x.split("|")[0] for x in a15_blocks]) == ctx)
        tgt = list(np.asarray(D["targets"]).astype(str))
        for j, i in enumerate(sel):
            t = tgt.index(a15_blocks[i].split("|")[1])
            # `nL` is a function of the picked control libraries only, so the two arms share
            # it.  Here the profile is known EXACTLY (it is a closed form of the effect), so
            # `nL` is pinned by intersecting `T_g == floor(nL p_g)` over all 18,533 genes
            # rather than reconstructed: each gene gives `nL in [T_g/p_g, (T_g+1)/p_g)`, and
            # the intersection being non-empty is itself proof that this profile is the one
            # the artifact was rendered from.
            p_off = profile_from_effect(eff_off[t], b)
            nz = p_off > 0
            lo = float(np.max(cd_tot[i][nz] / p_off[nz]))
            hi = float(np.min((cd_tot[i][nz] + 1.0) / p_off[nz]))
            if not lo < hi:
                raise SystemExit(f"nL constraint set EMPTY on {a15_blocks[i]} — the effect "
                                 "matrix or b_ctrl is not the one cd4f034 was rendered from")
            nL = 0.5 * (lo + hi)
            T_off = np.floor(nL * p_off).astype(np.int64)
            T_on = np.floor(nL * profile_from_effect(eff_on[t], b)).astype(np.int64)
            rows_ka.append(float(np.mean(T_off == cd_tot[i])))
            rows_h1.append(magnitudes(T_off, T_on, np.ones(b.size, bool)))
    ka.append(("KA-E  effects->totals reproduces cd4f034's totals",
               f"{np.mean(rows_ka):.8f}", "1.00000000",
               float(np.min(rows_ka)) == 1.0))
    mags["h1union"] = {k: float(np.mean([x[k] for x in rows_h1])) for k in
                       ("frac_moved", "relL1", "TV", "max_abs_d")}
    mags["h1union"]["digest_same"] = int(sum(x["digest_same"] for x in rows_h1))
    mags["h1union"]["blocks"] = len(rows_h1)

    # ---- the rest, from the two Phase-1 measurements ------------------------------------
    M = json.loads((WORK / "magnitude230.json").read_text())["summary"]
    E = json.loads((WORK / "emit230.json").read_text())["summary"]
    for a in ("a8", "a10", "a21", "a22", "rc_hybrid", "rc_callset"):
        mags[a] = dict(frac_moved=0.0, relL1=0.0, TV=0.0, max_abs_d=0.0,
                       digest_same=900, blocks=900)
    mags["a23"] = dict(M["a23"], digest_same=0, blocks=900)
    mags["ham"] = dict(M["ham"], digest_same=0, blocks=900)
    ka.append(("KA-F  emit230: every emission law kept the totals",
               f"{sum(E[l]['same'] for l in E)}/{sum(E[l]['n'] for l in E)}", "24/24",
               all(E[l]["same"] == E[l]["n"] for l in E)))

    ok = all(k[3] for k in ka)
    print(f"{'known-answer row':<52}{'got':>20}{'want':>20}")
    for name, got, want, good in ka:
        print(f"  [{'PASS' if good else 'FAIL'}] {name:<46}{got:>20}{want:>20}")
    if not ok:
        raise SystemExit("known-answer failure — STOP")

    # ---- the table the fit is run on ----------------------------------------------------
    print(f"\n=== EVERY SCORED ARM WITH A KNOWN EMISSION LAW ===")
    print(f"{'arm':<12}{'vs':<12}{'digest':>10}{'frac moved':>12}{'relL1':>11}{'TV':>11}"
          f"{'lib CV':>9}{'board d_reach':>15}")
    tab = []
    for name, inc, dr, what in ARMS:
        m = mags[name]
        tab.append(dict(arm=name, incumbent=inc, d_reach=dr, what=what,
                        relL1=m["relL1"], TV=m["TV"], frac_moved=m["frac_moved"],
                        digest_same=m["digest_same"], lib_cv=CELL_CHANNEL[name],
                        frac_resampled=FRAC_RESAMPLED.get(name, 0.0)))
        print(f"{name:<12}{inc:<12}{m['digest_same']:>6}/900{m['frac_moved']:>12.4f}"
              f"{m['relL1']:>11.3e}{m['TV']:>11.3e}{CELL_CHANNEL[name]:>9.4f}"
              f"{dr:>+15.6f}")
    hm = mags["ham"]
    print(f"{'ham (208 A)':<12}{'cd4f034':<12}{0:>6}/900{hm['frac_moved']:>12.4f}"
          f"{hm['relL1']:>11.3e}{hm['TV']:>11.3e}{0.0187:>9.4f}{'?':>15}")

    # ---- STEP vs LINEAR vs SATURATING ---------------------------------------------------
    y = np.array([t["d_reach"] for t in tab])
    axes = {"relL1": np.array([t["relL1"] for t in tab]),
            "TV": np.array([t["TV"] for t in tab]),
            "frac_moved": np.array([t["frac_moved"] for t in tab])}
    print(f"\n=== STEP vs LINEAR vs SATURATING, n = {len(y)} board points ===")
    print(f"{'axis':<12}{'law':<12}{'R^2':>9}{'AIC':>9}{'slope':>13}{'toll at ham':>14}"
          f"{'LOO RMSE':>11}")
    fits = {}
    for ax, m in axes.items():
        for law, X in (("STEP", np.column_stack([np.ones_like(m), (m > 0).astype(float)])),
                       ("LINEAR", np.column_stack([np.ones_like(m), m])),
                       ("SATURATE", np.column_stack(
                           [np.ones_like(m), 1.0 - np.exp(-m / max(m[m > 0].min(), 1e-12))]))):
            beta, *_ = np.linalg.lstsq(X, y, rcond=None)
            res = y - X @ beta
            ss = float(res @ res)
            r2 = 1.0 - ss / float(((y - y.mean()) ** 2).sum())
            aic = len(y) * np.log(max(ss, 1e-30) / len(y)) + 2 * X.shape[1]
            loo = []
            for k in range(len(y)):
                keep = np.arange(len(y)) != k
                bk, *_ = np.linalg.lstsq(X[keep], y[keep], rcond=None)
                loo.append(y[k] - X[k] @ bk)
            hx = {"relL1": hm["relL1"], "TV": hm["TV"], "frac_moved": hm["frac_moved"]}[ax]
            hv = (1.0 if law == "STEP" else hx if law == "LINEAR"
                  else 1.0 - np.exp(-hx / max(m[m > 0].min(), 1e-12)))
            fits[f"{ax}|{law}"] = dict(r2=r2, aic=float(aic), beta=beta.tolist(),
                                       toll_at_ham=float(beta[1] * hv),
                                       loo_rmse=float(np.sqrt(np.mean(np.square(loo)))))
            print(f"{ax:<12}{law:<12}{r2:>9.4f}{aic:>9.2f}{beta[1]:>+13.4f}"
                  f"{beta[1]*hv:>+14.6f}{np.sqrt(np.mean(np.square(loo))):>11.6f}")

    # ---- falsifier 1: retrodict a held-out calibration arm ------------------------------
    print("\n=== FALSIFIER 1 — leave out a calibration arm, then predict it ===")
    print(f"{'held out':<12}{'law':<20}{'predicted':>12}{'realised':>12}{'error':>12}")
    f1 = {}
    for held in ("a21", "a22", "a23", "cd4f034"):
        k = [t["arm"] for t in tab].index(held)
        keep = np.arange(len(y)) != k
        for ax, law in (("relL1", "STEP"), ("relL1", "LINEAR")):
            m = axes[ax]
            X = (np.column_stack([np.ones_like(m), (m > 0).astype(float)]) if law == "STEP"
                 else np.column_stack([np.ones_like(m), m]))
            bk, *_ = np.linalg.lstsq(X[keep], y[keep], rcond=None)
            pred = float(X[k] @ bk)
            f1[f"{held}|{law}"] = dict(pred=pred, real=float(y[k]), err=pred - float(y[k]))
            print(f"{held:<12}{ax + ' ' + law:<20}{pred:>+12.6f}{y[k]:>+12.6f}"
                  f"{pred - y[k]:>+12.6f}")

    # ---- falsifier 2: the zero-magnitude cloud ------------------------------------------
    z = y[axes["relL1"] == 0.0]
    nz = y[axes["relL1"] > 0.0]
    print(f"\n=== FALSIFIER 2 — the spread WITHIN each digest class ===")
    print(f"  digest-IDENTICAL arms (n={z.size}): d_reach in "
          f"[{z.min():+.6f}, {z.max():+.6f}], mean {z.mean():+.6f}, sd {z.std(ddof=1):.6f}")
    print(f"  digest-DIFFERING  arms (n={nz.size}): d_reach in "
          f"[{nz.min():+.6f}, {nz.max():+.6f}], mean {nz.mean():+.6f}, "
          f"sd {nz.std(ddof=1):.6f}")
    print(f"  rule 13's charge is {-0.0146:+.6f}; it is calibrated on a21+a22, BOTH of which "
          f"are digest-IDENTICAL.")
    sep = dict(zero_mean=float(z.mean()), zero_sd=float(z.std(ddof=1)),
               nonzero_mean=float(nz.mean()), nonzero_sd=float(nz.std(ddof=1)),
               zero_n=int(z.size), nonzero_n=int(nz.size))

    json.dump(dict(table=tab, ham=hm, fits=fits, falsifier1=f1, falsifier2=sep,
                   ka=[dict(name=n, got=g, want=w, ok=bool(o)) for n, g, w, o in ka]),
              open(WORK / "fit230.json", "w"), indent=1)
    print("\nwrote work/fit230.json")


if __name__ == "__main__":
    main()
