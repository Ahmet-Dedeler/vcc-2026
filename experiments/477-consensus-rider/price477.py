#!/usr/bin/env python3
"""477 (third session, at the coordinator's request) — 480's donor-consensus rider priced analytically through 490's exact
reconstruction of the production K-class emitter (recon_h1.py), on 435's three coarse-MDE H1 runs.

Rider (480 §4.6/§7): each carrier's intended arithmetic mean A_pg (the emitter's own frame, after kf, prior and the coarse
MDE, at the recorded doses) becomes A_pg * 2^(k * X_pg), with X = cons (mixed: the panel donors' same-target consensus
log2 ratio, pc 1 CPM, mean over K562-GWPS / Orion HCT116 / Orion HEK293T that carry the target, 0 elsewhere) or
X = s * |cons| (magnitude only, s = the sign of our delivered lfc). The DELIVERED 400-cell mean becomes
m_del + (n_car / 400) * A_pg * (2^(k X) - 1): the ballast cells are carried over exactly as delivered (490's price475
arithmetic). Calls, p-values and the reach order are held at the delivered ones; the lfc and therefore the signs read by
fid / reach / nmae move. That is exact for the lfc channel and ignores class capacity (426/435: adds above the MDE's dose
leaked calls), which only an emitter run can price.

Known answers (abort): recon_h1's per-target known answers (kf, mde_fu/fd, mde_n_up/dn, flip, n_block == the run's diag),
and the delivered tables reproduce the harness's metrics JSON (|d| <= 1e-9). k = 0 must return a delta of exactly 0.

    uv run python experiments/477-consensus-rider/price477.py      (< ~2 GB)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "experiments/490-de-emitter-review"))
sys.path.insert(0, str(ROOT / "experiments/480-nmae-instrument"))
import recon_h1 as RC  # noqa: E402
import score_de  # noqa: E402
from common480 import load_donor_rows  # noqa: E402

SP = {"fid": 0.30, "reach": 0.90, "nmae": -0.60, "jac": 0.37}      # 475's slopes (pds/mse are 0 by construction)
ESS = {l.strip() for l in open(ROOT / "data/donors/k562_essential/perts.txt") if l.strip()}
DON = ["k562_gwps", "orion_hct116", "orion_hek293t"]
ARMS = [("mixed", 0.1), ("mixed", 0.2), ("mixed", 0.3), ("magnitude", 0.1), ("magnitude", 0.2), ("mixed", 0.0)]


def run_seed(seed):
    tag = f"base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse"
    diag = json.loads((RC.H1 / f"diag-{tag}.json").read_text())
    st = np.load(RC.H1 / f"stats-{tag}.npz", allow_pickle=False)
    zr = np.load(RC.REF, allow_pickle=False)
    r = {k: zr[k] for k in ("genes", "targets", "keep", "ref_lin", "emit_ref_lin", "libraries", "b_ctrl")}
    zr.close()
    keep = np.asarray(r["keep"], dtype=bool)
    ref_emit = np.asarray(r["emit_ref_lin"], dtype=np.float64)
    ref_score = np.asarray(r["ref_lin"], dtype=np.float64)
    with RC.R.model_env("uniform"):
        b = RC.R.load_module("builder477", ROOT / "scripts/make_submission_pm_cd4.py")
    with np.load(RC.EFF, allow_pickle=False) as z:
        eff = RC.R.map_axis(np.asarray(z["cap15"], dtype=np.float64), z["genes"], r["genes"])
    prior_g = RC.prior_signs_full(r["genes"])[keep]
    cnt = ref_score * float(np.median(r["libraries"])) / 1e6
    madd = np.where(cnt >= RC.MDE_MIN, np.minimum(RC.MDE_C / np.sqrt(np.maximum(cnt, 1e-3)), RC.MDE_CAP), 0.0)
    T = r["targets"].astype(str)
    G = np.asarray(r["genes"])[keep].astype(str)
    rows = np.stack([load_donor_rows(d, T, G, pc_cpm=1.0)[0] for d in DON])
    with np.errstate(invalid="ignore"):
        cons = np.nan_to_num(np.nanmean(rows, 0))
    pick_rng, _ = b.make_rngs(seed)
    lfc_del = st["lfc"]
    A_orig, n_car = [], []
    for i, t in enumerate(r["targets"]):
        picks = pick_rng.choice(len(r["libraries"]), 400, replace=False)
        o = RC.reconstruct(i, diag[i], eff[i], r, keep, ref_emit, ref_score, prior_g, madd, b, picks, fix=False)
        for name, (got, want) in o["ka"].items():
            if want is not None and abs(float(got) - float(want)) > (1e-9 if name.startswith(("kf", "mde_f")) else 0):
                raise AssertionError(f"seed {seed} {t}: known answer {name} {got} != {want}")
        A_orig.append(o["A_int"]); n_car.append(o["n_car"])
    A_orig = np.asarray(A_orig)
    n_car = np.asarray(n_car, float)
    m_del = (ref_score + 1e-9) * np.power(2.0, lfc_del) - 1e-9
    raw_o, per_o = score_de.score(T, G, st["p"], st["q"], lfc_del)
    want = json.loads((RC.H1 / f"metrics-{tag}.json").read_text())["raw"]
    ka = max(abs(raw_o[m] - want[m]) for m in raw_o)
    if ka > 1e-9:
        raise AssertionError(f"seed {seed}: delivered tables do not reproduce the harness ({ka})")
    s = np.sign(lfc_del)
    s = np.where(s == 0, np.sign(cons), s)
    rng = np.random.default_rng(477)
    res = {"cons_coverage": float(np.mean(cons != 0)), "ka_harness_max_abs": float(ka)}
    for kind, k in ARMS:
        X = cons if kind == "mixed" else s * np.abs(cons)
        m_new = m_del + (n_car[:, None] / 400.0) * A_orig * (np.power(2.0, k * X) - 1.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            lfc_new = np.log2((np.maximum(m_new, 0.0) + 1e-9) / (ref_score + 1e-9))
        lfc_new = np.where(np.isfinite(lfc_del), lfc_new, lfc_del)
        raw_n, per_n = score_de.score(T, G, st["p"], st["q"], lfc_new)
        flips = float(np.mean(np.sign(lfc_new) != np.sign(lfc_del)))
        arm = {"flip_frac": flips}
        for sub in ("all", "nonessential"):
            idx = [q for q, t in enumerate(T) if sub == "all" or t not in ESS]
            D = {m: np.array([per_n[m].get(T[q], np.nan) - per_o[m].get(T[q], np.nan) for q in idx]) for m in SP}
            mem = {m: float(np.nanmean(D[m])) for m in SP}
            ov = lambda ii: sum(np.nanmean(D[m][ii]) / SP[m] for m in SP) / 6
            bo = np.array([ov(rng.integers(0, len(idx), len(idx))) for _ in range(3000)])
            arm[sub] = {**{m: round(v, 6) for m, v in mem.items()}, "overall_equiv": round(float(ov(np.arange(len(idx)))), 6),
                        "ci95": [round(float(np.quantile(bo, .025)), 6), round(float(np.quantile(bo, .975)), 6)],
                        "p_le_0": round(float((bo <= 0).mean()), 3), "n_targets": len(idx)}
        if k == 0.0:
            assert all(abs(arm["all"][m]) < 1e-12 for m in SP), f"k=0 known answer failed: {arm['all']}"
        res[f"{kind}_k{k:g}"] = arm
        print(f"seed {seed} {kind:9s} k {k:.1f}: " + " ".join(f"{m} {arm['all'][m]:+.5f}" for m in SP)
              + f" | overall-equiv {arm['all']['overall_equiv']:+.5f} {arm['all']['ci95']} P<=0 {arm['all']['p_le_0']}"
              + f" | non-ess {arm['nonessential']['overall_equiv']:+.5f} | flips {flips:.4f}", flush=True)
    return res


def main():
    t0 = time.monotonic()
    out = {}
    for seed in (0, 1, 2):
        out[f"seed{seed}"] = run_seed(seed)
    out["seconds"] = time.monotonic() - t0
    (HERE / "work" / "price477.json").write_text(json.dumps(out, indent=1))
    print(f"done {out['seconds']:.0f}s")


if __name__ == "__main__":
    main()
