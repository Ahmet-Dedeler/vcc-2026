#!/usr/bin/env python3
"""494 (third session) — per-target q80 protection by design-strength tercile, priced analytically through 490's exact
reconstruction on 435's three coarse-MDE H1 seeds (calls / p / order held; cell-eval2's own DE functions).

q80 (395) keeps the DESIGN sign on each target's strongest 20 % of block genes (|E0| above the block's 80th percentile);
the rest take the CLEAN4 prior. For design-weak targets those design signs are mostly noise (467b: the weak third wants
CLEAN4). Arms (protection quantile for weak / mid / strong design-strength tercile; strength = the design row's non-panel
L2 norm, as 492d's strength):
  base  0.8 / 0.8 / 0.8  (shipped; known answer: == delivered, Δ 0)
  W     1.0 / 0.8 / 0.8  (no protection on the weak third)
  W7    1.0 / 0.8 / 0.7  (and more protection on the strong third)
  S7    0.8 / 0.8 / 0.7  (more protection on the strong third only)
  R     0.7 / 0.8 / 1.0  (CONTROL: the reverse)
The protection change re-derives the kf shrink and the MDE balance inside recon_h1.reconstruct (its module-level
PROTECT_Q is set per target). Exploratory; H1 prices DE ordering changes at ~0.4x on the board and sign-prior reach not at
all (docs/H1-BOARD-CALIBRATION.md) — read the members, not a board number.

    uv run python experiments/494-pertarget-protection/price494.py
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
import recon_h1 as RC  # noqa: E402
import score_de  # noqa: E402

SP = {"fid": 0.30, "reach": 0.90, "nmae": -0.60, "jac": 0.37}
ARMS = {"base": (0.8, 0.8, 0.8), "W": (1.0, 0.8, 0.8), "W7": (1.0, 0.8, 0.7), "S7": (0.8, 0.8, 0.7), "R": (0.7, 0.8, 1.0)}


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
        b = RC.R.load_module("builder494", ROOT / "scripts/make_submission_pm_cd4.py")
    with np.load(RC.EFF, allow_pickle=False) as z:
        eff = RC.R.map_axis(np.asarray(z["cap15"], dtype=np.float64), z["genes"], r["genes"])
    prior_g = RC.prior_signs_full(r["genes"])[keep]
    cnt = ref_score * float(np.median(r["libraries"])) / 1e6
    madd = np.where(cnt >= RC.MDE_MIN, np.minimum(RC.MDE_C / np.sqrt(np.maximum(cnt, 1e-3)), RC.MDE_CAP), 0.0)
    T = r["targets"].astype(str)
    G = np.asarray(r["genes"])[keep].astype(str)
    panel = np.isin(np.asarray(r["genes"]).astype(str), T)
    strength = np.linalg.norm(np.where(panel[None, :], 0.0, eff), axis=1)
    ter = np.digitize(strength, np.quantile(strength, [1 / 3, 2 / 3]))
    lfc_del = st["lfc"]
    m_del = (ref_score + 1e-9) * np.power(2.0, lfc_del) - 1e-9
    raw_o, per_o = score_de.score(T, G, st["p"], st["q"], lfc_del)
    want = json.loads((RC.H1 / f"metrics-{tag}.json").read_text())["raw"]
    assert max(abs(raw_o[m] - want[m]) for m in raw_o) <= 1e-9, "harness tables"
    res = {"tercile_counts": np.bincount(ter, minlength=3).tolist()}
    rng = np.random.default_rng(494)
    for arm, qs in ARMS.items():
        pick_rng, _ = b.make_rngs(seed)
        lfc_new = lfc_del.copy()
        for i, t in enumerate(r["targets"]):
            picks = pick_rng.choice(len(r["libraries"]), 400, replace=False)
            RC.PROTECT_Q = 0.8
            o = RC.reconstruct(i, diag[i], eff[i], r, keep, ref_emit, ref_score, prior_g, madd, b, picks, fix=False)
            if arm == "base":
                for name, (got, want_) in o["ka"].items():
                    if want_ is not None and abs(float(got) - float(want_)) > (1e-9 if name.startswith(("kf", "mde_f")) else 0):
                        raise AssertionError(f"seed {seed} {t}: known answer {name} {got} != {want_}")
            RC.PROTECT_Q = float(qs[ter[i]])
            f = RC.reconstruct(i, diag[i], eff[i], r, keep, ref_emit, ref_score, prior_g, madd, b, picks, fix=False)
            RC.PROTECT_Q = 0.8
            m_new = m_del[i] + (o["n_car"] / 400.0) * (f["A_int"] - o["A_int"])
            lfc_new[i] = np.log2((np.maximum(m_new, 0.0) + 1e-9) / (ref_score + 1e-9))
        lfc_new = np.where(np.isfinite(lfc_del), lfc_new, lfc_del)
        raw_n, per_n = score_de.score(T, G, st["p"], st["q"], lfc_new)
        D = {m: np.array([per_n[m].get(t, np.nan) - per_o[m].get(t, np.nan) for t in T]) for m in SP}
        if arm == "base":
            assert all(np.nanmax(np.abs(D[m])) == 0.0 for m in SP), "KA: the base arm must reproduce the delivered tables"
        ov = lambda ii: sum(np.nanmean(D[m][ii]) / SP[m] for m in SP) / 6
        bo = np.array([ov(rng.integers(0, len(T), len(T))) for _ in range(2000)])
        res[arm] = {**{m: round(float(np.nanmean(D[m])), 6) for m in SP}, "overall_equiv": round(float(ov(np.arange(len(T)))), 6),
                    "ci95": [round(float(np.quantile(bo, .025)), 6), round(float(np.quantile(bo, .975)), 6)],
                    "by_tercile_reach": [round(float(np.nanmean(D["reach"][ter == k])), 5) for k in range(3)]}
        print(f"seed {seed} {arm:5s} {qs}: " + " ".join(f"{m} {res[arm][m]:+.5f}" for m in SP)
              + f" | ov {res[arm]['overall_equiv']:+.5f} {res[arm]['ci95']} | reach by tercile {res[arm]['by_tercile_reach']}", flush=True)
    return res


def main():
    t0 = time.monotonic()
    out = {f"seed{s}": run_seed(s) for s in (0, 1, 2)}
    out["seconds"] = time.monotonic() - t0
    (HERE / "work" / "price494.json").write_text(json.dumps(out, indent=1))
    print(f"done {out['seconds']:.0f}s")


if __name__ == "__main__":
    main()
