#!/usr/bin/env python3
"""490 X-475: an independent, exact pricing of 484's q80 patch (475's arm) through recon_h1.py, before 475 reports.

For each of 435's coarse-MDE seeds: reconstruct the emitter's intended carrier arithmetic means with the production mask
(`orig`, known answer = the run's own diag) and with 484's patch (`fix`: protected genes out of `side`, kf re-solved,
the MDE's fu/fd and coarse split re-derived). The patched DELIVERED mean over all 400 cells is the delivered one plus
(n_car / 400) x (A_fix - A_orig): the ballast cells are carried over exactly as delivered. Calls, p-values and the
reach order are held at the delivered ones. That is exact for the lfc channel; it ignores the second-order change of
the |z| realisation when the per-class split moves (475's own run measures that; its preregistration expects
|delta| <= 0.003 on the ordering members). Readout mirrors experiments/475-q80-protect/summary475.py: members, the
overall-equiv with its slopes, all 85 targets and the 67 non-essential ones, target bootstrap 3,000 draws.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import recon_h1 as RC  # noqa: E402
import score_de  # noqa: E402

ROOT = HERE.parents[1]
SP = {"fid": 0.30, "reach": 0.90, "nmae": -0.60, "jac": 0.37}      # 475's slopes (pds is 0 by construction)
ESS = {l.strip() for l in open(ROOT / "data/donors/k562_essential/perts.txt") if l.strip()}


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
        b = RC.R.load_module("builder490p", ROOT / "scripts/make_submission_pm_cd4.py")
    with np.load(RC.EFF, allow_pickle=False) as z:
        eff = RC.R.map_axis(np.asarray(z["cap15"], dtype=np.float64), z["genes"], r["genes"])
    prior_g = RC.prior_signs_full(r["genes"])[keep]
    cnt = ref_score * float(np.median(r["libraries"])) / 1e6
    madd = np.where(cnt >= RC.MDE_MIN, np.minimum(RC.MDE_C / np.sqrt(np.maximum(cnt, 1e-3)), RC.MDE_CAP), 0.0)
    pick_rng, _ = b.make_rngs(seed)
    lfc_del = st["lfc"]
    lfc_fix = lfc_del.copy()
    kf_o, kf_f, lam = [], [], []
    for i, t in enumerate(r["targets"]):
        picks = pick_rng.choice(len(r["libraries"]), 400, replace=False)
        o = RC.reconstruct(i, diag[i], eff[i], r, keep, ref_emit, ref_score, prior_g, madd, b, picks, fix=False)
        for name, (got, want) in o["ka"].items():
            if want is not None and abs(float(got) - float(want)) > (1e-9 if name.startswith(("kf", "mde_f")) else 0):
                raise AssertionError(f"seed {seed} {t}: known answer {name} {got} != {want}")
        f = RC.reconstruct(i, diag[i], eff[i], r, keep, ref_emit, ref_score, prior_g, madd, b, picks, fix=True)
        m_del = (ref_score + 1e-9) * np.power(2.0, lfc_del[i]) - 1e-9
        m_fix = m_del + (o["n_car"] / 400.0) * (f["A_int"] - o["A_int"])
        lfc_fix[i] = np.log2((np.maximum(m_fix, 0.0) + 1e-9) / (ref_score + 1e-9))
        kf_o.append(o["kf"]); kf_f.append(f["kf"]); lam.append(o["lam"])
    T = r["targets"].astype(str)
    G = np.asarray(r["genes"])[keep].astype(str)
    raw_o, per_o = score_de.score(T, G, st["p"], st["q"], lfc_del)
    raw_f, per_f = score_de.score(T, G, st["p"], st["q"], lfc_fix)
    want = json.loads((RC.H1 / f"metrics-{tag}.json").read_text())["raw"]
    ka = max(abs(raw_o[m] - want[m]) for m in raw_o)
    if ka > 1e-9:
        raise AssertionError(f"seed {seed}: delivered tables do not reproduce the harness ({ka})")
    rng = np.random.default_rng(475)
    res = {"kf_orig_median": float(np.median(kf_o)), "kf_fix_median": float(np.median(kf_f)),
           "prior_lambda_lt1_blocks": int((np.array(lam) < 1).sum())}
    for sub in ("all", "nonessential"):
        idx = [k for k, t in enumerate(T) if sub == "all" or t not in ESS]
        D = {m: np.array([per_f[m].get(T[k], np.nan) - per_o[m].get(T[k], np.nan) for k in idx]) for m in SP}
        mem = {m: float(np.nanmean(D[m])) for m in SP}
        ov = lambda ii: sum(np.nanmean(D[m][ii]) / SP[m] for m in SP) / 6
        bo = np.array([ov(rng.integers(0, len(idx), len(idx))) for _ in range(3000)])
        res[sub] = {**{m: round(v, 6) for m, v in mem.items()}, "overall_equiv": round(float(ov(np.arange(len(idx)))), 6),
                    "ci95": [round(float(np.quantile(bo, .025)), 6), round(float(np.quantile(bo, .975)), 6)],
                    "p_le_0": round(float((bo <= 0).mean()), 3), "n_targets": len(idx)}
    return res


def main():
    t0 = time.monotonic()
    out = {}
    for seed in (0, 1, 2):
        out[f"seed{seed}"] = run_seed(seed)
        print(seed, json.dumps(out[f"seed{seed}"]), flush=True)
    out["seconds"] = time.monotonic() - t0
    (HERE / "work" / "price475.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
