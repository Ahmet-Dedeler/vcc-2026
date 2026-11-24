#!/usr/bin/env python3
"""553 — DIAGNOSTIC (emission side only; prices nothing): why the production emitter's per-cell placement fails on a
surface. Runs run496.main's per-target setup (the same calls, the surface's inputs via run553's patches, the bulk from the
cached key) and calls 410b's ballast block per target inside try/except, with the emitter's `place_k` instrumented to
record each target's class layout (levels Lk, cells per class nk, the cell totals), the sub-gate mass per class, the
residual the placement leaves (resid = cell total - placed counts) and whether the repair resolved it.
A failed target leaves the shared emission stream in a different state for the targets after it, so per-target rows
after the first failure are illustrative, not the production stream. Nothing reads any truth.

usage: uv run python experiments/553-r7b-harness/diag553.py <surface> <MODE> <seed> [--classes q25,q50,q75|q40,q50,q60|h1|h1ratio|r7bratio|spread<r>]
       [--exclude-targets A,B]   (passed to run553: the targets are not emitted, so the stream equals run553's)
writes work/diag553_<surface>_<MODE>_s<seed>[_<classes>][-x<targets>].json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
surface, mode, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])
cls_opt = sys.argv[sys.argv.index("--classes") + 1] if "--classes" in sys.argv else "q25,q50,q75"
excl = sys.argv[sys.argv.index("--exclude-targets") + 1] if "--exclude-targets" in sys.argv else ""
sys.argv = ["run553.py", surface, mode, "base", str(seed), "--bulk", "ps", "--mde-c", "0.1", "--mde-cap", "0.1",
            "--mde-min-count", "0.3", "--mde-coarse", "--bulk-cache", "use", "--tag", "-diag"]
if excl:
    sys.argv += ["--exclude-targets", excl]
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import run553 as G  # noqa: E402  (module-level setup only; its main() is not run)

M, K = G.M, G.KPROD
if cls_opt != "q25,q50,q75":
    libs = np.asarray(G.R_DICT["libraries"], dtype=np.float64)
    if cls_opt == "h1":
        M.H.CLASSES3 = ((1 / 3, 16078.75), (1 / 3, 20190.5), (1 / 3, 25330.5))
    elif cls_opt.startswith("spread"):
        # geometric classes around the surface's median with q75/q25 = r: (med / sqrt r, med, med * sqrt r)
        r_ = float(cls_opt[len("spread"):])
        med = float(np.quantile(libs, 0.5))
        M.H.CLASSES3 = ((1 / 3, med / np.sqrt(r_)), (1 / 3, med), (1 / 3, med * np.sqrt(r_)))
    elif cls_opt in ("h1ratio", "r7bratio"):
        # the surface's own median library x another context's relative class levels (q25/q50, 1, q75/q50)
        rel = {"h1ratio": (16078.75 / 20190.5, 1.0, 25330.5 / 20190.5),
               "r7bratio": (2815.5 / 11418.0, 1.0, 20378.5 / 11418.0)}[cls_opt]
        med = float(np.quantile(libs, 0.5))
        M.H.CLASSES3 = tuple((1 / 3, med * r) for r in rel)
    else:
        qs = [float(x[1:]) / 100 for x in cls_opt.split(",")]
        q = np.quantile(libs, qs)
        M.H.CLASSES3 = tuple((1 / 3, float(v)) for v in q)
REC = {"surface": surface, "mode": mode, "seed": seed, "classes_option": cls_opt, "classes": [list(c) for c in M.H.CLASSES3],
       "targets": []}
_place = K.place_k
LAST = {}


def place_k(dsg, gidx, sub_block, T_cells, cell_class, Lk, nk, pass2_rng):
    B, dg = _place(dsg, gidx, sub_block, T_cells, cell_class, Lk, nk, pass2_rng)
    resid = T_cells.astype(np.int64) - B.sum(1)
    sub = np.asarray(sub_block).sum(1)
    LAST.update(Lk=[int(v) for v in Lk], nk=[int(v) for v in nk],
                T_min_by_class=[int(T_cells[cell_class == k].min()) for k in range(len(nk))],
                sub_gate_mass_by_class=[int(sub[cell_class == k].sum()) for k in range(len(nk))],
                resid_nonzero=int((resid != 0).sum()), resid_abs_sum=int(np.abs(resid).sum()),
                resid_by_class=[int(resid[cell_class == k].sum()) for k in range(len(nk))], **{k: int(v) for k, v in dg.items()})
    return B, dg


K.place_k = place_k


def main():
    t0 = time.monotonic()
    b = M.H.builder389()
    eff, genes, targets = M.H.native_cap15(b)
    eff_bulk = M.ps_bulk(eff, genes, targets)
    r, cc, emit_c = M.E.load_reference()
    cnt = np.asarray(cc.ref_lin, dtype=np.float64) * float(np.median(r["libraries"])) / 1e6
    M.MDE_ADD = np.where(cnt >= 0.3, np.minimum(0.1 / np.sqrt(np.maximum(cnt, 1e-3)), 0.1), 0.0)
    M.MDE_SCOPE, M.MDE_COARSE = "all", True
    M.H.PRIOR_FULL, _ = M.H.prior_signs_full(r["genes"])
    eff_de = M.R.map_axis(eff, genes, r["genes"])
    eff_bk = M.R.map_axis(eff_bulk, genes, r["genes"])
    pick_rng, emit_rng = b.make_rngs(seed)
    n_fail = 0
    for i, target in enumerate(r["targets"]):
        prof_de = b.profile_from_effect(eff_de[i], r["b_ctrl"], min_cpm=10.)
        prof_bk = b.profile_from_effect(eff_bk[i], r["b_ctrl"], min_cpm=10.)
        picks = pick_rng.choice(len(r["libraries"]), 400, replace=False)
        spread = b.emission_spread_vector(eff_de[i], r["b_ctrl"], min_cpm=10., tail_spread=1., keep_top=8000,
                                          keep_top_frac=None, keep_spread=.25, rank_score=None, spread_map="block",
                                          spread_lo=0., spread_hi=.60)
        gate = b.emission_rank_axis(eff_de[i], r["b_ctrl"], min_cpm=10.)[0]
        LAST.clear()
        row = {"i": i, "target": str(target), "picked_lib_min": float(r["libraries"][picks].min()),
               "picked_lib_median": float(np.median(r["libraries"][picks]))}
        try:
            _, dg = M.ballast_block(prof_de, prof_bk, r["libraries"][picks], emit_rng, spread,
                                    np.where(gate, eff_de[i], 0.), i, emit_c, 16)
            row.update(ok=True, ballast_K=int(dg["ballast_K"]),
                       # the emitter's call-set record (prescribed / delivered / missed calls of this block)
                       **{k: (None if dg.get(k) is None else float(dg[k])) for k in
                          ("call_set_prescribed", "call_set_delivered", "call_set_missed", "call_set_unsuppressible")})
        except AssertionError as e:
            n_fail += 1
            row.update(ok=False, error=str(e))
        if K.LAST_PRIOR_REC:
            pr = K.LAST_PRIOR_REC.pop()
            row["prior_lambda"] = float(pr.get("prior_lambda", float("nan")))
            row["mde_lambda"] = float(pr["mde_lambda"]) if "mde_lambda" in pr else None
        row["placement"] = dict(LAST)
        REC["targets"].append(row)
        if (i + 1) % 20 == 0:
            print(f"{i + 1}/{len(r['targets'])} failures so far {n_fail} ({time.monotonic() - t0:.0f}s)", flush=True)
    REC.update(n_targets=len(REC["targets"]), n_fail=n_fail, seconds=time.monotonic() - t0)
    out = HERE / "work" / (f"diag553_{surface}_{mode}_s{seed}{'' if cls_opt == 'q25,q50,q75' else '_' + cls_opt.replace(',', '-')}"
                           f"{'-x' + excl.replace(',', '+') if excl else ''}.json")
    out.write_text(json.dumps(REC, indent=1))
    fails = [t for t in REC["targets"] if not t["ok"]]
    lam = np.array([t.get("prior_lambda", np.nan) for t in REC["targets"] if t["ok"]], dtype=float)
    mlam = np.array([np.nan if t.get("mde_lambda") is None else t["mde_lambda"] for t in REC["targets"] if t["ok"]], dtype=float)
    miss = np.array([t.get("call_set_missed", np.nan) for t in REC["targets"] if t["ok"]], dtype=float)
    REC.update(prior_lambda_quantiles=[float(np.nanquantile(lam, q)) for q in (0, .25, .5, .75, 1)] if lam.size else None,
               n_prior_lambda_zero=int((lam == 0).sum()), n_mde_lambda_missing=int(np.isnan(mlam).sum()),
               mde_lambda_q10_median=([float(np.nanquantile(mlam, .1)), float(np.nanmedian(mlam))]
                                      if np.isfinite(mlam).any() else None),
               call_set_missed_q10_median_q90=([float(np.nanquantile(miss, q)) for q in (.1, .5, .9)]
                                               if np.isfinite(miss).any() else None))
    out.write_text(json.dumps(REC, indent=1))
    print(f"steering among placed targets: prior lambda q0/q25/q50/q75/q100 {REC['prior_lambda_quantiles']}; "
          f"lambda 0 on {REC['n_prior_lambda_zero']}; no MDE dose on {REC['n_mde_lambda_missing']}; MDE lambda q10/median "
          f"{REC['mde_lambda_q10_median']}; call-set missed q10/median/q90 {REC['call_set_missed_q10_median_q90']}")
    print(f"{surface} {mode} s{seed} classes {cls_opt} {REC['classes']}: {n_fail}/{len(REC['targets'])} targets fail the "
          f"per-cell placement")
    for t in fails[:5]:
        print(json.dumps({k: t[k] for k in ("i", "target", "picked_lib_min", "picked_lib_median")}), json.dumps(t["placement"]))


if __name__ == "__main__":
    main()
