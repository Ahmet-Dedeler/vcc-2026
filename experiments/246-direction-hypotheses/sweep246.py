"""246 Phase 2 — one PLOLO fold per process: known-answer rows, then the arms.

    uv run python sweep246.py --fold k562_gwps          (then orion_hct116, orion_hek293t, kolf2)

KNOWN-ANSWER ROWS (abort on failure, tolerance 1e-9 on cos_agg / a23_cos / pds):
  KA-1  `f0`            == resid218 `f0`                     (whitened base, nothing added)
  KA-2  `self0.1_bare`  == resid218 `self0.1`                (099's self-residual, bare base)
  KA-3  `cd4only`       == compose218 `cd40.34+none0`
  KA-4  `STACK`         == compose218 `cd40.34+h1PAIRCARRY0.28`   (the shipped stack;
        on the three cancer folds this is 218's +7.6803 % macro / +0.004887 pds row)
  KA-5  `A_none@0` (deflation mode 'none', dose 0) is BIT-IDENTICAL to `STACK`.

ARMS, all on top of STACK unless the name says otherwise:
  A_<mode><k>@f     self-residual (099/218) with the shared subspace removed, dose f
  A_<mode><k>_SHUF@f   target-shuffled (derangement) control
  A_rand<k>@f       random-subspace control (energy removal without direction selection)
  B_cd4_<mode><k>   STACK with the CD4 residual deflated at the SHIPPED solved fraction
  B_h1_<mode><k>    STACK with the H1 union deflated at the shipped solved fraction
  B_both_<mode><k>  both
  C_partial[_eq]@g  CD4's one-/two-condition rows (production discards them) as a third
                    residual, carrying-row sizing; _SHUF and the row-matched SELF control
"""
from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import lib246 as L  # noqa: E402

TOL = 1e-9
DOSES_A = (0.05, 0.10, 0.20)
MODES_A = [("none", 0), ("centre", 0), ("svd", 1), ("svd", 2), ("svd", 5), ("svd", 10),
           ("svd", 20), ("op", 2), ("op", 10)]
SHUF_MODES = [("none", 0), ("centre", 0), ("svd", 5), ("svd", 10)]
RAND_K = (1, 2, 5, 10, 20)
MODES_B_CD4 = [("centre", 0), ("svd", 1), ("svd", 2), ("svd", 5), ("svd", 10)]
MODES_B_H1 = [("centre", 0), ("svd", 1), ("svd", 2), ("svd", 5)]
MODES_B_BOTH = [("svd", 2), ("svd", 5)]
DOSES_C = (0.17, 0.34)


def peak_gib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold", required=True)
    args = ap.parse_args()
    held = args.fold
    t0 = time.time()

    ref_resid = {(r["fold"], r["arm"]): r
                 for r in json.loads((L.E218 / "work/resid218.json").read_text())}
    ref_comp = {(r["fold"], r["arm"]): r
                for r in json.loads((L.E218 / "work/compose218.json").read_text())}

    st = L.stage(held)
    fold, targets = st["fold"], st["targets"]
    print(f"[{held}] staged: {len(targets)} targets, cd4 frac {st['frac_cd4']:.5f}, "
          f"union frac {st['frac_x']:.5f}  ({time.time()-t0:.0f}s, peak {peak_gib():.2f} GiB)",
          flush=True)

    rows: list[dict] = []
    per_target: dict[str, np.ndarray] = {}

    def run(name: str, parts, extra: dict | None = None, ka=None) -> dict:
        t1 = time.time()
        final = L.compose(st, parts)
        e = L.score(st, final)
        per_target[f"{name}|pds"] = np.asarray(e["_pds_per_target_a23"])
        per_target[f"{name}|cos"] = np.asarray(e["_cos_per_target_a23"])
        row = dict(fold=held, arm=name, n_targets=len(targets),
                   **{k: v for k, v in e.items() if not k.startswith("_")},
                   audit=extra or {}, seconds=round(time.time() - t1, 1))
        rows.append(row)
        msg = (f"[{held}/{name}] cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} "
               f"pds {e['pds']:.6f} q_w {e['a15_sgn_w']:.6f} ({row['seconds']}s)")
        if ka is not None:
            src, key = ka
            r = src.get((held, key))
            assert r is not None, f"no reference row {key} for {held}"
            d = max(abs(r["cos_agg"] - e["cos_agg"]), abs(r["a23_cos"] - e["a23_cos"]),
                    abs(r["pds"] - e["pds"]))
            assert d < TOL, f"KA FAILED {name} vs {key}: max|d| {d:.3e}"
            msg += f"   KA PASS vs {key} (max|d| {d:.1e})"
            row["ka"] = dict(ref=key, max_abs_diff=float(d))
        print(msg, flush=True)
        return row

    # ------------------------------------------------------------------ known answers
    run("f0", [], ka=(ref_resid, "f0"))
    self_res = L.self_residual(st["cd4"], st["pre"])
    fs, _ = L.solve_fraction(st["whitened"], self_res, 0.10)
    run("self0.1_bare", [(self_res, fs)], dict(fraction=float(fs)),
        ka=(ref_resid, "self0.1"))
    run("cd4only", [(st["cd4"], st["frac_cd4"])], ka=(ref_comp, "cd40.34+none0"))
    stack_row = run("STACK", L.stack_parts(st), dict(frac_cd4=st["frac_cd4"],
                                                    frac_x=st["frac_x"]),
                    ka=(ref_comp, "cd40.34+h1PAIRCARRY0.28"))
    stack_final = L.compose(st, L.stack_parts(st))
    d0, aud0 = L.deflate(self_res, "none")
    k0_final = L.compose(st, L.stack_parts(st) + [(d0, 0.0)])
    assert np.array_equal(k0_final, stack_final), "KA-5 FAILED: k=0 is not the identity"
    print(f"[{held}] KA-5 PASS: mode none / dose 0 is bit-identical to STACK", flush=True)
    del k0_final

    # ------------------------------------------------------------------ family A
    carried_self = L._carried(self_res)
    for mode, k in MODES_A:
        d, aud = L.deflate(self_res, mode, k, op_top=st["op_top"], wcols=fold["wcols"])
        aud.update(L.residual_norms(d))
        for f in DOSES_A:
            fr, per_unit = L.solve_fraction(st["whitened"], d, f)
            assert fr > 0, f"VACUOUS deflated self-residual {mode}{k}"
            run(f"A_{mode}{k if k else ''}@{f:g}", L.stack_parts(st) + [(d, fr)],
                dict(family="A", dose=f, fraction=float(fr),
                     per_unit=float(per_unit), **aud))
        if (mode, k) in SHUF_MODES:
            ds = L.shuffled(d, seed=L.SHUF_SEED)
            for f in DOSES_A:
                fr, _ = L.solve_fraction(st["whitened"], ds, f)
                run(f"A_{mode}{k if k else ''}_SHUF@{f:g}", L.stack_parts(st) + [(ds, fr)],
                    dict(family="A_ctrl", control="shuffled", dose=f,
                         fraction=float(fr), **aud))
    for k in RAND_K:
        for seed_i in ((0, 1) if k == 5 else (0,)):
            d, aud = L.deflate(self_res, "rand", k, seed=L.RAND_SEED + seed_i)
            aud.update(L.residual_norms(d))
            for f in (DOSES_A if k in (5, 10) else (0.10,)):
                fr, _ = L.solve_fraction(st["whitened"], d, f)
                run(f"A_rand{k}{'' if seed_i == 0 else 'b'}@{f:g}",
                    L.stack_parts(st) + [(d, fr)],
                    dict(family="A_ctrl", control="random_subspace",
                         seed=seed_i, dose=f, fraction=float(fr), **aud))

    # ------------------------------------------------------------------ family B
    for mode, k in MODES_B_CD4:
        d, aud = L.deflate(st["cd4"], mode, k)
        aud.update(L.residual_norms(d))
        run(f"B_cd4_{mode}{k if k else ''}", [(d, st["frac_cd4"]), (st["union"], st["frac_x"])],
            dict(family="B", target="cd4", sizing="fixed", **aud))
    for seed_i in (0, 1):
        d, aud = L.deflate(st["cd4"], "rand", 5, seed=L.RAND_SEED + seed_i)
        run(f"B_cd4_rand5{'' if seed_i == 0 else 'b'}",
            [(d, st["frac_cd4"]), (st["union"], st["frac_x"])],
            dict(family="B_ctrl", control="random_subspace", target="cd4",
                 seed=seed_i, sizing="fixed", **aud))
    # the best-looking CD4 deflation with the fraction RE-SOLVED on the deflated object
    d, aud = L.deflate(st["cd4"], "svd", 5)
    fr, _ = L.solve_fraction(st["whitened"], d, L.CD4_F)
    run("B_cd4_svd5_resolved", [(d, fr), (st["union"], st["frac_x"])],
        dict(family="B", target="cd4", sizing="resolved",
             fraction=float(fr), **aud))
    for mode, k in MODES_B_H1:
        d, aud = L.deflate(st["union"], mode, k)
        aud.update(L.residual_norms(d))
        run(f"B_h1_{mode}{k if k else ''}", [(st["cd4"], st["frac_cd4"]), (d, st["frac_x"])],
            dict(family="B", target="h1", sizing="fixed", **aud))
    d, aud = L.deflate(st["union"], "rand", 2)
    run("B_h1_rand2", [(st["cd4"], st["frac_cd4"]), (d, st["frac_x"])],
        dict(family="B_ctrl", control="random_subspace", target="h1",
             sizing="fixed", **aud))
    for mode, k in MODES_B_BOTH:
        dc, ac = L.deflate(st["cd4"], mode, k)
        dh, ah = L.deflate(st["union"], mode, k)
        run(f"B_both_{mode}{k}", [(dc, st["frac_cd4"]), (dh, st["frac_x"])],
            dict(family="B", target="both", mode=mode, k=k, sizing="fixed",
                 cd4=ac, h1=ah))

    # ------------------------------------------------------------------ family C
    for eq in (False, True):
        part = L.cd4_partial(targets, equalise=eq)
        info = part.condition_details[0]
        tag = "C_partial_eq" if eq else "C_partial"
        if info["n_partial"] == 0:
            print(f"[{held}] {tag}: no partial rows on this fold", flush=True)
            continue
        for g in DOSES_C:
            fr, _ = L.solve_fraction_carrying(st["whitened"], part, g)
            assert fr > 0, "VACUOUS partial-CD4 residual"
            run(f"{tag}@{g:g}", L.stack_parts(st) + [(part, fr)],
                dict(family="C", equalise=eq, dose=g, fraction=float(fr), **info,
                     **L.residual_norms(part)))
        if not eq:
            ps = L.shuffled(part, seed=L.SHUF_SEED)
            ctrl = L.restrict_rows(self_res, part)
            for g in DOSES_C:
                fr, _ = L.solve_fraction_carrying(st["whitened"], ps, g)
                run(f"{tag}_SHUF@{g:g}", L.stack_parts(st) + [(ps, fr)],
                    dict(family="C_ctrl", control="shuffled", dose=g, fraction=float(fr),
                         **info))
                fr, _ = L.solve_fraction_carrying(st["whitened"], ctrl, g)
                run(f"{tag}_ROWSELF@{g:g}", L.stack_parts(st) + [(ctrl, fr)],
                    dict(family="C_ctrl", control="row_matched_self", dose=g,
                         fraction=float(fr), **info))

    # ------------------------------------------------------------------ persist
    meta = dict(fold=held, targets=list(targets), n_targets=len(targets),
                carried_self=carried_self.tolist(),
                carried_cd4=L._carried(st["cd4"]).tolist(),
                carried_union=L._carried(st["union"]).tolist(),
                frac_cd4=st["frac_cd4"], frac_x=st["frac_x"],
                op_top_evals=st["op_top_evals"].tolist(),
                seconds=round(time.time() - t0, 1), peak_rss_gib=peak_gib())
    (L.WORK / f"sweep246-{held}.json").write_text(json.dumps(dict(meta=meta, rows=rows),
                                                             indent=1))
    np.savez_compressed(L.WORK / f"pertarget246-{held}.npz", **per_target)
    print(f"[{held}] {len(rows)} arms in {time.time()-t0:.0f}s, peak RSS {peak_gib():.2f} GiB",
          flush=True)


if __name__ == "__main__":
    main()
