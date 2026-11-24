"""222 KNOWN-ANSWER ROW 3 — reproduce 218's headline through the PRODUCTION BUILDER's code.

218 measured the arm with its own residual builder (`resid218.multi_source_residual`,
`compose218.union_rows`, `compose218.solve_fraction_carrying`).  The thing that will actually
be shipped is `scripts/make_submission_pm_cd4.py`'s `prepare_h1_union` /
`solve_carrying_row_fraction`.  Those are two independent implementations of the same object,
so running 162's PLOLO instrument with **the production ones substituted in** is a real
known-answer row: if the union the builder makes differs from the union 218 measured -- a
different gene axis, a different panel exclusion, a different normalisation, a different
average over the 7 shared targets -- these numbers will not come back.

Reproduced, against `experiments/218-independent-donors/work/`:

  KA-3a  `h1PAIRCARRY` at `g = 0.28`, three cancer-line folds: relative MACRO cosine
         **+7.680%** and `d pds` **+0.004887**   (`diag218.json`, `rel_macro_nokolf` /
         `d_pds_nokolf`).
  KA-3b  per-fold `cos_agg` / `a23_cos` / `pds` for all four folds and all six doses
         (`compose218.json`).
  KA-3c  the NON-CARRYING-ROW bound: `gse295214_h1` alone, doses 0.17 / 0.28 / 0.45, four
         folds -- the largest `rel_out` excursion must be <= 0.018 pp (`rowsplit218.json`).
  KA-3d  the target-shuffled control must be strictly negative at every dose
         (`diag218.json`, `h1PAIRCARRY_SHUF`).

Nothing here writes outside this experiment's `work/`.  The fold caches were copied in from
218 so `build_fold` is a pure cache hit.

    uv run python repro222.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E218 = ROOT / "experiments/218-independent-donors"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(E218))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"

from plolo import (E062, FOLDS, OFFICIAL, PANEL, PROD_OPS, Operator,  # noqa: E402
                   build_fold, cd4_for, deliver, solve_fraction)
from resid218 import (ALPHA, GAMMA, LAM, score_effects, shuffled,  # noqa: E402
                      source_residual)
from compose218 import merged  # noqa: E402

CD4_F = 0.34
G = (0.02, 0.05, 0.10, 0.17, 0.28, 0.45)
ROWSPLIT_G = (0.17, 0.28, 0.45)
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
TOL = 1e-9


def load_builder():
    """The production builder, imported by PATH -- rule 11, and the object under test."""
    path = ROOT / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("repro222_builder", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["repro222_builder"] = module
    spec.loader.exec_module(module)
    return module, path


def per_target_cos(P, T):
    a = np.linalg.norm(P, axis=1)
    b = np.linalg.norm(T, axis=1)
    out = np.zeros(P.shape[0])
    ok = (a > 0) & (b > 0)
    out[ok] = (P[ok] * T[ok]).sum(1) / (a[ok] * b[ok])
    return out


def main() -> None:
    t0 = time.time()
    MOD, builder_path = load_builder()
    print(f"builder {builder_path}\n  sha256 {MOD.file_sha256(builder_path)}", flush=True)

    ref_rows = {(r["fold"], r["arm"]): r
                for r in json.loads((E218 / "work/compose218.json").read_text())}
    ref_macro = json.loads((E218 / "work/diag218.json").read_text())["macro"]
    ref_split = json.loads((E218 / "work/rowsplit218.json").read_text())

    out: dict[str, object] = {"builder_sha256": MOD.file_sha256(builder_path),
                              "arms": [], "rowsplit": [], "failures": []}

    def fail(msg: str) -> None:
        print(f"  *** {msg}", flush=True)
        out["failures"].append(msg)

    for held in FOLDS:
        fold = build_fold(held)
        targets = fold["targets"]
        scored = fold["scored"]
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, scored]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(targets)
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        whitened = op.apply(pre * ALPHA, fold["wcols"], GAMMA, LAM)
        frac_cd4, _ = solve_fraction(whitened, cd4, CD4_F)

        # ---- THE PRODUCTION OBJECTS, not 218's -----------------------------------------
        union, details = MOD.prepare_h1_union(targets, OFFICIAL, PANEL, max_rss_gib=40.0)
        carry_union = np.asarray(union.target_condition_count) > 0
        print(f"[{held}] {len(targets)} targets; production union carries "
              f"{int(carry_union.sum())} "
              f"({int(np.sum(np.asarray(union.target_condition_count) == 2))} by both); "
              f"per-source covered "
              f"{ {d['condition']: d['covered_targets'] for d in details} }", flush=True)

        # 218 built the same union from its own code.  Compare the ARRAYS, not just the
        # endpoints: this is the sharpest statement that the two implementations agree.
        ref_union = MOD.union_residual([source_residual(n, targets, "gene_norm")
                                        for n in ("gse295214_h1", "h1_2025")])
        d_eff = float(np.abs(np.asarray(union.mean_effect)
                             - np.asarray(ref_union.mean_effect)).max())
        d_vote = int(np.abs(np.asarray(union.gene_votes, dtype=int)
                            - np.asarray(ref_union.gene_votes, dtype=int)).max())
        print(f"    production union vs 218's builder: max|d effect| {d_eff:.3e}, "
              f"max|d votes| {d_vote}", flush=True)
        if d_eff > 1e-12 or d_vote != 0:
            fail(f"{held}: production union != 218's union "
                 f"(d_eff {d_eff:.3e}, d_votes {d_vote})")

        base_final, _ = E062.add_post_whiten_residual(
            whitened, merged(cd4, frac_cd4, None, 0.0), 1.0, fold["wcols"])
        base_final[:, fold["panel_positions"]] = 0.0
        base = score_effects(base_final, fold, truth)
        ref = ref_rows.get((held, f"cd4{CD4_F:g}+none0"))
        if ref is not None:
            for key in ("cos_agg", "a23_cos", "pds"):
                delta = abs(ref[key] - base[key])
                if delta > TOL:
                    fail(f"{held}: baseline {key} differs from 218 by {delta:.3e}")
        del base_final

        # the delivered per-target cosine split (KA-3c) uses context A only, as 218 did
        b_ctrl = plolo.challenge_bctrl("A")
        expression = np.expm1(b_ctrl)
        gate = (1e6 * expression / expression.sum() > plolo.MIN_CPM) & (b_ctrl > 0)

        def cos_of(residual, frac_x):
            final, _ = E062.add_post_whiten_residual(
                whitened, merged(cd4, frac_cd4, residual, frac_x), 1.0, fold["wcols"])
            final[:, fold["panel_positions"]] = 0.0
            return per_target_cos(deliver(final, b_ctrl, gate, "minnorm")[:, scored], truth)

        base_cos = cos_of(None, 0.0)

        single = source_residual("gse295214_h1", targets, "gene_norm")
        carry_single = np.asarray(single.target_condition_count) > 0

        for name, residual, carry in (("h1PAIRCARRY", union, carry_union),
                                      ("h1PAIRCARRY_SHUF", shuffled(union), carry_union)):
            for g in (G if name == "h1PAIRCARRY" else ROWSPLIT_G):
                frac_x, frame = MOD.solve_carrying_row_fraction(whitened, residual, g)
                assert frac_x > 0.0 and frame["carrying_rows"] > 0, "vacuous arm"
                final, preserve = E062.add_post_whiten_residual(
                    whitened, merged(cd4, frac_cd4, residual, frac_x), 1.0, fold["wcols"])
                final[:, fold["panel_positions"]] = 0.0
                e = score_effects(final, fold, truth)
                row = dict(fold=held, family=name, g=g, frac_x=float(frac_x),
                           carrying=int(frame["carrying_rows"]),
                           shipped_median_rule_fraction=float(
                               frame["shipped_median_rule_would_solve"]),
                           cos_agg=e["cos_agg"], a23_cos=e["a23_cos"], pds=e["pds"],
                           base_cos_agg=base["cos_agg"], base_a23_cos=base["a23_cos"],
                           base_pds=base["pds"])
                out["arms"].append(row)
                ref = ref_rows.get((held, f"cd4{CD4_F:g}+{name}{g:g}"))
                tag = ""
                if ref is not None:
                    worst = max(abs(ref[k] - e[k]) for k in ("cos_agg", "a23_cos", "pds"))
                    tag = f"  vs218 max|d| {worst:.2e}"
                    if worst > TOL:
                        fail(f"{held}/{name}@{g}: differs from 218 by {worst:.3e}")
                print(f"    {name:18s} g={g:<5g} fx {frac_x:.5f}  "
                      f"macro {100*(e['a23_cos']/base['a23_cos']-1):+7.3f}%  "
                      f"agg {100*(e['cos_agg']/base['cos_agg']-1):+7.3f}%  "
                      f"dpds {e['pds']-base['pds']:+.6f}{tag}", flush=True)
                del final

        # ---- KA-3c: the non-carrying-row bound, 218's own instrument -------------------
        for tag, residual, carry in (("gse295214_h1", single, carry_single),
                                     ("union", union, carry_union)):
            for g in ROWSPLIT_G:
                frac_x, _ = MOD.solve_carrying_row_fraction(whitened, residual, g)
                c = cos_of(residual, frac_x)
                rel_in = float(c[carry].mean() / base_cos[carry].mean() - 1)
                rel_out = float(c[~carry].mean() / base_cos[~carry].mean() - 1)
                out["rowsplit"].append(dict(fold=held, source=tag, g=g,
                                            n_carry=int(carry.sum()),
                                            rel_in=rel_in, rel_out=rel_out,
                                            rel_macro=float(c.mean()/base_cos.mean()-1)))
                note = ""
                if tag == "gse295214_h1":
                    ref = [r for r in ref_split
                           if r["fold"] == held and r["g"] == g and r["arm"] == "real"]
                    if ref:
                        note = (f"  vs218 d rel_in {abs(ref[0]['rel_in']-rel_in):.2e} "
                                f"d rel_out {abs(ref[0]['rel_out']-rel_out):.2e}")
                        if abs(ref[0]["rel_in"] - rel_in) > 1e-9:
                            fail(f"{held}/rowsplit {tag}@{g}: rel_in differs from 218")
                print(f"    rowsplit {tag:14s} g={g:<5g} n_carry {int(carry.sum()):3d}  "
                      f"carrying {100*rel_in:+7.3f}%   NOT-carrying {100*rel_out:+8.5f}%"
                      f"{note}", flush=True)
        del op, whitened, union, single

    # ------------------------------------------------------------------ the headline row
    print("\n" + "=" * 100)
    arms = out["arms"]
    summary = {}
    for name in ("h1PAIRCARRY", "h1PAIRCARRY_SHUF"):
        for g in G:
            sel = [r for r in arms if r["family"] == name and r["g"] == g]
            cancer = [r for r in sel if r["fold"] in CANCER]
            if not cancer:
                continue
            summary[f"{name}@{g:g}"] = dict(
                n_folds=len(sel), n_cancer=len(cancer),
                rel_macro_nokolf=float(np.mean(
                    [r["a23_cos"] / r["base_a23_cos"] - 1 for r in cancer])),
                rel_agg_nokolf=float(np.mean(
                    [r["cos_agg"] / r["base_cos_agg"] - 1 for r in cancer])),
                d_pds_nokolf=float(np.mean([r["pds"] - r["base_pds"] for r in cancer])),
            )
    print(f"{'arm':26s} {'rel MACRO cos':>15s} {'218':>12s} {'rel agg cos':>13s} "
          f"{'d pds':>11s} {'218':>12s}")
    for key, v in summary.items():
        r = ref_macro.get(key.replace("@", "@"))
        rm = r["rel_macro_nokolf"] if r else float("nan")
        rp = r["d_pds_nokolf"] if r else float("nan")
        print(f"{key:26s} {100*v['rel_macro_nokolf']:+15.4f} {100*rm:+12.4f} "
              f"{100*v['rel_agg_nokolf']:+13.4f} {v['d_pds_nokolf']:+11.6f} {rp:+12.6f}")
    out["summary"] = summary

    head = summary["h1PAIRCARRY@0.28"]
    ka = [
        ("KA-3a macro cosine, 3 cancer folds, g=0.28",
         head["rel_macro_nokolf"], 0.07680281456197664, 2e-6),
        ("KA-3a d pds, 3 cancer folds, g=0.28",
         head["d_pds_nokolf"], 0.004886737506532087, 2e-6),
    ]
    print()
    for label, got, want, tol in ka:
        ok = abs(got - want) <= tol
        print(f"{'PASS' if ok else 'FAIL'}  {label:52s} got {got:+.9f}  "
              f"want {want:+.9f}  |d| {abs(got-want):.2e}")
        if not ok:
            fail(f"{label}: got {got:.9f}, want {want:.9f}")

    real_out = [abs(r["rel_out"]) for r in out["rowsplit"]
                if r["source"] == "gse295214_h1"]
    union_out = [abs(r["rel_out"]) for r in out["rowsplit"] if r["source"] == "union"]
    print(f"\nKA-3c largest NON-CARRYING excursion, gse295214_h1 (218's instrument): "
          f"{100*max(real_out):.5f} pp  (218 reported 0.018 pp)")
    print(f"      largest NON-CARRYING excursion, the shipped UNION:              "
          f"{100*max(union_out):.5f} pp")
    ok3c = max(real_out) <= 1.8e-4
    print(f"{'PASS' if ok3c else 'FAIL'}  KA-3c non-carrying rows bounded at 0.018 pp")
    if not ok3c:
        fail("KA-3c: non-carrying rows moved more than 0.018 pp")

    shuf = [v for k, v in summary.items() if k.startswith("h1PAIRCARRY_SHUF")]
    ok3d = all(v["rel_macro_nokolf"] < 0 and v["d_pds_nokolf"] < 0 for v in shuf)
    print(f"{'PASS' if ok3d else 'FAIL'}  KA-3d shuffled control strictly negative at all "
          f"{len(shuf)} doses on macro cosine AND pds")
    if not ok3d:
        fail("KA-3d: the shuffled control is not strictly negative")

    (WORK / "repro222.json").write_text(json.dumps(out, indent=2))
    print(f"\n{len(out['failures'])} failures; wrote {WORK/'repro222.json'} "
          f"({time.time()-t0:.0f}s)")
    if out["failures"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
