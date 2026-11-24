"""229 — price what is actually BUILT, through the production builder's own union.

Three jobs, one instrument (162's PLOLO folds, 218's scoring, 225's rate card):

  KA-4  The union `scripts/make_submission_pm_cd4.py` builds for a given `VCC_UNION_SOURCES`
        must equal, to 0.0e+00, the union 225 priced (`resid218.source_residual` +
        `compose218.union_rows`) on every fold.  Two independent implementations; if they
        disagree the package is not the arm 225 measured and nothing after this reads.

  KA-5  225's own CUM4 row must come back: `base222 + stack4@g` on the three cancer folds,
        against `work/stack225d.json`.  That is 225's frame -- MARGINAL over the shipped
        `h1union` -- and it is NOT the frame a single-union build can ship (§ below).

  THE PRICE  Every construction a registry-driven build can actually emit, scored against
        BOTH baselines, all six members, three cancer folds, `kolf2` reported separately.

WHY THE THIRD JOB EXISTS.  The builder carries ONE extra residual with ONE dose.  225 priced
`stack4` as a SECOND residual riding on top of `h1union` (two residuals, doses 0.28 and 0.35).
A single-union build cannot express that: naming four sources REPLACES the H1 pair, and naming
six unions them at a single shared dose.  So the shipped arm is one of

    ARM4-g   cd4@0.34 + union(dld1, jiang5, igvf, myllia)@g              (H1 pair DROPPED)
    ARM6-g   cd4@0.34 + union(H1 pair + those four)@g                    (one shared dose)

and which of them is worth a slot is a measurement, not a reading.  Both are reported here
against `h1union` (the incumbent, board 0.2118310) and against `cd4f034`.

    uv run python price229.py
"""
from __future__ import annotations

import importlib.util
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
E218 = ROOT / "experiments/218-independent-donors"
E225 = ROOT / "experiments/225-source-sweep"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(E218))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"          # 229's OWN cache (copied from 225); writes nothing to 225

from plolo import (E062, FOLDS, OFFICIAL, PANEL, PROD_OPS, Operator,  # noqa: E402
                   build_fold, cd4_for, solve_fraction)
from resid218 import (ALPHA, GAMMA, LAM, score_effects,  # noqa: E402
                      source_residual)
from compose218 import merged, solve_fraction_carrying, union_rows  # noqa: E402

# --- 225's registry, replicated so `source_residual` can find the four banks --------------
import sources218  # noqa: E402

# MUTATE IN PLACE.  `resid218` does `from sources218 import SOURCES`, so it holds a reference
# to this dict object; rebinding `sources218.SOURCES` would leave `resid218`'s copy stale and
# the four banks unreachable.  (225's harness gets away with rebinding only because it does it
# BEFORE importing `resid218`.)
sources218.SOURCES["dld1_lowmoi"] = dict(
    path="experiments/223-silent-zero-audit/work/tables/dld1_lowmoi")
for _n in ("jiang_5line_noK562", "igvf_cardiomyocyte"):
    sources218.SOURCES[_n] = dict(path=f"experiments/225-source-sweep/work/tables/{_n}")

CD4_F = 0.34
H1_G = 0.28                       # 218/222's shipped dose for the H1 pair, PINNED
H1_PAIR = ("gse295214_h1", "h1_2025")
STACK4 = ("dld1_lowmoi", "jiang_5line_noK562", "igvf_cardiomyocyte", "myllia_a549")
STACK6 = H1_PAIR + STACK4
# INBOX-3's fallback and the two intermediate unions between it and STACK6.  Every one of
# these is a SINGLE union at a SINGLE dose, i.e. something a registry-driven build can emit.
STACK1 = ("dld1_lowmoi",)
STACK3 = H1_PAIR + ("dld1_lowmoi",)
STACK5 = H1_PAIR + ("dld1_lowmoi", "jiang_5line_noK562", "igvf_cardiomyocyte")
GRID = (0.10, 0.17, 0.28, 0.35, 0.45, 0.60)
BUILT_G = 0.35                    # 225's pre-registered dose.  NOT re-tuned here.
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")

# 218 §7.4's rate card == 223's == 225's.  All six members; never four at zero (§8.4).
RATE = {"mse": 0.003634, "reach": 0.002446, "nmae": 0.000980, "fid": 0.000847,
        "jac": 0.000100}
RATE_EX_PDS = sum(RATE.values())
PDS_SPAN = 0.452
GATE_SHAPE_ONLY = 0.0019


def price(d_cos_abs: float, d_pds_raw: float, k: float = 1.0) -> dict:
    """All six members, 225's `harness225.price` verbatim."""
    cos_gain = d_cos_abs * k
    five = RATE_EX_PDS * (cos_gain / 0.01)
    pds_local = (d_pds_raw / PDS_SPAN) / 6.0 * 1.011 * 2.0
    return dict(five_members=five, pds=pds_local, total=five + pds_local)


def peak_gib() -> float:
    return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1024**3


def load_builder():
    path = ROOT / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("price229_builder", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["price229_builder"] = module
    spec.loader.exec_module(module)
    return module, path


def production_union(MOD, names, targets):
    """The union the BUILD makes, for an arbitrary selection, without re-importing."""
    MOD.H1_UNION_SOURCES = tuple(
        (n, MOD.UNION_SOURCE_REGISTRY[n]["path"]) for n in names)
    MOD.verify_union_source_hashes()          # rule 11, enforced on every call
    union, details = MOD.prepare_h1_union(targets, OFFICIAL, PANEL, max_rss_gib=40.0)
    return union, details


def harness_union(names, targets):
    """The union 225 priced -- 218's own implementation, the independent one."""
    return union_rows([source_residual(n, targets, "gene_norm") for n in names])


def max_abs_delta(a, b) -> tuple[float, int]:
    da = float(np.max(np.abs(np.asarray(a.mean_effect) - np.asarray(b.mean_effect))))
    dv = int(np.max(np.abs(np.asarray(a.gene_votes, dtype=np.int64)
                           - np.asarray(b.gene_votes, dtype=np.int64))))
    return da, dv


def main() -> None:
    t0 = time.time()
    MOD, builder_path = load_builder()
    print(f"builder {builder_path}\n  sha256 {MOD.file_sha256(builder_path)}", flush=True)

    ref = {(r["fold"], r["tag"], r.get("variant", ""), r["g"]): r
           for r in json.loads((E225 / "work/stack225d.json").read_text())}

    out: dict[str, object] = {
        "builder_sha256": MOD.file_sha256(builder_path),
        "ka4": [], "ka5": [], "arms": [], "conc": {}, "failures": [],
    }

    def fail(msg: str) -> None:
        print(f"  *** FAIL {msg}", flush=True)
        out["failures"].append(msg)

    per_fold_conc: dict[str, dict] = {}

    for held in FOLDS:
        fold = build_fold(held)
        targets = list(fold["targets"])
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(targets)
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        whitened = op.apply(pre * ALPHA, fold["wcols"], GAMMA, LAM)
        del op
        frac_cd4, _ = solve_fraction(whitened, cd4, CD4_F)

        # ---- KA-4: the production union == 225's union, for every selection -------------
        unions = {}
        for tag, names in (("h1pair", H1_PAIR), ("stack1", STACK1), ("stack3", STACK3),
                           ("stack4", STACK4), ("stack5", STACK5), ("stack6", STACK6)):
            prod, _details = production_union(MOD, names, targets)
            ours = harness_union(names, targets)
            da, dv = max_abs_delta(prod, ours)
            row = dict(fold=held, selection=tag, names=list(names),
                       max_abs_delta_mean_effect=da, max_abs_delta_votes=dv,
                       carrying=int(np.sum(np.asarray(prod.target_condition_count) > 0)),
                       n_rows=len(targets))
            out["ka4"].append(row)
            print(f"[KA-4 {held:>14} {tag:>7}] |d effect| {da:.3e}  |d votes| {dv}  "
                  f"carrying {row['carrying']}/{len(targets)}", flush=True)
            if da != 0.0 or dv != 0:
                fail(f"KA-4 {held}/{tag}: production union != 225's union")
            unions[tag] = prod
            del ours

        frac_h1, _ = solve_fraction_carrying(whitened, unions["h1pair"], H1_G)
        assert frac_h1 > 0.0, "the H1 pair sized to zero -- broken instrument"

        base = {
            "cd4f034": merged(cd4, frac_cd4, None, 0.0),
            "h1union": merged(cd4, frac_cd4, unions["h1pair"], frac_h1),
        }

        def finish(res):
            f, _ = E062.add_post_whiten_residual(whitened, res, 1.0, fold["wcols"])
            f[:, fold["panel_positions"]] = 0.0
            return f

        scored_base = {}
        for name, res in base.items():
            f = finish(res)
            scored_base[name] = score_effects(f, fold, truth)
            del f
        print(f"[{held:>14}] cd4f034 cos {scored_base['cd4f034']['cos']:.6f} "
              f"pds {scored_base['cd4f034']['pds']:.6f} | h1union cos "
              f"{scored_base['h1union']['cos']:.6f} pds {scored_base['h1union']['pds']:.6f}",
              flush=True)

        # ---- the arms -------------------------------------------------------------------
        for arm_tag, sel, on_top in (("ARM1", "stack1", None),
                                     ("ARM3", "stack3", None),
                                     ("ARM4", "stack4", None),
                                     ("ARM5", "stack5", None),
                                     ("ARM6", "stack6", None),
                                     ("CUM1", "stack1", "h1union"),
                                     ("CUM4", "stack4", "h1union")):
            res = unions[sel]
            for g in GRID:
                fr, per_unit = solve_fraction_carrying(whitened, res, g)
                assert fr > 0.0, f"{arm_tag}@{g} sized to zero"
                if on_top is None:
                    m = merged(cd4, frac_cd4, res, fr)
                else:
                    m = merged(base[on_top], 1.0, res, fr)
                f = finish(m)
                e = score_effects(f, fold, truth)
                del f
                added = np.where(np.asarray(res.gene_votes) > 0,
                                 np.asarray(res.mean_effect), 0.0)
                norms = np.linalg.norm(added, axis=1)
                live = norms[norms > 0]
                row = dict(
                    fold=held, arm=arm_tag, selection=sel, on_top=on_top, g=g,
                    fraction=float(fr), per_unit=float(per_unit),
                    carrying_rows=int(live.size), n_rows=int(norms.size),
                    median_row_norm_carrying=float(np.median(live)),
                    median_row_norm_all_rows=float(np.median(norms)),
                    cos=e["cos"], pds=e["pds"], cos_agg=e["cos_agg"], qw=e["a15_sgn_w"],
                )
                for bname, b in scored_base.items():
                    row[f"d_cos_vs_{bname}"] = e["cos"] - b["cos"]
                    row[f"rel_cos_pct_vs_{bname}"] = (
                        100.0 * (e["cos"] - b["cos"]) / max(abs(b["cos"]), 1e-300))
                    row[f"d_pds_vs_{bname}"] = e["pds"] - b["pds"]
                    row[f"d_qw_vs_{bname}"] = e["a15_sgn_w"] - b["a15_sgn_w"]
                out["arms"].append(row)

                # KA-5: 225's CUM4 rows, marginal over base222 == our `h1union`
                if arm_tag in ("CUM1", "CUM4") and held in CANCER:
                    key = (held, arm_tag, "", g)
                    if key in ref:
                        r = ref[key]
                        dc = row["d_cos_vs_h1union"] - r["d_cos"]
                        dp = row["d_pds_vs_h1union"] - r["d_pds"]
                        out["ka5"].append(dict(fold=held, arm=arm_tag, g=g, d_cos_delta=dc,
                                               d_pds_delta=dp,
                                               ours_d_cos=row["d_cos_vs_h1union"],
                                               ref_d_cos=r["d_cos"],
                                               ours_d_pds=row["d_pds_vs_h1union"],
                                               ref_d_pds=r["d_pds"]))
                        ok = abs(dc) < 1e-12 and abs(dp) < 1e-12
                        print(f"[KA-5 {held:>14} {arm_tag} g={g:<5}] d(d_cos) {dc:+.3e} "
                              f"d(d_pds) {dp:+.3e}  {'PASS' if ok else 'FAIL'}", flush=True)
                        if not ok:
                            fail(f"KA-5 {held}@{g} {arm_tag}: 225's row did not reproduce")

                # concentration arrays at the BUILT dose, for ARM4/ARM6, cancer folds
                if g == BUILT_G and held in CANCER and arm_tag in ("ARM4", "ARM6"):
                    per_fold_conc.setdefault(held, {})[arm_tag] = dict(
                        targets=targets,
                        carry=(np.asarray(res.target_condition_count) > 0).tolist(),
                        d_pds=(np.asarray(e["_pds_per_target_a23"])
                               - np.asarray(scored_base["h1union"]
                                            ["_pds_per_target_a23"])).tolist(),
                        d_cos=(np.asarray(e["_cos_per_target_a23"])
                               - np.asarray(scored_base["h1union"]
                                            ["_cos_per_target_a23"])).tolist(),
                        d_pds_vs_cd4=(np.asarray(e["_pds_per_target_a23"])
                                      - np.asarray(scored_base["cd4f034"]
                                                   ["_pds_per_target_a23"])).tolist(),
                        d_cos_vs_cd4=(np.asarray(e["_cos_per_target_a23"])
                                      - np.asarray(scored_base["cd4f034"]
                                                   ["_cos_per_target_a23"])).tolist(),
                    )
                del m
        del whitened, unions, base, scored_base
        print(f"[{held}] done, {time.time() - t0:.0f}s, peak {peak_gib():.2f} GiB", flush=True)

    out["conc"] = per_fold_conc
    out["peak_rss_gib"] = peak_gib()
    out["seconds"] = round(time.time() - t0, 1)
    (WORK / "price229.json").write_text(json.dumps(out, indent=1))
    print(f"\n{len(out['failures'])} failures; peak {peak_gib():.2f} GiB; "
          f"{out['seconds']}s", flush=True)
    if out["failures"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
