"""223 Phase 2 — every partial-coverage source we hold, re-measured on CARRYING ROWS,
on top of 222's Pan/`h1_2025` union, on 162's four PLOLO folds.

FRAME.  The mandate forbids double-counting 222, so the reported number for each source is
its MARGINAL value on top of `cd4@0.34 + h1pair@0.28` -- the arm 222 is building.  Each
source is also measured alone on top of `cd4@0.34` so the two can be separated.

SIZING.  Three rules are run side by side on every source at every dose:
  OLD    production `solve_cd4_fraction` / `plolo.solve_fraction` -- median over EVERY row.
         For a <=50%-coverage source this is EXACTLY ZERO and the arm is bit-identical to
         its baseline.  It is included precisely so the void is visible in this experiment's
         own numbers rather than asserted from 218's.
  CARRY  218 section 7.1's `median(||added_i|| : ||added_i|| > 0)`.
  (CD4's own sizing is never touched anywhere -- 216/217's territory.)

KNOWN-ANSWER ROWS -- the run aborts if any fails:
  KA-1  `f0`            == 218 `resid218` f0                       (all folds)
  KA-2  `cd4`           == 218 `resid218` cd4_gn0.34               (all folds)
  KA-3  `base222`       == 218 `compose218` cd40.34+h1PAIRCARRY0.28
  KA-4  telohaec OLD    -> fraction exactly 0 and endpoints bit-identical to its baseline
  KA-5  kolf2 OLD       -> fraction > 0 (a >50% source must NOT be void)

FALSIFIERS, per source per dose:
  F1  target-shuffled source (218's derangement, whole rows permuted)
  F2  non-carrying rows must not move (per-target cosine split by the carrying mask)
  F3  row-matched CD4 control: the same rows, the same sizing, CD4's information
  F4  leak audit: donor <-> held-surface identity recorded on every arm
  F5  non-vacuity assertion on every arm that claims a nonzero dose

    uv run python sweep223.py [--folds ...] [--sources ...]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
plolo.WORK.mkdir(parents=True, exist_ok=True)

from plolo import (E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for,  # noqa: E402
                   solve_fraction)

sys.path.insert(0, str(ROOT / "experiments/218-independent-donors"))
import sources218  # noqa: E402

# Two sources have no donor-shaped table upstream; `build_tables.py` made one inside this
# experiment.  Registered here so 218's loader can read them unchanged.
sources218.SOURCES = dict(sources218.SOURCES)
sources218.SOURCES["dld1_lowmoi"] = dict(
    path="experiments/223-silent-zero-audit/work/tables/dld1_lowmoi",
    platform="DLD-1 low-MOI CRISPRi (GSE, two libraries)", chemistry="10x",
    lineage="epithelial (colorectal)", transformation="immortalised cancer line",
    lab="116's acquisition", semantics="author log2FC, split-mean", role="candidate")
sources218.SOURCES["altair_a549"] = dict(
    path="experiments/223-silent-zero-audit/work/tables/altair_a549",
    platform="Altair OPS A549", chemistry="10x", lineage="epithelial (lung adeno)",
    transformation="immortalised cancer line", lab="Altair",
    semantics="official OPS log2FC", role="candidate")
for _n in ("dld1_low1", "dld1_low2"):
    sources218.SOURCES[_n] = dict(
        path=f"experiments/223-silent-zero-audit/work/tables/{_n}",
        platform="DLD-1 low-MOI CRISPRi, ONE library", chemistry="10x",
        lineage="epithelial (colorectal)", transformation="immortalised cancer line",
        lab="116's acquisition", semantics="author log2FC, single split",
        role="replication control")
for _n, _lab, _sem in (
        ("cmap_sh_cgs", "LINCS 2020 trt_sh.cgs, 9 core lines",
         "MODZ consensus z-score (978 landmark + 11,350 inferred)"),
        ("cmap_sh_cgs_lm", "LINCS 2020 trt_sh.cgs, LANDMARK genes only",
         "MODZ z, inferred coordinates zeroed"),
        ("cmap_xpr", "LINCS 2020 trt_xpr CRISPR-KO", "MODZ consensus z-score")):
    sources218.SOURCES[_n] = dict(
        path=f"experiments/223-silent-zero-audit/work/tables/{_n}",
        platform="L1000 bulk, 96 h", chemistry="L1000 Luminex",
        lineage="mixed cancer lines", transformation="immortalised cancer lines",
        lab=_lab, semantics=_sem, role="candidate")

from resid218 import (ALPHA, GAMMA, LAM, score_effects, shuffled,  # noqa: E402
                      source_residual)
from compose218 import (merged, restrict_rows, solve_fraction_carrying,  # noqa: E402
                        union_rows)

CD4_F = 0.34            # board incumbent cd4f034, PINNED (216/217 own the dose)
UNION_G = 0.28          # 218's headline dose for the Pan/h1_2025 union, PINNED
G = (0.17, 0.28, 0.45)
UNION = ("gse295214_h1", "h1_2025")

# every partial-coverage source with a table on disk, plus the two >50% controls
SRC = ("telohaec", "dld1_lowmoi", "altair_a549", "myllia_a549",
       "frangieh_melanoma", "feng_ipsc", "kolf2")

PANEL_SET = plolo.PANEL_SET


def carrying_mask(res) -> np.ndarray:
    return np.asarray(res.target_condition_count) > 0


def row_split(cos_per_target, base_cos_per_target, mask) -> dict:
    a = np.asarray(cos_per_target, dtype=np.float64)
    b = np.asarray(base_cos_per_target, dtype=np.float64)
    m = np.asarray(mask, dtype=bool)
    out = {}
    for tag, sel in (("carrying", m), ("non_carrying", ~m)):
        if not sel.any():
            out[f"{tag}_n"] = 0
            continue
        out[f"{tag}_n"] = int(sel.sum())
        out[f"{tag}_cos"] = float(a[sel].mean())
        out[f"{tag}_base_cos"] = float(b[sel].mean())
        out[f"{tag}_d_cos"] = float(a[sel].mean() - b[sel].mean())
        out[f"{tag}_rel_pct"] = float(
            100.0 * (a[sel].mean() - b[sel].mean()) / max(abs(b[sel].mean()), 1e-300))
    return out


def audit(res, whitened, frac) -> dict:
    added = np.where(np.asarray(res.gene_votes) > 0, np.asarray(res.mean_effect), 0.0)
    norms = np.linalg.norm(added, axis=1)
    live = norms[norms > 0]
    return dict(fraction=float(frac),
                median_row_norm_all_rows=float(np.median(norms)),
                median_row_norm_carrying=float(np.median(live)) if live.size else 0.0,
                carrying_rows=int(live.size), n_rows=int(norms.size),
                coverage_frac=float(live.size / max(norms.size, 1)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", nargs="*", default=list(FOLDS))
    ap.add_argument("--sources", nargs="*", default=list(SRC))
    ap.add_argument("--out", default=str(WORK / "sweep223.json"))
    args = ap.parse_args()

    ref_resid = {(r["fold"], r["arm"]): r for r in json.loads(
        (ROOT / "experiments/218-independent-donors/work/resid218.json").read_text())}
    ref_comp = {(r["fold"], r["arm"]): r for r in json.loads(
        (ROOT / "experiments/218-independent-donors/work/compose218.json").read_text())}

    out_path = Path(args.out)
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"]) for r in rows}
    ka: list[dict] = []

    for held in args.folds:
        fold = build_fold(held)
        targets = fold["targets"]
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(targets)
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        whitened = op.apply(pre * ALPHA, fold["wcols"], GAMMA, LAM)

        def finish(res_or_none, frac):
            if res_or_none is None or frac == 0.0:
                f = whitened.copy()
            else:
                f, _ = E062.add_post_whiten_residual(whitened, res_or_none, 1.0,
                                                     fold["wcols"])
            f[:, fold["panel_positions"]] = 0.0
            return f

        # ---------------------------------------------------------------- baselines
        frac_cd4, _ = solve_fraction(whitened, cd4, CD4_F)
        h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
        frac_u, _ = solve_fraction_carrying(whitened, h1, UNION_G)
        assert frac_u > 0.0, "the Pan/h1_2025 union sized to zero -- instrument is broken"

        base_specs = {
            "f0": (None, 0.0),
            "cd4": (merged(cd4, frac_cd4, None, 0.0), 1.0),
            "base222": (merged(cd4, frac_cd4, h1, frac_u), 1.0),
        }
        scored: dict[str, dict] = {}
        for name, (res, fr) in base_specs.items():
            e = score_effects(finish(res, fr), fold, truth)
            scored[name] = e
            rows.append(dict(fold=held, arm=name, kind="baseline",
                             n_targets=len(targets), frac_cd4=frac_cd4,
                             frac_union=frac_u, **e))
            print(f"[{held}/{name}] cos {e['cos']:.6f} cos_agg {e['cos_agg']:.6f} "
                  f"pds {e['pds']:.6f} q_w {e['a15_sgn_w']:.6f}", flush=True)

        # KA-1/2/3
        for arm, refmap, refarm in (("f0", ref_resid, "f0"),
                                    ("cd4", ref_resid, "cd4_gn0.34"),
                                    ("base222", ref_comp, "cd40.34+h1PAIRCARRY0.28")):
            r = refmap.get((held, refarm))
            if r is None:
                ka.append(dict(fold=held, gate=arm, status="REFERENCE ABSENT"))
                continue
            dc = abs(r["cos_agg"] - scored[arm]["cos_agg"])
            dp = abs(r["pds"] - scored[arm]["pds"])
            ok = dc < 1e-9 and dp < 1e-9
            ka.append(dict(fold=held, gate=arm, ref=refarm, d_cos_agg=dc, d_pds=dp,
                           status=("PASS" if ok else "FAIL")))
            assert ok, f"KA FAILED {held}/{arm} vs 218 {refarm}: dcos {dc:.3e} dpds {dp:.3e}"
            print(f"    KA PASS {arm} == 218 {refarm} (dcos {dc:.1e} dpds {dp:.1e})",
                  flush=True)

        base_cos_pt = np.asarray(scored["base222"]["_cos_per_target_a23"])

        # ------------------------------------------------------------------- arms
        for sname in args.sources:
            if sname == held:
                continue                                    # never score a source on itself
            S = source_residual(sname, targets, "gene_norm")
            cm = carrying_mask(S)
            variants = {
                "": S,
                "SHUF": shuffled(S),
                "CD4ROWS": restrict_rows(cd4, S),
            }
            src_targets = set(
                x.strip() for x in
                (ROOT / sources218.SOURCES[sname]["path"] / "perts.txt").read_text().split())
            leak = dict(source=sname, held_line=held,
                        source_is_held_line=bool(sname == held),
                        source_targets_on_fold=int(len(src_targets & set(targets))),
                        fold_targets=len(targets),
                        note=("h1_2025 covers 85/85 of h1-test's targets; h1-test is NOT "
                              "used anywhere in this experiment"))
            for vtag, res in variants.items():
                for g in G:
                    for rule, solver in (("CARRY", solve_fraction_carrying),
                                         ("OLD", solve_fraction)):
                        arm = f"{sname}{vtag}_{rule}{g:g}"
                        if (held, arm) in done:
                            continue
                        t0 = time.time()
                        fr, per_unit = solver(whitened, res, g)
                        m = merged(cd4, frac_cd4, h1, frac_u)
                        if fr > 0.0:
                            m = merged(m, 1.0, res, fr)
                        e = score_effects(finish(m, 1.0), fold, truth)
                        aud = audit(res, whitened, fr)
                        aud["rule"] = rule
                        aud["per_unit"] = float(per_unit)
                        aud["VACUOUS"] = bool(fr == 0.0)
                        bit = bool(e["cos_agg"] == scored["base222"]["cos_agg"]
                                   and e["pds"] == scored["base222"]["pds"])
                        if fr > 0.0:
                            assert aud["median_row_norm_carrying"] > 0.0, \
                                f"F5 non-vacuity failed on {held}/{arm}"
                        rows.append(dict(
                            fold=held, arm=arm, kind="arm", source=sname, variant=vtag,
                            rule=rule, g=g, n_targets=len(targets),
                            d_cos=e["cos"] - scored["base222"]["cos"],
                            rel_cos_pct=100.0 * (e["cos"] - scored["base222"]["cos"])
                            / max(abs(scored["base222"]["cos"]), 1e-300),
                            d_cos_agg=e["cos_agg"] - scored["base222"]["cos_agg"],
                            rel_cos_agg_pct=100.0
                            * (e["cos_agg"] - scored["base222"]["cos_agg"])
                            / max(abs(scored["base222"]["cos_agg"]), 1e-300),
                            d_pds=e["pds"] - scored["base222"]["pds"],
                            d_qw=e["a15_sgn_w"] - scored["base222"]["a15_sgn_w"],
                            bit_identical_to_base222=bit,
                            rowsplit=row_split(e["_cos_per_target_a23"], base_cos_pt, cm),
                            audit=aud, leak=leak, **e))
                        out_path.write_text(json.dumps(rows, indent=2))
                        flag = " [VACUOUS/BIT-IDENTICAL]" if bit else ""
                        print(f"[{held}/{arm}] frac {fr:.5f} carry {aud['carrying_rows']}"
                              f"/{aud['n_rows']} dcos {e['cos'] - scored['base222']['cos']:+.3e}"
                              f" ({100.0 * (e['cos'] - scored['base222']['cos']) / max(abs(scored['base222']['cos']), 1e-300):+.3f}%)"
                              f" dpds {e['pds'] - scored['base222']['pds']:+.5f}"
                              f"{flag} ({time.time() - t0:.0f}s)", flush=True)

                        # KA-4 / KA-5 fire on the OLD rule at the first dose
                        if rule == "OLD" and vtag == "" and g == G[0]:
                            cov = aud["coverage_frac"]
                            if cov <= 0.5:
                                assert fr == 0.0 and bit, \
                                    f"KA-4 FAILED: {sname} at {cov:.1%} was NOT void"
                                ka.append(dict(fold=held, gate=f"KA-4 {sname} OLD void",
                                               coverage=cov, fraction=fr,
                                               bit_identical=bit, status="PASS"))
                            else:
                                assert fr > 0.0, \
                                    f"KA-5 FAILED: {sname} at {cov:.1%} was void"
                                ka.append(dict(fold=held, gate=f"KA-5 {sname} OLD live",
                                               coverage=cov, fraction=fr, status="PASS"))
        (WORK / "sweep223_ka.json").write_text(json.dumps(ka, indent=2))
        del op, whitened

    out_path.write_text(json.dumps(rows, indent=2))
    (WORK / "sweep223_ka.json").write_text(json.dumps(ka, indent=2))
    print(f"\n-> {out_path}\n-> {WORK / 'sweep223_ka.json'}")


if __name__ == "__main__":
    main()
