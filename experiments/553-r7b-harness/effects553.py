#!/usr/bin/env python3
"""553 — the carriers' donor effects for any target panel by 323's construction (run323.py stage -> operator -> compose),
cap15 arm only (VCC_NORMREL_EXP=-1.5, VCC_NORMREL_MODE=cap, under 323's model_env: ALPHA 0.982, CD4_TARGET_F 0.34,
CD4_DEFLATE_K 5, WHITEN_GAMMA 0.65, WHITEN_LAM 0). The production builder and 170's lean loaders are imported, never
copied; run323's helpers (builder, model_env, axis_columns, axis_digest, verify_axis) are imported read-only.

  stage   : the three production donors (k562_gwps, orion_hct116, orion_hek293t) -> MultiDonorModel (lean load), rows trimmed
            to the panel + A/B/C's targets (323's rule; a pure memory move), PerturbMatch (061) K562 replacement, the panel's
            frozen CD4 banks (extract_cd4.py's compact extraction for THIS panel), build_v9_pre_whiten_effects with the
            operator columns = union(donor coverage) minus the panel -> stage553.npz
  operator: scripts/build_union_whitener.py exactly as run323.operator() calls it (--basis-prune essential, the essential
            list, seed 0, --chunk-rows 512, --max-rss-gib 18), on a bundle whose pert_counts.csv is THIS panel, so the
            panel's genes are off the axis (DEF-RUNBOOK section 0.1)
  compose : load_production_operator (axis digest verified against the panel), apply_v9_whitening, the frozen CD4 residual
            (k 5 deflation re-fit on the panel's own CD4 rows) -> effects553.npz (cap15) + compose553.json (final sha256)

KA-E (h1): `stage` + `compose` on h1-test's 85 targets with 323's operator and 273's H1 CD4 banks must reproduce 323's
cap15 final sha256 (compose.json 565d57ab...) and 323's staged cap15 pre-whitening matrix bit for bit.

usage: uv run python experiments/553-r7b-harness/effects553.py stage|compose <h1|r7b> [--ka]
       uv run python experiments/553-r7b-harness/effects553.py operator r7b
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/323-current-family-validation"))
import lib553 as S  # noqa: E402
import run323 as R  # noqa: E402  (read-only: builder, model_env, load_module, axis_columns, axis_digest, verify_axis)

E323 = ROOT / "experiments/323-current-family-validation"
ARM = "cap15"


def panel_paths(surface: str) -> dict:
    """genes (official axis), targets (sorted panel), the CD4 bank dir, the operator, the output files."""
    genes = np.asarray(pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"], dtype=str)   # '<U', not object
    if surface == "h1":
        with np.load(ROOT / "experiments/170-gamma-mse/work/pre.npz", allow_pickle=False) as z:
            assert np.array_equal(z["genes"].astype(str), genes)
            targets = z["h1_targets"].astype(str)
        out = HERE / "work/h1ka"
        return dict(genes=genes, targets=targets, cd4=ROOT / "experiments/273-flex-reference-instrument/work/cd4/cd4_compact_h1",
                    operator=E323 / "work/whitener_h1_pruned.npz", stage=out / "stage553_h1.npz",
                    effects=out / "effects553_h1.npz", compose=out / "compose553_h1.json", bundle=E323 / "work/bundle")
    spec = S.surface_spec(surface)
    targets = np.asarray(S.surface_targets(spec), dtype=str)
    b_genes = np.asarray(pd.read_csv(Path(spec["bundle"]) / "gene_names.csv")["gene_name"], dtype=str)
    assert np.array_equal(b_genes, genes), "the surface's gene axis is not the official order"
    return dict(genes=genes, targets=targets, cd4=Path(spec["cd4"]), operator=Path(spec["operator"]),
                stage=Path(spec["op_columns"]), effects=Path(spec["effects"]), compose=Path(spec["compose"]),
                bundle=Path(spec["bundle"]))


def stage(surface: str, ka: bool):
    P = panel_paths(surface)
    genes, targets = P["genes"], P["targets"]
    S.assert_clean_sources([P["cd4"], *(ROOT / f"data/donors/{n}" for n in R.DONORS), ROOT / "data/external/perturbmatch"])
    t0 = time.time()
    dgenes = [(ROOT / f"data/donors/{n}/genes.txt").read_text().splitlines() for n in R.DONORS]
    cols = R.axis_columns(genes, dgenes, targets)
    b = R.builder()
    with R.model_env():
        lean = R.load_module("lean553", ROOT / "experiments/170-gamma-mse/stage.py")
    pm = b.EXP061.prepare_perturbmatch(ROOT / "data/external/perturbmatch", tuple(genes), chunk_rows=128, max_rss_gib=14)
    _, cd4 = b.prepare_frozen_cd4(P["cd4"], tuple(targets), tuple(genes), tuple(targets), chunk_rows=128, max_rss_gib=14)
    complete = int(np.count_nonzero(cd4.target_condition_count == 3))
    with lean.lean_norm():
        donors = {n: lean.load_donor_lean(b, ROOT / f"data/donors/{n}", genes) for n in R.DONORS}
        model = b.MultiDonorModel(donors, genes, gene_norm=False, eb=False)
    abc = pd.read_csv(ROOT / "data/challenge/pert_counts.csv").target_gene.tolist()
    lean.trim_model_rows(model, list(targets) + abc)
    del donors
    weights = b.donor_weights(model)
    m = R.builder(ARM)
    handled = set()
    pre = m.build_v9_pre_whiten_effects(model, targets, alpha=m.ALPHA, beta=m.BETA, weights=weights, pm=pm,
                                        operator_columns=cols, pm_handled_out=handled)
    cand, aud = m.apply_perturbmatch_k562_replacement(pre, model, pm, targets, alpha=m.ALPHA, weights=weights,
                                                      already_blended=handled)
    if not np.array_equal(cols, m.production_operator_columns(model, targets)):
        raise AssertionError("metadata and production operator axes disagree")
    coverage = {n: int(sum(t in model.pert_index[n] for t in targets)) for n in model.names}
    per_target_donors = {str(t): [n for n in model.names if t in model.pert_index[n]] for t in targets}
    rec = dict(surface=surface, n_targets=len(targets), arm=ARM, axis_columns=len(cols),
               axis_sha256=R.axis_digest(genes[cols]), panel_sha256=R.axis_digest(targets),
               cd4_dir=str(P["cd4"].relative_to(ROOT)), cd4_complete_targets=complete,
               cd4_condition_count_hist={int(k): int(v) for k, v in zip(*np.unique(cd4.target_condition_count, return_counts=True))},
               donor_coverage=coverage, donors_per_target_hist={int(k): int(v) for k, v in zip(
                   *np.unique([len(v) for v in per_target_donors.values()], return_counts=True))},
               per_target_donors=per_target_donors,
               baseline_pre_sha256=m.array_sha256(pre), candidate_pre_sha256=m.array_sha256(cand), pm=aud,
               seconds=round(time.time() - t0, 1), peak_rss_gb=S.rss_gb())
    if ka:
        st = json.loads((E323 / "work/stage.json").read_text())["arms"][ARM]
        with np.load(E323 / "work/stage.npz", allow_pickle=False) as z:
            same_stage = bool(np.array_equal(z[ARM], cand)) and bool(np.array_equal(z["columns"], cols))
            same_cd4 = all(bool(np.array_equal(z[k], v)) for k, v in (("cd4_mean", cd4.mean_effect),
                                                                        ("cd4_votes", cd4.gene_votes),
                                                                        ("cd4_counts", cd4.target_condition_count)))
        rec["ka"] = dict(baseline_pre_equal=rec["baseline_pre_sha256"] == st["baseline_pre_sha256"],
                         candidate_pre_equal=rec["candidate_pre_sha256"] == st["candidate_pre_sha256"],
                         stage_npz_cap15_and_columns_equal=same_stage, cd4_equal=same_cd4)
        print(f"KA-E stage: {rec['ka']}", flush=True)
        if not all(rec["ka"].values()):
            S.write_json(P["stage"].with_suffix(".json"), rec)
            raise SystemExit("KA-E stage FAILED")
    P["stage"].parent.mkdir(parents=True, exist_ok=True)
    np.savez(P["stage"], genes=genes, targets=targets, columns=cols, **{ARM: cand}, cd4_mean=cd4.mean_effect,
             cd4_votes=cd4.gene_votes, cd4_counts=cd4.target_condition_count)
    rec["stage_sha256"] = S.sha256(P["stage"])
    S.write_json(P["stage"].with_suffix(".json"), rec)
    print(json.dumps({k: v for k, v in rec.items() if k not in ("per_target_donors", "pm")}), flush=True)


def operator(surface: str):
    """build_union_whitener.py as run323.operator() runs it, on the surface's bundle (pert_counts.csv = the panel)."""
    assert surface != "h1", "the H1 operator is 323's (a46fd323...); it is not refit here"
    P = panel_paths(surface)
    if P["operator"].exists():
        raise FileExistsError(P["operator"])
    P["operator"].parent.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, str(ROOT / "scripts/build_union_whitener.py"),
           "--bundle", str(P["bundle"]), "--donor-root", str(ROOT / "data/donors"),
           "--donors", ",".join(R.DONORS), "--basis-prune", "essential",
           "--essential-list", str(ROOT / "data/donors/k562_essential/perts.txt"),
           "--out", str(P["operator"]), "--scratch", str(P["operator"].parent / "operator-scratch"),
           "--chunk-rows", "512", "--max-rss-gib", "18"]
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)
    meta = json.loads(P["operator"].with_suffix(".json").read_text())
    got = S.sha256(P["operator"])
    assert meta.get("sha256") == got, "operator metadata sha differs from the file"
    print(f"operator {P['operator'].relative_to(ROOT)} sha256 {got}; axis {meta['axis_sha256']} ({meta['n_cols']} cols)")


def compose(surface: str, ka: bool):
    P = panel_paths(surface)
    t0 = time.time()
    st = json.loads(P["stage"].with_suffix(".json").read_text())
    if S.sha256(P["stage"]) != st["stage_sha256"]:
        raise ValueError("stage changed")
    z = np.load(P["stage"], allow_pickle=False)
    genes, targets, cols = z["genes"], z["targets"], z["columns"]
    b = R.builder()
    model = SimpleNamespace(names=list(R.DONORS), target_genes=genes,
                            tcols={n: np.flatnonzero(np.isin(genes, (ROOT / f"data/donors/{n}/genes.txt").read_text().splitlines()))
                                   for n in R.DONORS})
    op_sha = S.sha256(P["operator"])
    op = b.load_production_operator(model, targets, P["operator"], P["operator"].with_suffix(".json"), expected_sha256=op_sha)
    R.verify_axis(genes, cols, targets, op.axis_sha256)
    if not np.array_equal(cols, op.columns):
        raise AssertionError("staged columns differ from the loaded operator")
    cd4 = b.EXP062.CD4Prepared(mode="gene_norm", mean_effect=z["cd4_mean"], gene_votes=z["cd4_votes"],
                               target_condition_count=z["cd4_counts"], condition_details=[])
    panel_positions = np.flatnonzero(np.isin(genes, targets))
    probe = np.zeros((2, len(genes)))
    probe[:, panel_positions] = np.arange(1, len(panel_positions) + 1)
    got, _ = b.apply_v9_whitening(probe, op.whitener, cols)
    scored = np.flatnonzero(~np.isin(genes, targets))
    if np.any(got[:, scored] != 0):
        raise AssertionError("excluded coordinates spread into scored coordinates")
    m = R.builder(ARM)
    w, wa = m.apply_v9_whitening(z[ARM], op.whitener, cols)
    final, ca, ua = m.apply_frozen_cd4_residual(w, cd4, cols, panel_positions)
    if ua is not None:
        raise AssertionError("H1 union unexpectedly active")
    audit = {"whitening": wa, "cd4": ca, "h1_union": "OFF", "final_sha256": m.array_sha256(final)}
    rec = dict(surface=surface, gates_pass=True, operator=str(P["operator"].relative_to(ROOT)), operator_sha256=op_sha,
               axis_sha256=op.axis_sha256, target_only_probe_scored_max_abs=0.0, arms={ARM: audit},
               row_norm_median=float(np.median(np.linalg.norm(final, axis=1))), seconds=round(time.time() - t0, 1),
               peak_rss_gb=S.rss_gb())
    if ka:
        want = json.loads((E323 / "work/compose.json").read_text())["arms"][ARM]["final_sha256"]
        with np.load(E323 / "work/effects.npz", allow_pickle=False) as ze:
            same = bool(np.array_equal(ze[ARM], final))
        rec["ka"] = dict(final_sha256_equal=audit["final_sha256"] == want, effects_npz_cap15_equal=same, want=want)
        print(f"KA-E compose: {rec['ka']}", flush=True)
        if not (rec["ka"]["final_sha256_equal"] and same):
            S.write_json(P["compose"], rec)
            raise SystemExit("KA-E compose FAILED")
    np.savez(P["effects"], genes=genes, targets=targets, **{ARM: final})
    rec["effects_sha256"] = S.sha256(P["effects"])
    S.write_json(P["compose"], rec)
    print(json.dumps({k: v for k, v in rec.items() if k != "arms"} | {"final_sha256": audit["final_sha256"]}), flush=True)


def main():
    step, surface = sys.argv[1], sys.argv[2]
    ka = "--ka" in sys.argv
    {"stage": lambda: stage(surface, ka), "operator": lambda: operator(surface),
     "compose": lambda: compose(surface, ka)}[step]()


if __name__ == "__main__":
    main()
