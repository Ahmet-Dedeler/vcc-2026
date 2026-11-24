#!/usr/bin/env python3
"""618: run323's faithful production-key pipeline (stage / operator / compose), parameterised by PANEL.

  ka : h1-test's 85 targets with 323's own inputs (273's cd4_compact_h1). KNOWN ANSWER: effects.npz cap15 must equal
       323/work/effects.npz cap15 bit for bit (ka618.py), with 323's own stage KAs (historical H1 pre-whitening, 319's
       A/B/C hashes) still asserted.
  tv : h1-trainval's 126 targets (disjoint from h1-test's). CD4 is a fresh 047 extract_cd4.py run on this panel, and the
       operator is refit on this panel's bundle (pruned of its own targets).
Usage: uv run python run618.py <ka|tv> <bundle|cd4|freeze|stage|operator|compose>
Differences from run323 (all mechanical): WORK = work/<PANEL>; targets / CD4 / operator path by PANEL; a `bundle` and
a `cd4` step before `freeze`; 323's "77 complete CD4 targets" and historical-pre KAs apply to `ka` only.
Every output stays in this experiment. The production builder is imported, not copied or edited.
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PANEL = sys.argv.pop(1)
assert PANEL in ("ka", "tv"), PANEL
WORK = HERE / "work" / PANEL
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
CD4 = (ROOT / "experiments/273-flex-reference-instrument/work/cd4/cd4_compact_h1" if PANEL == "ka"
       else WORK / "cd4_compact")
CD4_SOURCE = ROOT / "experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad"
TV = ROOT / "data/harness/h1-trainval.h5ad"
PRE = ROOT / "experiments/170-gamma-mse/work/pre.npz"
REAL = ROOT / "data/harness/h1-test.h5ad"
OLD_OPERATOR = ROOT / "data/donors/whitener_union_pruned.npz"
OPERATOR = WORK / ("whitener_h1_pruned.npz" if PANEL == "ka" else "whitener_tv_pruned.npz")
ARMS = {
    "uniform": (None, "plain", None),
    "plain15": (-1.5, "plain", None),
    "plain30": (-3.0, "plain", None),
    "cap15": (-1.5, "cap", None),
    "parity15": (-1.5, "plain", "blend"),
}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    tmp.replace(path)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@contextlib.contextmanager
def model_env(arm="uniform"):
    old = {k: v for k, v in os.environ.items() if k.startswith("VCC_")}
    try:
        for k in old:
            del os.environ[k]
        os.environ.update(VCC_ALPHA="0.982", VCC_CD4_TARGET_F="0.34",
                          VCC_CD4_DEFLATE_K="5", VCC_WHITEN_GAMMA="0.65",
                          VCC_WHITEN_LAM="0")
        exp, mode, pm = ARMS[arm]
        if exp is not None:
            os.environ.update(VCC_NORMREL_EXP=str(exp), VCC_NORMREL_MODE=mode)
        if pm is not None:
            os.environ["VCC_NORMREL_PM"] = pm
        yield
    finally:
        for k in list(os.environ):
            if k.startswith("VCC_"):
                del os.environ[k]
        os.environ.update(old)


def builder(arm="uniform"):
    with model_env(arm):
        b = load_module(f"builder323_{arm}", ROOT / "scripts/make_submission_pm_cd4.py")
    assert b.H1_UNION_G is None and b.MASS_FIX is None
    return b


def axis_columns(genes, donor_genes, panel):
    union = set().union(*(set(g) for g in donor_genes))
    return np.asarray([i for i, g in enumerate(genes)
                       if g in union and g not in set(panel)], dtype=np.int64)


def axis_digest(genes):
    return hashlib.sha256("\n".join(map(str, genes)).encode()).hexdigest()


def verify_axis(genes, columns, panel, expected_digest):
    columns = np.asarray(columns, dtype=np.int64)
    if len(set(columns.tolist())) != len(columns):
        raise ValueError("duplicate operator coordinate")
    if np.any(columns < 0) or np.any(columns >= len(genes)):
        raise ValueError("operator coordinate outside gene axis")
    symbols = np.asarray(genes)[columns]
    if set(symbols) & set(panel):
        raise ValueError("excluded target coordinate entered whitening")
    if axis_digest(symbols) != expected_digest:
        raise ValueError("operator axis digest mismatch")


def map_axis(effects, source, target):
    if len(set(source)) != len(source) or len(set(target)) != len(target):
        raise ValueError("duplicate gene symbol in effect mapping")
    pos = {g: i for i, g in enumerate(source)}
    out = np.zeros((len(effects), len(target)), dtype=np.float64)
    for j, g in enumerate(target):
        if g in pos:
            out[:, j] = effects[:, pos[g]]
    return out


def derangement(n, seed=20260909):
    if n < 2:
        raise ValueError("derangement needs at least two rows")
    rng = np.random.default_rng(seed)
    while True:
        p = rng.permutation(n)
        if np.all(p != np.arange(n)):
            return p


def axes():
    # Only strings and index arrays from the staging archive, never its effects.
    with np.load(PRE, allow_pickle=False) as z:
        genes = z["genes"].astype(str)
        targets = z["h1_targets"].astype(str)
        h1_genes = z["h1_genes"].astype(str)
        old_cols = z["operator_columns"].astype(np.int64)
    if PANEL == "tv":
        targets = np.array([l.strip() for l in (WORK / "targets.txt").read_text().splitlines() if l.strip()])
    dgenes = [(ROOT / f"data/donors/{n}/genes.txt").read_text().splitlines()
              for n in DONORS]
    cols = axis_columns(genes, dgenes, targets)
    return genes, targets, h1_genes, old_cols, cols


def bundle():
    """618: the panel's bundle (gene axis = PRE's official 18,533; targets) BEFORE the CD4 extraction and the freeze."""
    if PANEL == "tv":
        import anndata as ad
        a = ad.read_h5ad(TV, backed="r")
        lab = sorted(set(a.obs["target_gene"].astype(str)) - {"non-targeting"})
        a.file.close()
        (WORK / "targets.txt").write_text("\n".join(lab) + "\n")
    genes, targets, _, _, _ = axes()
    b = WORK / "bundle"
    b.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"gene_name": genes}).to_csv(b / "gene_names.csv", index=False)
    pd.DataFrame({"target_gene": targets}).to_csv(b / "pert_counts.csv", index=False)
    print(f"bundle {PANEL}: {len(genes)} genes, {len(targets)} targets", flush=True)


def cd4():
    """618: 047's extractor (the D/E/F cd4-extract stage) on this panel. For `ka` it writes a CHECK copy that ka618.py
    compares with 273's cd4_compact_h1 array for array; for `tv` it writes the compact the pipeline uses."""
    out = WORK / ("cd4_ka_check" if PANEL == "ka" else "cd4_compact")
    subprocess.run([sys.executable, str(ROOT / "experiments/047-cd4-direct/extract_cd4.py"), "--input", str(CD4_SOURCE),
                    "--output", str(out), "--panel-csv", str(WORK / "bundle/pert_counts.csv"),
                    "--genes-csv", str(ROOT / "data/challenge/gene_names.csv"), "--overwrite",
                    "--min-targets", "0", "--min-union-targets", "0"], check=True)


def freeze():
    manifest_path = WORK / "inputs.json"
    if manifest_path.exists():
        raise FileExistsError("inputs already frozen; do not overwrite")
    files = [
        Path(__file__), HERE / "PLAN.md", WORK / "bundle/gene_names.csv", WORK / "bundle/pert_counts.csv",
        ROOT / "scripts/make_submission_pm_cd4.py",
        ROOT / "scripts/build_union_whitener.py",
        ROOT / "scripts/build_submission.py",
        ROOT / "scripts/make_submission.py",
        ROOT / "src/vcc2026/multi_donor.py",
        ROOT / "src/vcc2026/donor_transfer.py",
        ROOT / "src/vcc2026/whiten.py",
        ROOT / "src/vcc2026/rc_emit.py",
        ROOT / "experiments/061-perturbmatch-hct-count/build_candidate.py",
        ROOT / "experiments/062-cd4-post-whiten/run.py",
        ROOT / "experiments/170-gamma-mse/stage.py", PRE, REAL, *([TV] if PANEL == "tv" else []),
        OLD_OPERATOR, OLD_OPERATOR.with_suffix(".json"),
        ROOT / "data/donors/k562_essential/perts.txt",
        ROOT / "data/challenge/gene_names.csv",
        ROOT / "data/challenge/pert_counts.csv",
        ROOT / "submissions/319-normrel15-cap.h5ad.meta.json",
        ROOT / "experiments/278-depth-inventory/work/depth278.csv",
        ROOT / "experiments/285-reliability-benchmark/work/bench285.csv",
    ]
    for folder in [*(ROOT / f"data/donors/{n}" for n in DONORS),
                   ROOT / "data/external/perturbmatch",
                   *(CD4 / f"cd4_{c}" for c in ("rest", "stim8hr", "stim48hr"))]:
        files.extend(folder / n for n in ("delta.npy", "b_ctrl.npy", "genes.txt",
                                         "perts.txt", "meta.json") if (folder / n).exists())
    hashes = {}
    for f in files:
        hashes[str(f.relative_to(ROOT))] = sha(f)
        print(f"pin {f.relative_to(ROOT)}", flush=True)
    genes, targets, h1_genes, old_cols, cols = axes()
    info = {
        "classification": "retrospective development validation; not confirmation",
        "inputs_sha256": hashes,
        "panel": PANEL, "training_donors": list(DONORS), "h1_union": "OFF",
        "h1_targets": targets.tolist(), "h1_genes": len(h1_genes),
        "historical_axis_columns": len(old_cols),
        "historical_axis_h1_target_overlap": sorted(set(genes[old_cols]) & set(targets)),
        "required_axis_columns": len(cols),
        "required_axis_sha256": axis_digest(genes[cols]),
        "required_panel_sha256": axis_digest(targets),
    }
    write_json(manifest_path, info)
    print(json.dumps({k: v for k, v in info.items() if k not in
                      ("inputs_sha256", "h1_targets")}, indent=2), flush=True)


def verify_inputs():
    frozen = json.loads((WORK / "inputs.json").read_text())
    for name, expected in frozen["inputs_sha256"].items():
        if sha(ROOT / name) != expected:
            raise ValueError(f"frozen input drift: {name}")
    # A commit is a prerequisite, not a claim that all later code is unchanged.
    commit = subprocess.check_output(
        ["git", "log", "-1", "--format=%H", "--", str(WORK / "inputs.json")],
        cwd=ROOT, text=True).strip()
    if not commit:
        raise ValueError("commit the preregistration and inputs before measurements")
    return frozen


def stage():
    frozen = verify_inputs()
    genes, targets, _, _, cols = axes()
    b = builder()
    with model_env():
        lean = load_module("lean323", ROOT / "experiments/170-gamma-mse/stage.py")
    pm = b.EXP061.prepare_perturbmatch(
        ROOT / "data/external/perturbmatch", tuple(genes),
        chunk_rows=128, max_rss_gib=14)
    _, cd4 = b.prepare_frozen_cd4(
        CD4, tuple(targets), tuple(genes), tuple(targets),
        chunk_rows=128, max_rss_gib=14)
    complete = int(np.count_nonzero(cd4.target_condition_count == 3))
    if PANEL == "ka" and complete != 77:
        raise AssertionError(f"CD4 complete coverage {complete}, expected 77")
    print(f"CD4 complete coverage {complete} of {len(targets)} targets", flush=True)
    with lean.lean_norm():
        donors = {n: lean.load_donor_lean(b, ROOT / f"data/donors/{n}", genes)
                  for n in DONORS}
        model = b.MultiDonorModel(donors, genes, gene_norm=False, eb=False)
    abc = pd.read_csv(ROOT / "data/challenge/pert_counts.csv").target_gene.tolist()
    lean.trim_model_rows(model, list(targets) + abc)
    del donors
    weights = b.donor_weights(model)
    staged, audits = {}, {}
    for arm in ARMS:
        m = builder(arm)
        handled = set()
        pre = m.build_v9_pre_whiten_effects(
            model, targets, alpha=m.ALPHA, beta=m.BETA, weights=weights,
            pm=pm, operator_columns=cols, pm_handled_out=handled)
        cand, aud = m.apply_perturbmatch_k562_replacement(
            pre, model, pm, targets, alpha=m.ALPHA, weights=weights,
            already_blended=handled)
        staged[arm] = cand
        audits[arm] = {"baseline_pre_sha256": m.array_sha256(pre),
                       "candidate_pre_sha256": m.array_sha256(cand),
                       "pm": aud}
        print(f"staged {arm}", flush=True)
    diff = None
    if PANEL == "ka":
        with np.load(PRE, allow_pickle=False) as z:
            diff = float(np.max(np.abs(staged["uniform"] - z["cand_h1"])))
        if diff > 1e-10:
            raise AssertionError(f"historical H1 pre-whitening mismatch {diff}")
    m = builder("cap15")
    abc_cols = m.production_operator_columns(model, abc)
    pre = m.build_v9_pre_whiten_effects(model, abc, alpha=m.ALPHA, beta=m.BETA,
                                      weights=weights)
    cand, _ = m.apply_perturbmatch_k562_replacement(
        pre, model, pm, abc, alpha=m.ALPHA, weights=weights)
    expected = json.loads((ROOT / "submissions/319-normrel15-cap.h5ad.meta.json").read_text())
    for key, value in (("baseline_pre_sha256", pre), ("candidate_pre_sha256", cand)):
        if m.array_sha256(value) != expected["effects"][key]:
            raise AssertionError(f"319 known-answer hash mismatch: {key}")
    if not np.array_equal(cols, m.production_operator_columns(model, targets)):
        raise AssertionError("metadata and production operator axes disagree")
    np.savez(WORK / "stage.npz", genes=genes, targets=targets,
             columns=cols, **staged, cd4_mean=cd4.mean_effect,
             cd4_votes=cd4.gene_votes, cd4_counts=cd4.target_condition_count)
    write_json(WORK / "stage.json", {
        "gates_pass": True, "h1_pre_max_abs": diff,
        "abc_cap15_hashes_match_319": True, "cd4_complete_targets": complete,
        "held_context_sources_loaded": [], "arms": audits,
        "axis_sha256": axis_digest(genes[cols]),
        "stage_sha256": sha(WORK / "stage.npz"),
        "historical_axis_h1_target_overlap": frozen["historical_axis_h1_target_overlap"],
        "abc_columns": len(abc_cols),
    })


def operator():
    verify_inputs()
    if OPERATOR.exists():
        raise FileExistsError(OPERATOR)
    subprocess.run([
        sys.executable, str(ROOT / "scripts/build_union_whitener.py"),
        "--bundle", str(WORK / "bundle"), "--donor-root", str(ROOT / "data/donors"),
        "--donors", ",".join(DONORS), "--basis-prune", "essential",
        "--essential-list", str(ROOT / "data/donors/k562_essential/perts.txt"),
        "--out", str(OPERATOR), "--scratch", str(WORK / "operator-scratch"),
        "--chunk-rows", "512", "--max-rss-gib", "18",
    ], check=True)


def compose():
    verify_inputs()
    stage_info = json.loads((WORK / "stage.json").read_text())
    if not stage_info["gates_pass"] or sha(WORK / "stage.npz") != stage_info["stage_sha256"]:
        raise ValueError("stage failed or changed")
    z = np.load(WORK / "stage.npz", allow_pickle=False)
    genes, targets, cols = z["genes"], z["targets"], z["columns"]
    b = builder()
    model = SimpleNamespace(
        names=list(DONORS), target_genes=genes,
        tcols={n: np.flatnonzero(np.isin(genes, (
            ROOT / f"data/donors/{n}/genes.txt").read_text().splitlines()))
               for n in DONORS})
    old_rejected = False
    try:
        b.load_production_operator(model, targets, OLD_OPERATOR,
                                   OLD_OPERATOR.with_suffix(".json"),
                                   expected_sha256=sha(OLD_OPERATOR))
    except ValueError as exc:
        old_rejected = True
        old_error = str(exc)
    if not old_rejected:
        raise AssertionError("old A/B/C operator accepted on H1")
    op = b.load_production_operator(model, targets, OPERATOR,
                                   OPERATOR.with_suffix(".json"),
                                   expected_sha256=sha(OPERATOR))
    verify_axis(genes, cols, targets, op.axis_sha256)
    if not np.array_equal(cols, op.columns):
        raise AssertionError("staged columns differ from loaded operator")
    cd4 = b.EXP062.CD4Prepared(mode="gene_norm", mean_effect=z["cd4_mean"],
                              gene_votes=z["cd4_votes"],
                              target_condition_count=z["cd4_counts"],
                              condition_details=[])
    panel_positions = np.flatnonzero(np.isin(genes, targets))
    # A target-coordinate-only row must not enter any whitened/scored coordinate.
    probe = np.zeros((2, len(genes)))
    probe[:, panel_positions] = np.arange(1, len(panel_positions) + 1)
    got, _ = b.apply_v9_whitening(probe, op.whitener, cols)
    scored = np.flatnonzero(~np.isin(genes, targets))
    if np.any(got[:, scored] != 0):
        raise AssertionError("excluded coordinates spread into scored coordinates")
    effects, audits = {}, {}
    for arm in ARMS:
        m = builder(arm)
        w, wa = m.apply_v9_whitening(z[arm], op.whitener, cols)
        final, ca, ua = m.apply_frozen_cd4_residual(
            w, cd4, cols, panel_positions)
        if ua is not None:
            raise AssertionError("H1 union unexpectedly active")
        effects[arm] = final
        audits[arm] = {"whitening": wa, "cd4": ca, "h1_union": "OFF",
                       "final_sha256": m.array_sha256(final)}
    permutation = derangement(len(targets))
    effects["shuffled_cap15"] = effects["cap15"][permutation]
    np.savez(WORK / "effects.npz", genes=genes, targets=targets,
             **effects, permutation=permutation)
    write_json(WORK / "compose.json", {
        "gates_pass": True, "old_operator_rejected": old_error,
        "new_operator_sha256": sha(OPERATOR), "new_axis_sha256": op.axis_sha256,
        "target_only_probe_scored_max_abs": 0.0,
        "arms": audits, "effects_sha256": sha(WORK / "effects.npz"),
    })


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", choices=["bundle", "cd4", "freeze", "stage", "operator", "compose"])
    args = ap.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    globals()[args.mode]()


if __name__ == "__main__":
    main()
