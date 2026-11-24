"""One pre-registered cap-dose comparison using the frozen 323 instrument."""

from __future__ import annotations

import argparse
import contextlib
import gc
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BASE = ROOT / "experiments/323-current-family-validation"
sys.path.insert(0, str(BASE))
import run323 as R
import evaluate323 as E

ARMS = {"cap15": (-1.5, "cap", None), "cap30": (-3.0, "cap", None)}


def inputs():
    files = [HERE / n for n in
             ("PREREGISTRATION.md", "run326.py", "report326.py", "test326.py")]
    files += [BASE / n for n in
              ("inputs.json", "evaluation-inputs.json", "anchor-inputs.json",
               "run323.py", "evaluate323.py", "report323.py")]
    files += [BASE / "work" / n for n in
              ("stage.npz", "stage.json", "effects.npz", "compose.json",
               "reference.npz", "reference.json", "reference_de.parquet",
               "h1-anchors.json", "whitener_h1_pruned.npz",
               "whitener_h1_pruned.json", "metrics-known.json")]
    for arm in ("cap15", "plain15", "plain30"):
        files += [BASE / "work" / f"{prefix}-{arm}.{extension}" for prefix, extension in
                  (("metrics", "json"), ("stats", "npz"), ("per-target", "csv"))]
    files += [ROOT / "submissions/319d3-cap30.h5ad.meta.json"]
    return files


def freeze():
    out = HERE / "inputs.json"
    if out.exists():
        raise FileExistsError(out)
    R.write_json(out, {"inputs_sha256": {str(p.relative_to(ROOT)): R.sha(p) for p in inputs()},
                       "arms": ARMS, "classification": "reused H1 development diagnostic"})


def verify():
    manifest = json.loads((HERE / "inputs.json").read_text())
    if set(manifest["inputs_sha256"]) != {str(p.relative_to(ROOT)) for p in inputs()}:
        raise ValueError("Input allowlist changed")
    for name, expected in manifest["inputs_sha256"].items():
        if R.sha(ROOT / name) != expected:
            raise ValueError(f"326 frozen input drift: {name}")
    tracked = [str(p.relative_to(ROOT)) for p in inputs() if p.parent == HERE]
    tracked += [str((HERE / "inputs.json").relative_to(ROOT))]
    for p in tracked:
        subprocess.run(["git", "ls-files", "--error-unmatch", p], cwd=ROOT,
                       check=True, stdout=subprocess.DEVNULL)
    if subprocess.check_output(["git", "status", "--porcelain", "--", *tracked],
                               cwd=ROOT, text=True):
        raise ValueError("Commit protocol, implementation and hashes before measurement")
    E.verify()


@contextlib.contextmanager
def configured_arms():
    """Extend the imported builder's settings only in this process."""
    old = R.ARMS.copy()
    R.ARMS.update(ARMS)
    try:
        yield
    finally:
        R.ARMS.clear()
        R.ARMS.update(old)


def require_same(name, observed, expected):
    if observed != expected:
        raise AssertionError(f"{name}: {observed} != {expected}")


def compose():
    verify()
    if (WORK / "effects.npz").exists():
        raise FileExistsError("Do not overwrite a composed prediction")
    WORK.mkdir(parents=True, exist_ok=True)
    with np.load(BASE / "work/stage.npz", allow_pickle=False) as stage:
        genes, targets, cols = stage["genes"], stage["targets"], stage["columns"]
        cd4_mean, cd4_votes, cd4_counts = (
            stage["cd4_mean"], stage["cd4_votes"], stage["cd4_counts"])
        old_cap15 = stage["cap15"]
    stage_info = json.loads((BASE / "work/stage.json").read_text())
    compose_info = json.loads((BASE / "work/compose.json").read_text())
    require_same("source stage gate", stage_info["gates_pass"], True)
    require_same("source composition gate", compose_info["gates_pass"], True)
    with configured_arms():
        b = R.builder("cap15")
        with R.model_env():
            lean = R.load_module("lean326", ROOT / "experiments/170-gamma-mse/stage.py")
        pm = b.EXP061.prepare_perturbmatch(
            ROOT / "data/external/perturbmatch", tuple(genes),
            chunk_rows=128, max_rss_gib=14)
        with lean.lean_norm():
            donors = {n: lean.load_donor_lean(b, ROOT / f"data/donors/{n}", genes)
                      for n in R.DONORS}
            model = b.MultiDonorModel(donors, genes, gene_norm=False, eb=False)
        abc = pd.read_csv(ROOT / "data/challenge/pert_counts.csv").target_gene.tolist()
        lean.trim_model_rows(model, list(targets) + abc)
        del donors
        weights = b.donor_weights(model)
        op = b.load_production_operator(
            model, targets, R.OPERATOR, R.OPERATOR.with_suffix(".json"),
            expected_sha256=compose_info["new_operator_sha256"])
        R.verify_axis(genes, cols, targets, op.axis_sha256)
        if not np.array_equal(cols, op.columns):
            raise AssertionError("operator columns changed")
        cd4 = b.EXP062.CD4Prepared(
            mode="gene_norm", mean_effect=cd4_mean, gene_votes=cd4_votes,
            target_condition_count=cd4_counts, condition_details=[])
        require_same("CD4 complete targets", int(np.count_nonzero(cd4_counts == 3)), 77)
        positions = np.flatnonzero(np.isin(genes, targets))
        effects, audit = {}, {}
        for name in ARMS:
            m = R.builder(name)
            handled = set()
            pre = m.build_v9_pre_whiten_effects(
                model, targets, alpha=m.ALPHA, beta=m.BETA, weights=weights,
                pm=pm, operator_columns=cols, pm_handled_out=handled)
            cand, pm_audit = m.apply_perturbmatch_k562_replacement(
                pre, model, pm, targets, alpha=m.ALPHA, weights=weights,
                already_blended=handled)
            if name == "cap15":
                for key, value in (("baseline_pre_sha256", pre), ("candidate_pre_sha256", cand)):
                    require_same(f"cap15 {key}", m.array_sha256(value),
                                 stage_info["arms"]["cap15"][key])
                if not np.array_equal(cand, old_cap15):
                    raise AssertionError("cap15 staged array differs")
            w, wa = m.apply_v9_whitening(cand, op.whitener, cols)
            final, ca, ua = m.apply_frozen_cd4_residual(w, cd4, cols, positions)
            if ua is not None or m.H1_UNION_G is not None or m.MASS_FIX is not None:
                raise AssertionError("H1 union/mass-fix configuration drift")
            final_sha = m.array_sha256(final)
            if name == "cap15":
                require_same("cap15 final effect SHA", final_sha,
                             compose_info["arms"]["cap15"]["final_sha256"])
            effects[name] = final
            audit[name] = {"pre_sha256": m.array_sha256(pre),
                           "pm_corrected_sha256": m.array_sha256(cand),
                           "final_sha256": final_sha,
                           "pm": pm_audit, "whitening": wa, "cd4": ca}
            print(f"composed {name}", flush=True)

        # Source-side binding to the actual unsubmitted A/B/C package metadata.
        m = R.builder("cap30")
        pre = m.build_v9_pre_whiten_effects(
            model, abc, alpha=m.ALPHA, beta=m.BETA, weights=weights)
        cand, _ = m.apply_perturbmatch_k562_replacement(
            pre, model, pm, abc, alpha=m.ALPHA, weights=weights)
        expected = json.loads((ROOT / "submissions/319d3-cap30.h5ad.meta.json").read_text())
        require_same("319d3 exponent", expected["config"]["normrel_exp"], -3.0)
        abc_hashes = {}
        for key, value in (("baseline_pre_sha256", pre), ("candidate_pre_sha256", cand)):
            abc_hashes[key] = m.array_sha256(value)
            require_same(f"319d3 {key}", abc_hashes[key], expected["effects"][key])
    permutation = R.derangement(len(targets), seed=20260909)
    effects["shuffled_cap30"] = effects["cap30"][permutation]
    np.savez(WORK / "effects.npz", genes=genes, targets=targets,
             **effects, permutation=permutation)
    R.write_json(WORK / "compose.json", {
        "gates_pass": True, "arms": audit, "abc_cap30_matches_319d3": abc_hashes,
        "operator_sha256": R.sha(R.OPERATOR), "h1_union": "OFF",
        "held_h1_sources_loaded": [], "effects_sha256": R.sha(WORK / "effects.npz")})


@contextlib.contextmanager
def evaluator_output(work):
    """Redirect outputs only; the scoring implementation remains hash-pinned."""
    original = E.R
    link = work / "reference_de.parquet"
    expected = BASE / "work/reference_de.parquet"
    if link.is_symlink():
        require_same("reference DE link", link.resolve(), expected.resolve())
    elif link.exists():
        raise FileExistsError(link)
    else:
        link.symlink_to(expected)
    E.R = SimpleNamespace(**{**vars(original), "WORK": work})
    try:
        yield
    finally:
        E.R = original


def evaluate(name):
    verify()
    comp = json.loads((WORK / "compose.json").read_text())
    require_same("composition gates", comp["gates_pass"], True)
    require_same("effect archive SHA", R.sha(WORK / "effects.npz"), comp["effects_sha256"])
    if (WORK / f"metrics-{name}.json").exists():
        raise FileExistsError(f"refusing to rescore {name}")
    if name != "cap15":
        require_same("cap15 known answer", json.loads(
            (WORK / "known-answer.json").read_text())["pass"], True)
    r, c, emit_c = E.load_reference()
    with np.load(WORK / "effects.npz", allow_pickle=False) as z:
        if not np.array_equal(z["targets"], r["targets"]):
            raise ValueError("target order mismatch")
        effect = R.map_axis(z[name], z["genes"], r["genes"])
    with evaluator_output(WORK):
        result = E.emit_and_score(name, r, c, emit_c, effect)
    if name == "cap15":
        expected = json.loads((BASE / "work/metrics-cap15.json").read_text())
        differences = {m: abs(result["raw"][m] - expected["raw"][m]) for m in E.MEMBERS}
        digests = {k: result[k] == expected[k] for k in
                   ("emitted_count_digest_sha256", "panel_total_digest_sha256")}
        passed = all(digests.values()) and all(
            differences[m] <= E.TOLERANCE[m] for m in E.MEMBERS)
        R.write_json(WORK / "known-answer.json",
                     {"pass": passed, "absolute_differences": differences,
                      "tolerances": E.TOLERANCE, "digests_equal": digests})
        if not passed:
            raise AssertionError("cap15 emitted-count known-answer failed")
    gc.collect()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("freeze", "compose", "cap15", "cap30",
                                        "shuffled_cap30"))
    args = parser.parse_args()
    if args.mode in ("freeze", "compose"):
        globals()[args.mode]()
    else:
        evaluate(args.mode)


if __name__ == "__main__":
    main()
