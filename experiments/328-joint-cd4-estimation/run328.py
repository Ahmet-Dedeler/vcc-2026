"""Metadata registration, bounded source extraction, and held-target mixture test."""

import argparse
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import h5py
import numpy as np
import pandas as pd

import mixture328 as M

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
SOURCE = ROOT / "experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad"
COMPACT = ROOT / "experiments/273-flex-reference-instrument/work/cd4/cd4_compact_h1"
BASE = ROOT / "experiments/323-current-family-validation"
CONDITIONS = ("Rest", "Stim8hr", "Stim48hr")
PACKAGES = ("numpy", "scipy", "h5py", "anndata")
sys.path.insert(0, str(BASE))
import run323 as R


def text_column(group, key):
    obj = group[key]
    if isinstance(obj, h5py.Group):
        categories = obj["categories"].asstr()[:]
        codes = obj["codes"][:]
        return np.array([categories[c] if c >= 0 else "" for c in codes], dtype=object)
    return obj.asstr()[:]


def digest_text(value):
    return hashlib.sha256(value.encode()).hexdigest()


def select():
    """Only identifiers and schema; no selected response or standard-error values."""
    destination = HERE / "selection.json"
    if destination.exists():
        raise FileExistsError(destination)
    h1 = json.loads((BASE / "inputs.json").read_text())["h1_targets"]
    abc = pd.read_csv(ROOT / "data/challenge/pert_counts.csv").target_gene.tolist()
    excluded = {str(t).strip().casefold() for t in h1+abc}
    gene_map = pd.read_csv(COMPACT / "gene_map.csv")
    with h5py.File(SOURCE, "r") as f:
        names = text_column(f["obs"], "target_contrast_gene_name")
        ids = text_column(f["obs"], "target_contrast")
        conditions = text_column(f["obs"], "culture_condition")
        var_names = text_column(f["var"], "gene_name")
        var_ids = text_column(f["var"], "gene_ids")
        excluded_ids = {str(i).split(".")[0] for i, n in zip(var_ids, var_names)
                        if str(n).strip().casefold() in excluded}
        shape = list(f["layers/log_fc"].shape)
        if shape != list(f["layers/lfcSE"].shape):
            raise ValueError("effect and SE schema differs")
    groups = {}
    for row, (name, identifier, condition) in enumerate(zip(names, ids, conditions)):
        name, identifier = str(name).strip(), str(identifier).strip().split(".")[0]
        if (not name or not identifier or condition not in CONDITIONS
                or name.casefold() in excluded or identifier in excluded_ids):
            continue
        groups.setdefault(identifier, []).append((name, str(condition), row))
    eligible = {}
    for identifier, rows in groups.items():
        if (len(rows) == 3 and len({n for n, _, _ in rows}) == 1
                and {c for _, c, _ in rows} == set(CONDITIONS)):
            eligible[identifier] = {
                "id": identifier, "name": rows[0][0],
                "rows": {c: row for _, c, row in rows}}
    ordered = sorted(eligible, key=lambda i: digest_text("328|"+i))
    if len(ordered) < 1250:
        raise ValueError("fewer than 1250 unambiguous non-panel targets")
    targets = [{**eligible[i], "split": "fit" if j < 1000 else "development"}
               for j, i in enumerate(ordered[:1250])]
    genes = [r for r in gene_map.to_dict("records")
             if str(r["gene_name"]).casefold() not in excluded]
    genes = sorted(genes, key=lambda r: digest_text("328-gene|"+r["gene_name"]))[:256]
    if len(genes) != 256 or len({r["source_var_index"] for r in genes}) != 256:
        raise ValueError("gene selection failed")
    R.write_json(destination, {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "response_values_read": False, "source_layer_shape": shape,
        "eligible_nonpanel_targets": len(ordered), "targets": targets,
        "genes": [{"name": r["gene_name"], "column": int(r["source_var_index"])} for r in genes],
        "excluded_symbols": sorted(excluded), "excluded_ids": sorted(excluded_ids),
        "conditions": CONDITIONS})
    print(f"Registered {len(targets)} targets and {len(genes)} genes, metadata only", flush=True)


def input_paths():
    files = [HERE / n for n in ("PREREGISTRATION.md", "selection.json", "run328.py",
                                "mixture328.py", "test328.py")]
    files += [SOURCE, ROOT / "uv.lock", BASE / "inputs.json", BASE / "run323.py",
              COMPACT / "gene_map.csv", COMPACT / "target_row_map.csv",
              ROOT / "data/challenge/pert_counts.csv",
              ROOT / "experiments/047-cd4-direct/extract_cd4.py"]
    for condition in ("rest", "stim8hr", "stim48hr"):
        files += [COMPACT / ("cd4_"+condition) / n for n in
                  ("perts.txt", "genes.txt", "log_fc.npy", "lfc_se.npy", "delta.npy")]
    return files


def freeze():
    path = HERE / "inputs.json"
    if path.exists():
        raise FileExistsError(path)
    R.write_json(path, {
        "inputs_sha256": {str(p.relative_to(ROOT)): R.sha(p) for p in input_paths()},
        "packages": {p: version(p) for p in PACKAGES}})
    print("Frozen input manifest", flush=True)


def verify():
    manifest = json.loads((HERE / "inputs.json").read_text())
    if set(manifest["inputs_sha256"]) != {str(p.relative_to(ROOT)) for p in input_paths()}:
        raise ValueError("input allowlist changed")
    for name, expected in manifest["inputs_sha256"].items():
        if R.sha(ROOT / name) != expected:
            raise ValueError(f"source pilot input drift: {name}")
    if manifest["packages"] != {p: version(p) for p in PACKAGES}:
        raise ValueError("runtime package changed")
    own = [str(p.relative_to(ROOT)) for p in input_paths() if p.parent == HERE]
    own += [str((HERE / "inputs.json").relative_to(ROOT))]
    for path in own:
        subprocess.run(["git", "ls-files", "--error-unmatch", path],
                       cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    if subprocess.check_output(["git", "status", "--porcelain", "--", *own], cwd=ROOT, text=True):
        raise ValueError("commit code, protocol and selection before response access")


def known_answers(source):
    X = R.load_module("extract328", ROOT / "experiments/047-cd4-direct/extract_cd4.py")
    mapping = pd.read_csv(COMPACT / "target_row_map.csv")
    columns = pd.read_csv(COMPACT / "gene_map.csv").source_var_index.to_numpy()
    result = {}
    for condition in CONDITIONS:
        bank = COMPACT / ("cd4_"+condition.lower())
        target = (bank / "perts.txt").read_text().splitlines()[0]
        rows = mapping[(mapping.condition == condition) & (mapping.target_gene == target)
                       & (mapping.selected == True)]
        if len(rows) != 1:
            raise ValueError("ambiguous known-answer row")
        row = int(rows.iloc[0].source_row)
        measured = {}
        for layer, file in (("log_fc", "log_fc.npy"), ("lfcSE", "lfc_se.npy")):
            observed = X.read_matrix_subset(source["layers"][layer], np.array([row]), columns,
                                             row_block=1, max_rss_gib=3)[0]
            expected = np.load(bank / file, mmap_mode="r")[0]
            if not np.array_equal(observed, expected, equal_nan=True):
                raise AssertionError(f"compact {condition}/{layer} not reproduced")
            measured[layer] = observed
        beta = measured["log_fc"]
        delta = np.where(np.isfinite(beta), beta*X.LN2, 0.).astype(np.float32)
        if not np.array_equal(delta, np.load(bank / "delta.npy", mmap_mode="r")[0]):
            raise AssertionError(f"compact {condition}/delta not reproduced")
        result[condition] = {"target": target, "source_row": row,
                             "log_fc_exact": True, "lfcSE_exact": True, "delta_exact": True}
    return result


def extract():
    verify()
    if (WORK / "source.npz").exists():
        raise FileExistsError("source already extracted")
    selection = json.loads((HERE / "selection.json").read_text())
    columns = np.array([g["column"] for g in selection["genes"]])
    order = np.argsort(columns)
    reverse = np.argsort(order)
    targets = selection["targets"]
    y, se = np.empty((1250, 256, 3)), np.empty((1250, 256, 3))
    with h5py.File(SOURCE, "r") as f:
        ka = known_answers(f)
        for i, target in enumerate(targets):
            for c, condition in enumerate(CONDITIONS):
                row = target["rows"][condition]
                y[i, :, c] = f["layers/log_fc"][row, columns[order]][reverse]
                se[i, :, c] = f["layers/lfcSE"][row, columns[order]][reverse]
            if (i+1) % 250 == 0:
                print(f"source triples extracted: {i+1}/1250 targets", flush=True)
    valid = np.all(np.isfinite(y) & np.isfinite(se) & (se > 0), axis=2)
    counts = valid.sum(1)
    enough = counts >= 64
    support = {"fit_targets_64plus": int(enough[:1000].sum()),
               "development_targets_64plus": int(enough[1000:].sum()),
               "complete_triples": int(valid.sum()), "per_target_complete": counts.tolist()}
    support_pass = support["fit_targets_64plus"] >= 800 and support["development_targets_64plus"] >= 200
    WORK.mkdir(parents=True, exist_ok=True)
    np.savez(WORK / "source.npz", y=y, se=se, valid=valid, enough=enough,
             targets=np.array([t["id"] for t in targets]))
    R.write_json(WORK / "source.json", {
        "known_answers": ka, "support": support, "support_pass": support_pass,
        "source_archive_sha256": R.sha(WORK / "source.npz"),
        "selection_sha256": R.sha(HERE / "selection.json"),
        "response_scope": "selected CD4 source triples and three fixed compact known-answer rows",
        "reserve_access": False, "held_h1_outcomes_read": False})
    if not support_pass:
        raise ValueError("insufficient support; no replacement selection is allowed")
    print(json.dumps({k: v for k, v in support.items() if k != "per_target_complete"}), flush=True)


def target_means(values, valid):
    counts = valid.sum(1)
    if np.any(counts < 64):
        raise ValueError("ineligible target in evaluation")
    sums = np.zeros(len(valid))
    target_rows = np.nonzero(valid)[0]
    np.add.at(sums, target_rows, values)
    return sums/counts


def fits():
    verify()
    if (WORK / "results.json").exists():
        raise FileExistsError("source test is terminal")
    metadata = json.loads((WORK / "source.json").read_text())
    if not metadata["support_pass"] or R.sha(WORK / "source.npz") != metadata["source_archive_sha256"]:
        raise AssertionError("extraction gate/hash")
    with np.load(WORK / "source.npz", allow_pickle=False) as z:
        enough = z["enough"]
        y = z["y"][:1000][enough[:1000]]
        se = z["se"][:1000][enough[:1000]]
        dy = z["y"][1000:][enough[1000:]]
        ds = z["se"][1000:][enough[1000:]]
        dev_targets = z["targets"][1000:][enough[1000:]]
    valid = np.all(np.isfinite(y) & np.isfinite(se) & (se > 0), axis=2)
    dev_valid = np.all(np.isfinite(dy) & np.isfinite(ds) & (ds > 0), axis=2)
    models, per_target = {}, {}
    for name in ("joint", "diagonal", "shuffled_joint"):
        spec = M.components(diagonal=name == "diagonal")
        fy, fs = M.shuffled_targets(y, se) if name == "shuffled_joint" else (y, se)
        fvalid = np.all(np.isfinite(fy) & np.isfinite(fs) & (fs > 0), axis=2)
        if not np.array_equal(valid, fvalid):
            raise AssertionError("shuffled fit support changed")
        print(f"fitting {name}: {valid.sum()} triples, {len(spec)} components", flush=True)
        logpdf = M.log_densities(fy[valid], fs[valid], spec)
        result = M.fit(logpdf)
        del logpdf
        dev_logpdf = M.log_densities(dy[dev_valid], ds[dev_valid], spec)
        predicted = M.predictive(dev_logpdf, result["weights"])
        del dev_logpdf
        per_target[name] = target_means(predicted, dev_valid)
        models[name] = {**result, "components": spec,
                        "held_target_mean_log_density": float(per_target[name].mean())}
        R.write_json(WORK / f"fit-{name}.json", models[name])
    rng = np.random.default_rng(20260909)
    draws = rng.integers(0, len(dev_targets), size=(10000, len(dev_targets)))
    comparisons = {}
    for other in ("diagonal", "shuffled_joint"):
        difference = per_target["joint"]-per_target[other]
        boot = difference[draws].mean(1)
        comparisons["joint-minus-"+other] = {
            "mean_log_density_delta": float(difference.mean()),
            "ci95": np.quantile(boot, [.025, .975]).tolist(),
            "fraction_targets_favor_joint": float(np.mean(difference > 0))}
    gates = {
        "all_fits_converged": all(m["converged"] for m in models.values()),
        "joint_beats_diagonal_interval": comparisons["joint-minus-diagonal"]["ci95"][0] > 0,
        "joint_beats_shuffled_interval": comparisons["joint-minus-shuffled_joint"]["ci95"][0] > 0,
        "at_least_60_percent_targets": comparisons["joint-minus-diagonal"]["fraction_targets_favor_joint"] >= .6,
    }
    result = {
        "status": "EARNS_H1_COUNT_STAGE_REGISTRATION" if all(gates.values()) else "SOURCE_PILOT_NO_GO",
        "gates": gates, "comparisons": comparisons,
        "fit_targets": len(y), "development_targets": len(dy), "genes": y.shape[1],
        "fit_triples": int(valid.sum()), "development_triples": int(dev_valid.sum()),
        "bootstrap": {"draws": 10000, "seed": 20260909, "unit": "target",
                      "fixed": ["genes", "conditions", "source study", "models"]},
        "input_manifest_sha256": R.sha(HERE / "inputs.json"),
        "source_archive_sha256": R.sha(WORK / "source.npz"),
        "models": {k: {f: v[f] for f in ("converged", "iterations", "held_target_mean_log_density")}
                   for k, v in models.items()},
        "not_measured": ["latent biological accuracy", "H1 count metrics", "A/B/C", "independent study"],
        "reserve_access": False, "submission": False,
    }
    R.write_json(WORK / "results.json", result)
    pd.DataFrame({"target": dev_targets, **per_target}).to_csv(WORK / "per-target.csv", index=False)
    print(json.dumps(result, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("select", "freeze", "extract", "fits"))
    mode = parser.parse_args().mode
    started, timer = datetime.now(timezone.utc).isoformat(), time.monotonic()
    globals()[mode]()
    if mode in ("extract", "fits"):
        R.write_json(WORK / ("execution-"+mode+".json"), {
            "status": "complete", "started_utc": started,
            "finished_utc": datetime.now(timezone.utc).isoformat(),
            "seconds": time.monotonic()-timer,
            "process_max_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})


if __name__ == "__main__":
    main()
