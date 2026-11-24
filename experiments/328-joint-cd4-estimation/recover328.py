"""Continue the identical unfinished EM trajectory, preserving v1 evidence."""

import argparse
from datetime import datetime, timezone
import json
import resource
import subprocess
import time

import numpy as np
import pandas as pd
from scipy.special import logsumexp

import run328 as P

R, M = P.R, P.M
OUT = P.WORK / "recovery"


def continue_em(logpdf, initial):
    weights = np.asarray(initial["weights"], float)
    trace = list(initial["mean_loglik_trace"])
    initial_value = float(M.predictive(logpdf, weights).mean())
    if abs(initial_value-trace[-1]) > 1e-10:
        raise AssertionError("v1 training likelihood not reproduced")
    if initial["converged"]:
        return dict(initial)
    original_iterations = int(initial["iterations"])
    converged = False
    for iteration in range(original_iterations+1, 2001):
        with np.errstate(divide="ignore"):
            log_joint = logpdf+np.log(weights)[None, :]
        normalizer = logsumexp(log_joint, axis=1)
        weights = np.exp(log_joint-normalizer[:, None]).mean(0)
        weights /= weights.sum()
        current = float(M.predictive(logpdf, weights).mean())
        improvement = current-trace[-1]
        if not np.isfinite(current) or improvement < -1e-9:
            raise ArithmeticError("nonfinite or decreasing likelihood")
        trace.append(current)
        if iteration % 25 == 0:
            print(f"continued EM {iteration}: mean={current:.8f}, change={improvement:.3g}", flush=True)
        if iteration-original_iterations >= 10 and improvement < 1e-6:
            converged = True
            break
    return {**initial, "weights": weights.tolist(), "iterations": iteration,
            "converged": converged, "mean_loglik_trace": trace}


def files():
    return [P.HERE / n for n in ("NUMERICAL-RECOVERY.md", "recover328.py", "test_recovery328.py",
                                 "inputs.json", "selection.json")] + [
        P.WORK / n for n in ("source.npz", "source.json", "results.json", "per-target.csv",
                             "fit-joint.json", "fit-diagonal.json", "fit-shuffled_joint.json")]


def freeze():
    path = P.HERE / "recovery-inputs.json"
    if path.exists():
        raise FileExistsError(path)
    R.write_json(path, {"sha256": {str(p.relative_to(P.ROOT)): R.sha(p) for p in files()}})


def verify():
    P.verify()
    manifest = json.loads((P.HERE / "recovery-inputs.json").read_text())
    if set(manifest["sha256"]) != {str(p.relative_to(P.ROOT)) for p in files()}:
        raise ValueError("recovery input allowlist drift")
    for path, expected in manifest["sha256"].items():
        if R.sha(P.ROOT / path) != expected:
            raise ValueError("recovery input drift: "+path)
    own = [str((P.HERE / n).relative_to(P.ROOT)) for n in (
        "NUMERICAL-RECOVERY.md", "recover328.py", "test_recovery328.py", "recovery-inputs.json")]
    for path in own:
        subprocess.run(["git", "ls-files", "--error-unmatch", path], check=True,
                       cwd=P.ROOT, stdout=subprocess.DEVNULL)
    if subprocess.check_output(["git", "status", "--porcelain", "--", *own], cwd=P.ROOT, text=True):
        raise ValueError("commit recovery freeze before continuation")


def run():
    verify()
    if OUT.exists():
        raise FileExistsError("do not overwrite a recovery")
    started, timer = datetime.now(timezone.utc).isoformat(), time.monotonic()
    old_rows = pd.read_csv(P.WORK / "per-target.csv")
    with np.load(P.WORK / "source.npz", allow_pickle=False) as z:
        enough = z["enough"]
        y, se = z["y"][:1000][enough[:1000]], z["se"][:1000][enough[:1000]]
        dy, ds = z["y"][1000:][enough[1000:]], z["se"][1000:][enough[1000:]]
        targets = z["targets"][1000:][enough[1000:]]
    if not np.array_equal(old_rows.target.to_numpy(), targets):
        raise ValueError("v1 target order changed")
    valid = np.all(np.isfinite(y) & np.isfinite(se) & (se > 0), axis=2)
    dev_valid = np.all(np.isfinite(dy) & np.isfinite(ds) & (ds > 0), axis=2)
    models, rows, known = {}, {}, {}
    for name in ("joint", "diagonal", "shuffled_joint"):
        initial = json.loads((P.WORK / f"fit-{name}.json").read_text())
        specifications = M.components(diagonal=name == "diagonal")
        if specifications != initial["components"]:
            raise AssertionError("component grid changed")
        fy, fs = M.shuffled_targets(y, se) if name == "shuffled_joint" else (y, se)
        logpdf = M.log_densities(fy[valid], fs[valid], specifications)
        train_error = abs(float(M.predictive(logpdf, initial["weights"]).mean())
                          - initial["mean_loglik_trace"][-1])
        dev_logpdf = M.log_densities(dy[dev_valid], ds[dev_valid], specifications)
        old_pred = P.target_means(M.predictive(dev_logpdf, initial["weights"]), dev_valid)
        dev_error = float(np.max(np.abs(old_pred-old_rows[name].to_numpy())))
        if train_error > 1e-10 or dev_error > 1e-10:
            raise AssertionError(f"v1 reproduction failed: {name}")
        known[name] = {"training_endpoint_abs_error": train_error,
                       "held_target_max_abs_error": dev_error}
        result = continue_em(logpdf, initial)
        del logpdf
        if initial["converged"] and result["weights"] != initial["weights"]:
            raise AssertionError("changed an already-converged control")
        rows[name] = P.target_means(M.predictive(dev_logpdf, result["weights"]), dev_valid)
        del dev_logpdf
        models[name] = {**result, "held_target_mean_log_density": float(rows[name].mean())}
    rng = np.random.default_rng(20260909)
    draws = rng.integers(0, len(targets), (10000, len(targets)))
    comparisons = {}
    for other in ("diagonal", "shuffled_joint"):
        delta = rows["joint"]-rows[other]
        comparisons["joint-minus-"+other] = {
            "mean_log_density_delta": float(delta.mean()),
            "ci95": np.quantile(delta[draws].mean(1), [.025, .975]).tolist(),
            "fraction_targets_favor_joint": float(np.mean(delta > 0))}
    gates = {
        "all_fits_converged": all(v["converged"] for v in models.values()),
        "joint_beats_diagonal_interval": comparisons["joint-minus-diagonal"]["ci95"][0] > 0,
        "joint_beats_shuffled_interval": comparisons["joint-minus-shuffled_joint"]["ci95"][0] > 0,
        "at_least_60_percent_targets": comparisons["joint-minus-diagonal"]["fraction_targets_favor_joint"] >= .6}
    result = {
        "status": "EARNS_H1_COUNT_STAGE_REGISTRATION" if all(gates.values()) else "RECOVERY_NO_GO",
        "classification": "post-result numerical amendment; repeated source development, not confirmation",
        "gates": gates, "comparisons": comparisons, "known_answer": known,
        "models": {n: {k: v[k] for k in ("iterations", "converged", "held_target_mean_log_density")}
                   for n, v in models.items()},
        "bootstrap": {"draws": 10000, "seed": 20260909, "unit": "target"},
        "original_result_sha256": R.sha(P.WORK / "results.json"),
        "recovery_manifest_sha256": R.sha(P.HERE / "recovery-inputs.json"),
        "execution": {"started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(),
                      "seconds": time.monotonic()-timer,
                      "process_max_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
        "held_h1_outcomes_read": False, "submission": False,
    }
    OUT.mkdir()
    for name, model in models.items():
        R.write_json(OUT / f"fit-{name}.json", model)
    pd.DataFrame({"target": targets, **rows}).to_csv(OUT / "per-target.csv", index=False)
    R.write_json(OUT / "results.json", result)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("freeze", "run"))
    globals()[parser.parse_args().mode]()
