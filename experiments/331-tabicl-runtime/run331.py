"""A bounded CPU-only test with synthetic data and a hash-verified public model."""

import argparse
from datetime import datetime, timezone
import gc
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import resource
import subprocess
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
CHECKPOINT = WORK / "tabicl-regressor-v2-20260212.ckpt"
CHECKPOINT_SHA = "0db9cb538f114e79026bf08f45f41ad8dd7ad2de2aaca9a5ca8cd3bd9748ae7a"
CHECKPOINT_SIZE = 114324594
MODEL_REVISION = "4dcd344ece2c00be9e831fdd35bed57b5ad83e19"


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix+".tmp")
    temp.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False)+"\n")
    temp.replace(path)


def inputs():
    return [HERE / n for n in ("PREREGISTRATION.md", "requirements.in",
                               "requirements.lock", "run331.py", "test331.py")] + [ROOT / "uv.lock"]


def freeze():
    path = HERE / "inputs.json"
    if path.exists():
        raise FileExistsError(path)
    write(path, {"sha256": {str(p.relative_to(ROOT)): sha(p) for p in inputs()},
                 "checkpoint_sha256": CHECKPOINT_SHA, "checkpoint_bytes": CHECKPOINT_SIZE,
                 "model_revision": MODEL_REVISION})


def verify():
    manifest = json.loads((HERE / "inputs.json").read_text())
    if set(manifest["sha256"]) != {str(p.relative_to(ROOT)) for p in inputs()}:
        raise ValueError("input allowlist drift")
    for path, expected in manifest["sha256"].items():
        if sha(ROOT / path) != expected:
            raise ValueError("331 input drift: "+path)
    own = [str(p.relative_to(ROOT)) for p in inputs() if p.parent == HERE]
    own.append(str((HERE / "inputs.json").relative_to(ROOT)))
    for path in own:
        subprocess.run(["git", "ls-files", "--error-unmatch", path], cwd=ROOT,
                       check=True, stdout=subprocess.DEVNULL)
    if subprocess.check_output(["git", "status", "--porcelain", "--", *own], cwd=ROOT, text=True):
        raise ValueError("commit 331 before execution")


def download():
    verify()
    if CHECKPOINT.exists():
        raise FileExistsError("checkpoint already exists; inspect, do not overwrite")
    WORK.mkdir(parents=True, exist_ok=True)
    url = ("https://huggingface.co/jingang/TabICL/resolve/"+MODEL_REVISION+
           "/tabicl-regressor-v2-20260212.ckpt")
    subprocess.run(["curl", "--fail", "--show-error", "--location", "--max-time", "600",
                    "--output", str(CHECKPOINT), url], check=True)
    if CHECKPOINT.stat().st_size != CHECKPOINT_SIZE or sha(CHECKPOINT) != CHECKPOINT_SHA:
        raise ValueError("public checkpoint size/hash mismatch; do not load")
    write(WORK / "download.json", {
        "verified": True, "bytes": CHECKPOINT_SIZE, "sha256": CHECKPOINT_SHA,
        "repository": "jingang/TabICL", "revision": MODEL_REVISION,
        "utc": datetime.now(timezone.utc).isoformat(), "paid_api": False})


def synthetic(n, seed):
    import numpy as np
    rng = np.random.default_rng(seed)
    x = rng.uniform(size=(n+256, 12)).astype(np.float32)
    y = (10*np.sin(np.pi*x[:, 0]*x[:, 1]) + 20*(x[:, 2]-.5)**2 +
         10*x[:, 3] + 5*x[:, 4] + rng.normal(0, .1, len(x))).astype(np.float32)
    return x[:n], y[:n], x[n:], y[n:]


def test_size(n):
    import numpy as np
    from tabicl import TabICLRegressor

    x, y, xt, yt = synthetic(n, 20260909+n)
    config = {
        "n_estimators": 1, "batch_size": 1, "norm_methods": ["none"],
        "kv_cache": False, "device": "cpu", "use_amp": False, "use_fa3": False,
        "model_path": str(CHECKPOINT), "allow_auto_download": False,
        "random_state": 20260909, "n_jobs": 4, "verbose": False}
    started = time.monotonic()
    model = TabICLRegressor(**config)
    fit_start = time.monotonic()
    model.fit(x, y)
    fit_seconds = time.monotonic()-fit_start
    predict_start = time.monotonic()
    output = model.predict(xt, output_type=["mean", "median"])
    predict_seconds = time.monotonic()-predict_start
    again_start = time.monotonic()
    repeated = model.predict(xt)
    repeat_seconds = time.monotonic()-again_start
    if any(np.asarray(v).shape != (256,) or not np.isfinite(v).all() for v in output.values()):
        raise ArithmeticError("invalid prediction")
    if repeated.shape != (256,) or not np.isfinite(repeated).all():
        raise ArithmeticError("invalid repeat prediction")
    error = float(np.max(np.abs(repeated-output["mean"])))
    del model
    gc.collect()
    shuffled = TabICLRegressor(**config)
    perm = np.random.default_rng(20260909).permutation(n)
    shuffle_start = time.monotonic()
    shuffled.fit(x, y[perm])
    shuffled_pred = shuffled.predict(xt)
    shuffle_seconds = time.monotonic()-shuffle_start
    if shuffled_pred.shape != (256,) or not np.isfinite(shuffled_pred).all():
        raise ArithmeticError("invalid shuffled prediction")
    null_mse = float(np.mean((yt-float(y.mean()))**2))
    mse = {k: float(np.mean((yt-v)**2)) for k, v in output.items()}
    shuffled_mse = float(np.mean((yt-shuffled_pred)**2))
    elapsed = time.monotonic()-started
    result = {
        "n_train": n, "n_test": 256, "features": 12, "repeat_max_abs_difference": error,
        "mse": mse, "training_mean_null_mse": null_mse,
        "r2_vs_training_mean": {k: 1-v/null_mse for k, v in mse.items()},
        "shuffled_mse": shuffled_mse,
        "seconds": {"fit": fit_seconds, "predict": predict_seconds, "repeat": repeat_seconds,
                    "shuffled_fit_predict": shuffle_seconds, "total": elapsed},
        "gates": {"repeatable": error <= 1e-6,
                  "beats_mean_null": mse["mean"] < null_mse,
                  "beats_shuffled_labels": mse["mean"] < shuffled_mse,
                  "within_15_minutes": elapsed <= 900},
        "process_max_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "synthetic_input_sha256": hashlib.sha256(
            b"".join(np.ascontiguousarray(v).tobytes() for v in (x, y, xt, yt))).hexdigest()}
    write(WORK / f"synthetic-{n}.json", result)
    print(json.dumps(result, indent=2), flush=True)
    del shuffled
    gc.collect()
    return result


def run():
    verify()
    if (WORK / "results.json").exists():
        raise FileExistsError("do not rerun a terminal probe")
    if CHECKPOINT.stat().st_size != CHECKPOINT_SIZE or sha(CHECKPOINT) != CHECKPOINT_SHA:
        raise ValueError("checkpoint drift")
    if version("tabicl") != "2.2.0":
        raise ValueError("unexpected tabicl version")
    import torch
    torch.set_num_threads(4)
    torch.set_num_interop_threads(1)
    start, timer = datetime.now(timezone.utc).isoformat(), time.monotonic()
    rows = []
    for n in (512, 2048):
        row = test_size(n)
        rows.append(row)
        if not all(row["gates"].values()):
            break
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    passed = len(rows) == 2 and all(all(r["gates"].values()) for r in rows) and peak <= 12*1024**3
    result = {
        "status": "ENGINEERING_FEASIBLE_ONLY" if passed else "ENGINEERING_NO_GO",
        "classification": "synthetic-only runtime/sanity test; no biological evidence",
        "workloads": rows, "device": "cpu", "mps_tested": False,
        "packages": {p: version(p) for p in ("tabicl", "torch", "numpy", "scipy",
                                           "scikit-learn", "huggingface-hub")},
        "started_utc": start, "finished_utc": datetime.now(timezone.utc).isoformat(),
        "seconds": time.monotonic()-timer, "process_max_rss_bytes_macos": peak,
        "checkpoint_sha256": CHECKPOINT_SHA, "inputs_sha256": sha(HERE / "inputs.json"),
        "biological_data_read": False, "submission": False}
    write(WORK / "results.json", result)
    print(json.dumps(result, indent=2), flush=True)
    verify()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("freeze", "download", "run"))
    args = parser.parse_args()
    globals()[args.mode]()
