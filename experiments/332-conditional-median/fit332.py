"""Reads training labels and query FEATURES only; no outer evaluation file."""

import gc
from importlib.metadata import version
import json
import resource
import sys
import time

import numpy as np
import torch
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from tabicl import TabICLRegressor

import core332 as C
from run332 import WORK, RUNTIME, CHECKPOINT, CHECKPOINT_SHA, sha, write


def predict_chunks(model, x, median_and_mean=False):
    keys = ["mean", "median"] if median_and_mean else ["median"]
    output = {k: [] for k in keys}
    for offset in range(0, len(x), 256):
        pred = model.predict(x[offset:offset+256], output_type=keys)
        pred = C.prediction_columns(pred, keys, min(256, len(x)-offset))
        for key in keys:
            output[key].append(pred[key])
    return {k: np.concatenate(v) for k, v in output.items()}


def main(line, block):
    if line not in C.LINES or line in C.RESERVED or block not in (0, 1):
        raise ValueError("unregistered split")
    torch.set_num_threads(4)
    torch.set_num_interop_threads(1)
    directory = WORK / line
    source, output_path = directory / f"inputs-{block}.npz", directory / f"prediction-{block}.npz"
    if output_path.exists():
        raise FileExistsError("completed split cannot be overwritten")
    prep = json.loads((directory / "prepared.json").read_text())
    assert sha(source) == prep["splits"][block]["inputs_sha256"]
    assert sha(CHECKPOINT) == CHECKPOINT_SHA
    assert version("tabicl") == "2.2.0"
    with np.load(source) as z:
        x, y, ys = z["x_train"], z["y_train"], z["y_shuffled"]
        xt, xw, indices = z["x_test"], z["x_wrong"], z["test_indices"]
    sx, st = C.standardize(x, xt)
    _, sw = C.standardize(x, xw)
    config = dict(n_estimators=1, batch_size=1, norm_methods=["none"], kv_cache=False,
                  device="cpu", use_amp=False, use_fa3=False, model_path=str(CHECKPOINT),
                  allow_auto_download=False, random_state=C.SEED, n_jobs=4, verbose=False)
    start = time.monotonic()
    outputs, timings = {}, {}
    for name, xx, yy, query in (("full", sx, y, st),
                                ("no_recipient", sx[:, :10], y, st[:, :10]),
                                ("shuffled_labels", sx, ys, st)):
        now = time.monotonic()
        model = TabICLRegressor(**config).fit(xx, yy)
        pred = predict_chunks(model, query, median_and_mean=name == "full")
        if name == "full":
            outputs["tabicl_median"], outputs["tabicl_mean"] = pred["median"], pred["mean"]
            outputs["wrong_recipient"] = predict_chunks(model, sw)["median"]
        else:
            outputs[name] = pred["median"]
        timings[name] = time.monotonic()-now
        print(line, block, name, f"{timings[name]:.1f}s", flush=True)
        del model
        gc.collect()
    for name, model in (("ridge", Ridge(alpha=10)),
                        ("absolute_boost", HistGradientBoostingRegressor(
                            loss="absolute_error", max_iter=100, learning_rate=.05,
                            max_leaf_nodes=15, min_samples_leaf=20, l2_regularization=1,
                            early_stopping=False, random_state=C.SEED))):
        now = time.monotonic()
        outputs[name] = model.fit(sx, y).predict(st)
        timings[name] = time.monotonic()-now
    if any(v.shape != (len(st),) or not np.isfinite(v).all() for v in outputs.values()):
        raise ArithmeticError("invalid fitted predictions")
    np.savez_compressed(output_path, test_indices=indices, **outputs)
    receipt = {"line": line, "block": block, "input_sha256": sha(source),
               "prediction_sha256": sha(output_path), "seconds": time.monotonic()-start,
               "timings": timings, "max_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
               "packages": {p: version(p) for p in ("tabicl", "torch", "numpy", "scikit-learn")},
               "checkpoint_sha256": CHECKPOINT_SHA, "outer_truth_read": False}
    write(directory / f"fit-{block}.json", receipt)
    print("completed", line, block, f"{receipt['seconds']:.1f}s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
