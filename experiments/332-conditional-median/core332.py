"""Pure split, feature and endpoint functions for experiment 332."""

import hashlib
import numpy as np

SEED = 332
N_TRAIN = 2048
N_GENES = 512
LINES = ("C4I", "HCT15", "HKGZCC", "KMRC20", "KYSE140", "KYSE450",
         "ONE58", "RPMI7951", "SKGII", "SNU8", "T3M5", "U343")
RESERVED = frozenset(("JVE127", "MG63", "SLR23", "UMRC3"))
FEATURES = ("mean", "median", "sd", "positive", "min", "max", "log_scale",
            "source_basal_mean", "source_basal_sd", "source_covariance",
            "recipient_basal", "basal_difference", "mean_x_difference",
            "mean_saturation")
ARMS = ("source_mean", "zero", "tabicl_median", "tabicl_mean", "no_recipient",
        "wrong_recipient", "shuffled_labels", "ridge", "absolute_boost")


def hash_order(names, prefix):
    return sorted(range(len(names)),
                  key=lambda i: hashlib.sha256((prefix + str(names[i])).encode()).hexdigest())


def target_blocks(names):
    ix = hash_order(names, "332-target|")
    return np.array(ix[::2]), np.array(ix[1::2])


def source_indices(n_lines, outer, recipient=None):
    return [i for i in range(n_lines) if i != outer and i != recipient]


def features(source, basal_sources, recipient_basal, source_rms):
    """source: S,T,G; source_rms: S,T measured on the full eligible axis."""
    scale = np.maximum(np.sqrt(np.mean(np.square(source_rms), axis=0)), 1e-3)
    x = source.astype(np.float64) / scale[None, :, None]
    m = x.mean(axis=0)
    bmean = basal_sources.mean(axis=0)
    bsd = basal_sources.std(axis=0)
    diff = recipient_basal - bmean
    cov = ((x - m[None]) * (basal_sources - bmean)[:, None]).mean(axis=0)
    broad = lambda v: np.broadcast_to(v, m.shape)
    out = np.stack((m, np.median(x, axis=0), x.std(axis=0),
                    (x > 0).mean(axis=0), x.min(axis=0), x.max(axis=0),
                    broad(np.log1p(scale)[:, None]), broad(bmean), broad(bsd), cov,
                    broad(recipient_basal), broad(diff), m * diff,
                    m / (1 + recipient_basal)), axis=-1)
    if not np.isfinite(out).all():
        raise ArithmeticError("nonfinite feature")
    return out.astype(np.float32), scale


def fold_features(response, controls, rms, outer):
    """Construct all inputs without reading the held response or its RMS."""
    train_lines = source_indices(len(response), outer)
    blocks = []
    for recipient in train_lines:
        donors = source_indices(len(response), outer, recipient)
        x, scale = features(response[donors], controls[donors], controls[recipient], rms[donors])
        y = response[recipient] / scale[:, None]
        blocks.append((recipient, x, y.astype(np.float32)))
    query, scale = features(response[train_lines], controls[train_lines],
                            controls[outer], rms[train_lines])
    wrong = (outer + 1) % len(response)
    xwrong, scale_wrong = features(response[train_lines], controls[train_lines],
                                  controls[wrong], rms[train_lines])
    np.testing.assert_array_equal(scale, scale_wrong)
    np.testing.assert_array_equal(query[..., :10], xwrong[..., :10])
    return blocks, query, xwrong, scale


def sample_training(blocks, targets, split_seed, n_train=N_TRAIN):
    """Balanced recipient sampling; complete target blocks stay out of training."""
    rng = np.random.default_rng(split_seed)
    xs, ys, shuffles, provenance = [], [], [], []
    count = len(targets) * blocks[0][1].shape[1]
    for i, (line, x, y) in enumerate(blocks):
        n = n_train // len(blocks) + (i < n_train % len(blocks))
        if n > count:
            raise ValueError("not enough distinct pairs")
        take = rng.choice(count, n, replace=False)
        xx = x[targets].reshape(count, -1)
        yy = y[targets].reshape(count)
        # One coherent cyclic target derangement per recipient, keeping genes aligned.
        yy_shuf = y[np.roll(targets, 1)].reshape(count)
        xs.append(xx[take]); ys.append(yy[take]); shuffles.append(yy_shuf[take])
        provenance.extend((line, int(targets[k // x.shape[1]]), int(k % x.shape[1]))
                          for k in take)
    return (np.concatenate(xs), np.concatenate(ys), np.concatenate(shuffles),
            np.array(provenance, dtype=np.int32))


def standardize(x, query):
    mean = x.mean(axis=0, dtype=np.float64)
    sd = np.maximum(x.std(axis=0, dtype=np.float64), 1e-6)
    return ((x - mean) / sd).astype(np.float32), ((query - mean) / sd).astype(np.float32)


def prediction_columns(pred, keys, n_rows):
    """TabICL returns an array even for a one-element output_type list."""
    if len(keys) == 1 and not isinstance(pred, dict):
        pred = {keys[0]: pred}
    if not isinstance(pred, dict) or set(pred) != set(keys):
        raise ValueError("unexpected prediction schema")
    output = {k: np.asarray(v) for k, v in pred.items()}
    if any(v.shape != (n_rows,) or not np.isfinite(v).all() for v in output.values()):
        raise ArithmeticError("nonfinite/misaligned model output")
    return output


def cosine(a, b):
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    return np.sum(a*b, axis=-1) / np.maximum(
        np.linalg.norm(a, axis=-1) * np.linalg.norm(b, axis=-1), 1e-12)


def retrieval(prediction, truth):
    p = prediction / np.maximum(np.linalg.norm(prediction, axis=1, keepdims=True), 1e-12)
    t = truth / np.maximum(np.linalg.norm(truth, axis=1, keepdims=True), 1e-12)
    similarities = p @ t.T
    own = np.diag(similarities)
    rank = ((similarities > own[:, None]).sum(axis=1)
            + .5 * ((similarities == own[:, None]).sum(axis=1) - 1))
    return 1 - rank / (len(truth) - 1)


def orthogonal(a, baseline):
    denom = np.maximum(np.square(baseline).sum(axis=1, keepdims=True), 1e-12)
    return a - baseline * np.sum(a * baseline, axis=1, keepdims=True) / denom


def endpoints(prediction, truth, baseline, scale):
    error = (prediction - truth) / scale[:, None]
    return {"absolute_error": np.mean(np.abs(error), axis=1),
            "squared_error": np.mean(np.square(error), axis=1),
            "cosine": cosine(prediction, truth),
            "retrieval": retrieval(prediction, truth),
            "orthogonal_cosine": cosine(orthogonal(prediction, baseline),
                                       orthogonal(truth, baseline))}


def paired_interval(delta, n=10000):
    """Paired target and crossed line/target bootstraps; no gene resampling."""
    rng = np.random.default_rng(SEED)
    nl, nt = delta.shape
    tw = rng.multinomial(nt, np.full(nt, 1/nt), size=n) / nt
    lw = rng.multinomial(nl, np.full(nl, 1/nl), size=n) / nl
    target = tw @ delta.mean(axis=0)
    crossed = np.einsum("bi,ij,bj->b", lw, delta, tw, optimize=True)
    lm = delta.mean(axis=1)
    return {"mean": float(delta.mean()),
            "target_ci95": np.quantile(target, [.025, .975]).tolist(),
            "crossed_ci95": np.quantile(crossed, [.025, .975]).tolist(),
            "positive_lines": int((lm > 0).sum()), "line_means": lm.tolist(),
            "minimum_leave_one_line_out": float(((lm.sum()-lm)/(nl-1)).min())}
