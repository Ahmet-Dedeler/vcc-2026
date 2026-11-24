"""Pure, held-recipient-safe primitives for the registered text pilot."""

import hashlib
import re
import numpy as np
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.ensemble import HistGradientBoostingRegressor

SEED = 3340909
N_TARGETS, N_GENES, N_TRAIN, DIMS = 1024, 512, 32768, 16
ARMS = ("source_mean", "zero", "diagonal", "source_only", "string",
        "text", "text_shuffled", "labels_shuffled")
FITTED = ARMS[3:]
SOURCE_FEATURES = ("mean", "median", "std", "positive_fraction", "mean_abs",
                   "min", "max", "basal_mean", "basal_std",
                   "raw_rms_mean", "raw_rms_std", "row_consensus_rms")


def order(names, salt):
    return sorted(range(len(names)), key=lambda i: hashlib.sha256(
        (salt + str(names[i])).encode()).hexdigest())


def rms_direction(a):
    a = np.asarray(a, dtype=np.float64)
    rms = np.sqrt(np.mean(a*a, axis=-1))
    return a / np.maximum(rms[..., None], 1e-12), rms


def unit(a):
    a = np.asarray(a, dtype=np.float64)
    return a / np.maximum(np.linalg.norm(a, axis=-1, keepdims=True), 1e-12)


def features(direction, basal, rms, donors, targets):
    """Nothing outside explicit donor indices can affect the output."""
    donors, targets = np.asarray(donors), np.asarray(targets)
    v = direction[donors][:, targets]
    mean = v.mean(0)
    shape = mean.shape
    broad = lambda x: np.broadcast_to(x, shape)
    b = basal[donors]
    r = rms[donors][:, targets]
    x = np.stack((mean, np.median(v, axis=0), v.std(0), (v > 0).mean(0),
                  np.abs(v).mean(0), v.min(0), v.max(0),
                  broad(b.mean(0)), broad(b.std(0)),
                  broad(r.mean(0)[:, None]), broad(r.std(0)[:, None]),
                  broad(np.sqrt(np.mean(mean**2, axis=1))[:, None])), axis=-1)
    assert x.shape == (*shape, len(SOURCE_FEATURES))
    assert np.isfinite(x).all()
    return x


def source_blocks(direction, basal, rms, outer, targets):
    """Five inner recipients; outer line absent even as a donor."""
    train = [i for i in range(len(direction)) if i != outer]
    return [(recipient, features(direction, basal, rms,
                                 [i for i in train if i != recipient], targets),
             direction[recipient, targets]) for recipient in train]


def query_blocks(direction, basal, rms, outer, targets):
    train = [i for i in range(len(direction)) if i != outer]
    return [features(direction, basal, rms,
                     [i for i in train if i != omit], targets) for omit in train]


def sample_indices(n_recipients, n_targets, n_genes, n, seed):
    rng = np.random.default_rng(seed)
    base, rem = divmod(n, n_recipients)
    rows = []
    for r in range(n_recipients):
        ix = rng.choice(n_targets*n_genes, base+(r < rem), replace=False)
        rows.append(np.column_stack((np.full(len(ix), r), ix//n_genes, ix % n_genes)))
    out = np.concatenate(rows).astype(np.int64)
    return out[rng.permutation(len(out))]


def diagonal(blocks):
    num = sum((x[..., 0]*y).sum(0) for _, x, y in blocks)
    den = sum((x[..., 0]**2).sum(0) for _, x, _ in blocks)
    prior = float(np.median(den))
    return (num+prior)/np.maximum(den+prior, 1e-12)


def description(mapping, key):
    value = mapping.get(key, "")
    if not isinstance(value, str):
        raise ValueError("text must be string")
    m = re.match(r"^Gene Symbol " + re.escape(key) + r"(?:\s+|$)", value)
    return value[m.end():].strip() if m else ""


def text_representation(docs, fit, seed):
    vectorizer = TfidfVectorizer(lowercase=True, stop_words="english",
                               ngram_range=(1, 2), sublinear_tf=True,
                               min_df=2, max_features=12000, norm="l2")
    a = vectorizer.fit_transform([docs[i] for i in fit])
    svd = TruncatedSVD(n_components=DIMS, n_iter=7, random_state=seed)
    svd.fit(a)
    emb = unit(svd.transform(vectorizer.transform(docs)))
    emb[np.array([not bool(t) for t in docs])] = 0
    return emb, vectorizer, svd


def string_representation(values, available, fit, seed):
    pca = PCA(n_components=DIMS, svd_solver="full", random_state=seed)
    pca.fit(values[fit])
    emb = unit(pca.transform(values))
    emb[~available] = 0
    return emb, pca


def pair_features(target_vectors, gene_vectors, t, g):
    a, b = target_vectors[t], gene_vectors[g]
    return np.column_stack((a, b, a*b, np.sum(a*b, axis=1)))


def derangement(n, seed):
    if n < 2:
        raise ValueError("need at least two rows")
    rng = np.random.default_rng(seed)
    cycle = rng.permutation(n)
    result = np.empty(n, dtype=np.int64)
    result[cycle] = np.roll(cycle, 1)
    assert np.all(result != np.arange(n))
    return result


def model(seed):
    return HistGradientBoostingRegressor(
        loss="squared_error", max_iter=150, learning_rate=.05,
        max_leaf_nodes=15, max_depth=4, min_samples_leaf=64,
        l2_regularization=10., early_stopping=False, random_state=seed)


def endpoints(prediction, truth):
    p, t = unit(prediction), unit(truth)
    s = p @ t.T
    own = np.diag(s).copy()
    beaten = (s < own[:, None]).sum(1)
    ties = (s == own[:, None]).sum(1)-1
    return {"absolute_error": np.abs(prediction-truth).mean(1),
            "squared_error": ((prediction-truth)**2).mean(1),
            "cosine": own,
            "retrieval": (beaten + .5*ties)/(len(truth)-1),
            "decoy_cosine": (s.sum(1)-own)/(len(truth)-1)}


def interval(delta, draws=10000):
    """Common target resamples, plus crossed recipient resamples."""
    rng = np.random.default_rng(SEED)
    l, t = delta.shape
    a, b = [], []
    for start in range(0, draws, 100):
        n = min(100, draws-start)
        ti = rng.integers(t, size=(n, t))
        li = rng.integers(l, size=(n, l))
        byline = delta[:, ti].mean(-1).T
        a.extend(byline.mean(1))
        b.extend(byline[np.arange(n)[:, None], li].mean(1))
    m = delta.mean(1)
    return {"mean": float(m.mean()),
            "target_ci95": np.quantile(a, [.025, .975]).tolist(),
            "crossed_ci95": np.quantile(b, [.025, .975]).tolist(),
            "line_means": m.tolist(), "positive_lines": int((m > 0).sum()),
            "minimum_leave_one_line_out": float(((m.sum()-m)/(l-1)).min())}
