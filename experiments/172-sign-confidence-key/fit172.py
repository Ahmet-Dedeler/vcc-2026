"""172 — fit `c(t, g)` = P(our transferred sign is correct) on the LODO folds, score on h1-test.

PROTOCOL, and the reason for its shape
--------------------------------------
Nothing may be selected on h1-test.  So every hyperparameter, every model class and every
top-heavy weighting is chosen by a LEAVE-ONE-FOLD-OUT loop *inside* the donor folds: train on
two donor lines, score the third with the same `reach` rule the metric uses.  Only after a
model class is fixed is it refitted on all three folds and evaluated ONCE on h1-test.

`reach` IS THE OBJECTIVE, AND IT IS NOT A RANK CORRELATION.  The brief's calibration is the
governing fact: `|H|` has pool Spearman +0.069 with sign-correctness and reaches 0.235, while a
synthetic key needs Spearman 0.316 to reach only 0.224.  `reach` reads the FRONT of the
ordering and nothing else -- `k*` is the deepest prefix at >=90% purity, and at a front accuracy
near 0.6 that prefix is a handful of genes.  So the models below are graded by fold `reach` and
by front precision, never by AUC or Spearman, and the top-heavy arms exist because a plain
logistic fit spends its capacity on the 99% of rows `reach` never reads.

MODEL CLASSES (constraint 5: small model first, then bigger)
------------------------------------------------------------
  M1  plain logistic on the 21 standardised features                    22 params (the BASELINE)
  M2  additive step-basis logistic, 24 quantile bins per feature       ~484 params
      -- AgingLab's scale (they publish "a 578-parameter head"), and inspectable: each
         feature's contribution is a monotone-ish step curve that can be printed.
  M3  M2 refitted with a top-heavy weighting (three variants, see TOPHEAVY)
  M4  HistGradientBoostingClassifier                                   ~10^4 params
  M5  small MLP (32, 16)                                              ~1.4x10^3 params

THE STEP BASIS.  Feature `f` is quantile-binned into `B` bins and encoded as
`Z_k = 1[bin(x) >= k]`, so the fitted level at bin `j` is the running sum of the coefficients.
An L2 penalty on those coefficients is therefore a penalty on the FIRST DIFFERENCES of the
response curve -- a smoothness prior, not a shrink-to-zero prior on the level.  That is what
makes a few-hundred-parameter binned model behave like a spline instead of like a lookup table,
and it is why the same code can be read as "what did the model learn about expression".

TIES (097 §10).  Bin edges are quantiles computed with midranks; a feature whose mass is
concentrated on one value (`concord`, `n_don`, `rel_dis` at 1.0) collapses to fewer distinct
bins rather than being handed an arbitrary ramp.  Duplicate edges are dropped explicitly.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
from scipy.stats import rankdata, spearmanr
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
WORK = os.path.join(HERE, "work")

from feats172 import FEATURES, Recipient, build, stack   # noqa: E402
from harness172 import H1Test, check                     # noqa: E402

PROD = ("k562_gwps", "orion_hct116", "orion_hek293t")
SOURCES_H1 = ("k562_gwps", "orion_hek293t")
NBINS = 24
PURITY = 0.90


# ---- the metric, on a fold ------------------------------------------------------------

def fold_reach(score: np.ndarray, t: np.ndarray, y: np.ndarray,
               rng: np.random.Generator | None = None) -> float:
    """`debench.Surface.reach`'s rule, applied to a fold's pool rows.

    R_p = k*/n_p with k* the DEEPEST prefix at >=90% purity ("deepest" is literal: purity is
    not monotone, so a prefix that dips and recovers counts at its deeper crossing).  Ties in
    `score` are broken by a fixed random permutation, which is the fold analogue of debench's
    gene-name lexsort: deterministic, and identical for every key so it cannot favour one.
    """
    rng = rng or np.random.default_rng(12345)
    jitter = rng.permutation(score.size)
    order = np.lexsort((jitter, -score))
    ts, ys = t[order], y[order]
    out = []
    for p in np.unique(ts):
        m = ts == p
        ok = ys[m].astype(np.float64)
        pur = np.cumsum(ok) / np.arange(1, ok.size + 1)
        deep = np.flatnonzero(pur >= PURITY)
        out.append((int(deep[-1]) + 1) / ok.size if deep.size else 0.0)
    return float(np.mean(out))


def fold_front(score: np.ndarray, t: np.ndarray, y: np.ndarray, k: int = 5) -> float:
    """Mean over targets of the accuracy of the top-`k` genes -- what `reach` actually reads."""
    order = np.lexsort((np.random.default_rng(12345).permutation(score.size), -score))
    ts, ys = t[order], y[order]
    out = []
    for p in np.unique(ts):
        ok = ys[ts == p]
        out.append(float(ok[:min(k, ok.size)].mean()))
    return float(np.mean(out))


# ---- the step basis -------------------------------------------------------------------

class StepBasis:
    def __init__(self, X: np.ndarray, nbins: int = NBINS):
        self.edges = []
        for j in range(X.shape[1]):
            q = np.quantile(X[:, j], np.linspace(0, 1, nbins + 1)[1:-1])
            self.edges.append(np.unique(q))
        self.n = int(sum(len(e) for e in self.edges))

    def __call__(self, X: np.ndarray) -> np.ndarray:
        cols = [np.greater.outer(X[:, j], self.edges[j]) for j in range(X.shape[1])
                if len(self.edges[j])]
        return np.concatenate(cols, axis=1).astype(np.float32)

    def blocks(self):
        out, o = [], 0
        for j, e in enumerate(self.edges):
            if len(e):
                out.append((FEATURES[j], o, o + len(e), e))
                o += len(e)
        return out


# ---- model classes --------------------------------------------------------------------

def _std(X, mu, sd):
    return (X - mu) / sd


class Model:
    """Fit on (X, y[, w]); predict a confidence in [0, 1]."""

    def __init__(self, kind: str, **kw):
        self.kind, self.kw = kind, kw

    def fit(self, X, y, w=None, t=None):
        self.mu = X.mean(0)
        self.sd = np.maximum(X.std(0), 1e-6)
        if self.kind == "logistic":
            self.clf = LogisticRegression(C=self.kw.get("C", 1.0), max_iter=2000)
            self.clf.fit(_std(X, self.mu, self.sd), y, sample_weight=w)
        elif self.kind == "step":
            self.basis = StepBasis(X, self.kw.get("nbins", NBINS))
            self.clf = LogisticRegression(C=self.kw.get("C", 0.05), max_iter=4000)
            self.clf.fit(self.basis(X), y, sample_weight=w)
        elif self.kind == "gbt":
            self.clf = HistGradientBoostingClassifier(
                max_iter=self.kw.get("max_iter", 200),
                learning_rate=self.kw.get("lr", 0.06),
                max_leaf_nodes=self.kw.get("leaves", 15),
                min_samples_leaf=self.kw.get("msl", 200),
                l2_regularization=self.kw.get("l2", 1.0), random_state=0)
            self.clf.fit(X, y, sample_weight=w)
        elif self.kind == "mlp":
            self.clf = MLPClassifier(hidden_layer_sizes=self.kw.get("hidden", (32, 16)),
                                     alpha=self.kw.get("alpha", 1e-2), max_iter=self.kw.get(
                                         "max_iter", 60), random_state=0, early_stopping=True)
            self.clf.fit(_std(X, self.mu, self.sd), y)     # MLPClassifier takes no weights
        else:
            raise ValueError(self.kind)
        return self

    def predict(self, X):
        if self.kind == "step":
            return self.clf.predict_proba(self.basis(X))[:, 1]
        if self.kind == "gbt":
            return self.clf.predict_proba(X)[:, 1]
        return self.clf.predict_proba(_std(X, self.mu, self.sd))[:, 1]

    @property
    def n_params(self) -> int:
        if self.kind == "gbt":
            return int(sum(len(p[0].nodes) for p in self.clf._predictors))
        if self.kind == "mlp":
            return int(sum(c.size for c in self.clf.coefs_) +
                       sum(b.size for b in self.clf.intercepts_))
        return int(self.clf.coef_.size + 1)


# ---- top-heavy weightings (constraint 4) ----------------------------------------------

def topheavy_weights(kind: str, X: np.ndarray, t: np.ndarray, score: np.ndarray | None):
    """Sample weights that concentrate the loss where `reach` actually reads."""
    if kind == "none":
        return None
    r = np.empty(score.size)
    for p in np.unique(t):                       # midrank within target, 1 = front
        m = t == p
        r[m] = rankdata(-score[m], method="average")
    if kind == "top25":                          # hard cut: the top quarter of each target
        n = np.zeros(t.size)
        for p in np.unique(t):
            n[t == p] = (t == p).sum()
        return (r <= np.maximum(1, 0.25 * n)).astype(np.float64)
    if kind == "invrank":                        # 1/rank: reach reads a handful of genes
        return 1.0 / r
    if kind == "exp10":                          # exp(-rank/10): sharper still
        return np.exp(-r / 10.0)
    raise ValueError(kind)


TOPHEAVY = ("none", "top25", "invrank", "exp10")


# ---- data -----------------------------------------------------------------------------

def load_folds() -> dict:
    F = {}
    for h in PROD:
        z = np.load(os.path.join(WORK, f"fold.{h}.npz"), allow_pickle=True)
        assert list(np.asarray(z["features"]).astype(str)) == FEATURES, "feature order drift"
        F[h] = dict(X=z["X"], y=z["y"], t=z["t"])
    return F


def h1_features(T: H1Test):
    rec = Recipient("h1test", T.genes, T.b_ctrl, T.perts, T.gate.copy())
    B = build(rec, SOURCES_H1, verbose=False)
    Tn, G = len(T.perts), len(T.genes)
    rt, rg = np.divmod(np.arange(Tn * G), G)
    return stack(B["F"], rt, rg), B


def key_from(model: Model, Xh: np.ndarray, shape) -> np.ndarray:
    return model.predict(Xh).reshape(shape)
