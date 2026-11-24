"""Shared machinery for experiment 093: does the emission `spread` channel truncate
our significant call set?

Experiment 092 established the *counterfactual*: keeping only the top ~850 genes of our
own call list raises `fid` 0.5313 -> 0.5906 (+0.0299 on avg6) with `reach` and `nmae`
exactly invariant. 092 achieved that by fiat -- it relabelled significance in a table.
Here the truncation has to be produced by the only channel that can produce it for real:
`emit_counts(mode="quasi")`'s per-gene `spread`.

Key structural fact that makes this cheap. In the quasi emitter the ONLY thing carried
from the profile into the emitted block is the per-gene panel total

    totals_g = floor(n * L * profile_p_g)

and that integer is preserved EXACTLY by the emitter whatever `spread` is. So a
re-emission needs nothing but `totals`, which we read straight off the existing
prediction as its per-group column sums. Feeding
`profile_p := (totals + 0.5) / (n * L)` back into the real `emit_counts` recovers the same
`totals` bit-exactly (floor(totals + 0.5) == totals for integers), so we drive the
PRODUCTION function unmodified rather than re-implementing its quasi branch.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import h5py
from scipy import stats
from scipy.stats import false_discovery_control

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "experiments", "092-call-set-truncation"))

from vcc2026.donor_transfer import emit_counts  # noqa: E402
from de_kernel import ctrl_tiesums, wilcoxon_z  # noqa: E402

SCRATCH = ("/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/"
           "76a42756-c5fb-4398-aea0-8513243571e9/scratchpad")
E092 = os.path.join(ROOT, "experiments", "092-call-set-truncation")
REAL = os.path.join(ROOT, "data", "harness", "orion-hct116.h5ad")
PRED = os.path.join(ROOT, "experiments", "031-somatic", "out_som_gamma", "runs",
                    "g0.65_q1.0_a1.08_spread0.25", "pred.h5ad")
CTRL, PERT_COL, EPS, ALPHA = "non-targeting", "target_gene", 1e-9, 0.05
N_EMIT = 400
PROD_SPREAD = 0.25
# `L` only enters the quasi emitter through `totals`; any value that keeps
# (totals + 0.5) / (n * L) in a well-conditioned float64 range reproduces `totals`.
EMIT_L = 100_000

# Official 2026 anchors (baseline, replicate), AGENTS.md section 4 midpoints.
ANCHORS = {
    "fid": (0.5135, 0.8135),
    "jac": (0.029, 0.399),
    "reach": (0.072, 0.968),
    "nmae": (1.0013, 0.400),
}


def scaled(name: str, raw: float) -> float:
    b, r = ANCHORS[name]
    return (raw - b) / (r - b)


# --------------------------------------------------------------------------- gate / DE

class Gate:
    """The reference-only CPM gate plus the sorted real-control block."""

    def __init__(self) -> None:
        g = np.load(os.path.join(E092, "gate.npz"), allow_pickle=True)
        self.keep = g["keep"]
        self.var = g["var"].astype(str)
        self.perts = g["perts"].astype(str)
        self.ctrl_mean_cpm = g["ctrl_mean_cpm"]
        self.keep_idx = np.flatnonzero(self.keep).astype(np.int64)
        self.genes = self.var[self.keep]
        self.ref_mean = self.ctrl_mean_cpm[self.keep]
        self.G = self.genes.size
        self.G_full = self.var.size
        self._ctrl_sorted = None
        self._ctrl_tie = None

    def ctrl(self):
        """The sorted control block (mmapped, so N workers share ONE resident copy) and
        its per-gene tie sums (cached, because computing them needs the block dense)."""
        if self._ctrl_sorted is None:
            self._ctrl_sorted = np.load(os.path.join(SCRATCH, "ctrl_log_sorted.npy"),
                                        mmap_mode="r")
            tie_path = os.path.join(SCRATCH, "ctrl_tiesums093.npy")
            if os.path.exists(tie_path):
                self._ctrl_tie = np.load(tie_path)
            else:
                self._ctrl_tie = ctrl_tiesums(np.ascontiguousarray(self._ctrl_sorted))
                np.save(tie_path, self._ctrl_tie)
        return self._ctrl_sorted, self._ctrl_tie

    def de(self, dense_counts: np.ndarray):
        """Raw counts (n, G_full) -> (p_value, p_adj, lfc) over gated genes.

        Byte-for-byte the path in experiments/092/step2_de.py, which reproduces
        cell-eval2 0.16.0 to 2.5e-6 on the mean metric.
        """
        ctrl_sorted, ctrl_tie = self.ctrl()
        L = dense_counts.sum(axis=1)
        scale = L / 1e6
        scale = scale + (scale == 0)
        cpm = dense_counts / scale[:, None]
        sub = cpm[:, self.keep_idx]
        mean_cpm = sub.mean(axis=0)
        lg = np.ascontiguousarray(np.log1p(sub).T)
        z = wilcoxon_z(ctrl_sorted, ctrl_tie, lg)
        pval = 2 * stats.distributions.norm.sf(np.abs(z))
        lfc = np.log2((mean_cpm + EPS) / (self.ref_mean + EPS))
        padj = false_discovery_control(pval, method="bh")
        padj = np.where(np.isnan(lfc), 1.0, padj)
        return pval, padj, lfc


# --------------------------------------------------------------------------- emission

def totals_to_profile(totals: np.ndarray, n: int = N_EMIT, L: int = EMIT_L) -> np.ndarray:
    """A profile that the real `emit_counts` turns back into exactly `totals`."""
    return (np.asarray(totals, dtype=np.float64) + 0.5) / (n * L)


def emit_block(totals: np.ndarray, spread: np.ndarray, rng,
               n: int = N_EMIT, L: int = EMIT_L) -> np.ndarray:
    """Re-emit one perturbation's cells at a per-gene `spread`. Returns dense (n, G)."""
    values, cols, indptr = emit_counts(
        totals_to_profile(totals, n, L),
        np.full(n, float(L)),
        0.0,
        rng,
        mode="quasi",
        spread=np.asarray(spread, dtype=np.float64),
    )
    out = np.zeros((n, totals.size), dtype=np.float64)
    for i in range(n):
        a, b = indptr[i], indptr[i + 1]
        out[i, cols[a:b]] = values[a:b]
    return out


# --------------------------------------------------------------------------- pred I/O

def pred_labels(path: str = PRED) -> np.ndarray:
    with h5py.File(path, "r") as f:
        o = f["obs"][PERT_COL]
        cats = np.array([x.decode() if isinstance(x, bytes) else str(x)
                         for x in o["categories"][:]])
        return cats[o["codes"][:]]


def read_rows(path: str, rows: np.ndarray) -> np.ndarray:
    """CSR rows -> dense (len(rows), n_genes) float64. `rows` must be sorted."""
    with h5py.File(path, "r") as f:
        Xg = f["X"]
        indptr = Xg["indptr"]
        nvar = int(Xg.attrs["shape"][1])
        lo, hi = int(rows[0]), int(rows[-1]) + 1
        ip = indptr[lo:hi + 1]
        s, e = int(ip[0]), int(ip[-1])
        data = Xg["data"][s:e].astype(np.float64)
        idx = Xg["indices"][s:e]
        ip = ip - s
        out = np.zeros((len(rows), nvar), dtype=np.float64)
        rowmap = {r: i for i, r in enumerate(rows.tolist())}
        for k in range(hi - lo):
            i = rowmap.get(lo + k)
            if i is None:
                continue
            a, b = ip[k], ip[k + 1]
            out[i, idx[a:b]] = data[a:b]
    return out


def build_totals_cache(out_path: str) -> None:
    """Per-perturbation emitted panel totals, read off the production prediction."""
    gate = Gate()
    lab = pred_labels()
    tot = np.zeros((gate.perts.size, gate.G_full), dtype=np.int64)
    for i, p in enumerate(gate.perts):
        rows = np.flatnonzero(lab == p)
        assert rows.size == N_EMIT, (p, rows.size)
        block = read_rows(PRED, rows)
        col = block.sum(axis=0)
        assert np.all(col == np.round(col)), p
        tot[i] = col.astype(np.int64)
        if i % 20 == 0:
            print(f"  totals {i + 1}/{gate.perts.size} {p}", flush=True)
    np.save(out_path, tot)


# --------------------------------------------------------------------------- ranking

def rank_key(totals: np.ndarray, gate: Gate) -> np.ndarray:
    """|predicted effect| on the gated axis, computed from the emitted totals alone.

    `log2(pred_cpm / ctrl_cpm)` is the emitted prediction's own log-fold change: a
    monotone function of the model's per-gene effect (`profile = expm1(b_ctrl + eff)`),
    computable inside the submission builder before a single cell is emitted, and --
    because the quasi emitter preserves `totals` exactly -- completely invariant to
    whatever `spread` we then choose. That last property is what makes it usable as the
    ranking key for the spread policy without circularity.
    """
    t = np.asarray(totals, dtype=np.float64)
    cpm = 1e6 * t / t.sum()
    return np.abs(np.log2((cpm[gate.keep_idx] + EPS) / (gate.ref_mean + EPS)))
