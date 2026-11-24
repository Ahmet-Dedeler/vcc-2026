"""517 — memory-bounded streaming of CSR h5ad files: group sums of raw counts / per-cell CPM (optionally log1p'd) by a
per-cell group label, with every per-chunk array preallocated (read_direct + np.take/np.multiply with out=, bincount).
macOS did not reuse freed, variably sized per-chunk allocations in the first versions of these scripts (6.7 GB RSS for a
16 s job); fully preallocated, the RSS stays at the buffers' size.

    it = CSRStream(path)                                  # opens the file, reads indptr, allocates the buffers
    lib = it.libraries()                                  # per-cell total counts (float64)
    S = it.group_sums(labels, n_groups, cols=None, cpm=False, lib=None, log1p=False)   # (n_groups, n_cols) float64
labels: per-cell int group (< 0 = skip); cols: gene -> output column (-1 = drop); binary: count detections.
"""
from __future__ import annotations

import h5py
import numpy as np


class CSRStream:
    def __init__(self, path, chunk_rows=2000):
        self.f = h5py.File(path, "r")
        self.ipd = self.f["X/indptr"][:].astype(np.int64)
        self.n = self.ipd.size - 1
        self.G = int(self.f["X"].attrs["shape"][1])
        self.CH = chunk_rows
        mx = int(max(self.ipd[min(s + chunk_rows, self.n)] - self.ipd[s] for s in range(0, self.n, chunk_rows)))
        self.mx = mx
        self.d32 = np.empty(mx, self.f["X/data"].dtype)
        self.i32 = np.empty(mx, self.f["X/indices"].dtype)
        self.w = np.empty(mx, np.float64)
        self.lab = np.empty(mx, np.int64)
        self.col = np.empty(mx, np.int64)

    def close(self):
        self.f.close()

    def _chunks(self):
        for s in range(0, self.n, self.CH):
            e = min(self.n, s + self.CH)
            a, b = int(self.ipd[s]), int(self.ipd[e])
            m = b - a
            self.f["X/data"].read_direct(self.d32, np.s_[a:b], np.s_[0:m])
            self.f["X/indices"].read_direct(self.i32, np.s_[a:b], np.s_[0:m])
            self.w[:m] = self.d32[:m]
            yield s, e, m, self.ipd[s:e + 1] - a

    def libraries(self):
        out = np.zeros(self.n)
        for s, e, m, ip in self._chunks():
            for r in range(e - s):
                out[s + r] = self.w[ip[r]:ip[r + 1]].sum()
        return out

    def group_sums(self, labels, n_groups, cols=None, cpm=False, lib=None, log1p=False, binary=False):
        if cpm and lib is None:
            lib = self.libraries()
        nc = self.G if cols is None else int(cols.max()) + 1
        colmap = (np.arange(self.G, dtype=np.int64) if cols is None else np.where(cols >= 0, cols, nc).astype(np.int64))
        out = np.zeros((n_groups + 1) * (nc + 1))
        for s, e, m, ip in self._chunks():
            for r in range(e - s):
                a_, b_ = ip[r], ip[r + 1]
                lb = int(labels[s + r])
                self.lab[a_:b_] = lb if lb >= 0 else n_groups
                if cpm:
                    self.w[a_:b_] *= 1e6 / max(float(lib[s + r]), 1.0)
            if log1p:
                np.log1p(self.w[:m], out=self.w[:m])
            if binary:                        # count of nonzero entries (detections) instead of their values
                self.w[:m] = 1.0
            np.take(colmap, self.i32[:m], out=self.col[:m], mode="clip")   # mode raise would buffer out
            np.multiply(self.lab[:m], nc + 1, out=self.lab[:m])
            np.add(self.lab[:m], self.col[:m], out=self.lab[:m])
            out += np.bincount(self.lab[:m], weights=self.w[:m], minlength=out.size)
        return out.reshape(n_groups + 1, nc + 1)[:n_groups, :nc]
