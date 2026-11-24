"""813: read the HDF5 *metadata* of a remote h5ad via HTTP range requests (no full download).

Usage: uv run python -I h5peek813.py <portal-download-url> <out.json> [max_mb]
Prints/records: n_obs, n_var, obs column names (+ categories size for categoricals), layer/X storage
(format, shape, nnz from indptr tail), and for a guide-call column, if present, a few categories.
Total bytes fetched are capped (default 60 MB) so this stays a metadata read.
"""
import json
import sys
import urllib.request

import h5py

URL, OUT = sys.argv[1], sys.argv[2]
CAP = int(float(sys.argv[3]) * 2**20) if len(sys.argv) > 3 else 60 * 2**20
BLK = 1 << 20


def resolve(url):
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.geturl(), int(r.headers["Content-Length"])


class RangeFile:
    def __init__(self, url):
        self.url, self.size = resolve(url)
        self.pos, self.cache, self.fetched = 0, {}, 0

    def _block(self, i):
        if i not in self.cache:
            if self.fetched + BLK > CAP:
                raise RuntimeError(f"byte cap {CAP} reached")
            a, b = i * BLK, min(self.size, (i + 1) * BLK) - 1
            req = urllib.request.Request(self.url, headers={"Range": f"bytes={a}-{b}"})
            with urllib.request.urlopen(req, timeout=120) as r:
                self.cache[i] = r.read()
            self.fetched += len(self.cache[i])
        return self.cache[i]

    def read(self, n=-1):
        if n < 0:
            n = self.size - self.pos
        out = bytearray()
        while n > 0 and self.pos < self.size:
            i, off = divmod(self.pos, BLK)
            chunk = self._block(i)[off:off + n]
            out += chunk
            self.pos += len(chunk)
            n -= len(chunk)
        return bytes(out)

    def readinto(self, b):
        d = self.read(len(b))
        b[:len(d)] = d
        return len(d)

    def seek(self, p, whence=0):
        self.pos = p if whence == 0 else (self.pos + p if whence == 1 else self.size + p)
        return self.pos

    def tell(self):
        return self.pos

    def seekable(self):
        return True

    def readable(self):
        return True


rf = RangeFile(URL)
rec = {"url_base": URL, "size_bytes": rf.size}
with h5py.File(rf, "r") as f:
    rec["top"] = list(f.keys())
    obs = f["obs"]
    rec["obs_cols"] = {}
    for k in obs.keys():
        o = obs[k]
        if isinstance(o, h5py.Group):
            cats = o.get("categories")
            rec["obs_cols"][k] = {"kind": "categorical", "n_categories": int(cats.shape[0]) if cats is not None else None,
                                  "n": int(o["codes"].shape[0]) if "codes" in o else None}
        else:
            rec["obs_cols"][k] = {"kind": str(o.dtype), "n": int(o.shape[0])}
    idx = obs.attrs.get("_index", "_index")
    rec["n_obs"] = int(obs[idx].shape[0]) if idx in obs else None
    var = f["var"]
    vidx = var.attrs.get("_index", "_index")
    rec["n_var"] = int(var[vidx].shape[0]) if vidx in var else None
    rec["var_cols"] = list(var.keys())
    X = f["X"]
    if isinstance(X, h5py.Group):
        rec["X"] = {"enc": X.attrs.get("encoding-type"), "shape": [int(s) for s in X.attrs.get("shape", [])],
                    "data_dtype": str(X["data"].dtype), "nnz": int(X["data"].shape[0])}
    else:
        rec["X"] = {"dense_shape": list(X.shape), "dtype": str(X.dtype)}
    for grp in ("layers", "obsm", "uns", "obsp"):
        if grp in f:
            rec[grp] = list(f[grp].keys())[:40]
rec["bytes_fetched"] = rf.fetched
json.dump(rec, open(OUT, "w"), indent=1, default=str)
print(json.dumps(rec, indent=1, default=str)[:4000])
