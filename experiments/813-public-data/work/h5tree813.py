"""813: list datasets (name, shape, dtype) of a remote HDF5 via range reads, byte-capped. Usage: uv run python h5tree813.py <url> <out.json> <cap_mb>"""
import json, sys, h5py
src = open(__file__.replace("h5tree813.py", "h5peek813.py")).read().split("rf = RangeFile(URL)")[0]
ns = {"__name__": "x"}; sys.argv = sys.argv[:4]; exec(src, ns)
rf = ns["RangeFile"](sys.argv[1]); out = {}
with h5py.File(rf, "r") as f:
    def v(name, o):
        if isinstance(o, h5py.Dataset): out[name] = [list(o.shape), str(o.dtype)]
    f.visititems(v)
out["_bytes_fetched"] = rf.fetched; out["_size"] = rf.size
json.dump(out, open(sys.argv[2], "w"), indent=1); print(json.dumps(out, indent=1)[:3000])
