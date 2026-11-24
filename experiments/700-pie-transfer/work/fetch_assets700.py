"""700: resolve every asset pinned by the PIE xdataset checkpoint config into PIE_DATA_ROOT.

Run with the 700 venv: PIE_DATA_ROOT=... .venv/bin/python work/fetch_assets700.py work/xd_config.yaml
Order: small assets first, the two large ones (Tahoe preprocessed, ncbi_text) last.
"""
import sys, time, yaml, logging
from pie.assets import resolve_asset

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
cfg = yaml.safe_load(open(sys.argv[1]))
data = cfg["data"]
jobs = []
for p in data["preprocessed_dirs"]:
    jobs.append((p, "preprocessed"))
for name, p in data["source_dirs"].items():
    jobs.append((p, "source"))
jobs.append((data["gene_text_dir"], "source"))
big = lambda j: ("tahoe" in j[0]) or ("ncbi_text" in j[0])
jobs = [j for j in jobs if not big(j)] + [j for j in jobs if big(j)]
for ref, kind in jobs:
    t = time.time()
    path = resolve_asset(ref, kind=kind)
    print(f"{time.strftime('%H:%M:%SZ', time.gmtime())} OK {kind} {ref} -> {path} ({time.time()-t:.0f}s)", flush=True)
print("ALL ASSETS RESOLVED", flush=True)
