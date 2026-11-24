"""413: DLD-1 GSE337988 sublib2 — panel overlap of the held Low1/Low2 arms vs the unheld Med/High arms (header only, range-read)."""
import gzip, json, re, urllib.request, zlib
from overlap import report, ROOT
def header_local(p):
    with gzip.open(p, "rt") as fh: return fh.readline().rstrip("\n").split(",")[1:]
def header_remote(url, nbytes=2_000_000):
    req = urllib.request.Request(url, headers={"Range": f"bytes=0-{nbytes}"})
    raw = urllib.request.urlopen(req, timeout=120).read()
    txt = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(raw).decode()
    return txt.split("\n", 1)[0].split(",")[1:]
tgt = lambda cols: sorted({re.split(r"_(P\d+|ENST)", c)[0] for c in cols})
out = {}
for arm in ["Low1", "Low2"]:
    cols = header_local(ROOT / f"data/external/dld1-lowmoi/GSE337988_sublib2_de_matrices_lfc_matrix_{arm}.csv.gz")
    out[arm] = report(f"DLD-1 {arm} (held)", tgt(cols)) | {"n_columns": len(cols)}
for arm in ["Med", "High"]:
    cols = header_remote(f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE337nnn/GSE337988/suppl/GSE337988_sublib2_de_matrices_lfc_matrix_{arm}.csv.gz")
    out[arm] = report(f"DLD-1 {arm} (unheld)", tgt(cols)) | {"n_columns": len(cols)}
for k, v in out.items(): print(k, v["n_columns"], v["n_unique"], v["panel_exact"], v["panel_resolved"])
print("Med/High panel hits not in Low1:", sorted(set(out["Med"]["panel_hits"]) | set(out["High"]["panel_hits"]) - set(out["Low1"]["panel_hits"])))
json.dump(out, open("dld1_arms_overlap.json", "w"), indent=1)
