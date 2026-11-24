"""413: parallel version of biorxiv_crawl.py (same keyword filter), cursor range split across workers."""
import json, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor
from biorxiv_kw import KW
start, end, total, nw = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
def fetch(cursor):
    url = f"https://api.biorxiv.org/details/biorxiv/{start}/{end}/{cursor}"
    for attempt in range(6):
        try: return json.load(urllib.request.urlopen(url, timeout=180)).get("collection", [])
        except Exception: time.sleep(5 * (attempt + 1))
    print("FAILED", cursor, file=sys.stderr, flush=True); return None
def worker(cursors):
    hits, n, failed = [], 0, []
    for c in cursors:
        col = fetch(c)
        if col is None: failed.append(c); continue
        n += len(col)
        hits += [{k: p.get(k) for k in ["doi", "title", "date", "version", "category", "authors", "abstract"]}
                 for p in col if KW.search(p.get("title", "") + " " + p.get("abstract", ""))]
    return hits, n, failed
cursors = list(range(0, total, 30))
chunks = [cursors[i::nw] for i in range(nw)]
with ThreadPoolExecutor(nw) as ex: res = list(ex.map(worker, chunks))
hits = [h for r in res for h in r[0]]; n = sum(r[1] for r in res); failed = sorted(c for r in res for c in r[2])
json.dump({"n_scanned": n, "failed_cursors": failed, "hits": hits}, open(f"biorxiv_hits_{start}_{end}.json", "w"), indent=1)
print("done scanned", n, "hits", len(hits), "failed", failed, file=sys.stderr)
