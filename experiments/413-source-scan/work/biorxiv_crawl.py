"""413: crawl bioRxiv details API over a date window and keep preprints whose title/abstract match
perturbation-transcriptomics keywords. Output biorxiv_hits_<window>.json."""
import json, re, sys, time, urllib.request
start, end = sys.argv[1], sys.argv[2]
KW = re.compile(r"perturb-?seq|crop-?seq|crispri|crispra|crispr interference|crispr activation|single-cell crispr|"
                r"pooled crispr|genome-(wide|scale) (crispr|perturb)|perturbation atlas|perturbation screen|"
                r"transcriptom\w* (screen|readout)|l1000|virtual cell", re.I)
hits, cursor, total = [], 0, None
while True:
    url = f"https://api.biorxiv.org/details/biorxiv/{start}/{end}/{cursor}"
    for attempt in range(5):
        try:
            d = json.load(urllib.request.urlopen(url, timeout=120)); break
        except Exception as e:
            time.sleep(5 * (attempt + 1))
    else:
        print("FAILED at", cursor, file=sys.stderr); break
    msg = d["messages"][0]
    total = int(msg.get("total", msg.get("count", 0)))
    col = d.get("collection", [])
    for p in col:
        txt = (p.get("title", "") + " " + p.get("abstract", ""))
        if KW.search(txt):
            hits.append({k: p.get(k) for k in ["doi", "title", "date", "version", "category", "authors", "abstract"]})
    cursor += len(col)
    if cursor % 2000 < 100: print(cursor, "/", total, "hits", len(hits), file=sys.stderr, flush=True)
    if not col or cursor >= total: break
json.dump(hits, open(f"biorxiv_hits_{start}_{end}.json", "w"), indent=1)
print("done", cursor, total, len(hits), file=sys.stderr)
