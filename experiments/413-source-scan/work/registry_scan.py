"""413: figshare / Zenodo / BioStudies keyword scans for recent perturbation-transcriptomics deposits."""
import json, urllib.request, urllib.parse, time
UA = {"User-Agent": "Mozilla/5.0 (vcc-413 research scan)"}
TERMS = ["Perturb-seq", "CROP-seq", "CRISPRi single-cell", "single-cell CRISPR screen", "genome-wide Perturb-seq", "CRISPR screen scRNA-seq", "perturbation atlas"]
out = {"figshare": {}, "zenodo": {}, "biostudies": {}}
for t in TERMS:
    # figshare
    try:
        body = json.dumps({"search_for": t, "published_since": "2026-05-01", "page_size": 100, "item_type": 3}).encode()
        req = urllib.request.Request("https://api.figshare.com/v2/articles/search", data=body, headers={**UA, "Content-Type": "application/json"})
        for a in json.load(urllib.request.urlopen(req, timeout=60)):
            out["figshare"][a["id"]] = {"title": a["title"], "date": a.get("published_date"), "url": a.get("url_public_html"), "doi": a.get("doi")}
    except Exception as e: print("figshare ERR", t, e)
    # zenodo
    try:
        url = "https://zenodo.org/api/records?" + urllib.parse.urlencode({"q": f'"{t}"', "sort": "newest", "size": 50, "type": "dataset"})
        d = json.load(urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60))
        for h in d["hits"]["hits"]:
            if h.get("created", "") >= "2026-05-01":
                out["zenodo"][h["id"]] = {"title": h["metadata"]["title"], "date": h["created"][:10], "url": h["links"].get("self_html"),
                                         "size_GB": round(sum(f["size"] for f in h.get("files", [])) / 1e9, 2)}
    except Exception as e: print("zenodo ERR", t, e)
    # biostudies / arrayexpress
    try:
        url = "https://www.ebi.ac.uk/biostudies/api/v1/search?" + urllib.parse.urlencode({"query": t, "sortBy": "release_date", "sortOrder": "descending", "pageSize": 50})
        d = json.load(urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60))
        for h in d.get("hits", []):
            if (h.get("release_date") or "") >= "2026-05-01":
                out["biostudies"][h["accession"]] = {"title": h.get("title"), "date": h.get("release_date"), "type": h.get("type")}
    except Exception as e: print("biostudies ERR", t, e)
    time.sleep(1)
json.dump(out, open("registry_scan.json", "w"), indent=1)
for reg, d in out.items():
    print("=====", reg, len(d))
    for k, v in sorted(d.items(), key=lambda kv: kv[1].get("date") or "", reverse=True):
        print(" ", (v.get("date") or "")[:10], k, (v.get("title") or "")[:130], v.get("size_GB", ""))
