import json, urllib.request
for q in ["DepMap Public", "DepMap 25Q", "DepMap 26Q1", "Cancer Dependency Map DepMap Public"]:
    body = json.dumps({"search_for": q, "page_size": 50, "order": "published_date",
                       "order_direction": "desc"}).encode()
    req = urllib.request.Request("https://api.figshare.com/v2/articles/search", data=body,
                                 headers={"Content-Type": "application/json"})
    d = json.load(urllib.request.urlopen(req, timeout=60))
    print("QUERY:", q, "->", len(d), "hits")
    for a in d:
        t = a.get("title", "")
        tl = t.lower()
        if "depmap" in tl and ("public" in tl or "q1" in tl or "q2" in tl or "q3" in tl or "q4" in tl):
            print("  ", a.get("id"), "|", t, "|", a.get("published_date"))
