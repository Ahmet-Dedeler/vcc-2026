"""413: Hugging Face dataset scan. Lists datasets matching perturbation keywords, keeps those modified since 2026-06-01."""
import json, urllib.request, urllib.parse, time
TERMS = ["perturb", "perturb-seq", "perturbseq", "crispri", "crispr", "crop-seq", "virtual cell", "virtual-cell",
         "single cell perturbation", "genetic perturbation", "knockdown", "knockout", "l1000", "lincs", "xpr",
         "scperturb", "replogle", "gwps", "genome-wide", "billion", "atlas", "hela", "jurkat", "cal33", "cancer cell line",
         "vcc", "tahoe", "x-atlas", "pisces", "perturbation atlas", "sc-perturb", "screen"]
seen = {}
for t in TERMS:
    url = "https://huggingface.co/api/datasets?" + urllib.parse.urlencode({"search": t, "limit": 500, "sort": "lastModified", "direction": -1, "full": "false"})
    try:
        data = json.load(urllib.request.urlopen(url, timeout=60))
    except Exception as e:
        print("ERR", t, e); continue
    for d in data:
        seen.setdefault(d["id"], {"lastModified": d.get("lastModified"), "createdAt": d.get("createdAt"), "downloads": d.get("downloads"), "terms": []})["terms"].append(t)
    print(t, len(data))
    time.sleep(0.3)
json.dump(seen, open("hf_scan.json", "w"), indent=1)
recent = {k: v for k, v in seen.items() if (v["lastModified"] or "") >= "2026-06-01"}
print("total", len(seen), "recent", len(recent))
for k, v in sorted(recent.items(), key=lambda kv: kv[1]["lastModified"], reverse=True):
    print(v["lastModified"][:10], k, v["terms"][:4])
