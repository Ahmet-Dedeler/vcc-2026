"""413: inspect candidate HF datasets: file list (sizes) + card head."""
import json, sys, urllib.request
ids = sys.argv[1:]
for i in ids:
    print("=" * 100); print(i)
    try:
        d = json.load(urllib.request.urlopen(f"https://huggingface.co/api/datasets/{i}?blobs=true", timeout=60))
    except Exception as e:
        print("ERR", e); continue
    print("lastModified", d.get("lastModified"), "created", d.get("createdAt"), "gated", d.get("gated"), "license", [t for t in d.get("tags", []) if t.startswith("license")])
    sib = d.get("siblings", [])
    tot = sum(s.get("size") or 0 for s in sib)
    print(f"{len(sib)} files, {tot/1e9:.2f} GB")
    for s in sib[:25]:
        print(f"   {s.get('size',0)/1e6:10.1f} MB  {s['rfilename']}")
    try:
        card = urllib.request.urlopen(f"https://huggingface.co/datasets/{i}/raw/main/README.md", timeout=60).read().decode()[:1800]
        print(card)
    except Exception as e:
        print("no card", e)
