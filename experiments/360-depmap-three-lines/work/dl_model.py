import json, urllib.request, hashlib
d = json.load(urllib.request.urlopen("https://api.figshare.com/v2/articles/27993248", timeout=60))
urls = {f["name"]: f["download_url"] for f in d["files"]}
json.dump(urls, open("~/Code/virtual-cell-challenge-2026/experiments/360-depmap-three-lines/work/fs_urls.json","w"), indent=1)
print("Model.csv:", urls["Model.csv"][:120])
data = urllib.request.urlopen(urls["Model.csv"], timeout=120).read()
open("~/Code/virtual-cell-challenge-2026/experiments/360-depmap-three-lines/work/Model.csv","wb").write(data)
print("bytes:", len(data), "sha256:", hashlib.sha256(data).hexdigest())
