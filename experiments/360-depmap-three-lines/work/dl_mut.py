import json, urllib.request, time
W = "~/Code/virtual-cell-challenge-2026/experiments/360-depmap-three-lines/work/"
urls = json.load(open(W+"fs_urls.json"))
KEEP = {b"ACH-000995",b"ACH-001086",b"ACH-000518"}
t0=time.time()
resp = urllib.request.urlopen(urls["OmicsSomaticMutations.csv"], timeout=300)
buf=b""
while buf.count(b"\n")<2: buf += resp.read(1<<16)
first, second, rest = buf.split(b"\n",2)
print("header:", first.decode()[:400])
print("row1:", second.decode()[:300])
# find ModelID column index
hcols = first.decode().split(",")
print("n cols:", len(hcols))
for i,c in enumerate(hcols):
    if "model" in c.lower() or "depmap" in c.lower() or "line" in c.lower(): print(" idcol?", i, c)
