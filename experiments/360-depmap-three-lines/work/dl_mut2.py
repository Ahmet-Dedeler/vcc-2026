import json, urllib.request, time, hashlib
W = "~/Code/virtual-cell-challenge-2026/experiments/360-depmap-three-lines/work/"
urls = json.load(open(W+"fs_urls.json"))
KEEP = {b"ACH-000995",b"ACH-001086",b"ACH-000518"}
t0=time.time()
resp = urllib.request.urlopen(urls["OmicsSomaticMutations.csv"], timeout=300)
out = open(W+"Mutations_3lines.csv","wb")
nrows=nkept=0
carry=b""
first=True
nbytes=0
while True:
    blk = resp.read(1<<20)
    if not blk: break
    nbytes+=len(blk)
    carry+=blk
    *lines, carry = carry.split(b"\n")
    for ln in lines:
        if first: out.write(ln+b"\n"); first=False; continue
        if not ln: continue
        nrows+=1
        # ModelID is col 67; cheap check: line contains one of the ids
        if ln.find(b"ACH-000995")>0 or ln.find(b"ACH-001086")>0 or ln.find(b"ACH-000518")>0:
            nkept+=1; out.write(ln+b"\n")
if carry.strip(): 
    nrows+=1
    if carry.find(b"ACH-000995")>0 or carry.find(b"ACH-001086")>0 or carry.find(b"ACH-000518")>0:
        nkept+=1; out.write(carry+b"\n")
out.close()
print(f"scanned {nrows} mut rows ({nbytes/1e6:.0f} MB) in {time.time()-t0:.0f}s, kept {nkept}")
