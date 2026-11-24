import json, urllib.request, hashlib, csv, sys, time
W = "~/Code/virtual-cell-challenge-2026/experiments/360-depmap-three-lines/work/"
urls = json.load(open(W+"fs_urls.json"))
KEEP = {"ACH-000995","ACH-001086","ACH-000518","ACH-000551","ACH-000971","ACH-000739",
        "ACH-000997","ACH-000823","ACH-000865","ACH-000348","ACH-001652","ACH-000036",
        "ACH-001334","ACH-001992","ACH-002042","ACH-000250","ACH-000460","ACH-001852"}
def stream_keep_rows(fname, outname, id_col="ModelID"):
    url = urls[fname]
    t0=time.time()
    resp = urllib.request.urlopen(url, timeout=300)
    # read header
    buf = b""
    while b"\n" not in buf:
        buf += resp.read(1<<16)
    header, rest = buf.split(b"\n", 1)
    hcols = header.decode().split(",")
    assert hcols[0] in ("", id_col), hcols[:3]
    rows = {}
    chunk = rest
    done = False
    nrows = 0
    leftover = b""
    import io
    def feed(data, last=False):
        nonlocal nrows
        lines = data.split(b"\n")
        return lines
    # simpler: stream lines via TextIOWrapper over resp with leftover prepended
    import itertools
    stream = itertools.chain([rest], iter(lambda: resp.read(1<<20), b""))
    carry = b""
    for blk in stream:
        carry += blk
        *lines, carry = carry.split(b"\n")
        for ln in lines:
            if not ln: continue
            nrows += 1
            # fast id check: row starts with ACH-xxxxxx,
            try: rid = ln.split(b",",1)[0].decode()
            except Exception: continue
            if rid in KEEP:
                rows[rid] = ln
        if len(rows)==len(KEEP) and fname!="OmicsSomaticMutations.csv":
            break
    # for mutations (long format, many rows per line) keep streaming all
    if fname=="OmicsSomaticMutations.csv":
        # rows dict keyed differently; redo: collect all matching lines
        pass
    with open(W+outname,"wb") as f:
        f.write(header+b"\n")
        for k in sorted(rows): f.write(rows[k]+b"\n")
    nb = sum(len(v) for v in rows.values())
    print(f"{fname}: scanned {nrows} rows in {time.time()-t0:.0f}s, kept {len(rows)}/{len(KEEP)} ids", flush=True)
    return header.decode()
for fname, out in [("CRISPRGeneEffect.csv","CRISPRGeneEffect_subset.csv"),
                   ("OmicsExpressionProteinCodingGenesTPMLogp1.csv","Expression_subset.csv")]:
    stream_keep_rows(fname, out)
