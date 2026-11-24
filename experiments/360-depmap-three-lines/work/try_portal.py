import urllib.request
urls = [
    "https://depmap.org/portal/api/download/files",
    "https://depmap.org/portal/download/api/download?file=Model.csv&release=DepMap+Public+26Q1",
]
for u in urls:
    try:
        req = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126.0 Safari/537.36"})
        r = urllib.request.urlopen(req, timeout=30)
        body = r.read(300)
        print(u, "->", r.status, body[:150])
    except Exception as e:
        print(u, "-> ERROR", repr(e)[:200])
