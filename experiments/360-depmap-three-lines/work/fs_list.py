import json, urllib.request
d = json.load(urllib.request.urlopen("https://api.figshare.com/v2/articles/27993248", timeout=60))
print("title:", d.get("title"))
print("license:", json.dumps(d.get("license"))[:500])
print("defined_type:", d.get("defined_type"), "| doi:", d.get("doi"), "| published:", d.get("published_date"))
print("n files:", len(d.get("files", [])))
want = ("Model.csv", "CRISPRGeneEffect.csv", "OmicsExpressionProteinCodingGenesTPMLogp1",
        "OmicsSomaticMutations", "CCLE_expression", "sample_info")
for f in d.get("files", []):
    n = f.get("name", "")
    if n.startswith(want):
        print(f'{f.get("size",0)/1e6:9.1f} MB  {n}')
