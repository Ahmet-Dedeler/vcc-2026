"""413: GEO E-utilities scan of human CRISPR-perturbation expression series by publication date window.
Writes geo_scan.json with GSE accession, title, n samples, date, summary snippet."""
import json, time, urllib.request, urllib.parse, sys
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
QUERIES = {
 "perturbseq": '(Perturb-seq OR perturbseq OR "CROP-seq" OR "CRISPR droplet sequencing" OR "single-cell CRISPR screen" OR "pooled CRISPR screen" AND single-cell)',
 "crispri_sc": '(CRISPRi OR "CRISPR interference" OR dCas9-KRAB OR "CRISPR activation" OR CRISPRa) AND ("single cell" OR scRNA-seq OR "single-cell")',
 "crispr_screen_rna": '("CRISPR screen" OR "CRISPR knockout screen" OR "genome-wide CRISPR") AND (transcriptom* OR "RNA-seq")',
 "hela_crispr": '(HeLa) AND (CRISPR OR CRISPRi OR sgRNA OR knockout)',
 "hnscc_crispr": '(CAL33 OR CAL-33 OR CAL27 OR CAL-27 OR FaDu OR Detroit562 OR "head and neck squamous" OR HNSCC) AND (CRISPR OR CRISPRi OR sgRNA)',
 "jurkat_crispr": '(Jurkat) AND (CRISPR OR CRISPRi OR sgRNA OR Perturb-seq)',
}
def esearch(term, mindate, maxdate):
    q = f'({term}) AND "Homo sapiens"[Organism] AND gse[ETYP]'
    url = E + "esearch.fcgi?" + urllib.parse.urlencode({"db": "gds", "term": q, "retmax": 2000, "retmode": "json",
                                                       "datetype": "pdat", "mindate": mindate, "maxdate": maxdate})
    return json.load(urllib.request.urlopen(url, timeout=120))["esearchresult"]["idlist"]
def esummary(ids):
    out = []
    for i in range(0, len(ids), 150):
        url = E + "esummary.fcgi?" + urllib.parse.urlencode({"db": "gds", "id": ",".join(ids[i:i+150]), "retmode": "json"})
        r = json.load(urllib.request.urlopen(url, timeout=120))["result"]
        out += [r[k] for k in r.get("uids", [])]
        time.sleep(0.4)
    return out
mindate, maxdate = sys.argv[1], sys.argv[2]
res = {}
for name, term in QUERIES.items():
    ids = esearch(term, mindate, maxdate); time.sleep(0.4)
    for s in esummary(ids):
        acc = s.get("accession")
        d = res.setdefault(acc, {"acc": acc, "title": s.get("title"), "n_samples": s.get("n_samples"), "pdat": s.get("pdat"),
                                 "gpl": s.get("gpl"), "summary": (s.get("summary") or "")[:600], "suppfile": s.get("suppfile"), "queries": []})
        d["queries"].append(name)
    print(name, len(ids), file=sys.stderr)
json.dump(res, open(f"geo_scan_{mindate.replace('/','')}_{maxdate.replace('/','')}.json", "w"), indent=1)
for d in sorted(res.values(), key=lambda d: -int(d["n_samples"] or 0)):
    print(f'{d["acc"]:11s} {d["pdat"]} n={d["n_samples"]:>4} {",".join(q[:6] for q in d["queries"]):30s} {d["title"][:120]}')
