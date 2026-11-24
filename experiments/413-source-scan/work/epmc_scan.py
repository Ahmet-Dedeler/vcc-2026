"""413: Europe PMC scan (papers + preprints) for perturbation-transcriptomics studies in a date window."""
import json, sys, urllib.request, urllib.parse
start, end = sys.argv[1], sys.argv[2]
Q = {
 "perturbseq": '(TITLE_ABS:"Perturb-seq" OR TITLE_ABS:"CROP-seq" OR TITLE_ABS:"single-cell CRISPR screen" OR TITLE_ABS:"single-cell CRISPR screens" OR TITLE_ABS:"perturbation atlas")',
 "gw_crispri": '(TITLE_ABS:"genome-wide" OR TITLE_ABS:"genome-scale") AND (TITLE_ABS:"CRISPRi" OR TITLE_ABS:"Perturb-seq" OR TITLE_ABS:"transcriptomic readout")',
 "lines": '(TITLE_ABS:"Perturb-seq" OR TITLE_ABS:"CRISPRi" OR TITLE_ABS:"CROP-seq") AND (TITLE_ABS:"HeLa" OR TITLE_ABS:"Jurkat" OR TITLE_ABS:"CAL27" OR TITLE_ABS:"CAL33" OR TITLE_ABS:"squamous" OR TITLE_ABS:"cervical")',
 "multi_line": '(TITLE_ABS:"Perturb-seq" OR TITLE_ABS:"CRISPRi") AND (TITLE_ABS:"cell lines" OR TITLE_ABS:"cellular contexts")',
}
out = {}
for name, q in Q.items():
    full = f'({q}) AND FIRST_PDATE:[{start} TO {end}]'
    cursor = "*"
    while True:
        url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(
            {"query": full, "format": "json", "pageSize": 1000, "cursorMark": cursor, "resultType": "lite"})
        d = json.load(urllib.request.urlopen(url, timeout=120))
        if "resultList" not in d: print("ERR", name, str(d)[:300], file=sys.stderr); break
        for r in d["resultList"]["result"]:
            k = r.get("doi") or r.get("id")
            e = out.setdefault(k, {"title": r.get("title"), "date": r.get("firstPublicationDate"), "src": r.get("source"),
                                   "journal": r.get("journalTitle") or r.get("bookOrReportDetails", {}).get("publisher"), "doi": r.get("doi"), "q": []})
            e["q"].append(name)
        nxt = d.get("nextCursorMark")
        if not nxt or nxt == cursor or not d["resultList"]["result"]: break
        cursor = nxt
    print(name, d.get("hitCount"), file=sys.stderr)
json.dump(out, open(f"epmc_{start}_{end}.json", "w"), indent=1)
for k, e in sorted(out.items(), key=lambda kv: kv[1]["date"] or "", reverse=True):
    print(e["date"], e["src"], ",".join(e["q"]), "|", (e["title"] or "")[:150], "|", e["doi"])
