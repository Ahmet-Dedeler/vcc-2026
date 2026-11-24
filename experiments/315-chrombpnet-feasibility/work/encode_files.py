"""Gate 2 detail -- what a ChromBPNet ENCODE4 annotation actually contains, and how big it is.

Metadata only. Nothing is downloaded.
"""
from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path

WORK = Path(__file__).resolve().parent
API = "https://www.encodeproject.org"


def get(path: str, params: dict | None = None) -> dict:
    url = f"{API}{path}"
    if params:
        url += "?" + urllib.parse.urlencode(params, doseq=True)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=90) as fh:
        return json.loads(fh.read().decode())


cov = json.loads((WORK / "gate2_coverage.json").read_text())
interesting = [r for r in cov if r["n_chrombpnet_annotations"] > 0]

out = []
for r in interesting:
    for acc in r["chrombpnet_accessions"]:
        a = get(f"/annotations/{acc}/")
        files = []
        for f in a.get("files", []):
            fo = f if isinstance(f, dict) else get(f)
            files.append({
                "accession": fo.get("accession"), "file_format": fo.get("file_format"),
                "output_type": fo.get("output_type"), "file_size": fo.get("file_size"),
                "assembly": fo.get("assembly"), "href": fo.get("href"),
            })
        rec = {
            "line": r["line"], "annotation": acc,
            "biosample": a.get("biosample_ontology", {}).get("term_name"),
            "description": a.get("description"),
            "annotation_type": a.get("annotation_type"),
            "lab": a.get("lab", {}).get("title"),
            "award_project": a.get("award", {}).get("project"),
            "status": a.get("status"),
            "n_files": len(files),
            "total_bytes": sum(f["file_size"] or 0 for f in files),
            "files": files,
        }
        out.append(rec)
        print(f"{r['line']:16s} {acc}  {rec['biosample']:20s} files={len(files):2d} "
              f"{rec['total_bytes'] / 1e6:8.1f} MB  {(rec.get('description') or '')[:70]}")

(WORK / "gate2_files.json").write_text(json.dumps(out, indent=2) + "\n")

# the one HEK293 DNase hit -- what is it?
print("\n--- HEK293* DNase/ATAC experiments ---")
for assay in ("DNase-seq", "ATAC-seq"):
    r = get("/search/", {"type": "Experiment", "assay_title": assay, "status": "released",
                         "searchTerm": "HEK293", "limit": "all", "format": "json",
                         "field": ["accession", "biosample_ontology.term_name", "assay_title",
                                   "description"]})
    for e in r["@graph"]:
        print(f"  {assay:10s} {e['accession']}  "
              f"{e.get('biosample_ontology', {}).get('term_name')}  {e.get('description', '')[:60]}")

print(f"\nwrote {WORK / 'gate2_files.json'}")
