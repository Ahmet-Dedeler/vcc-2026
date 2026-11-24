"""Gate 2 — ENCODE coverage for the lines that matter.

Two catalogues, queried separately:
  (a) ChromBPNet ENCODE4 model annotations (lab = Anshul Kundaje, Stanford)
  (b) MEASURED DNase-seq / ATAC-seq experiments (the cheaper substitute route)

Public REST API, no auth, JSON only. Nothing is downloaded here except metadata.
"""
from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

WORK = Path(__file__).resolve().parent
API = "https://www.encodeproject.org"

# Lines of interest, grouped by why we care.
LINES = {
    # Gate 1's ranked identification of the three held contexts (041 / def_identity --verify)
    "A=Jurkat": ["Jurkat", "Jurkat, Clone E6-1", "Jurkat clone E61"],
    "B=HeLa": ["HeLa-S3", "HeLa"],
    "C=CAL33": ["CAL33", "CAL 33", "CAL-33"],
    # the three donor lines the incumbent actually trains on
    "donor K562": ["K562"],
    "donor HCT116": ["HCT116"],
    "donor HEK293T": ["HEK293T", "HEK293"],
    # Broad-16's 12 WORKING lines -- the only surface where Gate 3 could be run
    "b16 C4I": ["C4I", "C-4 I", "C-4I"],
    "b16 HCT15": ["HCT15", "HCT-15"],
    "b16 HKGZCC": ["HKGZ-CC", "HKGZCC"],
    "b16 KMRC20": ["KMRC-20", "KMRC20"],
    "b16 KYSE140": ["KYSE-140", "KYSE140"],
    "b16 KYSE450": ["KYSE-450", "KYSE450"],
    "b16 ONE58": ["ONE58", "ONE-58"],
    "b16 RPMI7951": ["RPMI7951", "RPMI-7951"],
    "b16 SKGII": ["SKG-II", "SKGII"],
    "b16 SNU8": ["SNU8", "SNU-8"],
    "b16 T3M5": ["T3M-5", "T3M5"],
    "b16 U343": ["U343", "U-343"],
    # Broad-16's 4 RESERVED lines -- metadata only, never read, never used
    "b16res JVE127": ["JVE-127", "JVE127"],
    "b16res MG63": ["MG63", "MG-63"],
    "b16res SLR23": ["SLR23", "SLR-23"],
    "b16res UMRC3": ["UMRC3", "UMRC-3"],
}


def get(path: str, params: dict) -> dict:
    q = urllib.parse.urlencode(params, doseq=True)
    req = urllib.request.Request(f"{API}{path}?{q}", headers={"Accept": "application/json"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=90) as fh:
                return json.loads(fh.read().decode())
        except Exception as exc:  # noqa: BLE001
            if attempt == 2:
                raise
            print(f"  retry {attempt + 1} after {exc}")
            time.sleep(3)
    raise RuntimeError


def norm(s: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", s.upper())


# ---------------------------------------------------------------- (a) ChromBPNet annotations
ann = json.loads((WORK / "chrombpnet_annotations.json").read_text())["@graph"]
by_bio: dict[str, list[dict]] = {}
for a in ann:
    tn = a.get("biosample_ontology", {}).get("term_name", "<none>")
    by_bio.setdefault(tn, []).append(a)
all_bios = sorted(by_bio)

# ------------------------------------------- (b) measured DNase-seq / ATAC-seq per biosample
# One query for the whole catalogue of released human DNase/ATAC experiments; then match names.
acc_exp = {}
for assay in ("DNase-seq", "ATAC-seq"):
    r = get("/search/", {
        "type": "Experiment", "assay_title": assay, "status": "released",
        "replicates.library.biosample.donor.organism.scientific_name": "Homo sapiens",
        "limit": "all", "format": "json",
        "field": ["accession", "biosample_ontology.term_name",
                  "biosample_ontology.classification", "assay_title"],
    })
    for e in r["@graph"]:
        tn = e.get("biosample_ontology", {}).get("term_name", "<none>")
        acc_exp.setdefault((assay, tn), []).append(e["accession"])
    print(f"[encode] {assay}: {r['total']} released human experiments, "
          f"{len({k[1] for k in acc_exp if k[0] == assay})} biosamples")

rows = []
for label, aliases in LINES.items():
    qs = {norm(a) for a in aliases}
    cb_hits = [b for b in all_bios if norm(b) in qs]
    if not cb_hits:  # substring fallback, so "Jurkat, Clone E6-1" is found from "Jurkat"
        cb_hits = [b for b in all_bios if any(q and q in norm(b) for q in qs)]
    dn = sorted({a for (asy, tn), accs in acc_exp.items() if asy == "DNase-seq"
                 and (norm(tn) in qs or any(q and q in norm(tn) for q in qs)) for a in accs})
    at = sorted({a for (asy, tn), accs in acc_exp.items() if asy == "ATAC-seq"
                 and (norm(tn) in qs or any(q and q in norm(tn) for q in qs)) for a in accs})
    n_models = sum(len(by_bio[b]) for b in cb_hits)
    rows.append({
        "line": label, "chrombpnet_biosamples": cb_hits, "n_chrombpnet_annotations": n_models,
        "chrombpnet_accessions": [a["accession"] for b in cb_hits for a in by_bio[b]][:8],
        "dnase_experiments": dn, "atac_experiments": at,
    })
    print(f"{label:16s} chrombpnet={n_models:3d} {cb_hits!s:28s} DNase={len(dn):3d} ATAC={len(at):3d}")

(WORK / "gate2_coverage.json").write_text(json.dumps(rows, indent=2) + "\n")
print(f"\nwrote {WORK / 'gate2_coverage.json'}")
