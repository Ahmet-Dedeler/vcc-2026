"""Pick the surface for the substitute test: ENCODE DNase-seq cell lines that also have a CCLE
basal expression profile.

Frozen selection rule (written before any accessibility value was computed):
  * ENCODE Experiment, assay_title = DNase-seq, status = released, Homo sapiens,
    biosample_ontology.classification = "cell line";
  * the biosample term_name, normalised to [A-Z0-9], must equal the leading token of exactly one
    CCLE sample id in `ccle_cache.npz` (CCLE ids are `<LINE>_<TISSUE>`);
  * the experiment must expose a GRCh38 `bed narrowPeak` peak file in a released analysis.
Ambiguous matches (a normalised name hitting 0 or >1 CCLE line) are dropped, not guessed.
"""
from __future__ import annotations

import json
import re
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
WORK = Path(__file__).resolve().parent
API = "https://www.encodeproject.org"
PEAK_TYPES = ("pseudoreplicated peaks", "replicated peaks", "IDR thresholded peaks", "peaks")


def get(path: str, params: dict | None = None) -> dict:
    url = f"{API}{path}" + ("?" + urllib.parse.urlencode(params, doseq=True) if params else "")
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as fh:
        return json.loads(fh.read().decode())


def norm(s: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", s.upper())


# ---- CCLE line names ---------------------------------------------------------------------
z = np.load(ROOT / "data/external/ccle_broad_2025/ccle_cache.npz", allow_pickle=False)
ccle_samples = [str(s) for s in z["samples"]]
ccle_by_norm: dict[str, list[str]] = defaultdict(list)
for s in ccle_samples:
    ccle_by_norm[norm(s.split("_")[0])].append(s)
print(f"[ccle] {len(ccle_samples)} samples, {len(ccle_by_norm)} distinct normalised line names")

# ---- ENCODE DNase-seq on cell lines ------------------------------------------------------
r = get("/search/", {
    "type": "Experiment", "assay_title": "DNase-seq", "status": "released",
    "replicates.library.biosample.donor.organism.scientific_name": "Homo sapiens",
    "biosample_ontology.classification": "cell line",
    "limit": "all", "format": "json",
    "field": ["accession", "biosample_ontology.term_name", "date_released",
              "files.accession", "files.file_format", "files.file_format_type",
              "files.output_type", "files.assembly", "files.status", "files.file_size",
              "files.href", "files.preferred_default"],
})
print(f"[encode] {r['total']} released human DNase-seq experiments on cell lines")

# Two hand-declared aliases, and only two: the ENCODE term names for Gate 1's contexts A and B
# do not normalise onto CCLE's ids ("Jurkat, Clone E6-1" -> JURKATCLONEE61; "HeLa-S3" -> HELAS3).
# Both are the same line in both catalogues; declared here so the held contexts are not dropped
# by a string rule. No other alias is added.
ALIASES = {"JURKATCLONEE61": "JURKAT", "HELAS3": "HELA"}

cands: dict[str, list[dict]] = defaultdict(list)
for e in r["@graph"]:
    tn = e.get("biosample_ontology", {}).get("term_name", "")
    key = ALIASES.get(norm(tn), norm(tn))
    hits = ccle_by_norm.get(key, [])
    if len(hits) != 1:
        continue
    peaks = [f for f in e.get("files", [])
             if f.get("file_format") == "bed" and f.get("file_format_type") == "narrowPeak"
             and f.get("assembly") == "GRCh38" and f.get("status") == "released"
             and f.get("output_type") in PEAK_TYPES]
    if not peaks:
        continue
    # deterministic pick: preferred_default first, then the narrowest output_type, then accession
    peaks.sort(key=lambda f: (not f.get("preferred_default", False),
                              PEAK_TYPES.index(f["output_type"]), f["accession"]))
    cands[tn].append({"experiment": e["accession"], "ccle": hits[0],
                      "date_released": e.get("date_released"), "peak": peaks[0],
                      "n_peak_files": len(peaks)})

print(f"[match] {len(cands)} ENCODE cell lines matched to exactly one CCLE line "
      f"and carry GRCh38 narrowPeak")
rows = []
for tn, exps in sorted(cands.items(), key=lambda kv: -len(kv[1])):
    exps.sort(key=lambda x: (x["date_released"] or "", x["experiment"]), reverse=True)
    rows.append({"encode_term_name": tn, "ccle_sample": exps[0]["ccle"],
                 "n_experiments": len(exps), "experiments": exps})
    print(f"  {tn:26s} -> {exps[0]['ccle']:44s} n_exp={len(exps):2d} "
          f"peak={exps[0]['peak']['accession']} "
          f"{(exps[0]['peak']['file_size'] or 0)/1e6:6.2f} MB {exps[0]['peak']['output_type']}")

(WORK / "line_selection.json").write_text(json.dumps(rows, indent=2) + "\n")
tot = sum((x["experiments"][0]["peak"]["file_size"] or 0) for x in rows)
print(f"\n{len(rows)} lines; one peak file each = {tot/1e6:.0f} MB; "
      f"two each (where available) = "
      f"{sum((e['peak']['file_size'] or 0) for x in rows for e in x['experiments'][:2])/1e6:.0f} MB")
print(f"wrote {WORK / 'line_selection.json'}")
