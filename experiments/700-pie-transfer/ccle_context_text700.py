#!/usr/bin/env python3
"""700 — precompute PIE context-text embeddings for every CCLE line, so the D/E/F day needs no embedding API.

For each DepMap ID in data/external/ccle_broad_2025/data_clinical_sample.txt (the def_identity reference): find its
Cellosaurus accession (batched API search on the DepMap cross-reference), fetch the record (cached under
work/cellosaurus/), render PIE's context text with pie.sources.text.contexts.build_context_entry under the runtime key =
the accession (PIE's Tahoe convention: "Context Name: CVCL_xxxx"), and embed with text-embedding-3-large through the Hack
Club AI proxy (the jurkat re-embed matched PIE's stored vector at cos 0.99997).
Writes data/external/pie/local_sources/context_text_ccle/{embeddings.npy, descriptions.json, meta.json} (PIE source
format) and work/ccle_cellosaurus_map.tsv (depmap_id, name, accession). Resumable: cached records and embeddings are reused.
usage: .venv/bin/python ccle_context_text700.py
"""
from __future__ import annotations

import csv
import json
import subprocess
import time
from pathlib import Path

import numpy as np
import requests

from pie.sources.text.contexts import build_context_entry

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W = HERE / "work"
CACHE = W / "cellosaurus"
ECACHE = W / "ccle_emb"
OUT = ROOT / "data/external/pie/local_sources/context_text_ccle"
BASE_SRC = ROOT / "data/external/pie/local_sources/context_text"
API = "https://api.cellosaurus.org"


def depmap_ids():
    rows = []
    with open(ROOT / "data/external/ccle_broad_2025/data_clinical_sample.txt") as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    for r in csv.DictReader(lines, delimiter="\t"):
        did = (r.get("DepMap ID") or r.get("DEPMAP_ID") or "").strip()
        if did.startswith("ACH-"):
            rows.append((did, (r.get("Cell Line Name") or r.get("CELL_LINE_NAME") or "").strip()))
    return sorted(set(rows))


def map_accessions(ids):
    out = {}
    for i in range(0, len(ids), 40):
        chunk = [d for d, _ in ids[i:i + 40]]
        q = " OR ".join(f"dr:{d}" for d in chunk)
        r = requests.get(f"{API}/search/cell-line", params={"q": q, "fields": "ac,dr", "format": "tsv", "rows": 200}, timeout=60)
        for ln in r.text.splitlines():
            if ln.startswith("#") or ln.startswith("ac\t") or not ln.strip():
                continue
            ac, dr = ln.split("\t", 1)
            for x in dr.split("|"):
                x = x.strip()
                if x.startswith("DepMap; "):
                    d = x.split("; ", 1)[1].strip()
                    if d in chunk and d not in out:
                        out[d] = ac
        time.sleep(0.2)
    return out


def record(ac):
    p = CACHE / f"{ac}.json"
    if not p.exists():
        r = requests.get(f"{API}/cell-line/{ac}", params={"format": "json"}, timeout=60)
        r.raise_for_status()
        p.write_text(r.text)
        time.sleep(0.1)
    return json.loads(p.read_text())


def embed(ac, text):
    p = ECACHE / f"{ac}.npy"
    if p.exists():
        return np.load(p)
    tf = ECACHE / f"{ac}.txt"
    tf.write_text(text)
    res = subprocess.run(["hcai", "embed", "-m", "openai/text-embedding-3-large", "--file", str(tf), "-o", "json"],
                         capture_output=True, text=True, timeout=120)
    v = np.asarray(json.loads(res.stdout)["data"][0]["embedding"], dtype=np.float32)
    assert v.shape == (3072,)
    np.save(p, v)
    return v


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    ECACHE.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    ids = depmap_ids()
    print("depmap ids", len(ids), flush=True)
    acc = map_accessions(ids)
    print("mapped", len(acc), flush=True)
    with open(W / "ccle_cellosaurus_map.tsv", "w") as f:
        f.write("depmap_id\tname\taccession\n")
        for d, n in ids:
            f.write(f"{d}\t{n}\t{acc.get(d, '')}\n")
    meta = json.loads((BASE_SRC / "meta.json").read_text())
    desc = json.loads((BASE_SRC / "descriptions.json").read_text())
    E = [np.load(BASE_SRC / "embeddings.npy")]
    keys = list(meta["keys"])
    have = set(keys)
    n_new = 0
    for d, n in ids:
        ac = acc.get(d)
        if not ac or ac in have:
            continue
        try:
            text = build_context_entry(ac, record(ac))["embedding_text"]
            v = embed(ac, text)
        except Exception as ex:  # one bad record must not stop the batch
            print("skip", d, ac, repr(ex)[:200], flush=True)
            continue
        E.append(v[None]); keys.append(ac); desc[ac] = text; have.add(ac); n_new += 1
        if n_new % 100 == 0:
            print("embedded", n_new, flush=True)
    meta["keys"] = keys
    meta["provenance"]["params"]["700_ccle"] = (f"{n_new} CCLE lines (DepMap 25Q3 via cBioPortal) rendered from Cellosaurus "
                                                "56.0 and embedded with text-embedding-3-large via the Hack Club AI proxy")
    np.save(OUT / "embeddings.npy", np.concatenate(E).astype(np.float32))
    (OUT / "descriptions.json").write_text(json.dumps(desc, indent=1))
    (OUT / "meta.json").write_text(json.dumps(meta, indent=1))
    print("wrote", len(keys), "contexts to", OUT, flush=True)


if __name__ == "__main__":
    main()
