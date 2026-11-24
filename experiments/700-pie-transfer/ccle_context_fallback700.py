#!/usr/bin/env python3
"""700: second pass over the CCLE lines PIE's context renderer rejected (Cellosaurus ages such as '>50Y', '4-5Y',
'25-26Y' that `_humanize_age` does not parse), so no D/E/F line is left without a context text.

Rendering is PIE's own `build_context_entry`, with one change for these records only: an age `_humanize_age` cannot parse
is carried through verbatim (e.g. "Age: 25-26Y") instead of raising. Embedding is text-embedding-3-large via Hack Club
AI, as in ccle_context_text700.py. Appends to data/external/pie/local_sources/context_text_ccle in place (PIE source format).
usage: .venv/bin/python ccle_context_fallback700.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

import pie.sources.text.contexts as C
from ccle_context_text700 import CACHE, ECACHE, OUT, embed, record  # noqa: E402  (same dir, same caches)

HERE = Path(__file__).resolve().parent
_orig = C._humanize_age


def _lenient(raw):
    try:
        return _orig(raw)
    except ValueError:
        return C._clean(raw)


def main():
    skipped = []
    for ln in (HERE / "work/ccle_context_text700.log").read_text().splitlines():
        m = re.match(r"skip (ACH-\d+) (CVCL_\w+) ", ln)
        if m:
            skipped.append(m.group(2))
    meta = json.loads((OUT / "meta.json").read_text())
    desc = json.loads((OUT / "descriptions.json").read_text())
    E = [np.load(OUT / "embeddings.npy")]
    keys = list(meta["keys"])
    have = set(keys)
    C._humanize_age = _lenient
    added = []
    for ac in sorted(set(skipped) - have):
        text = C.build_context_entry(ac, record(ac))["embedding_text"]
        v = embed(ac, text)
        E.append(v[None]); keys.append(ac); desc[ac] = text; added.append(ac)
    meta["keys"] = keys
    meta["provenance"]["params"]["700_ccle_fallback"] = (f"{len(added)} lines whose Cellosaurus age PIE's parser rejects, "
                                                         "rendered with the age carried verbatim")
    np.save(OUT / "embeddings.npy", np.concatenate(E).astype(np.float32))
    (OUT / "descriptions.json").write_text(json.dumps(desc, indent=1))
    (OUT / "meta.json").write_text(json.dumps(meta, indent=1))
    print(f"skipped in pass 1: {len(set(skipped))}; added now: {len(added)}; total contexts {len(keys)}")


if __name__ == "__main__":
    main()
