"""252 — build the SWAPPED operator for the identity probes: a byte copy of 062's
`held_orion_hct116.npz` with ONE eigenvector entry moved by one float32 ULP (re-saved with
`np.savez`, uncompressed, as the fitter saved the original), plus a byte copy of its metadata
json under the new name.

Every metadata check `E064.load_frozen_operator` runs (axis, donors, gamma, n_cols,
wcols_sha256, eigenvector shape) passes on this file, because none of them reads the
eigenbasis; and the pre-252 sha check passes too, because the caller hands it the file's own
hash.  So this is the file that lets finding #4 (the tautology) and #5 (`plolo.Operator`
accepts anything of the right shape) be shown FIRING, not argued.

    uv run python mk_tampered252.py        # -> work/probe/held_orion_hct116_TAMPERED.{npz,json}
"""
from __future__ import annotations

import hashlib
import json
import shutil
import struct
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SRC = ROOT / "experiments/062-cd4-post-whiten/cloud-output/whiteners/held_orion_hct116.npz"
OUT = HERE / "work/probe/held_orion_hct116_TAMPERED.npz"
FLIP_ELEMENT = 12_345_678          # which float32 of `evecs` (row-major) gets its low byte flipped


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(1 << 22):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    import numpy as np
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with np.load(SRC) as z:
        evals = np.asarray(z["evals"])
        evecs = np.asarray(z["evecs"])              # float32 (n, n), as stored
    assert evecs.dtype == np.float32 and evecs.flags.c_contiguous
    flat = evecs.reshape(-1)
    before = float(flat[FLIP_ELEMENT])
    flat[FLIP_ELEMENT] = np.nextafter(flat[FLIP_ELEMENT], np.float32(np.inf))   # one ULP
    after = float(flat[FLIP_ELEMENT])
    # np.savez, uncompressed, exactly as the fitter wrote the original; the zip CRC is recomputed
    np.savez(OUT, evals=evals, evecs=evecs)
    shutil.copyfile(SRC.with_suffix(".json"), OUT.with_suffix(".json"))
    rec = dict(source=str(SRC.relative_to(ROOT)), source_sha256=sha256_file(SRC),
               tampered=str(OUT.relative_to(ROOT)), tampered_sha256=sha256_file(OUT),
               flipped_element=FLIP_ELEMENT, value_before=before, value_after=after,
               change="one ULP of one float32 eigenvector entry",
               bytes=OUT.stat().st_size, source_bytes=SRC.stat().st_size)
    (HERE / "work/probe/tampered.json").write_text(json.dumps(rec, indent=2))
    print(json.dumps(rec, indent=2))


if __name__ == "__main__":
    main()
