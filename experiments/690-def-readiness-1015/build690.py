"""690 — two mock D/E/F bundles for the October dress rehearsal on the CURRENT tree (2026-10-15).

Both bundles carry the SAME three foreign control files (595's bundle X recipe, re-used through its own readers):
  D = HCT116  (X-Atlas/Orion NTC, FiCS; in CCLE and in the same-line registry: strong bank orion_hct116)
  E = H1      (Arc VCC 2025 10x Flex controls from data/harness/h1-test.h5ad; absent from CCLE -> abstain path)
  F = HEK293T (Orion NTC; absent from CCLE; 590's false-HAP1 case, blocked by the 0.65 floor)
They differ only in the panel:
  new/  = ef-bundle-r's 300 targets (238's essential-free Arc-like draw): the "new panel" day (CD4 gate, whitener refit,
          LODO key refit on a new panel)
  same/ = A/B/C's 300 targets (data/challenge/pert_counts.csv): the "same panel" day that Arc's "same sgRNAs in all six
          lines" makes plausible. On it the day's whitener fit reproduces A/B/C's operator e2f65ecf byte for byte, and
          def_pipeline aborts by design (453 attempt 1). This bundle rehearses the branch that handles it.
The context files are written once (under new/) and APFS-cloned into same/ (`cp -c`), so the two rehearsals read
byte-identical controls.

Usage: uv run python experiments/690-def-readiness-1015/build690.py <out dir>
Writes <out>/{new,same}/{context_D,E,F.h5ad, gene_names.csv, pert_counts.csv, manifest.json} and <out>/build690.json.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("b595", ROOT / "experiments/595-foreign-control-rehearsal/build595.py")
b595 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b595)  # 595's readers; its main() is guarded and does not run

PANELS = {
    "new": ROOT / "experiments/238-essential-free-panel/work/ef-bundle-r/pert_counts.csv",
    "same": ROOT / "data/challenge/pert_counts.csv",
}
CONTROLS = {"D": ("HCT116", lambda: b595.orion("HCT116")),
            "E": ("H1", b595.h1),
            "F": ("HEK293T", lambda: b595.orion("HEK293T"))}


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 24), b""):
            h.update(blk)
    return h.hexdigest()


def main() -> None:
    out = Path(sys.argv[1]).resolve()
    template = json.loads((ROOT / "experiments/180-def-readiness/work/def-bundle/manifest.json").read_text())
    record = {"seed": b595.SEED, "official_axis": "data/challenge/gene_names.csv",
              "panels": {k: {"path": str(v.relative_to(ROOT)), "sha256": sha(v)} for k, v in PANELS.items()},
              "controls": {}, "bundles": {}}
    for name, panel in PANELS.items():
        d = out / name
        d.mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / "data/challenge/gene_names.csv", d / "gene_names.csv")
        shutil.copy(panel, d / "pert_counts.csv")
        man = dict(template)
        man["partition"] = "synthetic-def"
        man["panel_id"] = f"vcc2026-mock-690-{name}"
        (d / "manifest.json").write_text(json.dumps(man, indent=2) + "\n")
    first = out / "new"
    for ctx, (line, fn) in CONTROLS.items():
        t0 = time.time()
        b595.log(f"context {ctx} = {line}")
        m, src = fn()
        rec = b595.write_context(m, ctx, first / f"context_{ctx}.h5ad")
        rec.update(line=line, **src, seconds=round(time.time() - t0, 1))
        record["controls"][ctx] = rec
        b595.log(f"  wrote {first / f'context_{ctx}.h5ad'}: {rec['cells']} cells, nnz {rec['nnz']}, "
                 f"q75/q25 {rec['q75_q25']:.3f}, median lib {rec['lib_median']:.0f}")
        del m
        # byte-identical controls in the same-panel bundle (APFS clone: no extra disk)
        subprocess.run(["cp", "-c", str(first / f"context_{ctx}.h5ad"), str(out / "same" / f"context_{ctx}.h5ad")],
                       check=True)
    for name in PANELS:
        record["bundles"][name] = {p.name: sha(p) for p in sorted((out / name).iterdir())}
    (out / "build690.json").write_text(json.dumps(record, indent=1) + "\n")
    b595.log(f"wrote {out / 'build690.json'}")


if __name__ == "__main__":
    main()
