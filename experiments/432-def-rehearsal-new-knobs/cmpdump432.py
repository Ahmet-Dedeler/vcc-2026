"""432: compare the three D/E/F bulk dumps (VCC_BULK_DUMP) against each other.

  dump_stack  H1 x2 (VCC_H1_UNION_G=0.56) + CD4 x1.5 written the D/E/F way (VCC_CD4_TARGET_G = 1.5 x 0.3443655...)
  dump_lit    H1 x2 + 417's LITERAL A/B/C recipe for CD4 x1.5 (VCC_CD4_TARGET_F=0.51)
  dump_h1     H1 x2 alone

Claim under test: on a non-A/B/C panel the driver's DEF_ENV sets VCC_CD4_SIZING=carrying, which sizes CD4 by
VCC_CD4_TARGET_G and only RECORDS VCC_CD4_TARGET_F, so 417's recipe silently delivers H1 x2 alone:
dump_lit == dump_h1 bit for bit, and dump_stack differs from both.
usage: uv run python cmpdump432.py
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

W = Path(__file__).resolve().parent / "work"
FILES = {"stack": W / "dump_stack/dump_h1x2_cd4x15.npz", "lit": W / "dump_lit/dump_h1x2_cd4f051.npz",
         "h1": W / "dump_h1/dump_h1x2.npz"}


def load(p: Path):
    with np.load(p, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def main() -> None:
    d = {k: load(p) for k, p in FILES.items() if p.exists()}
    out = {}
    for k, z in d.items():
        env = json.loads(str(z["env"]))
        out[k] = dict(sha256=hashlib.sha256(FILES[k].read_bytes()).hexdigest(), bytes=FILES[k].stat().st_size,
                      shape=list(z["bulk"].shape), perts_first=str(z["perts"][0]),
                      bulk_frob=float(np.linalg.norm(z["bulk"])), carriers_frob=float(np.linalg.norm(z["carriers"])),
                      env={e: env.get(e) for e in ("VCC_H1_UNION_G", "VCC_CD4_TARGET_F", "VCC_CD4_TARGET_G",
                                                   "VCC_CD4_SIZING", "VCC_PERT_SCALE", "VCC_BULK_DUMP")})
        print(k, json.dumps(out[k]))
    pairs = [("lit", "h1"), ("stack", "lit"), ("stack", "h1")]
    for a, b in pairs:
        if a in d and b in d:
            for arr in ("bulk", "carriers"):
                x, y = d[a][arr], d[b][arr]
                same = bool(np.array_equal(x, y))
                rel = float(np.linalg.norm(x - y) / max(np.linalg.norm(y), 1e-300))
                out[f"{a}_vs_{b}_{arr}"] = dict(bit_identical=same, rel_frob_diff=rel)
                print(f"{a} vs {b} [{arr}]: bit-identical {same}  |x-y|/|y| {rel:.6e}")
    (W / "cmpdump432.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
