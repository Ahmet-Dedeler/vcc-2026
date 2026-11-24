#!/usr/bin/env python3
"""453's slot gates, checked by the coordinator after chain453 ENDs; writes work/gates453.log for arm0928.sh.
PASS requires: def_day converged (or shipped its nearest safe build) with every final rho <= 1; the zero-gene audit at 0;
the relabelled meta's config differs from 446b's ONLY in the preregistered components; chain453's novelty rc 0 and a
package sha that matches the file.   usage: uv run python experiments/453-abc-dress-rehearsal/gate453.py [--write]"""
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
X = ROOT / "experiments/453-abc-dress-rehearsal"; W = X / "work"; RUN = W / "run"
PKG = ROOT / "submissions/453dr-october-path-abc-q1.75.vcc"
H5 = ROOT / "submissions/453dr-october-path-abc-q1.75.h5ad"
REF = ROOT / "submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
ALLOWED = {"bulk_ctx_scale", "bulk_load", "bulk_load_sha256", "contexts", "expect_contexts", "cd4_deflate_k",
           "cd4_deflate_basis", "cd4_sizing", "cd4_target_g", "pert_scale", "pert_scale_sha256", "sameline"}
out, ok = [], True
dd = (RUN / "def_day.log").read_text()
conv = re.findall(r"CONVERGED at (\S+) after", dd) or re.findall(r"NOT converged .*nearest safe build: (\S+)", dd)
fin = re.findall(r"final: (\S+) factors (\S+) ->", dd)
if not fin:
    print("def_day has no final line"); sys.exit(1)
final_name, factors = fin[-1]
nt = json.loads((RUN / f"{final_name}.ntclevel.json").read_text())
rho = {c: round(v["rho"], 4) for c, v in nt["per_context"].items()}
cliff_ok = all(v <= 1.0 for v in rho.values())
ok &= cliff_ok and bool(conv)
out.append(f"[453dr] {'CLIFF OK' if cliff_ok else 'CLIFF FAIL'} {json.dumps(rho)} final {final_name} factors {factors} "
           f"({'converged' if 'CONVERGED at' in dd else 'nearest safe build'})")
zl = [l for l in dd.splitlines() if "ZERO block total" in l]   # audit428's line: "...ZERO block total (target excluded): 0 in 0 blocks; ..."
zero_ok = bool(zl) and re.search(r"ZERO block total \(target excluded\): 0 in 0 blocks", zl[-1]) is not None
out.append(f"zero-gene audit: {'0 (OK)' if zero_ok else 'NOT 0 / unparsed'} :: {zl[-1][:160] if zl else '<no line>'}")
ok &= zero_ok
meta = json.loads(Path(str(H5) + ".meta.json").read_text()); ref = json.loads(REF.read_text())
c, r = meta["config"], ref["config"]
diff = sorted(k for k in set(c) | set(r) if c.get(k) != r.get(k))
extra = [k for k in diff if k not in ALLOWED]
out.append(f"[453dr] {'PASS' if not extra else 'FAIL'}: config diff vs 446b {diff}" + (f"  UNEXPECTED {extra}" if extra else ""))
ok &= not extra
out.append(f"relabel: {meta.get('relabel453', {}).get('categories_after')} cells {meta.get('relabel453', {}).get('cells')}")
cl = (W / "chain453.log").read_text()
nov = re.findall(r"novelty rc=(\d+)", cl); pk = re.findall(r"package (\S+) sha256 ([0-9a-f]{64})", cl)
h = hashlib.sha256()
with open(PKG, "rb") as f:
    for chunk in iter(lambda: f.read(1 << 24), b""):
        h.update(chunk)
sha_ok = bool(pk) and pk[-1][1] == h.hexdigest()
ok &= bool(nov) and nov[-1] == "0" and sha_ok
out.append(f"novelty rc={nov[-1] if nov else 'missing'}")
out.append(f"package sha256: {h.hexdigest()}" + ("" if sha_ok else "  (MISMATCH with chain453.log)"))
print("\n".join(out)); print("ALL GATES", "PASS" if ok else "FAIL")
if "--write" in sys.argv:
    if not ok:
        print("not writing gates453.log (a gate failed)"); sys.exit(1)
    (W / "gates453.log").write_text("\n".join(out) + "\n")
    print("wrote", W / "gates453.log")
