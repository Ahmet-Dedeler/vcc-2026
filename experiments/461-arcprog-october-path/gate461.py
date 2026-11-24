#!/usr/bin/env python3
"""461dr's slot gates (the October path + Arc's denoised programme for every context, vs 453dr), checked after chain461
ENDs; writes work/gates461.log. PASS requires: def_day converged with every final rho <= 1; the programme source was
h1_2025_full for every context; the zero-gene audit at 0; the relabelled meta's config differs from 453dr's ONLY in
`sameline` (and, flagged not failed, `bulk_ctx_scale`, which should be identical because the bulk is untouched);
chain461's novelty rc 0 and a package sha that matches the file.
usage: uv run python experiments/461-arcprog-october-path/gate461.py [--write]"""
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
X = ROOT / "experiments/461-arcprog-october-path"; W = X / "work"; RUN = W / "run"
PKG = ROOT / "submissions/461dr-october-path-arcprog-abc-q1.75.vcc"
H5 = ROOT / "submissions/461dr-october-path-arcprog-abc-q1.75.h5ad"
REF = ROOT / "submissions/453dr-october-path-abc-q1.75.h5ad.meta.json"
ALLOWED = {"sameline", "bulk_ctx_scale"}
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
out.append(f"[461dr] {'CLIFF OK' if cliff_ok else 'CLIFF FAIL'} {json.dumps(rho)} final {final_name} factors {factors} "
           f"({'converged' if 'CONVERGED at' in dd else 'nearest safe build'})")
import numpy as np
sl = np.load(RUN / "sameline.npz", allow_pickle=False)
prov = json.loads(str(sl["provenance"]))
src = {c: v.get("programme_source", {}).get("bank") for c, v in prov.get("contexts", {}).items()}
prog_ok = sorted(src) == ["D", "E", "F"] and set(src.values()) == {"h1_2025_full"} and "LODO key ON" not in dd
mode = [l for l in dd.splitlines() if "same-line mode ON" in l]
prog_ok &= all(bool(sl[f"programme_ok_{c}"]) for c in "DEF") and bool(mode) and mode[-1].rstrip().endswith(":auto:pt")
out.append(f"programme source: {json.dumps(src)} {'OK' if prog_ok else 'NOT h1_2025_full on D/E/F (or LODO on / noprog)'}")
ok &= prog_ok
zl = [l for l in dd.splitlines() if "ZERO block total" in l]
zero_ok = bool(zl) and re.search(r"ZERO block total \(target excluded\): 0 in 0 blocks", zl[-1]) is not None
out.append(f"zero-gene audit: {'0 (OK)' if zero_ok else 'NOT 0 / unparsed'} :: {zl[-1][:160] if zl else '<no line>'}")
ok &= zero_ok
meta = json.loads(Path(str(H5) + ".meta.json").read_text()); ref = json.loads(REF.read_text())
c, r = dict(meta["config"]), dict(ref["config"])
diff = sorted(k for k in set(c) | set(r) if c.get(k) != r.get(k))
# the stack dump is rewritten per run and records its env (run path; 453's also recorded q 1.25, before 455 moved the
# driver default): its file sha may differ, but its ARRAYS (bulk, carriers, perts, genes) must be identical
if "bulk_load_sha256" in diff:
    import numpy as np
    da = np.load(RUN / "dump_h1x2_cd4x1.npz", allow_pickle=True)
    db = np.load(ROOT / "experiments/453-abc-dress-rehearsal/work/run/dump_h1x2_cd4x1.npz", allow_pickle=True)
    arrs = ("bulk", "carriers", "perts", "genes")
    same = all(np.array_equal(da[k], db[k]) for k in arrs)
    if same:
        diff.remove("bulk_load_sha256"); out.append(f"bulk_load: identical dump arrays {arrs} (the file differs only in its env record)")
    else:
        out.append("bulk_load: DUMP ARRAYS DIFFER from 453dr's")
extra = [k for k in diff if k not in ALLOWED]
out.append(f"[461dr] {'PASS' if not extra else 'FAIL'}: config diff vs 453dr {diff}" + (f"  UNEXPECTED {extra}" if extra else ""))
if "bulk_ctx_scale" in diff:
    out.append(f"NOTE bulk_ctx_scale differs from 453dr's ({c.get('bulk_ctx_scale')} vs {r.get('bulk_ctx_scale')}): the bulk "
               "should be untouched -- explain before submitting")
ok &= not extra and "sameline" in diff
out.append(f"relabel: {meta.get('relabel453', {}).get('categories_after')} cells {meta.get('relabel453', {}).get('cells')}")
cl = (W / "chain461.log").read_text()
nov = re.findall(r"novelty rc=(\d+)", cl); pk = re.findall(r"package (\S+) sha256 ([0-9a-f]{64})", cl)
h = hashlib.sha256()
with open(PKG, "rb") as f:
    for chunk in iter(lambda: f.read(1 << 24), b""):
        h.update(chunk)
sha_ok = bool(pk) and pk[-1][1] == h.hexdigest()
ok &= bool(nov) and nov[-1] == "0" and sha_ok
out.append(f"novelty rc={nov[-1] if nov else 'missing'}")
out.append(f"package sha256: {h.hexdigest()}" + ("" if sha_ok else "  (MISMATCH with chain461.log)"))
print("\n".join(out)); print("ALL GATES", "PASS" if ok else "FAIL")
if "--write" in sys.argv:
    if not ok:
        print("not writing gates461.log (a gate failed)"); sys.exit(1)
    (W / "gates461.log").write_text("\n".join(out) + "\n")
    print("wrote", W / "gates461.log")
