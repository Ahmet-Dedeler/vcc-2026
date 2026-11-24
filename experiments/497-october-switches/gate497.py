#!/usr/bin/env python3
"""497dr's slot gates (the October path + the LODO key at the percentile cap + S_P2 per-target sharpening, vs 460dr),
checked after chain497 ENDs; writes work/gates497.log. PASS requires: def_day converged (or shipped its nearest safe build)
with every final rho <= 1; the LODO key was ON for the dump under the percentile cap (pin id lodo_lo0.5_pct_def, hi read
from its sidecar); S_P2 in the final meta; the zero-gene audit at 0; the relabelled meta's config differs from 460dr's
ONLY in the dump, the re-level factors and the sharpening form (`sameline` compared by content: its path is per run);
chain497's novelty rc 0 and a package sha that matches the file.
usage: uv run python experiments/497-october-switches/gate497.py [--write]"""
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
X = ROOT / "experiments/497-october-switches"; W = X / "work"; RUN = W / "run"
PKG = ROOT / "submissions/497dr-october-path-lodopct-spt-abc.vcc"
H5 = ROOT / "submissions/497dr-october-path-lodopct-spt-abc.h5ad"
REF = ROOT / "submissions/460dr-october-path-lodo-abc-q1.75.h5ad.meta.json"
LODO_ID = "lodo_lo0.5_pct_def"
ALLOWED = {"bulk_ctx_scale", "bulk_load", "bulk_load_sha256", "bulk_sharpen_q", "bulk_sharpen_pt"}
out, ok = [], True
dd = (RUN / "def_day.log").read_text()
pb_path = RUN / "postbuild.json"
pb = json.loads(pb_path.read_text()) if pb_path.is_file() else {}
failed = [c.get("check") for c in pb.get("checks", []) if not c.get("pass")]
if not pb or failed:
    print("postbuild.json missing or failing:", failed or "no file")
    sys.exit(1)
conv = re.findall(r"CONVERGED at (\S+) after", dd) or re.findall(r"NOT converged .*nearest safe build: (\S+)", dd)
fin = re.findall(r"final: (\S+) factors (\S+) ->", dd)
if not fin:
    print("def_day has no final line"); sys.exit(1)
final_name, factors = fin[-1]
nt = json.loads((RUN / f"{final_name}.ntclevel.json").read_text())
rho = {c: round(v["rho"], 4) for c, v in nt["per_context"].items()}
cliff_ok = all(v <= 1.0 for v in rho.values())
ok &= cliff_ok and bool(conv)
out.append(f"[497dr] {'CLIFF OK' if cliff_ok else 'CLIFF FAIL'} {json.dumps(rho)} final {final_name} factors {factors} "
           f"({'converged' if 'CONVERGED at' in dd else 'nearest safe build'})")
key = RUN / "pert_scale_lodo_def_lo0.5_pct.npz"
ksha = hashlib.sha256(key.read_bytes()).hexdigest() if key.exists() else ""
lodo_ok = bool(ksha) and f"LODO key ON for the stack dump: {LODO_ID} {ksha[:16]}" in dd
side = json.loads(key.with_suffix(".json").read_text()) if key.exists() else {}
lodo_ok &= side.get("hi_pct") == 96.64697522769474
out.append(f"LODO key: {'ON, ' + LODO_ID + ' ' + ksha[:16] + ' cap ' + str(side.get('hi')) + ' (p' + str(side.get('hi_pct')) + '), at cap ' + str(side.get('frac_at_cap')) if lodo_ok else 'NOT ON / wrong artifact'}")
pt = json.loads(Path(str(H5) + ".meta.json").read_text())["config"].get("bulk_sharpen_pt") if Path(str(H5) + ".meta.json").exists() else None
spt_ok = isinstance(pt, dict) and pt.get("mode") == "stercile" and pt.get("q") == [2.5, 2.0, 1.75]
out.append(f"S_P2: {'ON ' + json.dumps({k: pt[k] for k in ('mode', 'q', 'counts')}) if spt_ok else 'NOT ON'}")
ok &= spt_ok
ok &= lodo_ok
audit_path = RUN / "audit428.txt"
audit_lines = audit_path.read_text().splitlines() if audit_path.is_file() else []
zero_ok = bool(audit_lines) and re.search(r"ZERO block total \(target excluded\): 0 in 0 blocks", audit_lines[0]) is not None
out.append(f"zero-gene audit: {'0 (OK)' if zero_ok else 'NOT 0 / unparsed'} :: {audit_lines[0][:160] if audit_lines else '<no audit428.txt>'}")
ok &= zero_ok
meta = json.loads(Path(str(H5) + ".meta.json").read_text()); ref = json.loads(REF.read_text())
c, r = dict(meta["config"]), dict(ref["config"])
diff = sorted(k for k in set(c) | set(r) if c.get(k) != r.get(k))
# the same-line artifact is rewritten per work dir (its provenance names the identity file's path): only its path and
# file sha may differ (457's rule on the recorded summary), and its ARRAYS must be identical (checked on the files)
if "sameline" in diff and isinstance(c.get("sameline"), dict) and isinstance(r.get("sameline"), dict):
    strip = lambda d: {k: v for k, v in d.items() if k not in ("path", "sha256")}
    same = strip(c["sameline"]) == strip(r["sameline"])
    import numpy as np
    za, zb = np.load(c["sameline"]["path"], allow_pickle=False), np.load(r["sameline"]["path"], allow_pickle=False)
    keys = (set(za.files) | set(zb.files)) - {"provenance"}
    arrays_same = all(k in za.files and k in zb.files and np.array_equal(za[k], zb[k], equal_nan=za[k].dtype.kind == "f")
                      for k in keys)
    if same and arrays_same:
        diff.remove("sameline"); out.append(f"sameline: identical content ({len(keys)} arrays; per-run path/sha only)")
    else:
        out.append(f"sameline: CONTENT DIFFERS (summary equal {same}, arrays equal {arrays_same})")
extra = [k for k in diff if k not in ALLOWED]
out.append(f"[497dr] {'PASS' if not extra else 'FAIL'}: config diff vs 460dr {diff}" + (f"  UNEXPECTED {extra}" if extra else ""))
ok &= not extra and "bulk_load_sha256" in diff
out.append(f"relabel: {meta.get('relabel453', {}).get('categories_after')} cells {meta.get('relabel453', {}).get('cells')}")
cl = (W / "chain497.log").read_text()
nov = re.findall(r"novelty rc=(\d+)", cl); pk = re.findall(r"package (\S+) sha256 ([0-9a-f]{64})", cl)
h = hashlib.sha256()
with open(PKG, "rb") as f:
    for chunk in iter(lambda: f.read(1 << 24), b""):
        h.update(chunk)
sha_ok = bool(pk) and pk[-1][1] == h.hexdigest()
ok &= bool(nov) and nov[-1] == "0" and sha_ok
out.append(f"novelty rc={nov[-1] if nov else 'missing'}")
out.append(f"package sha256: {h.hexdigest()}" + ("" if sha_ok else "  (MISMATCH with chain497.log)"))
print("\n".join(out)); print("ALL GATES", "PASS" if ok else "FAIL")
if "--write" in sys.argv:
    if not ok:
        print("not writing gates497.log (a gate failed)"); sys.exit(1)
    (W / "gates497.log").write_text("\n".join(out) + "\n")
    print("wrote", W / "gates497.log")
