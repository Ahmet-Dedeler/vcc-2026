#!/usr/bin/env python3
"""469dr's gates (the October path + Arc's denoised programme + the strength gate q33, vs 461dr), checked after chain469
ENDs; writes work/gates469.log with --write. PASS requires:
  - def_day converged with every final rho <= 1, and its re-level factors;
  - the programme source h1_2025_full on D/E/F, flags auto:pt, and the prior gate q33 logged by def_day AND by the builds;
  - the zero-gene audit at 0;
  - the relabelled meta's config differs from 461dr's ONLY in `de_prior_gate` (q 0.33, tau 0.2, D/E/F audited) and in
    `sameline`'s artifact path (the same-line artifact's sha256 and every per-context record must be 461dr's).
    `bulk_ctx_scale` is flagged, not failed: the gate does not touch the bulk, so the re-level should retrace 461dr's;
  - the stack dump's arrays identical to 461dr's;
  - chain469's novelty rc 0 and a package sha that matches the file.
usage: uv run python experiments/469-gated-october-path/gate469.py [--write]"""
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
X = ROOT / "experiments/469-gated-october-path"; W = X / "work"; RUN = W / "run"
PKG = ROOT / "submissions/469dr-october-path-arcprog-gate33-abc-q1.75.vcc"
H5 = ROOT / "submissions/469dr-october-path-arcprog-gate33-abc-q1.75.h5ad"
REF = ROOT / "submissions/461dr-october-path-arcprog-abc-q1.75.h5ad.meta.json"
REF_DUMP = ROOT / "experiments/461-arcprog-october-path/work/run/dump_h1x2_cd4x1.npz"
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
out.append(f"[469dr] {'CLIFF OK' if cliff_ok else 'CLIFF FAIL'} {json.dumps(rho)} final {final_name} factors {factors} "
           f"({'converged' if 'CONVERGED at' in dd else 'nearest safe build'})")
sl = np.load(RUN / "sameline.npz", allow_pickle=False)
prov = json.loads(str(sl["provenance"]))
src = {c: v.get("programme_source", {}).get("bank") for c, v in prov.get("contexts", {}).items()}
prog_ok = sorted(src) == ["D", "E", "F"] and set(src.values()) == {"h1_2025_full"} and "LODO key ON" not in dd
mode = [l for l in dd.splitlines() if "same-line mode ON" in l]
prog_ok &= all(bool(sl[f"programme_ok_{c}"]) for c in "DEF") and bool(mode) and mode[-1].rstrip().endswith(":auto:pt")
out.append(f"programme source: {json.dumps(src)} {'OK' if prog_ok else 'NOT h1_2025_full on D/E/F (or LODO on / noprog)'}")
ok &= prog_ok
gate_log = [l for l in dd.splitlines() if "prior gate q33" in l]
built = sorted({l.strip() for p in (RUN / "logs").glob("build*.log") for l in p.read_text().splitlines() if l.startswith("[prior_gate]")})
gate_ok = bool(gate_log) and {l.split(":")[0] for l in built} == {"[prior_gate] D", "[prior_gate] E", "[prior_gate] F"}
out.append(f"prior gate: def_day {'logged' if gate_log else 'NOT logged'}; builds: {' | '.join(built)[:400] or '<none>'}")
ok &= gate_ok
zl = [l for l in dd.splitlines() if "ZERO block total" in l]
zero_ok = bool(zl) and re.search(r"ZERO block total \(target excluded\): 0 in 0 blocks", zl[-1]) is not None
out.append(f"zero-gene audit: {'0 (OK)' if zero_ok else 'NOT 0 / unparsed'} :: {zl[-1][:160] if zl else '<no line>'}")
ok &= zero_ok
meta = json.loads(Path(str(H5) + ".meta.json").read_text()); ref = json.loads(REF.read_text())
c, r = dict(meta["config"]), dict(ref["config"])
diff = sorted(k for k in set(c) | set(r) if c.get(k) != r.get(k))
if "sameline" in diff:
    # the artifact's sha differs by construction: its provenance records the run's own identity.json path. Accept iff every
    # array is equal (NaN-aware: the programmes carry NaN for uncovered genes), the provenance differs only in
    # `identity_source`, and the config records differ only in `path` / `sha256`.
    a = {k: v for k, v in c["sameline"].items() if k not in ("path", "sha256")}
    b = {k: v for k, v in r["sameline"].items() if k not in ("path", "sha256")}
    na = np.load(RUN / "sameline.npz", allow_pickle=False)
    nb = np.load(ROOT / "experiments/461-arcprog-october-path/work/run/sameline.npz", allow_pickle=False)
    arr_ok = sorted(na.files) == sorted(nb.files) and all(
        k == "provenance" or (na[k].shape == nb[k].shape and np.array_equal(na[k], nb[k], equal_nan=na[k].dtype.kind == "f"))
        for k in na.files)
    pa = json.loads(str(na["provenance"])); pb = json.loads(str(nb["provenance"]))
    pa.pop("identity_source", None); pb.pop("identity_source", None)
    if a == b and arr_ok and pa == pb:
        diff.remove("sameline")
        out.append("sameline: identical to 461dr's (every array NaN-aware; provenance apart from identity_source; config "
                   f"apart from path/sha) -- artifact sha {c['sameline']['sha256'][:16]} vs {r['sameline']['sha256'][:16]}")
    else:
        out.append(f"sameline: DIFFERS from 461dr's (config equal {a == b}, arrays equal {arr_ok}, provenance equal {pa == pb})")
if "bulk_load_sha256" in diff:
    da = np.load(RUN / "dump_h1x2_cd4x1.npz", allow_pickle=True); db = np.load(REF_DUMP, allow_pickle=True)
    arrs = ("bulk", "carriers", "perts", "genes")
    if all(np.array_equal(da[k], db[k]) for k in arrs):
        diff.remove("bulk_load_sha256"); out.append(f"bulk_load: identical dump arrays {arrs} (the file differs only in its env record)")
    else:
        out.append("bulk_load: DUMP ARRAYS DIFFER from 461dr's")
pg = c.get("de_prior_gate") or {}
pg_ok = pg.get("q") == 0.33 and pg.get("tau") == 0.2 and pg.get("weak") == "clean4" and sorted(pg.get("contexts", {})) == ["D", "E", "F"]
out.append("de_prior_gate: " + json.dumps({k: pg.get(k) for k in ("q", "tau", "weak")}) + " "
           + json.dumps({k: {kk: v[kk] for kk in ("threshold", "strong_programme", "weak_clean4", "weak_kept_pt")}
                         for k, v in pg.get("contexts", {}).items()}) + ("" if pg_ok else "  NOT AS EXPECTED"))
ok &= pg_ok
extra = [k for k in diff if k not in ("de_prior_gate", "bulk_ctx_scale")]
out.append(f"[469dr] {'PASS' if not extra and 'de_prior_gate' in diff else 'FAIL'}: config diff vs 461dr {diff}"
           + (f"  UNEXPECTED {extra}" if extra else ""))
if "bulk_ctx_scale" in diff:
    out.append(f"NOTE bulk_ctx_scale differs from 461dr's ({c.get('bulk_ctx_scale')} vs {r.get('bulk_ctx_scale')}): the bulk "
               "should be untouched -- explain before submitting")
ok &= not extra and "de_prior_gate" in diff
out.append(f"relabel: {meta.get('relabel453', {}).get('categories_after')} cells {meta.get('relabel453', {}).get('cells')}")
cl = (W / "chain469.log").read_text()
nov = re.findall(r"novelty rc=(\d+)", cl); pk = re.findall(r"package (\S+) sha256 ([0-9a-f]{64})", cl)
h = hashlib.sha256()
with open(PKG, "rb") as f:
    for chunk in iter(lambda: f.read(1 << 24), b""):
        h.update(chunk)
sha_ok = bool(pk) and pk[-1][1] == h.hexdigest()
ok &= bool(nov) and nov[-1] == "0" and sha_ok
out.append(f"novelty rc={nov[-1] if nov else 'missing'}")
out.append(f"package sha256: {h.hexdigest()}" + ("" if sha_ok else "  (MISMATCH with chain469.log)"))
print("\n".join(out)); print("ALL GATES", "PASS" if ok else "FAIL")
if "--write" in sys.argv:
    if not ok:
        print("not writing gates469.log (a gate failed)"); sys.exit(1)
    (W / "gates469.log").write_text("\n".join(out) + "\n")
    print("wrote", W / "gates469.log")
