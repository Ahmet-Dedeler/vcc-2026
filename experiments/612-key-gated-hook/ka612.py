#!/usr/bin/env python3
"""612 KA-c checker: 496a's env (446bj2's + VCC_DE_MDE_SIGN=clean4:0.05) WITH VCC_DE_MDE_SIGN_KEYGATE=0.5 must differ from
496a's artifact ONLY through the key gate. Pass conditions:
  1. config diff vs 496a's meta is exactly {de_mde_sign_keygate} (plus de_sign_emitter_sha256 if the emitter moved);
  2. de_mde_sign_keygate.q == 0.5, and de_mde_sign is unchanged ({"source": "clean4", "tau": 0.05});
  3. the output sha256 differs from 496a's (e21ff752...): the gate kept the delivered sign on strong-key genes;
  4. the driver exited 0 (run_ka_527.sh's log), checked by the caller.
usage: uv run python experiments/612-key-gated-hook/ka612.py <KA612c meta.json> [<496a meta.json>]"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / ("submissions/496a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-"
              "mde010coarse-zr5-as05.h5ad.meta.json")
a = json.load(open(sys.argv[1]))
b = json.load(open(sys.argv[2] if len(sys.argv) > 2 else REF))
ca, cb = a["config"], b["config"]
diff = sorted(k for k in set(ca) | set(cb) if ca.get(k) != cb.get(k))
ok1 = set(diff) <= {"de_mde_sign_keygate", "de_sign_emitter_sha256"} and "de_mde_sign_keygate" in diff
ok2 = (ca.get("de_mde_sign_keygate", {}).get("q") == 0.5
       and {k: v for k, v in ca.get("de_mde_sign", {}).items() if k != "contexts"} == {"source": "clean4", "tau": 0.05})
ok3 = a["output"]["sha256"] != b["output"]["sha256"]
print(f"config diff vs 496a: {diff} -> {'OK' if ok1 else 'UNEXPECTED'}")
print(f"de_mde_sign_keygate: {ca.get('de_mde_sign_keygate')}; de_mde_sign: {ca.get('de_mde_sign')} -> {'OK' if ok2 else 'MISMATCH'}")
print(f"output {a['output']['sha256'][:16]} vs 496a {b['output']['sha256'][:16]} -> {'differs (OK)' if ok3 else 'IDENTICAL (FAIL)'}")
print("KA612c PASS" if (ok1 and ok2 and ok3) else "KA612c FAIL")
sys.exit(0 if (ok1 and ok2 and ok3) else 1)
