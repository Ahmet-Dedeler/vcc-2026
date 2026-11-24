#!/usr/bin/env python3
"""701 deliverability (PREREGISTRATION, "Deliverability"): from build701's VCC_RC_AUDIT_DUMP (work/ku30.rc_audit.jsonl),
prescribed vs delivered call-set sizes per (context, target). If delivered calls grow by less than half of the prescribed
+800 over K 2,600 (i.e. mean delivered < 3,000), the arm under-delivers, and a NULL / NEG board read means "capacity-bound",
not "size does not pay".
usage: uv run python experiments/701-kup-callset/deliver701.py [audit.jsonl]
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
p = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "work/ku30.rc_audit.jsonl"
rows = [json.loads(l) for l in p.read_text().splitlines() if l.strip()]
out = {}
for c in sorted({r["context"] for r in rows}):
    R = [r for r in rows if r["context"] == c and "call_set_prescribed" in r]
    pr = np.array([r["call_set_prescribed"] for r in R], float)
    de = np.array([r["call_set_delivered"] for r in R], float)
    out[c] = {"n": len(R), "prescribed_mean": float(pr.mean()), "delivered_mean": float(de.mean()),
              "delivered_q10": float(np.quantile(de, .1)), "increment_delivered_frac": float((de.mean() - 2600) / 800)}
allde = np.mean([v["delivered_mean"] for v in out.values()])
verdict = "DELIVERS" if allde >= 3000 else "UNDER-DELIVERS (a NULL/NEG read is capacity-bound)"
print(json.dumps({"per_context": out, "delivered_mean_all": float(allde), "verdict": verdict}, indent=1))
