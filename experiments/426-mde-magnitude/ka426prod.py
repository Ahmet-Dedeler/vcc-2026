"""426 production known answers (the patched src/vcc2026/rc_emit_k.py through run426.py, RUN426_EMITTER=prod):
KA1: knob off  == 416's base_s0_bps (the shipped 410b construction), every stats array bit-equal;
KA2: knob on (c 0.1, cap 0.1, min 0.3) == k_emit426's base_s0_bps_mde0.1cap0.1min0.3, every stats array bit-equal."""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
pairs = {"KA1": (HERE / "work/h1/stats-base_s0_bps_prodem.npz", HERE.parent / "416-de-carrier-menu/work/h1/stats-base_s0_bps.npz"),
         "KA2": (HERE / "work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3.npz", HERE / "work/h1/stats-base_s0_bps_mde0.1cap0.1min0.3.npz")}
out, ok_all = {}, True
for name, (a_, b_) in pairs.items():
    a, b = np.load(a_), np.load(b_)
    eq = {k: bool(np.array_equal(a[k], b[k], equal_nan=(k != "targets"))) for k in ("targets", "sums", "jk", "p", "q", "lfc")}
    ok = all(eq.values())
    ok_all &= ok
    out[name] = dict(pass_=ok, arrays=eq, a=str(a_.relative_to(ROOT)), b=str(b_.relative_to(ROOT)))
    print(name, "PASS" if ok else "FAIL", eq)
out["emitter_sha256"] = hashlib.sha256((ROOT / "src/vcc2026/rc_emit_k.py").read_bytes()).hexdigest()
out["builder_sha256"] = hashlib.sha256((ROOT / "scripts/make_submission_pm_cd4.py").read_bytes()).hexdigest()
(HERE / "work/known-answer-prod.json").write_text(json.dumps(out, indent=1))
print("emitter", out["emitter_sha256"], "builder", out["builder_sha256"])
sys.exit(0 if ok_all else 1)
