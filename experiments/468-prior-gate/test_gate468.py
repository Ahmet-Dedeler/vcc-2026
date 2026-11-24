"""468 unit tests for VCC_DE_PRIOR_GATE's parse on a builder file (default: the production builder).

The build-time behaviour (which targets fall back to CLEAN4) is proven by the known-answer builds in chain468.sh and by the
gated arm's audit; this file covers the parse: unset is a no-op, q<NN> parses, everything else is refused.
usage: uv run python experiments/468-prior-gate/test_gate468.py [builder.py]"""
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILDER = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "scripts/make_submission_pm_cd4.py"
KA = ROOT / "experiments/458-arc-lab-programme/work/arch1d_programme_abc.npz"
SHA = __import__("hashlib").sha256(KA.read_bytes()).hexdigest()
assert SHA.startswith("6951dbc7"), SHA


def load(env):
    for k in [k for k in os.environ if k.startswith("VCC_")]:
        del os.environ[k]
    os.environ.update(env)
    name = f"b468_{abs(hash(tuple(sorted(env.items()))))}"
    spec = importlib.util.spec_from_file_location(name, BUILDER)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def refused(env):
    code = ("import importlib.util,sys; s=importlib.util.spec_from_file_location('x', sys.argv[1]); "
            "m=importlib.util.module_from_spec(s); sys.modules['x']=m; s.loader.exec_module(m)")
    e = {k: v for k, v in os.environ.items() if not k.startswith("VCC_")}
    e.update(env)
    r = subprocess.run([sys.executable, "-c", code, str(BUILDER)], env=e, capture_output=True, text=True)
    return r.returncode != 0, (r.stderr.strip().splitlines() or [""])[-1][:160]


ok = True
m0 = load({})
assert m0.DE_PRIOR_GATE is None and m0.DE_PRIOR_GATE_AUDIT == {}, "unset must leave DE_PRIOR_GATE None"
print("PASS unset -> DE_PRIOR_GATE None")
base = {"VCC_DE_SIGN_PRIOR": "clean4q80", "VCC_SAMELINE": f"{KA}:{SHA}:0"}
m1 = load(base)
assert m1.DE_PRIOR_GATE is None and m1.SAMELINE["prog"] is True, "SAMELINE alone must leave the gate off"
print("PASS VCC_SAMELINE alone -> DE_PRIOR_GATE None")
m2 = load({**base, "VCC_DE_PRIOR_GATE": "q33"})
assert m2.DE_PRIOR_GATE == {"q": 0.33, "tau": 0.2}, m2.DE_PRIOR_GATE
m3 = load({**base, "VCC_DE_PRIOR_GATE": " Q50 "})
assert m3.DE_PRIOR_GATE == {"q": 0.5, "tau": 0.2}, m3.DE_PRIOR_GATE
m4 = load({**base, "VCC_SAMELINE": f"{KA}:{SHA}:0:pt", "VCC_DE_PRIOR_GATE": "q1"})
assert m4.DE_PRIOR_GATE == {"q": 0.01, "tau": 0.2} and m4.SAMELINE["pt"] is True
print("PASS parse: q33, ' Q50 ' (trimmed, case-folded), q1 with pt")
for env, why in (({**base, "VCC_DE_PRIOR_GATE": "q0"}, "q0"),
                 ({**base, "VCC_DE_PRIOR_GATE": "q100"}, "q100"),
                 ({**base, "VCC_DE_PRIOR_GATE": "x33"}, "bad prefix"),
                 ({**base, "VCC_DE_PRIOR_GATE": "33"}, "no prefix"),
                 ({**base, "VCC_DE_PRIOR_GATE": "q3.3"}, "decimal"),
                 ({**base, "VCC_SAMELINE": f"{KA}:{SHA}:0:noprog", "VCC_DE_PRIOR_GATE": "q33"}, "noprog"),
                 ({"VCC_DE_SIGN_PRIOR": "clean4q80", "VCC_DE_PRIOR_GATE": "q33"}, "no VCC_SAMELINE")):
    r, msg = refused(env)
    ok &= r
    print(("PASS" if r else "FAIL") + f" refused: {why} ({msg})")
print("ALL PASS" if ok else "SOME FAILED")
sys.exit(0 if ok else 1)
