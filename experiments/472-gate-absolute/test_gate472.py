"""472 unit tests: `VCC_DE_PRIOR_GATE=le<N>` beside `q<NN>` (default: the production builder).
usage: uv run python experiments/472-gate-absolute/test_gate472.py [builder.py]"""
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILDER = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "scripts/make_submission_pm_cd4.py"
KA = ROOT / "experiments/458-arc-lab-programme/work/arch1d_programme_abc.npz"
SHA = __import__("hashlib").sha256(KA.read_bytes()).hexdigest()


def load(env):
    for k in [k for k in os.environ if k.startswith("VCC_")]:
        del os.environ[k]
    os.environ.update(env)
    name = f"b472_{abs(hash(tuple(sorted(env.items()))))}"
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
base = {"VCC_DE_SIGN_PRIOR": "clean4q80", "VCC_SAMELINE": f"{KA}:{SHA}:0"}
assert load({}).DE_PRIOR_GATE is None
for val, want in (("q33", {"q": 0.33, "tau": 0.2}), ("le1", {"le": 1, "tau": 0.2}), ("le0", {"le": 0, "tau": 0.2}),
                  ("LE12", {"le": 12, "tau": 0.2}), (" le3 ", {"le": 3, "tau": 0.2})):
    got = load({**base, "VCC_DE_PRIOR_GATE": val}).DE_PRIOR_GATE
    good = got == want
    ok &= good
    print(("PASS" if good else "FAIL") + f" parse {val!r} -> {got}")
for env, why in (({**base, "VCC_DE_PRIOR_GATE": "le"}, "le without a number"),
                 ({**base, "VCC_DE_PRIOR_GATE": "lex"}, "le with letters"),
                 ({**base, "VCC_DE_PRIOR_GATE": "le100"}, "le with 3 digits"),
                 ({**base, "VCC_DE_PRIOR_GATE": "le-1"}, "negative"),
                 ({**base, "VCC_DE_PRIOR_GATE": "lq1"}, "bad prefix"),
                 ({**base, "VCC_DE_PRIOR_GATE": "q0"}, "q0"),
                 ({**base, "VCC_DE_PRIOR_GATE": "q100"}, "q100"),
                 ({**base, "VCC_SAMELINE": f"{KA}:{SHA}:0:noprog", "VCC_DE_PRIOR_GATE": "le1"}, "le1 with noprog"),
                 ({"VCC_DE_SIGN_PRIOR": "clean4q80", "VCC_DE_PRIOR_GATE": "le1"}, "le1 without VCC_SAMELINE")):
    r, msg = refused(env)
    ok &= r
    print(("PASS" if r else "FAIL") + f" refused: {why} ({msg})")
print("ALL PASS" if ok else "SOME FAILED")
sys.exit(0 if ok else 1)
