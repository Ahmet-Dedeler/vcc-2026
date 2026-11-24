"""437: T6 guard parse tests (import only). TEST_BUILDER=<path> tests a dry copy."""
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILDER = Path(os.environ.get("TEST_BUILDER", ROOT / "scripts/make_submission_pm_cd4.py"))
CODE = (f"import importlib.util as u, sys; s=u.spec_from_file_location('b', {str(BUILDER)!r}); m=u.module_from_spec(s); "
        "sys.modules['b']=m; s.loader.exec_module(m); print('OK', m.CD4_SIZING, m.CD4_TARGET_F, m.CD4_TARGET_G)")


def run(env):
    e = {k: v for k, v in os.environ.items() if not k.startswith("VCC_")}
    e.update(env)
    return subprocess.run([sys.executable, "-c", CODE], env=e, capture_output=True, text=True, cwd=str(ROOT))


G = "0.3443655278649641"
cases = [
    (dict(), True, "unset: A/B/C legacy path"),
    (dict(VCC_CD4_TARGET_F="0.51"), True, "A/B/C: TARGET_F 0.51 without carrying (the 417 A/B/C dump recipe) still works"),
    (dict(VCC_CD4_SIZING="carrying", VCC_CD4_TARGET_G=G), True, "D/E/F: carrying + G"),
    (dict(VCC_CD4_SIZING="carrying", VCC_CD4_TARGET_G=G, VCC_CD4_TARGET_F="0.34"), True, "D/E/F: + the driver's F 0.34 (inert)"),
    (dict(VCC_CD4_SIZING="carrying", VCC_CD4_TARGET_G="0.5165482917974462", VCC_CD4_TARGET_F="0.34"), True, "D/E/F: CD4 x1.5 the right way (G)"),
    (dict(VCC_CD4_SIZING="carrying", VCC_CD4_TARGET_G=G, VCC_CD4_TARGET_F="0.51"), False, "D/E/F: CD4 x1.5 via F -> REFUSE (T6)"),
]
fails = 0
for env, ok_want, what in cases:
    r = run(env)
    ok = r.returncode == 0 and "OK" in r.stdout
    good = ok == ok_want and (ok_want or "432 T6" in (r.stderr + r.stdout))
    fails += not good
    print(f"{'PASS' if good else 'FAIL'}  {what}: rc {r.returncode} {(r.stdout.strip() or r.stderr.strip()[-160:])[:160]}")
print(f"{len(cases) - fails}/{len(cases)} PASS")
sys.exit(1 if fails else 0)
