"""517: riders_match on the 517-patched def_day_helpers.py COPY: every form of the requested VCC_DE_MDE_SIGN string against
metas carrying the builder's config dict (as patch_builder517.py writes it), including mismatches and 502's forms.
usage: uv run python experiments/517-flex-programme/test_defday517.py <ddh517.py>"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

H = Path(sys.argv[1])
SHA = "6951dbc7498e1f2aeabf405d5cbe94e7b695a7d3c0848a3a9460daa52d87a5ba"


def run(cfg_sign, request, mag="-"):
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
        json.dump({"config": ({} if cfg_sign is None else {"de_mde_sign": cfg_sign})}, f)
    return subprocess.run([sys.executable, str(H), "riders_match", f.name, mag, request], capture_output=True, text=True).returncode


cases = [
    ({"source": "clean4", "tau": 0.05}, "clean4:0.05", 0, "502 form"),
    ({"source": "clean4", "tau": 0.05}, "clean4:0.1", 1, "502 form, tau mismatch"),
    (None, "-", 0, "absent, not requested"),
    ({"source": "clean4", "tau": 0.05}, "-", 1, "present, not requested"),
    ({"source": "programme", "tau": 0.05, "tau_c4": 0.05, "contexts": {}}, "programme:0.05", 0, "programme, default tau_c4"),
    ({"source": "programme", "tau": 0.05, "tau_c4": 0.0, "contexts": {}}, "programme:0.05:0", 0, "programme, tau_c4 0"),
    ({"source": "programme", "tau": 0.05, "tau_c4": 0.05, "contexts": {}}, "programme:0.05:0", 1, "programme, tau_c4 mismatch"),
    ({"source": "clean4", "tau": 0.05}, "programme:0.05", 1, "source mismatch"),
    ({"source": "file", "path": "/w/hookprog.npz", "sha256": SHA, "tau": 0.05, "tau_c4": 0.05, "le": 1, "contexts": {}},
     f"file:/w/hookprog.npz:{SHA}:0.05:0.05:le1", 0, "file with le1"),
    ({"source": "file", "path": "/w/hookprog.npz", "sha256": SHA, "tau": 0.05, "tau_c4": 0.05, "contexts": {}},
     f"file:/w/hookprog.npz:{SHA}:0.05", 0, "file, defaults"),
    ({"source": "file", "path": "/w/hookprog.npz", "sha256": SHA, "tau": 0.05, "tau_c4": 0.05, "contexts": {}},
     f"file:/w/hookprog.npz:{SHA}:0.05:0.05:le1", 1, "file, le requested but absent"),
    ({"source": "file", "path": "/w/hookprog.npz", "sha256": "0" * 64, "tau": 0.05, "tau_c4": 0.05, "contexts": {}},
     f"file:/w/hookprog.npz:{SHA}:0.05", 1, "file, sha mismatch"),
    ({"source": "file", "path": "/w/hookprog.npz", "sha256": SHA, "tau": 0.05, "tau_c4": 0.05}, "file:x", 1, "unparseable request"),
]
ok = 0
for cfg, req, want, why in cases:
    got = run(cfg, req)
    assert got == want, (why, got, want)
    ok += 1
    print(f"PASS {why}: rc {got}")
print(f"ALL {ok} PASS")
