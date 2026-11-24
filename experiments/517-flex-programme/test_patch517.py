"""517 unit tests for patch_builder517.py on a builder COPY (never the production file): the parse of VCC_DE_MDE_SIGN's
new `programme:<tau>[:<tau_c4>]` form, that `clean4:<tau>` parses exactly as 496's patched builder does, the refusals,
and that sameline_load carries the artifact's pdown_<c> (and nothing else changes in its record).
usage: uv run python experiments/517-flex-programme/test_patch517.py <bld517.py> <bld496.py>"""
import hashlib
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
B517, B496 = Path(sys.argv[1]), Path(sys.argv[2])
KA = ROOT / "experiments/458-arc-lab-programme/work/arch1d_programme_abc.npz"
SHA = hashlib.sha256(KA.read_bytes()).hexdigest()
assert SHA.startswith("6951dbc7"), SHA
MDE = {"VCC_DE_SIGN_PRIOR": "clean4q80", "VCC_DE_MDE": "0.1,0.1,0.3,coarse"}
SL = {"VCC_SAMELINE": f"{KA}:{SHA}:0"}


def load(path, env):
    for k in [k for k in os.environ if k.startswith("VCC_")]:
        del os.environ[k]
    os.environ.update(env)
    name = f"b517_{path.stem}_{abs(hash(tuple(sorted(env.items()))))}"
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def refused(path, env):
    code = ("import importlib.util,sys; s=importlib.util.spec_from_file_location('x', sys.argv[1]); "
            "m=importlib.util.module_from_spec(s); sys.modules['x']=m; s.loader.exec_module(m)")
    e = {k: v for k, v in os.environ.items() if not k.startswith("VCC_")}
    e.update(env)
    r = subprocess.run([sys.executable, "-c", code, str(path)], env=e, capture_output=True, text=True)
    return r.returncode != 0, (r.stderr.strip().splitlines() or [""])[-1][:160]


n_pass = 0
m0 = load(B517, {})
assert m0.DE_MDE_SIGN is None and m0.DE_MDE_SIGN_AUDIT == {}
n_pass += 1; print("PASS unset -> DE_MDE_SIGN None, no audit")
for env in ({**MDE, "VCC_DE_MDE_SIGN": "clean4:0.05"}, {**MDE, "VCC_DE_MDE_SIGN": " CLEAN4:0 "},
            {**MDE, **SL, "VCC_DE_MDE_SIGN": "clean4:0.1"}):
    a, b = load(B517, env), load(B496, env)
    assert a.DE_MDE_SIGN == b.DE_MDE_SIGN, (env, a.DE_MDE_SIGN, b.DE_MDE_SIGN)
    n_pass += 1; print(f"PASS {env['VCC_DE_MDE_SIGN']!r} parses exactly as 496's: {a.DE_MDE_SIGN}")
m = load(B517, {**MDE, **SL, "VCC_DE_MDE_SIGN": "programme:0.05"})
assert m.DE_MDE_SIGN == {"source": "programme", "tau": 0.05, "tau_c4": 0.05}, m.DE_MDE_SIGN
n_pass += 1; print("PASS programme:0.05 ->", m.DE_MDE_SIGN)
m = load(B517, {**MDE, **SL, "VCC_DE_MDE_SIGN": "programme:0.1:0", "VCC_DE_PRIOR_GATE": "le1"})
assert m.DE_MDE_SIGN == {"source": "programme", "tau": 0.1, "tau_c4": 0.0}, m.DE_MDE_SIGN
n_pass += 1; print("PASS programme:0.1:0 with the le1 gate ->", m.DE_MDE_SIGN)
m = load(B517, {**MDE, "VCC_DE_MDE_SIGN": f"file:{KA}:{SHA}:0.05"})
assert m.DE_MDE_SIGN == {"source": "file", "path": str(KA), "sha256": SHA, "tau": 0.05, "tau_c4": 0.05}, m.DE_MDE_SIGN
n_pass += 1; print("PASS file:<458d artifact>:<sha>:0.05 without VCC_SAMELINE ->", {k: v for k, v in m.DE_MDE_SIGN.items() if k != "path"})
m = load(B517, {**MDE, **SL, "VCC_DE_MDE_SIGN": f"file:{KA}:{SHA}:0.1:0"})
assert m.DE_MDE_SIGN["tau"] == 0.1 and m.DE_MDE_SIGN["tau_c4"] == 0.0, m.DE_MDE_SIGN
n_pass += 1; print("PASS file:...:0.1:0 with VCC_SAMELINE")
m = load(B517, {**MDE, "VCC_DE_MDE_SIGN": f"file:{KA}:{SHA}:0.05:0.05:le1"})
assert m.DE_MDE_SIGN.get("le") == 1 and m.DE_MDE_SIGN["tau"] == 0.05, m.DE_MDE_SIGN
n_pass += 1; print("PASS file:...:0.05:0.05:le1 (the hook's own gate) ->", {k: v for k, v in m.DE_MDE_SIGN.items() if k not in ("path", "sha256")})
for env, why in (({**MDE, "VCC_DE_MDE_SIGN": "programme:0.05"}, "no VCC_SAMELINE"),
                 ({**MDE, "VCC_DE_MDE_SIGN": f"file:{KA}:{SHA}:0.05:0.05:q33"}, "file: 6th field not le<N>"),
                 ({**MDE, **SL, "VCC_DE_PRIOR_GATE": "le1", "VCC_DE_MDE_SIGN": f"file:{KA}:{SHA}:0.05:0.05:le1"}, "file le<N> with VCC_DE_PRIOR_GATE"),
                 ({**MDE, "VCC_DE_MDE_SIGN": f"file:{KA}:{SHA[:20]}:0.05"}, "file: short sha"),
                 ({**MDE, "VCC_DE_MDE_SIGN": f"file:{KA}:{SHA}"}, "file: no tau"),
                 ({**MDE, "VCC_DE_MDE_SIGN": f"file:{KA}:{SHA}:0.6"}, "file: tau 0.6"),
                 ({**MDE, "VCC_DE_MDE_SIGN": f"file:{KA}:{SHA}:x"}, "file: bad tau"),
                 ({**MDE, "VCC_SAMELINE": f"{KA}:{SHA}:0:noprog", "VCC_DE_MDE_SIGN": "programme:0.05"}, "noprog"),
                 ({**MDE, **SL, "VCC_DE_MDE_SIGN": "programme:0.5"}, "tau 0.5"),
                 ({**MDE, **SL, "VCC_DE_MDE_SIGN": "programme:0.05:0.7"}, "tau_c4 0.7"),
                 ({**MDE, **SL, "VCC_DE_MDE_SIGN": "programme"}, "no tau"),
                 ({**MDE, **SL, "VCC_DE_MDE_SIGN": "programme:x"}, "bad tau"),
                 ({**MDE, **SL, "VCC_DE_MDE_SIGN": "programme:0.05:0.05:1"}, "4 parts"),
                 ({**MDE, **SL, "VCC_DE_MDE_SIGN": "clean4:0.05:0.05"}, "clean4 with 3 parts"),
                 ({**MDE, **SL, "VCC_DE_MDE_SIGN": "arc:0.05"}, "unknown source"),
                 ({"VCC_DE_SIGN_PRIOR": "clean4q80", **SL, "VCC_DE_MDE_SIGN": "programme:0.05"}, "no VCC_DE_MDE")):
    ok, msg = refused(B517, env)
    assert ok, (why, msg)
    n_pass += 1; print(f"PASS refused ({why}): {msg}")
# sameline_load: the record gains 'pdown' (== the artifact's pdown_<c>) and nothing else changes
import csv  # noqa: E402
genes = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
perts = [r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))]
a = load(B517, {**MDE, **SL}); b = load(B496, {**MDE, **SL})
ra, rb = a.sameline_load(genes, perts), b.sameline_load(genes, perts)
z = np.load(KA, allow_pickle=False)
for c in "ABC":
    assert set(ra[c]) - set(rb[c]) == {"pdown"}, set(ra[c]) ^ set(rb[c])
    assert all(np.array_equal(ra[c][k], rb[c][k]) for k in rb[c])
    assert np.array_equal(ra[c]["pdown"], np.asarray(z[f"pdown_{c}"], dtype=np.float64), equal_nan=True)
n_pass += 1; print("PASS sameline_load: + 'pdown' (== artifact pdown_<c>), every other key identical to 496's")
print(f"ALL {n_pass} PASS")
