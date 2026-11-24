"""451 unit tests for VCC_SAMELINE on a builder file (default: the production builder).
usage: uv run python experiments/451-sameline-mode/test_sameline451.py [builder.py]"""
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BUILDER = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "scripts/make_submission_pm_cd4.py"
KA = ROOT / "experiments/451-sameline-mode/work/sameline_ka_A.npz"


def load(env):
    for k in [k for k in os.environ if k.startswith("VCC_")]:
        del os.environ[k]
    os.environ.update(env)
    name = f"b451_{abs(hash(tuple(sorted(env.items()))))}"
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


sha = __import__("hashlib").sha256(KA.read_bytes()).hexdigest()
ok = True
m0 = load({})
assert m0.SAMELINE is None and m0.SAMELINE_AUDIT == {}, "unset must leave SAMELINE None"
print("PASS unset -> SAMELINE None")
m1 = load({"VCC_DE_SIGN_PRIOR": "clean4q80", "VCC_SAMELINE": f"{KA}:{sha}:0"})
assert m1.SAMELINE["w"] == 0.0 and m1.SAMELINE["pt"] is False and Path(m1.SAMELINE["path"]) == KA
m2 = load({"VCC_DE_SIGN_PRIOR": "clean4q80", "VCC_SAMELINE": f"experiments/451-sameline-mode/work/sameline_ka_A.npz:{sha}:1.5:pt"})
assert m2.SAMELINE["w"] == 1.5 and m2.SAMELINE["pt"] is True and m2.SAMELINE["prog"] is True and Path(m2.SAMELINE["path"]) == KA
m3 = load({"VCC_DE_SIGN_PRIOR": "clean4q80", "VCC_SAMELINE": f"{KA}:{sha}:1:pt:noprog"})
assert m3.SAMELINE["pt"] is True and m3.SAMELINE["prog"] is False
m4 = load({"VCC_DE_SIGN_PRIOR": "clean4q80", "VCC_SAMELINE": f"{KA}:{sha}:1:noprog"})
assert m4.SAMELINE["pt"] is False and m4.SAMELINE["prog"] is False
assert m1.SAMELINE["prog"] is True
print("PASS parse: absolute / ROOT-relative path, w, pt, noprog")
for env, why in (({"VCC_DE_SIGN_PRIOR": "clean4q80ja", "VCC_SAMELINE": f"{KA}:{sha}:1"}, "ja variant"),
                 ({"VCC_DE_SIGN_PRIOR": "clean4q80", "VCC_SAMELINE": f"{KA}:{sha}:11"}, "w > 10"),
                 ({"VCC_DE_SIGN_PRIOR": "clean4q80", "VCC_SAMELINE": f"{KA}:abc:1"}, "short sha"),
                 ({"VCC_DE_SIGN_PRIOR": "clean4q80", "VCC_SAMELINE": f"{KA}:{sha}:1:xx"}, "bad flag"),
                 ({"VCC_DE_SIGN_PRIOR": "clean4q80", "VCC_SAMELINE": f"{KA}:{sha}:1:pt:pt"}, "repeated flag")):
    r, msg = refused(env)
    ok &= r
    print(("PASS" if r else "FAIL") + f" refused: {why} ({msg})")
# loader: sha, gene axis, panel
import csv
genes = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
perts = [r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))]
rec = m1.sameline_load(genes, perts)
z = np.load(ROOT / "experiments/396-jurkat-tilt-A/work/prior_ctx_challenge.npz")
assert list(rec) == ["A"] and np.array_equal(rec["A"]["programme"], z["psign_A"]) and len(rec["A"]["targets"]) == 0
print("PASS loader: A only, programme == 396 psign_A, 0 covered targets")
assert bool(rec["A"]["programme_ok"]) is True   # absent key (the 09-25 KA artifact) defaults to True
EF = ROOT / "experiments/451-sameline-mode/work/sameline_efr_forced.npz"
if EF.exists():
    import csv as _csv
    eg = [r["gene_name"] for r in _csv.DictReader(open(ROOT / "experiments/238-essential-free-panel/work/ef-bundle-r/gene_names.csv"))]
    ep = [r["target_gene"] for r in _csv.DictReader(open(ROOT / "experiments/238-essential-free-panel/work/ef-bundle-r/pert_counts.csv"))]
    esha = __import__("hashlib").sha256(EF.read_bytes()).hexdigest()
    m5 = load({"VCC_DE_SIGN_PRIOR": "clean4q80", "VCC_SAMELINE": f"{EF}:{esha}:1:pt"})
    er = m5.sameline_load(eg, ep)
    okf = bool(er["D"]["programme_ok"]) and bool(er["E"]["programme_ok"]) and not bool(er["F"]["programme_ok"])
    ok &= okf
    print(("PASS" if okf else "FAIL") + f" programme_ok on the forced rehearsal artifact: D {bool(er['D']['programme_ok'])} "
          f"E {bool(er['E']['programme_ok'])} F(A549, stimulated banks) {bool(er['F']['programme_ok'])}; covered "
          f"{[len(er[c]['targets']) for c in 'DEF']}")
try:
    m1.sameline_load(genes, perts[::-1]); ok = False; print("FAIL loader accepted a different panel order")
except SystemExit as e:
    print(f"PASS loader refuses another panel ({e})")
# blend: synthetic
rng = np.random.default_rng(0)
G, P = 50, 6
M = rng.normal(size=(P, G)); panel = np.zeros(G, bool); panel[:6] = True
rows = rng.normal(size=(2, G)); rows[:, panel] = 0; rows /= np.linalg.norm(rows[:, ~panel], axis=1, keepdims=True)
rec_s = {"targets": np.array([1, 4]), "rows": rows.astype(np.float32), "signs": np.sign(rows).astype(np.int8)}
assert m1.sameline_blend(M, rec_s, panel, 0.0) is M
B = m1.sameline_blend(M, rec_s, panel, 1.0)
assert np.array_equal(B[[0, 2, 3, 5]], M[[0, 2, 3, 5]]) and np.array_equal(B[:, panel], M[:, panel])
for k, ti in enumerate([1, 4]):
    assert abs(np.linalg.norm(B[ti, ~panel]) - np.linalg.norm(M[ti, ~panel])) < 1e-10
    u = lambda v: v / np.linalg.norm(v)
    want = u(u(M[ti, ~panel]) + u(rows[k, ~panel].astype(np.float32).astype(np.float64)))
    assert np.allclose(u(B[ti, ~panel]), want)
BB = m1.sameline_blend(M, rec_s, panel, 1e6)
assert np.corrcoef(BB[1, ~panel], rows[0, ~panel])[0, 1] > 0.999
print("PASS blend: w 0 is the identity object; uncovered rows and panel columns untouched; norms kept; direction = unit sum")
# 451b (per-context w, only when the builder supports `auto`)
if "else float(_slr[\"w\"])" in BUILDER.read_text() or 'float(_slr["w"])' in BUILDER.read_text():
    AUTO = ROOT / "experiments/451-sameline-mode/work/sameline_efr_auto.npz"
    if AUTO.exists():
        asha = __import__("hashlib").sha256(AUTO.read_bytes()).hexdigest()
        m6 = load({"VCC_DE_SIGN_PRIOR": "clean4q80", "VCC_SAMELINE": f"{AUTO}:{asha}:auto:pt:noprog"})
        assert m6.SAMELINE["w"] is None and m6.SAMELINE["pt"] and not m6.SAMELINE["prog"]
        import csv as _csv2
        eg = [r["gene_name"] for r in _csv2.DictReader(open(ROOT / "experiments/238-essential-free-panel/work/ef-bundle-r/gene_names.csv"))]
        ep = [r["target_gene"] for r in _csv2.DictReader(open(ROOT / "experiments/238-essential-free-panel/work/ef-bundle-r/pert_counts.csv"))]
        er = m6.sameline_load(eg, ep)
        ws = {c: float(er[c]["w"]) for c in er}
        okw = ws == {"D": 2.0, "E": 2.0, "F": 0.5}
        ok &= okw
        print(("PASS" if okw else "FAIL") + f" 451b auto: per-context w from the artifact {ws}")
        r, msg = refused({"VCC_DE_SIGN_PRIOR": "clean4q80", "VCC_SAMELINE": f"{KA}:{sha}:auto"})
        # the KA artifact has no w_<c>: loading it under auto must refuse (at load time, inside the build)
        m7 = load({"VCC_DE_SIGN_PRIOR": "clean4q80", "VCC_SAMELINE": f"{KA}:{sha}:auto"})
        try:
            m7.sameline_load(genes, perts); ok = False; print("FAIL 451b: auto accepted an artifact without w_<c>")
        except SystemExit as e:
            print(f"PASS 451b: auto refuses an artifact without w_<c> ({str(e)[:70]})")
print("ALL PASS" if ok else "SOME FAILED")
sys.exit(0 if ok else 1)
