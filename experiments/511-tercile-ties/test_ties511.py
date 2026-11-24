"""511: the tie-rule fix changes nothing on A/B/C and fixes the tie cases.
Extracts `sharpen_q_by_tercile` from two builder sources (the old one and the 511-patched one) and compares:
  1. 446bj2's LODO dump and 495's cap-6 dump, tercile:2.5,2.0,1.75 and stercile:2.5,2.0,1.75: q vectors EQUAL;
  2. synthetic ties (501's cases): the new function gives the corrected splits.
usage: uv run python experiments/511-tercile-ties/test_ties511.py <old builder.py> <new builder.py>
"""
import ast
import sys

import numpy as np


def load_fn(path):
    src = open(path).read()
    tree = ast.parse(src)
    fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "sharpen_q_by_tercile")
    code = ast.get_source_segment(src, fn)
    ns = {"np": np}
    exec(code, ns)
    return ns["sharpen_q_by_tercile"]


old, new = load_fn(sys.argv[1]), load_fn(sys.argv[2])
bad = 0
PT = [{"mode": "tercile", "q": [2.5, 2.0, 1.75]}, {"mode": "stercile", "q": [2.5, 2.0, 1.75]}]
for dump in ("experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz", "experiments/495-lodo-cap6/work/dump_h1x2_cd4x15_lodo05hi6.npz"):
    z = np.load(dump, allow_pickle=True)
    genes, perts = [str(g) for g in z["genes"]], [str(p) for p in z["perts"]]
    for pt in PT:
        qa, aa = old(z["bulk"], z["carriers"], genes, perts, pt)
        qb, ab = new(z["bulk"], z["carriers"], genes, perts, pt)
        ok = np.array_equal(qa, qb)
        bad += not ok
        print(f"{'PASS' if ok else 'FAIL'}  {dump.split('/')[-1]} {pt['mode']}: old counts {aa['counts']} new {ab['counts']} (tie_to_mid {aa['tie_to_mid']} -> {ab['tie_to_mid']})")
# synthetic ties (501 section 6), 300 targets on a 50-gene axis
rng = np.random.default_rng(511)
genes = [f"g{i}" for i in range(50)]; perts = [f"p{i}" for i in range(300)]


def rows_with_norms(n):
    x = rng.standard_normal((len(n), 50)); x /= np.linalg.norm(x, axis=1, keepdims=True); return x * np.asarray(n)[:, None]


def run(fn, carriers_norm, bulk_norm, mode):
    c = rows_with_norms(carriers_norm); b = rows_with_norms(bulk_norm)
    q, a = fn(b, c, genes, perts, {"mode": mode, "q": [2.5, 2.0, 1.75]})
    return a["counts"]


base = np.linspace(1.0, 3.0, 300)
cases = []
cap = base.copy(); cap[150:] = cap[150]                 # 150 capped targets on the key ratio (carriers norm 1)
cases.append(("tercile cap of 150", np.ones(300), cap, "tercile", [100, 50, 150]))
low = base.copy(); low[:140] = low[0]                   # a 140-row tie at the bottom of the strength
cases.append(("stercile 140 tied weakest", low, np.ones(300), "stercile", [140, 60, 100]))
col = base.copy(); col[:220] = col[0]                   # both cuts collapse onto the tie
cases.append(("stercile both cuts collapsed", col, np.ones(300), "stercile", [220, 0, 80]))
cases.append(("control distinct", base, np.ones(300), "stercile", [100, 100, 100]))
for name, cn, bn, mode, want in cases:
    got_old, got_new = run(old, cn, bn, mode), run(new, cn, bn, mode)
    ok = got_new == want
    bad += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}: old {got_old} -> new {got_new} (want {want})")
print("ALL PASS" if not bad else f"{bad} FAILED"); sys.exit(1 if bad else 0)
