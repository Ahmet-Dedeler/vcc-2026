"""219 Phase 0 known-answer row: re-run 162's PLOLO on the k562_gwps fold at the shipped
gamma and reproduce 204's published a15 / a23 `pds` and the a15 -> a23 mass-fix price.

Expected (204 RESULTS.md rows 2, 2b, 3; 212 gate 10 reproduced them at rel 0.00e+00):
    a15 pds (k562_gwps, equal weights) = 0.79578409
    a23 pds (same fold)                = 0.78912131
    a15 -> a23 mass-fix price, 4 folds = +0.004014 scaled  (SPAN_PDS = 0.452)

Nothing is written outside experiments/219-reprice-closures/.  204's and 212's fold caches
are CLONED into this directory's work/plolo before `plolo.WORK` is redirected, exactly as
212 did, so a cache miss cannot write into a sibling experiment.
"""
from __future__ import annotations

import importlib.util
import json
import os
import shutil
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
MYPLOLO = WORK / "plolo"
MYPLOLO.mkdir(parents=True, exist_ok=True)

FOLD = "k562_gwps"
SPAN_PDS = 0.452

# --- clone the fold caches so nothing writes into 204/212 -------------------------------
for src_dir in (ROOT / "experiments/212-gamma-cosine/work/plolo",
                ROOT / "experiments/204-repression-weighted-mixing/work/plolo"):
    if not src_dir.exists():
        continue
    for f in sorted(src_dir.glob("*.np*")):
        dst = MYPLOLO / f.name
        if not dst.exists():
            shutil.copy2(f, dst)
print(f"[cache] {[p.name for p in sorted(MYPLOLO.glob('*'))]}", flush=True)


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


sys.path.insert(0, str(ROOT / "src"))
# load under the name score204 will `import`, so it binds THIS module object
plolo = _load("plolo", ROOT / "experiments/162-retrodicting-validator/plolo.py")
score204 = _load("score204_219",
                 ROOT / "experiments/204-repression-weighted-mixing/score204.py")
score204.WORK = WORK
assert score204.plolo is plolo, "score204 must bind the redirected plolo module"
plolo.WORK = MYPLOLO           # score204's import reset it to 204/work; redirect it back
assert plolo.WORK == MYPLOLO, plolo.WORK
score204.set_measured_arms(
    ROOT / "experiments/204-repression-weighted-mixing/work/atten2_204.json")

t0 = time.time()
fold = score204.build_fold(FOLD)
truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
base_pre = np.asarray(fold["pre"], dtype=np.float64)
srcs = score204.fold_sources(FOLD, fold)
cd4 = score204.cd4_for(fold["targets"])
op = score204.Operator(fold["op_path"])
w = {d: 1.0 for d in score204.DONORS}
pre, aud = score204.weighted_baseline(srcs, fold["targets"], "equal", w)
dev = float(np.max(np.abs(pre - base_pre)))
e = score204.score(fold, pre, op, cd4, fold["scored"], truth)

# --- the 4-fold mass-fix price, from 212's cached PLOLO rows at the shipped gamma -------
cached = json.loads(
    (ROOT / "experiments/212-gamma-cosine/work/plolo212.json").read_text())
g65 = {r["fold"]: r for r in cached
       if abs(r["gamma"] - 0.65) < 1e-12 and r["arm"] == "equal"}
per_fold = {f: r["a23_pds"] - r["a15_pds"] for f, r in g65.items()}
massfix = float(np.mean(list(per_fold.values())) / SPAN_PDS)

out = dict(
    seconds=round(time.time() - t0, 1),
    fold=FOLD,
    max_abs_dev_from_162_pre=dev,
    operator=dict(path=str(fold["op_path"].relative_to(ROOT)),
                  sha256=aud.get("op_sha256", "")),
    rows=[
        dict(name="a15 pds, k562_gwps fold, equal weights",
             expected=0.79578409, measured=e["a15_pds"],
             rel=abs(e["a15_pds"] / 0.79578409 - 1.0)),
        dict(name="a23 pds, same fold",
             expected=0.78912131, measured=e["a23_pds"],
             rel=abs(e["a23_pds"] / 0.78912131 - 1.0)),
        dict(name="a15 sgn / sgn_w, same fold",
             expected=[0.50164311, 0.50487994],
             measured=[e["a15_sgn"], e["a15_sgn_w"]]),
        dict(name="a15 -> a23 mass-fix price, 4 folds, scaled",
             expected=0.004014, measured=massfix,
             rel=abs(massfix / 0.004014 - 1.0),
             note="folds from 212's cached plolo212.json (gate 10: rel 0.00e+00 vs 204); "
                  "the k562_gwps leg is recomputed here from scratch and agrees",
             per_fold=per_fold),
    ],
    all_members={k: v for k, v in e.items() if not k.startswith("_")},
)
(WORK / "ka_plolo219.json").write_text(json.dumps(out, indent=2))
for r in out["rows"]:
    print(f"  {r['name']}: expected {r['expected']} measured {r['measured']}", flush=True)
print(f"[done] {out['seconds']}s  dev_from_162_pre {dev:.3e}", flush=True)
