#!/usr/bin/env python3
"""402 builder patch: knob values knn25o400ja / knn25o400jaq80 = 399's per-target artifact extended with context-A rows from Jurkat-only
votes (`psign_by_target_A`, default `psign_A`). Contexts B/C read `psign_by_target`/`psign` exactly as 399. Idempotent; --write to apply."""
import sys, hashlib, difflib, json
from pathlib import Path
P = Path("scripts/make_submission_pm_cd4.py"); src = P.read_text(); orig = src
if "knn25o400jaq80" in src: sys.exit("already applied")
SHA = json.load(open("experiments/402-ctxA-jurkat-knn/work/knn25o400_A_jurkat50_prior_challenge.json"))["artifact_sha256"]
def rep(old, new):
    global src
    assert src.count(old) == 1, (src.count(old), old[:90]); src = src.replace(old, new)
rep('DE_SIGN_PRIOR_VALUES = ("clean4", "clean4q80", "down", "clean4ja", "clean4q80ja", "knn25o400", "knn25o400q80")',
    'DE_SIGN_PRIOR_VALUES = ("clean4", "clean4q80", "down", "clean4ja", "clean4q80ja", "knn25o400", "knn25o400q80", "knn25o400ja", "knn25o400jaq80")')
rep('DE_SIGN_PROTECT_Q = 0.8 if DE_SIGN_PRIOR in ("clean4q80", "clean4q80ja", "knn25o400q80") else None',
    'DE_SIGN_PROTECT_Q = 0.8 if DE_SIGN_PRIOR in ("clean4q80", "clean4q80ja", "knn25o400q80", "knn25o400jaq80") else None')
rep('DE_SIGN_PRIOR_PER_TARGET = DE_SIGN_PRIOR in ("knn25o400", "knn25o400q80")\n',
    'DE_SIGN_PRIOR_PER_TARGET = DE_SIGN_PRIOR in ("knn25o400", "knn25o400q80", "knn25o400ja", "knn25o400jaq80")\n'
    '# `knn25o400ja` / `knn25o400jaq80` (402, 2026-09-16): 399\'s rows for B/C; context A reads rows voted by the Nadig-Jurkat bank alone\n'
    '# (same-line: 401 split-half 0.607 vs 0.581 for the six-bank family), shrink 50 toward the Jurkat marginal (`psign_by_target_A`,\n'
    '# default `psign_A`).\n'
    'DE_SIGN_PRIOR_PER_TARGET_CTX_A = DE_SIGN_PRIOR in ("knn25o400ja", "knn25o400jaq80")\n')
rep('elif DE_SIGN_PRIOR_PER_TARGET:\n    DE_SIGN_PRIOR_PATH = ROOT / "experiments/399-knn-tilt-arm/work/knn25o400_prior_challenge.npz"\n',
    'elif DE_SIGN_PRIOR_PER_TARGET_CTX_A:\n    DE_SIGN_PRIOR_PATH = ROOT / "experiments/402-ctxA-jurkat-knn/work/knn25o400_A_jurkat50_prior_challenge.npz"\n'
    f'    DE_SIGN_PRIOR_SHA256 = "{SHA}"\n'
    'elif DE_SIGN_PRIOR_PER_TARGET:\n    DE_SIGN_PRIOR_PATH = ROOT / "experiments/399-knn-tilt-arm/work/knn25o400_prior_challenge.npz"\n')
rep('        if DE_SIGN_PRIOR_PER_TARGET:\n            DE_SIGN_PRIOR_BY_TARGET.clear()\n            for t, row, ok in zip(map(str, z["targets"]), z["psign_by_target"], z["has_row"]):\n                if bool(ok):\n                    DE_SIGN_PRIOR_BY_TARGET[t] = np.asarray(row, dtype=np.float64)\n        key = "psign_A" if (DE_SIGN_PRIOR_CTX_A is not None and str(context) == "A") else "psign"\n',
    '        if DE_SIGN_PRIOR_PER_TARGET:\n            DE_SIGN_PRIOR_BY_TARGET.clear()\n'
    '            sfx = "_A" if (DE_SIGN_PRIOR_PER_TARGET_CTX_A and str(context) == "A") else ""\n'
    '            for t, row, ok in zip(map(str, z["targets"]), z["psign_by_target" + sfx], z["has_row" + sfx]):\n'
    '                if bool(ok):\n                    DE_SIGN_PRIOR_BY_TARGET[t] = np.asarray(row, dtype=np.float64)\n'
    '        key = "psign_A" if ((DE_SIGN_PRIOR_CTX_A is not None or DE_SIGN_PRIOR_PER_TARGET_CTX_A) and str(context) == "A") else "psign"\n')
rep('                    **({} if not DE_SIGN_PRIOR_PER_TARGET else {"de_sign_prior_targets_with_row": len(DE_SIGN_PRIOR_BY_TARGET)}),\n',
    '                    **({} if not DE_SIGN_PRIOR_PER_TARGET else {"de_sign_prior_targets_with_row": len(DE_SIGN_PRIOR_BY_TARGET)}),\n'
    '                    **({} if not DE_SIGN_PRIOR_PER_TARGET_CTX_A else {"de_sign_prior_context_A": "nadig_jurkat_knn25_shrink50"}),\n')
diff = "".join(difflib.unified_diff(orig.splitlines(True), src.splitlines(True), "a/builder", "b/builder", n=1)); print(diff)
if "--write" in sys.argv:
    P.write_text(src); print("WRITTEN sha256", hashlib.sha256(P.read_bytes()).hexdigest())
else: print("dry run (pass --write)")
