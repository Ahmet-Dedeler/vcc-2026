#!/usr/bin/env python3
"""Builder patch (2026-09-16): generalise VCC_DE_SIGN_PRIOR to per-context and per-target prior artifacts.
New knob values: clean4ja / clean4q80ja (396: CLEAN4 for B/C, Nadig-Jurkat tilt for A) and knn25o400 / knn25o400q80 (399: per-target
neighbour-conditioned tilt, CLEAN4 fallback). clean4 / clean4q80 / down / unset paths are unchanged. Idempotent; prints the diff;
writes only with --write. The knn artifact sha is read from its JSON sidecar."""
import sys, hashlib, difflib, json
from pathlib import Path
P = Path("scripts/make_submission_pm_cd4.py"); src = P.read_text(); orig = src
if "knn25o400q80" in src: sys.exit("already applied")
KNN_META = json.load(open("experiments/399-knn-tilt-arm/work/knn25o400_prior_challenge.json")); KNN_SHA = KNN_META["artifact_sha256"]
JA_SHA = json.load(open("experiments/396-jurkat-tilt-A/work/prior_ctx_challenge.json"))["artifact_sha256"]
def rep(old, new, count=1):
    global src
    assert src.count(old) == count, (src.count(old), old[:80])
    src = src.replace(old, new)
rep('if DE_SIGN_PRIOR is not None and DE_SIGN_PRIOR not in ("clean4", "clean4q80", "down"):\n    raise SystemExit(f"VCC_DE_SIGN_PRIOR must be \'clean4\', \'clean4q80\' or \'down\' when set: {DE_SIGN_PRIOR!r}")',
    'DE_SIGN_PRIOR_VALUES = ("clean4", "clean4q80", "down", "clean4ja", "clean4q80ja", "knn25o400", "knn25o400q80")\n'
    'if DE_SIGN_PRIOR is not None and DE_SIGN_PRIOR not in DE_SIGN_PRIOR_VALUES:\n    raise SystemExit(f"VCC_DE_SIGN_PRIOR must be one of {DE_SIGN_PRIOR_VALUES} when set: {DE_SIGN_PRIOR!r}")')
rep('DE_SIGN_PROTECT_Q = 0.8 if DE_SIGN_PRIOR == "clean4q80" else None',
    'DE_SIGN_PROTECT_Q = 0.8 if DE_SIGN_PRIOR in ("clean4q80", "clean4q80ja", "knn25o400q80") else None')
rep('DE_SIGN_PRIOR_PATH = ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz"\nDE_SIGN_PRIOR_SHA256 = "7e532f84a2425176a4582a2e42e458bc3fd25bc3267f890b27a530565df5936e"',
    '# `clean4ja` / `clean4q80ja` (396, 2026-09-16): CLEAN4 for B/C and, for context A only, the Nadig-Jurkat per-gene tilt (P_down over\n'
    '# its 2,270 perturbations; 370: 0.577 held-out on A\'s block vs cross-line ~0.55; covers 99.8 % of A\'s block). The artifact carries\n'
    '# `psign` (bit-identical to the pinned CLEAN4 file) and `psign_A`.\n'
    '# `knn25o400` / `knn25o400q80` (399, 2026-09-16): a PER-TARGET tilt -- for each panel target in K562-GWPS, P_down over the rows of its 25\n'
    '# nearest GWPS perturbations (Pearson on GWPS response rows) in the six-bank family (CLEAN4 + raw Orion HCT116/HEK293T), shrunk\n'
    '# toward the family marginal with prior weight 400; targets absent from GWPS use the family marginal sign, genes without a vote\n'
    '# the CLEAN4 sign. 397 (H1 analytic, shipped construction): block precision 0.534 vs CLEAN4 0.522, fid +0.0207 vs E0; 398 (leave-one-bank-out): +0.035 over the marginal tilt on Jurkat, +0.065 on K562-essential,\n'
    '# +0.009 on H1, every CI above zero; 398c/398d: the Orion rows add an epithelial marginal (+0.030 on RPE1) and m=400 keeps every\n'
    '# held-out line above its marginal (Jurkat +0.034, HepG2 +0.022, RPE1 +0.002). The artifact carries `psign` (the family marginal default), `psign_clean4`, `targets`, `psign_by_target`, `has_row`.\n'
    'DE_SIGN_PRIOR_CTX_A = "nadig_jurkat" if DE_SIGN_PRIOR in ("clean4ja", "clean4q80ja") else None\n'
    'DE_SIGN_PRIOR_PER_TARGET = DE_SIGN_PRIOR in ("knn25o400", "knn25o400q80")\n'
    'if DE_SIGN_PRIOR_CTX_A is not None:\n'
    '    DE_SIGN_PRIOR_PATH = ROOT / "experiments/396-jurkat-tilt-A/work/prior_ctx_challenge.npz"\n'
    f'    DE_SIGN_PRIOR_SHA256 = "{JA_SHA}"\n'
    'elif DE_SIGN_PRIOR_PER_TARGET:\n'
    '    DE_SIGN_PRIOR_PATH = ROOT / "experiments/399-knn-tilt-arm/work/knn25o400_prior_challenge.npz"\n'
    f'    DE_SIGN_PRIOR_SHA256 = "{KNN_SHA}"\n'
    'else:\n'
    '    DE_SIGN_PRIOR_PATH = ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz"\n'
    '    DE_SIGN_PRIOR_SHA256 = "7e532f84a2425176a4582a2e42e458bc3fd25bc3267f890b27a530565df5936e"')
rep('def de_sign_prior_full(genes) -> np.ndarray:\n    """The block sign per official gene (-1 down / +1 up) for `VCC_DE_SIGN_PRIOR`; verifies both pins."""',
    'DE_SIGN_PRIOR_BY_TARGET: dict[str, np.ndarray] = {}   # filled by de_sign_prior_full for the per-target knobs\n\n\n'
    'def de_sign_prior_full(genes, context=None) -> np.ndarray:\n'
    '    """The DEFAULT block sign per official gene (-1 down / +1 up) for `VCC_DE_SIGN_PRIOR`; verifies both pins.\n'
    '    Per-context knobs (396): context A reads `psign_A`, every other context `psign`. Per-target knobs (399): also fills\n'
    '    DE_SIGN_PRIOR_BY_TARGET[target] from `psign_by_target` for the targets that have a row; the default stays `psign`."""')
rep('        if list(map(str, z["genes"])) != list(map(str, genes)):\n            raise SystemExit("clean4 prior gene axis differs from the official gene order")\n        return np.asarray(z["psign"], dtype=np.float64)',
    '        if list(map(str, z["genes"])) != list(map(str, genes)):\n            raise SystemExit("clean4 prior gene axis differs from the official gene order")\n'
    '        if DE_SIGN_PRIOR_PER_TARGET:\n'
    '            DE_SIGN_PRIOR_BY_TARGET.clear()\n'
    '            for t, row, ok in zip(map(str, z["targets"]), z["psign_by_target"], z["has_row"]):\n'
    '                if bool(ok):\n'
    '                    DE_SIGN_PRIOR_BY_TARGET[t] = np.asarray(row, dtype=np.float64)\n'
    '        key = "psign_A" if (DE_SIGN_PRIOR_CTX_A is not None and str(context) == "A") else "psign"\n'
    '        return np.asarray(z[key], dtype=np.float64)')
rep('                de_sign_full = de_sign_prior_full(spec.genes)\n', '                de_sign_full = de_sign_prior_full(spec.genes, context)\n')
rep('                                "prior": np.ascontiguousarray(de_sign_full[rc_control.keep]),\n',
    '                                "prior": np.ascontiguousarray(\n'
    '                                    DE_SIGN_PRIOR_BY_TARGET.get(str(target_value), de_sign_full)[rc_control.keep]\n'
    '                                ),\n')
rep('                    "de_sign_prior_sha256": DE_SIGN_PRIOR_SHA256 if DE_SIGN_PRIOR != "down" else None,\n',
    '                    "de_sign_prior_sha256": DE_SIGN_PRIOR_SHA256 if DE_SIGN_PRIOR != "down" else None,\n'
    '                    **({} if DE_SIGN_PRIOR_CTX_A is None else {"de_sign_prior_context_A": DE_SIGN_PRIOR_CTX_A}),\n'
    '                    **({} if not DE_SIGN_PRIOR_PER_TARGET else {"de_sign_prior_targets_with_row": len(DE_SIGN_PRIOR_BY_TARGET)}),\n')
diff = "".join(difflib.unified_diff(orig.splitlines(True), src.splitlines(True), "a/builder", "b/builder", n=1)); print(diff)
if "--write" in sys.argv:
    P.write_text(src); print("WRITTEN sha256", hashlib.sha256(P.read_bytes()).hexdigest())
else: print("dry run (pass --write)")
