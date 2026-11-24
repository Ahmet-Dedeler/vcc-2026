"""200 — Rule 11: a path and a sha256 for every operator and frozen input this experiment loads.

AGENTS.md 11.1 rule 11: "Record a path and a sha, never a description -- and never infer an
operator from a filename."  Eight whitening operators exist on disk with near-identical names
and an evening was once lost to loading the wrong one, so this script is run FIRST and its
output is pasted into RESULTS.md verbatim.

It also does the rule-10 half that a config diff cannot: read the SHIPPED artifact's meta
(`submissions/a15-cd4f010.vcc.meta.json`, `submissions/a23-massfix.h5ad.meta.json`) and print
the config fields this experiment is going to ablate, so every "incumbent value" quoted later
comes from the artifact rather than from a memory of the command line.

Reads only.  No cell-eval2, no build, no submission, nothing written outside this directory.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")

# Every operator this experiment can load, plus the ones it must NOT confuse them with.
OPERATORS = [
    ("production union_pruned (a15/a23 ship this)", "data/donors/whitener_union_pruned.npz"),
    ("production union (a2-prune's `before`)", "data/donors/whitener_union.npz"),
    ("PLOLO held k562_gwps", "experiments/062-cd4-post-whiten/cloud-output/whiteners/held_k562_gwps.npz"),
    ("PLOLO held orion_hct116", "experiments/062-cd4-post-whiten/cloud-output/whiteners/held_orion_hct116.npz"),
    ("PLOLO held orion_hek293t", "experiments/062-cd4-post-whiten/cloud-output/whiteners/held_orion_hek293t.npz"),
    ("h1-test held-safe union (095)", "experiments/095-clean-harness/out-h1test/whitener_union.npz"),
    ("NOT USED: union_full", "data/donors/whitener_union_full.npz"),
    ("NOT USED: whitener_full", "data/donors/whitener_full.npz"),
    ("NOT USED: 095 out/ (orion-hct116 held)", "experiments/095-clean-harness/out/whitener_union.npz"),
    ("NOT USED: 095 out-essential/", "experiments/095-clean-harness/out-essential/whitener_union.npz"),
    ("NOT USED: 095 out-random/", "experiments/095-clean-harness/out-random/whitener_union.npz"),
    ("NOT USED: 096 out-alpha/", "experiments/096-mse-floor/out-alpha/whitener_union.npz"),
    ("NOT USED: 040", "experiments/040-union-gamma/out/whitener_union.npz"),
]

FROZEN = [
    ("challenge panel", "data/challenge/pert_counts.csv"),
    ("challenge gene axis", "data/challenge/gene_names.csv"),
    ("challenge context A", "data/challenge/context_A.h5ad"),
    ("challenge context B", "data/challenge/context_B.h5ad"),
    ("challenge context C", "data/challenge/context_C.h5ad"),
    ("donor k562_gwps delta", "data/donors/k562_gwps/delta.npy"),
    ("donor orion_hct116 delta", "data/donors/orion_hct116/delta.npy"),
    ("donor orion_hek293t delta", "data/donors/orion_hek293t/delta.npy"),
    ("donor kolf2 delta (PLOLO truth only)", "data/external/kolf2/delta.npy"),
    ("PerturbMatch", "data/external/perturbmatch/delta.npy"),
    ("h1-test harness", "data/harness/h1-test.h5ad"),
    ("h1-test prediction (095 build, 157/164/169 all read it)",
     "experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
    ("096 v11 challenge stats", "experiments/096-mse-floor/work/v11-stats.npz"),
    ("096 board solve", "experiments/096-mse-floor/work/board_solve.npz"),
    ("096 control stats", "experiments/096-mse-floor/work/ctrl-stats.npz"),
    ("096 C_p proxy", "experiments/096-mse-floor/work/cp-proxy-ntc.npz"),
    ("162 PLOLO arms (known-answer source)", "experiments/162-retrodicting-validator/work/arms.json"),
    ("169 h1-test recovered effect", "experiments/169-sharpening-six-member/work/h1test-effect.npz"),
    ("the shipped builder", "scripts/make_submission_pm_cd4.py"),
    ("PLOLO instrument", "experiments/162-retrodicting-validator/plolo.py"),
    ("062 effect pipeline", "experiments/062-cd4-post-whiten/run.py"),
    ("064 PM/operator pipeline", "experiments/064-pm-cd4-composition/run.py"),
]

METAS = [
    ("a15 as shipped", "submissions/a15-cd4f010.vcc.meta.json"),
    ("a23 mass fix", "submissions/a23-massfix.h5ad.meta.json"),
]

# The a15 config fields this experiment ablates, in the order RESULTS.md reports them.
ABLATED = (
    "whiten_axis", "whiten_gamma", "whiten_lambda", "whitener_id", "perturbmatch_blend",
    "cd4_target_f", "cd4_fraction", "cd4_mode", "alpha", "donor_weights", "min_cpm",
    "keep_top", "keep_spread", "tail_spread", "spread", "spread_map", "keep_key",
    "overdispersion", "dispersion_argument", "effect_scale", "shrink_key", "emission",
)


def sha256(path: str, block: int = 8 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(block)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def npz_eig(path: str):
    """Max eigenvalue and basis shape -- the fingerprint that separates the eight operators."""
    import numpy as np
    with np.load(path, mmap_mode="r") as z:
        ev = np.asarray(z["evals"])
        shape = z["evecs"].shape if "evecs" in z.files else None
    return float(ev.max()), int(ev.size), shape


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    rec: dict[str, object] = {}

    print("=" * 104)
    print("OPERATORS  (rule 11: eight exist with near-identical names; maxeig is the fingerprint)")
    print("=" * 104)
    print(f"{'what':<44} {'sha256[:12]':<14} {'maxeig':>9} {'n':>7}  path")
    ops = {}
    for label, rel in OPERATORS:
        p = os.path.join(ROOT, rel)
        if not os.path.exists(p):
            print(f"{label:<44} {'MISSING':<14} {'':>9} {'':>7}  {rel}")
            continue
        s = sha256(p)
        try:
            mx, n, _ = npz_eig(p)
        except Exception as exc:  # noqa: BLE001
            mx, n = float("nan"), -1
            print(f"  (eig read failed: {exc})")
        ops[rel] = dict(sha256=s, maxeig=mx, n=n)
        print(f"{label:<44} {s[:12]:<14} {mx:9.4f} {n:7d}  {rel}")
    rec["operators"] = ops

    print()
    print("=" * 104)
    print("FROZEN INPUTS")
    print("=" * 104)
    frozen = {}
    for label, rel in FROZEN:
        p = os.path.join(ROOT, rel)
        if not os.path.exists(p):
            print(f"{label:<52} {'MISSING':<14}  {rel}")
            continue
        s = sha256(p)
        frozen[rel] = dict(sha256=s, bytes=os.path.getsize(p))
        print(f"{label:<52} {s[:12]:<14} {os.path.getsize(p)/1e6:10.1f} MB  {rel}")
    rec["frozen"] = frozen

    print()
    print("=" * 104)
    print("THE SHIPPED ARTIFACT'S OWN META  (rule 10: diff the artifact, not the command)")
    print("=" * 104)
    metas = {}
    for label, rel in METAS:
        p = os.path.join(ROOT, rel)
        m = json.loads(open(p).read())
        metas[label] = dict(path=rel, sha256=sha256(p),
                            config={k: m["config"].get(k) for k in ABLATED},
                            operator_sha256=m["operator"]["sha256"],
                            operator_columns=m["operator"]["columns"],
                            whitening=m["effects"].get("whitening"),
                            cd4_fraction_solved=m["effects"]["cd4"].get("fraction"),
                            final_sha256=m["effects"]["final_sha256"],
                            post_whiten_sha256=m["effects"]["post_whiten_sha256"],
                            candidate_pre_sha256=m["effects"]["candidate_pre_sha256"],
                            baseline_pre_sha256=m["effects"]["baseline_pre_sha256"])
        print(f"\n--- {label}  ({rel})")
        print(f"  operator sha256 {m['operator']['sha256'][:16]}  columns {m['operator']['columns']}")
        print(f"  effects.whitening: {json.dumps(m['effects'].get('whitening'))}")
        print(f"  effects.cd4.fraction (solved): {m['effects']['cd4'].get('fraction')}")
        print(f"  baseline_pre {m['effects']['baseline_pre_sha256'][:12]}  "
              f"candidate_pre {m['effects']['candidate_pre_sha256'][:12]}  "
              f"post_whiten {m['effects']['post_whiten_sha256'][:12]}  "
              f"final {m['effects']['final_sha256'][:12]}")
        for k in ABLATED:
            print(f"    {k:<22} {json.dumps(m['config'].get(k))}")
    rec["metas"] = metas

    # a15 -> a23 config diff: rule 10 says EXACTLY the intended number of differences.
    a15 = metas["a15 as shipped"]["config"]
    a23 = metas["a23 mass fix"]["config"]
    diff = {k: (a15.get(k), a23.get(k)) for k in ABLATED if a15.get(k) != a23.get(k)}
    m15 = json.loads(open(os.path.join(ROOT, METAS[0][1])).read())["config"]
    m23 = json.loads(open(os.path.join(ROOT, METAS[1][1])).read())["config"]
    full = sorted(set(m15) | set(m23))
    full_diff = {k: (m15.get(k), m23.get(k)) for k in full if m15.get(k) != m23.get(k)}
    print(f"\n--- a15 -> a23 FULL config diff ({len(full_diff)} field(s)): "
          f"{json.dumps(full_diff)}")
    rec["a15_a23_config_diff"] = full_diff
    rec["a15_a23_ablated_diff"] = diff

    with open(os.path.join(WORK, "inputs200.json"), "w") as fh:
        json.dump(rec, fh, indent=1)
    print(f"\nwrote {os.path.join(WORK, 'inputs200.json')}")


if __name__ == "__main__":
    sys.exit(main())
