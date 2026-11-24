"""180 — the DE-side price of the knee, on h1-test, at the three per-context amplitudes.

Rule 12 says price all six or say which are unmeasured and why. `mse` is priced on the
challenge panel through 096's claim route (`knee180.py`) and `pds` is settled by an identity
there; `nmae`, `reach`, `fid` and `jac` are not computable on the challenge panel at all
(no reference DE table), so they are read on h1-test -- the surface 168d's standing rule
names as the sizing instrument.

**h1-test is ONE context.** A per-context knee does not exist on it, and its own `rho` is
1.0 (it sits ABOVE the knee), so it cannot price `mse` (170 section 4). What transfers is the
DIRECTION and rough size of each DE member's response to a pure amplitude change, which is
exactly what the knee is.

The knee is a pure per-context scale on the effect matrix, so each arm here is
`170/work/eff/h1_g0.65.npz` multiplied by `alpha / 0.982` and pushed through 170's own
`deliver` / `emit_and_score` / `Surface` stack, unchanged.

KNOWN-ANSWER ROW: `s = 1.0` must reproduce 170's `h1sweep.json` gamma-0.65 row exactly. That
row was produced by a different script on a different day; if this one does not land on it,
the reuse is broken and nothing below it means anything. 170's own `mse` gate (its
`known_answer_gate`, anchored on a real cell-eval2 `agg_results.csv`) runs first.

    uv run python h1alpha180.py
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
EXP170 = REPO / "experiments/170-gamma-mse"
WORK = HERE / "work"

os.environ.setdefault("VCC_ALPHA", "0.982")
os.environ.setdefault("VCC_CD4_TARGET_F", "0.10")
os.environ.setdefault("VCC_MASS_FIX", "minnorm")

sys.path.insert(0, str(EXP170))
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "experiments/096-mse-floor"))
sys.path.insert(0, str(REPO / "experiments/092-call-set-truncation"))
sys.path.insert(0, str(REPO / "experiments/153-de-bench"))

ANCHOR_170 = {
    "mse": 1.0159689545093225, "pds": 0.8753501400560225, "nmae": 0.959881122307442,
    "reach": 0.17751503650931924, "fid": 0.5081565598953068, "jac": 0.06424535489807914,
}
# The scaled-score spans (AGENTS.md section 4). `nmae` is lower-is-better.
SPAN = {"pds": 0.455, "mse": 0.9502, "fid": 0.305, "reach": 0.895, "jac": 0.365,
        "nmae": 0.60}
ALPHA_REF = 0.982


def main() -> int:
    os.chdir(REPO)
    import h1sweep
    from h1kit import NCELL, builder, deliver, emit_and_score, pds_exact, to_h1_axis, Surface170, DEKit
    from debench import Surface
    from decompose import terms
    from msekit import score

    h1sweep.known_answer_gate()
    B = builder()
    S = Surface170()
    pre = np.load(EXP170 / "work/pre.npz", allow_pickle=True)
    de = DEKit()
    real = de.real_table(S.targets)
    de.release_reference()

    knee = json.loads((WORK / "knee180.json").read_text())
    scales = {"anchor(a23)": 1.0}
    for context in ("A", "B", "C"):
        scales[f"knee_{context}"] = knee["per_context"][context]["knee_scale"]
    scales["knee_mean"] = float(np.mean([
        knee["per_context"][c]["knee_scale"] for c in ("A", "B", "C")
    ]))
    print(f"[h1] arms: " + "  ".join(f"{k}={v:.6f}" for k, v in scales.items()), flush=True)

    eff_prod = np.load(EXP170 / "work/eff/h1_g0.65.npz")["effect"]
    eff0, _ = to_h1_axis(eff_prod, pre["genes"], pre["h1_genes"])
    pos = {str(t): i for i, t in enumerate(pre["h1_targets"])}
    eff0 = eff0[np.asarray([pos[str(t)] for t in S.targets])]

    rows = []
    for seed in (0, 1):
        for name, scale in scales.items():
            if seed == 1 and name not in ("anchor(a23)", "knee_mean"):
                continue                      # the twin only has to size the noise floor
            t0 = time.time()
            profiles, dh = deliver(B, eff0 * scale, S.b_ctrl)
            sums, jk, (pp_, pq_, pl_) = emit_and_score(
                B, profiles, eff0 * scale, S.b_ctrl, S.libs, de, seed=seed)
            counts = np.full(S.targets.size, float(NCELL))
            mse, _, _ = score(S.targets, sums, counts, jk, S.rp, S.rs, S.rc, S.rj, S.genes)
            t = terms(S.targets, sums, counts, jk, S.rp, S.rs, S.rc, S.rj, S.genes)
            pds = pds_exact(S.targets, sums, S)
            path = WORK / f"de_h1_{name.replace('(', '_').replace(')', '')}_s{seed}.npz"
            np.savez_compressed(path, genes=de.genes, perts=S.targets, **real,
                                pred_p=pp_, pred_padj=pq_, pred_lfc=pl_)
            surf = Surface(str(path))
            n_real, n_pred = surf.counts()
            row = dict(
                arm=name, seed=seed, scale=scale, alpha=ALPHA_REF * scale,
                mse=float(mse), pds=float(pds), nmae=float(surf.nmae()),
                reach=float(surf.reach()), fid=float(surf.fid()), jac=float(surf.jac()),
                A=float(t["A"]), X=float(t["X"]), B=float(t["budget"]),
                claim=float(t["claim"]), rho=float(t["rho"]),
                n_pred=float(np.mean(n_pred)), n_real=float(np.mean(n_real)),
                sign_acc=float(surf.sign_accuracy()),
                row_norm=float(np.median(np.linalg.norm(dh, axis=1))),
                secs=round(time.time() - t0, 1),
            )
            rows.append(row)
            print(f"[h1] {name:<12} seed {seed} s {scale:.6f}  mse {row['mse']:.5f}  "
                  f"pds {row['pds']:.5f}  nmae {row['nmae']:.5f}  reach {row['reach']:.5f}  "
                  f"fid {row['fid']:.5f}  jac {row['jac']:.5f}  rho {row['rho']:.4f}  "
                  f"n_pred {row['n_pred']:.0f}  [{row['secs']:.0f}s]", flush=True)
            (WORK / "h1alpha180.json").write_text(json.dumps(rows, indent=1))

    base = [r for r in rows if r["arm"] == "anchor(a23)" and r["seed"] == 0][0]
    print("\n[gate] anchor vs 170's own gamma-0.65 row:")
    bad = []
    for key, want in ANCHOR_170.items():
        got = base[key]
        print(f"[gate]   {key:<6} {got!r} vs {want!r}  d {got - want:+.3e}")
        if abs(got - want) > 1e-12:
            bad.append(key)
    if bad:
        raise SystemExit(f"KNOWN-ANSWER FAIL on h1-test: {bad}")
    print("[gate] KNOWN-ANSWER PASS: this reproduces 170's h1-test row exactly\n")

    twin = {r["arm"]: r for r in rows if r["seed"] == 1}
    print("[h1] re-seeded twin (the emission noise floor on this surface):")
    for name, row in twin.items():
        other = [r for r in rows if r["arm"] == name and r["seed"] == 0][0]
        print("[h1]   " + name + "  " + "  ".join(
            f"d{k} {row[k] - other[k]:+.5f}" for k in
            ("pds", "nmae", "reach", "fid", "jac")))

    print("\n[h1] deltas against the anchor, raw and as scaled-score/6 (overall units):")
    for row in rows:
        if row["seed"] or row["arm"] == "anchor(a23)":
            continue
        parts = []
        for key in ("pds", "nmae", "reach", "fid", "jac"):
            raw = row[key] - base[key]
            sign = -1.0 if key == "nmae" else 1.0
            parts.append(f"{key} {raw:+.5f} ({sign * raw / SPAN[key] / 6:+.6f})")
        print(f"[h1]   {row['arm']:<12} " + "  ".join(parts), flush=True)
    print(f"\n[h1] peak RSS {B.peak_rss_gib():.2f} GiB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
