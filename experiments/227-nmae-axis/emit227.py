"""227 — the emission falsifier, measured.  Does `spread`/`overdispersion` move `nmae` at all?

Candidate 1 of the mandate: `emit_counts(mode="quasi", spread=...)` changes cell-to-cell variance
at fixed gene totals (034), and Jensen makes the ARITHMETIC mean of per-cell CPM a function of that
variance -- so in principle it could move `nmae` with the log-pseudobulk `mse` pinned.

The structural answer is that it cannot, and production's own source says so
(`src/vcc2026/rc_emit.py:38`, `scripts/make_submission_pm_cd4.py:3328`): the quasi law's per-gene
totals are `floor(N*L*profile)` -- independent of `spread` -- and every emitted cell carries the
same library total, so `mean_i (x_ig / T_i) = (Sum_i x_ig) / (N*L)` EXACTLY.  The arithmetic mean
of per-cell CPM is the design pseudobulk, whatever the spread does to the cells.

Structural arguments are how this project has been wrong before, so this measures it: the same
design emitted at four spreads and two dispersions, two seeds each, `nmae` read off each.  The
deterministic design profile is the zero-emission reference row.
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, os.path.join(REPO, "src"))
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
sys.path.insert(0, HERE)
from vcc2026.donor_transfer import emit_counts        # noqa: E402
from debench import Surface                            # noqa: E402
from phase1 import gate_of                             # noqa: E402

E169 = os.path.join(REPO, "experiments/169-sharpening-six-member")
TABLE = os.path.join(E169, "de169.q1.00_mf.npz")
EFFECT = os.path.join(E169, "work/h1test-effect.npz")
N_EMIT, EMIT_L = 400, 100_000
MIN_CPM, TAIL, KEEP_TOP, KEEP_SPREAD, PHI = 5.0, 1.0, 8000, 0.25, 4.0
N_TARGETS = int(os.environ.get("N_TARGETS", "25"))


def load_builder():
    spec = importlib.util.spec_from_file_location(
        "prod_builder_227", f"{REPO}/scripts/make_submission_pm_cd4.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_227"] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    t0 = time.time()
    os.makedirs(WORK, exist_ok=True)
    prod = load_builder()
    E = np.load(EFFECT, allow_pickle=True)
    H0 = E["H"].astype(np.float64)
    b = E["b"].astype(np.float64)
    a_gate = np.asarray(E["gate"], bool)
    ecpm = E["cpm"].astype(np.float64)
    egenes = np.asarray([str(g) for g in E["genes"]])
    etargets = np.asarray([str(t) for t in E["targets"]])
    esums = E["sums"].astype(np.float64)

    S = Surface(TABLE)
    G = gate_of(S)
    keep = np.array([g in set(S.genes) for g in egenes])
    kpos = {g: i for i, g in enumerate(S.genes)}
    kidx = np.array([kpos[g] for g in egenes[keep]])
    assert np.array_equal(egenes[keep][np.argsort(kidx)], S.genes)
    order = np.argsort(kidx)
    ref_lin = ecpm[keep][order]                      # control mean CPM on the table's gene axis
    T = S.real_lfc.astype(np.float64)
    epos = {t: i for i, t in enumerate(etargets)}
    tpos = {t: i for i, t in enumerate(S.perts)}

    rows = [(p, s) for p, s in G if S.perts[p] in epos][:N_TARGETS]
    print(f"{len(rows)} rows, {int(keep.sum())} scored genes  [{time.time()-t0:.0f}s]",
          flush=True)

    # ---- the deterministic design row (zero emission) -----------------------------------
    design = {}
    spreadvec = {}
    for p, s in rows:
        t = S.perts[p]
        k = epos[t]
        h = prod.mass_consistent_effect(H0[k], b, a_gate, mode="minnorm")
        lin = np.clip(np.expm1(b + h), 0.0, None)
        profile = lin / lin.sum() * (float(esums[k].sum()) / (N_EMIT * EMIT_L))
        totals = np.floor(N_EMIT * EMIT_L * profile).astype(np.int64)
        pcpm = 1e6 * totals / max(totals.sum(), 1)
        design[t] = (profile, pcpm[keep][order])
        # production's own spread vector (188/169 conventions)
        eff = np.abs(np.log2((pcpm + 1e-9) / (ecpm + 1e-9)))
        gated = (ecpm > MIN_CPM) & (b > 0)
        o = np.lexsort((np.arange(eff.size), -np.where(gated, eff, -np.inf)))
        sv = np.full(eff.size, float(TAIL))
        sv[o[:min(KEEP_TOP, int(gated.sum()))]] = float(KEEP_SPREAD)
        spreadvec[t] = sv

    def nmae_from(pred):
        v = []
        for p, s in rows:
            h = pred[S.perts[p]][s]
            v.append(np.abs(h - T[p, s]).sum() / np.abs(T[p, s]).sum())
        return float(np.mean(v))

    det = {t: np.log2((d[1] + 1e-9) / (ref_lin + 1e-9)) for t, d in design.items()}
    base = nmae_from(det)
    print(f"deterministic design profile:  nmae {base:.9f}  [{time.time()-t0:.0f}s]", flush=True)

    ARMS = [("production spread vector", None, PHI),
            ("uniform spread 0.25", 0.25, PHI),
            ("uniform spread 1.0", 1.0, PHI),
            ("uniform spread 4.0", 4.0, PHI),
            ("uniform spread 1.0, phi 1.0", 1.0, 1.0),
            ("uniform spread 1.0, phi 16.0", 1.0, 16.0)]
    res = []
    for name, sp, phi in ARMS:
        for seed in (90_000, 424_242):
            pred = {}
            cellsum = []
            for p, s in rows:
                t = S.perts[p]
                profile, pcpm = design[t]
                sv = spreadvec[t] if sp is None else np.full(profile.size, float(sp))
                rng = np.random.default_rng(seed + p)
                vals, cols, ptr = emit_counts(
                    profile, np.full(N_EMIT, EMIT_L), 0.1, rng,
                    mode="quasi", spread=sv, overdispersion=float(phi))
                B = np.zeros((N_EMIT, profile.size))
                for r in range(N_EMIT):
                    sl = slice(ptr[r], ptr[r + 1])
                    B[r, cols[sl]] = vals[sl]
                ss = B.sum(1)
                cellsum.append((float(ss.min()), float(ss.max())))
                ss[ss == 0] = 1.0
                m = (B[:, keep][:, order] * (1e6 / ss[:, None])).mean(0)
                pred[t] = np.log2((m + 1e-9) / (ref_lin + 1e-9))
                del B
            v = nmae_from(pred)
            lo = min(c[0] for c in cellsum)
            hi = max(c[1] for c in cellsum)
            res.append(dict(arm=name, spread=sp, phi=phi, seed=seed, nmae=v,
                            d_vs_design=v - base, cell_total_min=lo, cell_total_max=hi))
            print(f"  {name:30s} phi {phi:5.1f} seed {seed:7d}  nmae {v:.9f}  "
                  f"d {v-base:+.3e}  cell totals [{lo:.0f}, {hi:.0f}]  "
                  f"[{time.time()-t0:.0f}s]", flush=True)

    d = [abs(r["d_vs_design"]) for r in res]
    out = dict(n_rows=len(rows), nmae_design=base, arms=res,
               max_abs_deviation=float(max(d)),
               spread_range=float(max(r["nmae"] for r in res)
                                  - min(r["nmae"] for r in res)))
    json.dump(out, open(os.path.join(WORK, "emit227.json"), "w"), indent=1)
    print(f"\nmax |nmae(emitted) - nmae(design)| over "
          f"{len(res)} emissions = {out['max_abs_deviation']:.3e}")
    print(f"full range of nmae across all spreads/dispersions/seeds = "
          f"{out['spread_range']:.3e}")


if __name__ == "__main__":
    main()
