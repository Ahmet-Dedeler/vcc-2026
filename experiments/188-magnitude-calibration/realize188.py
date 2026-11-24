"""188 — h1-test emissions.  Two questions, one run.

Q1 (the one that matters): **is 169's h1-test `nmae` cost an artefact of the GATE WIDTH?**
169 priced `nmae` by sharpening the effect recovered from a v11-era h1-test artifact whose own
gate is `--min-cpm 5.0`, i.e. **10,783 of 18,080 genes**.  Production a15 ships `--min-cpm 10.0`
and its challenge gate is **9,135 / 8,668 / 9,127 of 18,533** (49%, against h1-test's 60%).  The
extra band is 5-10 CPM, exactly where `|h|` is smallest, so a power transform normalised back to
the row's L2 norm crushes it hardest and lifts the head hardest.  If the h1-test cost shrinks
when the transform is restricted to the CPM>10 sub-gate, then h1-test and HCT116 were never
measuring the same arm and 169's 56x surface disagreement has a mechanism.

The arm is defined so that ONLY the transform's support changes: the 5-10 CPM band is passed
through UNCHANGED (not zeroed, which would be a second variable), and the row L2 norm is
restored over the sub-gate alone.

Q2: price the fitted magnitude maps (`normfit188.py`) on the surface validated for `reach`.

CONVENTIONS ARE 169's, VERBATIM, so its five-seed baselines (`de169.q1.00_mf_s*.npz`) are the
matched controls and nothing has to be re-emitted for them.

KNOWN-ANSWER ROW: `KA_q1.15_mf` re-runs 169's own q = 1.15 arm at its own seed 90,000 and must
reproduce `de169.q1.15_mf.npz` bit-for-bit on every field.
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import anndata as ad
import numpy as np
from scipy.stats import false_discovery_control, norm

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
E169 = os.path.join(ROOT, "experiments/169-sharpening-six-member")
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/092-call-set-truncation")
sys.path.insert(0, E169)
sys.path.insert(0, HERE)
from vcc2026.donor_transfer import emit_counts       # noqa: E402
from de_kernel import ctrl_tiesums, wilcoxon_z       # noqa: E402
from sharp import sharpen                            # noqa: E402
from mapfit import Map                               # noqa: E402

REAL = f"{ROOT}/data/harness/h1-test.h5ad"
OUT = os.path.join(HERE, "de188")
CTRL, COL, MIN_CPM = "non-targeting", "target_gene", 5.0
N_EMIT, EMIT_L = 400, 100_000
TS = 50_000.0
TAIL, KEEP_TOP, KEEP_SPREAD, PHI = 1.0, 8000, 0.25, 4.0
SEEDS = (90_000, 7_000_000, 11_111, 22_222, 33_333)

# name, kind, param, seed, sub-gate CPM floor (None = the artifact's own gate)
ARMS = [
    ("KA_q1.15_mf",  "power", 1.15, 90_000,    None),
    ("g10_q1.00_s1", "power", 1.00, 90_000,    10.0),   # identity; a pure plumbing check
    ("g10_q1.15_s1", "power", 1.15, 90_000,    10.0),
    ("g10_q1.15_s2", "power", 1.15, 7_000_000, 10.0),
    ("g10_q1.15_s3", "power", 1.15, 11_111,    10.0),
    ("g10_q1.15_s4", "power", 1.15, 22_222,    10.0),
    ("g10_q1.15_s5", "power", 1.15, 33_333,    10.0),
    ("g20_q1.15_s1", "power", 1.15, 90_000,    20.0),
    ("g20_q1.15_s2", "power", 1.15, 7_000_000, 20.0),
    ("g20_q1.15_s3", "power", 1.15, 11_111,    20.0),
    ("map_med_s1",   "map",   "median|POOLED", 90_000,    None),
    ("map_bud_s1",   "map",   "budget|POOLED", 90_000,    None),
]
if os.environ.get("ONLY"):
    _only = set(os.environ["ONLY"].split(","))
    ARMS = [a for a in ARMS if a[0] in _only]


def sharpen_sub(h, q, full, sub):
    """`sign(h)|h|^q` on `sub` with the SUB-gate row norm restored; `full \\ sub` untouched."""
    out = np.where(full, h, 0.0)
    if q == 1.0:
        return out
    x = np.where(sub, out, 0.0)
    y = np.sign(x) * np.abs(x) ** q
    n0, n1 = np.linalg.norm(x), np.linalg.norm(y)
    y = y * (n0 / max(n1, 1e-300))
    out = np.where(sub, y, out)
    return out


def apply_norm_row(h, m, gate):
    g = np.asarray(gate, bool)
    x = np.where(g, h, 0.0)
    s = float(np.sqrt(max((x ** 2).sum() / max(g.sum(), 1), 1e-300)))
    y = np.where(g, np.sign(x) * s * m(np.abs(x) / s), 0.0)
    n0, n1 = np.linalg.norm(x), np.linalg.norm(y)
    return y * (n0 / max(n1, 1e-300))


def load_builder():
    spec = importlib.util.spec_from_file_location(
        "prod_builder_188", f"{ROOT}/scripts/make_submission_pm_cd4.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_188"] = mod
    spec.loader.exec_module(mod)
    assert mod.TS == TS
    return mod


def dense(A, rows):
    m = A.X[rows]
    return np.asarray(m.todense() if hasattr(m, "todense") else m, dtype=np.float64)


def main() -> None:
    t0 = time.time()
    E = np.load(os.path.join(E169, "work/h1test-effect.npz"), allow_pickle=True)
    H0 = E["H"].astype(np.float64)
    b = E["b"].astype(np.float64)
    a_gate = np.asarray(E["gate"], bool)
    ecpm = E["cpm"].astype(np.float64)
    etargets = np.asarray([str(t) for t in E["targets"]])
    esums = E["sums"].astype(np.float64)
    egenes = np.asarray([str(g) for g in E["genes"]])
    MAPS = {k: Map.from_dict(v) for k, v in
            json.loads(open(os.path.join(WORK, "maps_norm.json")).read()).items()}
    prod = load_builder()
    for floor in (10.0, 20.0):
        sub = a_gate & (ecpm > floor)
        print(f"artifact gate {int(a_gate.sum())} genes (min_cpm {float(E['min_cpm'])}); "
              f"sub-gate CPM>{floor}: {int(sub.sum())} "
              f"({100.0*sub.sum()/a_gate.sum():.1f}% of it)", flush=True)

    R = ad.read_h5ad(REAL)
    rlab = np.asarray(R.obs[COL].astype(str).values, dtype="U64")
    genes_full = np.asarray([str(g) for g in R.var_names])
    assert np.array_equal(genes_full, egenes)
    G_full = genes_full.size
    ci = np.flatnonzero(rlab == CTRL)

    cpm = np.zeros(G_full)
    for i in range(0, ci.size, 2000):
        m = dense(R, ci[i:i + 2000])
        s = m.sum(1)
        s[s == 0] = 1.0
        cpm += (m * (1e6 / s[:, None])).sum(0)
    cpm /= ci.size
    keep = cpm > MIN_CPM
    genes = genes_full[keep]
    G = int(keep.sum())
    print(f"harness: {ci.size} controls, {G}/{G_full} genes past the {MIN_CPM}-CPM gate "
          f"[{time.time()-t0:.0f}s]", flush=True)

    cs = np.empty((G, ci.size))
    for i in range(0, ci.size, 2000):
        m = dense(R, ci[i:i + 2000])
        s = m.sum(1)
        s[s == 0] = 1.0
        cs[:, i:i + m.shape[0]] = np.log1p(m[:, keep] * (1e6 / s[:, None])).T
    cs.sort(axis=1)
    cs = np.ascontiguousarray(cs)
    tie = ctrl_tiesums(cs)
    ref_lin = cpm[keep]
    b_cell = np.log1p(cpm)
    print(f"control block ready [{time.time()-t0:.0f}s]", flush=True)

    epos = {t: i for i, t in enumerate(etargets)}
    targets = [t for t in sorted(set(rlab) - {CTRL}) if t in epos]
    T = len(targets)
    print(f"{T} perturbations, arms: {[a[0] for a in ARMS]}", flush=True)

    real = {k: np.full((T, G), np.nan) for k in ("p", "q", "l")}
    out = {n: {k: np.full((T, G), np.nan) for k in ("p", "q", "l")} for n, *_ in ARMS}
    svec = {n: np.zeros((T, G), dtype=np.float32) for n, *_ in ARMS}
    npred = {n: [] for n, *_ in ARMS}

    for i, t in enumerate(targets):
        Br = dense(R, np.flatnonzero(rlab == t))
        s = Br.sum(1)
        s[s == 0] = 1.0
        Br = np.log1p(Br[:, keep] * (1e6 / s[:, None]))
        z = wilcoxon_z(cs, tie, np.ascontiguousarray(Br.T))
        real["p"][i] = 2 * norm.sf(np.abs(z))
        real["q"][i] = false_discovery_control(real["p"][i], method="bh")
        real["l"][i] = np.log2((np.expm1(Br).mean(0) + 1e-9) / (ref_lin + 1e-9))
        del Br

        k = epos[t]
        S = float(esums[k].sum())
        h = H0[k]

        for name, kind, param, seed, floor in ARMS:
            if kind == "power":
                if floor is None:
                    hq = sharpen(h, float(param), a_gate, mode="row")
                else:
                    hq = sharpen_sub(h, float(param), a_gate, a_gate & (ecpm > floor))
            else:
                hq = apply_norm_row(h, MAPS[str(param)], a_gate)
            hq = prod.mass_consistent_effect(hq, b, a_gate, mode="minnorm")
            lin = np.clip(np.expm1(b + hq), 0.0, None)
            profile = lin / lin.sum() * (S / (N_EMIT * EMIT_L))
            totals = np.floor(N_EMIT * EMIT_L * profile).astype(np.int64)
            pcpm = 1e6 * totals / max(totals.sum(), 1)
            eff = np.abs(np.log2((pcpm + 1e-9) / (cpm + 1e-9)))
            gated = (cpm > MIN_CPM) & (b_cell > 0)
            score = np.where(gated, eff, -np.inf)
            order = np.lexsort((np.arange(G_full), -score))
            kt = min(int(KEEP_TOP), int(gated.sum()))
            sv = np.full(G_full, float(TAIL))
            sv[order[:kt]] = float(KEEP_SPREAD)

            rng = np.random.default_rng(seed + i)
            vals, cols, ptr = emit_counts(
                profile, np.full(N_EMIT, EMIT_L), 0.1, rng,
                mode="quasi", spread=sv, overdispersion=float(PHI),
            )
            B = np.zeros((N_EMIT, G_full))
            for r in range(N_EMIT):
                sl = slice(ptr[r], ptr[r + 1])
                B[r, cols[sl]] = vals[sl]
            assert np.array_equal(B.sum(0).astype(np.int64), totals), name
            ss = B.sum(1)
            ss[ss == 0] = 1.0
            nb = np.log1p(B[:, keep] * (1e6 / ss[:, None]))
            del B
            z = wilcoxon_z(cs, tie, np.ascontiguousarray(nb.T))
            p = 2 * norm.sf(np.abs(z))
            qq = false_discovery_control(p, method="bh")
            out[name]["p"][i], out[name]["q"][i] = p, qq
            out[name]["l"][i] = np.log2((np.expm1(nb).mean(0) + 1e-9) / (ref_lin + 1e-9))
            svec[name][i] = sv[keep]
            npred[name].append(int((qq < 0.05).sum()))
            del nb
        if (i + 1) % 10 == 0:
            print(f"  {i+1}/{T} [{time.time()-t0:.0f}s]", flush=True)

    for name, *_ in ARMS:
        path = f"{OUT}.{name}.npz"
        np.savez_compressed(
            path, genes=genes, perts=np.array(targets),
            real_p=real["p"], real_padj=real["q"], real_lfc=real["l"],
            pred_p=out[name]["p"], pred_padj=out[name]["q"], pred_lfc=out[name]["l"],
            spread=svec[name])
        print(f"wrote {path}   median n_pred {np.median(npred[name]):.0f}", flush=True)

    # ---------------------------------------------------------------- known-answer row
    ka = f"{OUT}.KA_q1.15_mf.npz"
    ref = os.path.join(E169, "de169.q1.15_mf.npz")
    if os.path.exists(ka) and os.path.exists(ref):
        A, Bz = np.load(ka, allow_pickle=True), np.load(ref, allow_pickle=True)
        worst = {}
        for f in ("real_lfc", "real_padj", "pred_lfc", "pred_padj", "pred_p"):
            d = np.nanmax(np.abs(np.asarray(A[f], float) - np.asarray(Bz[f], float)))
            worst[f] = float(d)
        print("\nKNOWN-ANSWER: KA_q1.15_mf vs 169's de169.q1.15_mf.npz — max |diff| per field")
        for f, d in worst.items():
            print(f"  {f:<11} {d:.3e}  {'PASS' if d == 0.0 else 'FAIL (not bit-identical)'}")
    print(f"done [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
