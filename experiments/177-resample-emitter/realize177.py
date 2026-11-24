"""177 — realize the RESAMPLING emission law on h1-test and produce its DE tables.

WHY h1-test AND NOT THE ORION PAIR
----------------------------------
157 section 1: h1-test reproduces BOTH board raw `reach` points to within 0.014 (a8 0.1890
board / 0.2029 here; a15 0.1811 / 0.1828) while repaired HCT116 misses a15's by 2.6x.  It is
the only surface validated on this member, and 165 restricts that validation to EMISSION-LAW
changes at a fixed profile -- which is exactly what every arm here is.  Donors are `k562_gwps`
+ `orion_hek293t`; `data/donors/h1_2025` never enters.

WHAT CHANGED FROM THE COORDINATOR'S FIRST PASS (`work/realize177.log`)
---------------------------------------------------------------------
That pass rescaled every real control cell from its own ~20k library to the harness's nominal
`EMIT_L = 100_000` before resampling, while the emitted panel total stayed pinned at the
artifact's own ~20k/cell depth.  `emit_counts`'s exact-total repair then had to peel ~80% of
every gene's column, one count per cell per round, which is a maximally variance-DESTROYING
operation: it called 10,776 of 10,779 genes significant.  `signal177.py` reproduces that
geometry and reads 10,658/10,779, so the diagnosis is exact.

Three consequences, all implemented here:

  1. **`EMIT_L` was never the emitted depth.**  `profile` is scaled by `S/(N_EMIT*EMIT_L)`
     where `S` is the artifact's own panel total, so `totals = floor(N*L*profile) = floor(S*p)`
     regardless of `L`.  The artifact's `S/400` is **20,250 counts/cell** and the real h1-test
     controls carry a median of 20,190 -- **the shipped quasi arms are ALREADY at real depth**,
     and a `quasi20` arm is bit-identical to `base_mf` by construction (asserted below, not
     run as a separate row: an identity is not a control).  The informative depth control is
     the other direction: `quasi100`, the same profile emitted at a genuine 5x depth.
  2. The resample arms pass the sampled cells' OWN libraries and RAW counts, so
     `L = median(lib) ~ 20,190` and the repair only has the ~4% mean-vs-median residual to fix.
  3. `control_profile` is the POOLED control profile `expm1(b)/sum`, the space `b` and the
     rendered `profile` live in -- not the mean-per-cell-CPM profile, which differs from it by
     an sd of 0.020 in log2 and would make `fold != 1` at zero predicted effect.

KNOWN-ANSWER ROWS (rule 10's corollary):
  * `KA_a15_raw` -- 164's exact code path, no recovery, no mass fix.  Must reproduce 157's
    median `n_pred` 6035 (169/173 read 6034) and `reach` 0.18283 EXACTLY.
  * `base_mf`    -- 173's row: `reach` 0.18432, median `n_pred` 6066.
  * the recovery round trip `|base_mf - KA_a15_raw|` against the matched-seed spread.

No cell-eval2, no build, no submission.
"""
from __future__ import annotations

import hashlib
import importlib.util
import os
import sys
import time

import anndata as ad
import numpy as np
from scipy.stats import false_discovery_control, norm

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = f"{ROOT}/experiments/169-sharpening-six-member/work"   # reuse 169's recovered effect
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/092-call-set-truncation")
sys.path.insert(0, HERE)
from signal177 import fracpair  # noqa: E402
import vcc2026.donor_transfer as _dt  # noqa: E402


def _fast_repair_totals(block, totals, rng):
    """Vectorised equivalent of `donor_transfer._repair_totals`: same contract (each gene's
    column total forced to `totals`, single counts moved between cells, removals uniform over
    the cells that still hold one).  The original's one-count-at-a-time Python loop costs
    ~2 min per perturbation on a resample arm, where every gene carries a residual."""
    out = np.asarray(block, dtype=np.int64).copy()
    n = out.shape[0]
    resid = np.asarray(totals, dtype=np.int64) - out.sum(0)
    for j in np.flatnonzero(resid > 0):
        np.add.at(out[:, j], rng.integers(0, n, int(resid[j])), 1)
    for j in np.flatnonzero(resid < 0):
        d = int(-resid[j])
        col = out[:, j]
        while d > 0:
            nz = np.flatnonzero(col > 0)
            if nz.size == 0:
                break
            take = rng.choice(nz, size=min(d, nz.size), replace=False)
            col[take] -= 1
            d -= take.size
        out[:, j] = col
    return out


_dt._repair_totals = _fast_repair_totals
from vcc2026.donor_transfer import emit_counts  # noqa: E402
from de_kernel import ctrl_tiesums, wilcoxon_z  # noqa: E402

REAL = f"{ROOT}/data/harness/h1-test.h5ad"
OUT = os.path.join(HERE, "de177")
CTRL, COL, MIN_CPM = "non-targeting", "target_gene", 5.0
N_EMIT, EMIT_L = 400, 100_000
TS = 50_000.0
# a15's shipped emission law: --overdispersion 4.0 --tail-spread 1.0 --keep-top 8000
#                             --keep-spread 0.25
TAIL, KEEP_TOP, KEEP_SPREAD, PHI = 1.0, 8000, 0.25, 4.0
SCHED_REPS = 2                       # fracpair passes; `signal177.py` grid
CV_Q = 60.0                          # expression-compensating target quantile of the control CV
CPM_FLOOR = 15.0                     # INBOX-4 section 3b's matchable floor

# name, law, mass_fix, seed, use_recovery, depth, jitter
#   law: "quasi" | "rs" (all genes resampled) | "hy<K>" (top-K |effect| genes kept on the
#        parametric quasi law, the rest resampled -- `param_mask`'s designed use)
#   jitter: None | "sched" | "sched_perm" | "hkey"
ARMS = [
    ("KA_a15_raw",    "quasi",  False, 90_000,    False, "real", None),
    ("base_mf",       "quasi",  True,  90_000,    True,  "real", None),
    ("base_mf_b",     "quasi",  True,  7_000_000, True,  "real", None),
    ("quasi100",      "quasi",  True,  90_000,    True,  "x5",   None),
    ("rs_all",        "rs",     True,  90_000,    True,  "real", None),
    ("rs_all_b",      "rs",     True,  7_000_000, True,  "real", None),
    ("rs_sched",      "rs",     True,  90_000,    True,  "real", "sched"),
    ("rs_sched_b",    "rs",     True,  7_000_000, True,  "real", "sched"),
    ("rs_sched_perm", "rs",     True,  90_000,    True,  "real", "sched_perm"),
    ("rs_hkey",       "rs",     True,  90_000,    True,  "real", "hkey"),
    ("hy_kt2000",     "hy2000", True,  90_000,    True,  "real", None),
    ("hy_kt2000_b",   "hy2000", True,  7_000_000, True,  "real", None),
    ("hy_kt4000",     "hy4000", True,  90_000,    True,  "real", None),
    ("hy_kt6000",     "hy6000", True,  90_000,    True,  "real", None),
    ("hy_kt8000",     "hy8000", True,  90_000,    True,  "real", None),
    ("rs_hkey_b",     "rs",     True,  7_000_000, True,  "real", "hkey"),
]
if os.environ.get("ONLY"):
    _only = set(os.environ["ONLY"].split(","))
    ARMS = [a for a in ARMS if a[0] in _only]


def load_builder():
    """The PRODUCTION `mass_consistent_effect`, not a copy of it (156's lesson, a19's cost)."""
    path = f"{ROOT}/scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("prod_builder_177", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_177"] = mod
    spec.loader.exec_module(mod)
    assert mod.TS == TS, f"builder TS {mod.TS} != {TS}"
    return mod


def dense(A, rows):
    m = A.X[rows]
    return np.asarray(m.todense() if hasattr(m, "todense") else m, dtype=np.float64)


def main() -> None:
    t0 = time.time()
    E = np.load(os.path.join(WORK, "h1test-effect.npz"), allow_pickle=True)
    H0 = E["H"].astype(np.float64)
    b = E["b"].astype(np.float64)
    a_gate = E["gate"]                      # the ARTIFACT's own gate (measured: min_cpm 5.0)
    etargets = np.asarray([str(t) for t in E["targets"]])
    esums = E["sums"].astype(np.float64)
    egenes = np.asarray([str(g) for g in E["genes"]])
    print(f"recovered effect: {H0.shape[0]} targets x {H0.shape[1]} genes, "
          f"artifact gate {int(a_gate.sum())} genes at --min-cpm {float(E['min_cpm'])}",
          flush=True)
    prod = load_builder()
    R = ad.read_h5ad(REAL)
    rlab = np.asarray(R.obs[COL].astype(str).values, dtype="U64")
    genes_full = np.asarray([str(g) for g in R.var_names])
    assert np.array_equal(genes_full, egenes), "harness / effect gene axes differ"
    G_full = genes_full.size
    ci = np.flatnonzero(rlab == CTRL)

    cpm = np.zeros(G_full)
    cpm2 = np.zeros(G_full)
    ctot = np.zeros(G_full)
    clib = []
    for i in range(0, ci.size, 2000):
        m = dense(R, ci[i:i + 2000])
        s = m.sum(1)
        clib.append(s.copy())
        s[s == 0] = 1.0
        d = m * (1e6 / s[:, None])
        cpm += d.sum(0)
        cpm2 += (d * d).sum(0)
        ctot += m.sum(0)
    cpm /= ci.size
    cpm2 /= ci.size
    clib = np.concatenate(clib)
    cv_ctrl = np.sqrt(np.maximum(cpm2 - cpm ** 2, 0.0)) / np.maximum(cpm, 1e-12)
    keep = cpm > MIN_CPM                    # the SCORER's DE universe (mean-per-cell CPM > 5)
    lin_b = np.clip(np.expm1(b), 0.0, None)
    ctrl_rel = lin_b / lin_b.sum()          # POOLED control profile -- the space `b` lives in
    pooled_real = ctot / ctot.sum()
    dev = np.abs(np.log2((ctrl_rel[keep] + 1e-12) / (pooled_real[keep] + 1e-12)))
    assert dev.max() < 1e-3, f"b is not the harness pooled control profile (max {dev.max():.2e})"
    genes = genes_full[keep]
    G = int(keep.sum())
    print(f"harness: {ci.size} controls, {G}/{G_full} genes past the {MIN_CPM}-CPM gate; "
          f"control lib median {np.median(clib):.0f}; artifact depth "
          f"{np.median(esums.sum(1))/N_EMIT:.0f} counts/cell [{time.time()-t0:.0f}s]", flush=True)

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
    b_cell = np.log1p(cpm)                  # realize164's gate convention, kept verbatim
    print(f"control block ready [{time.time()-t0:.0f}s]", flush=True)

    # ---------------------------------------------------------------- the noise schedule
    # `signal177.py` measured that the ONLY kernel keeping the resample law's zero-effect null
    # is `fracpair`, and `null177b.py` that absolute/multiplicative kernels break it.
    gated_full = np.flatnonzero((cpm > MIN_CPM) & (b_cell > 0))
    ok = gated_full[cpm[gated_full] > CPM_FLOOR]
    T_cv = float(np.percentile(cv_ctrl[ok], CV_Q))
    sched = np.zeros(G_full)
    sched[ok] = np.clip(np.sqrt(np.maximum(T_cv ** 2 - cv_ctrl[ok] ** 2, 0.0)), 0.0, 1.0)
    perm_rng = np.random.default_rng(4242)
    sched_perm = np.zeros(G_full)
    sched_perm[gated_full] = perm_rng.permutation(sched[gated_full])
    print(f"schedule: T(cv q{CV_Q:.0f}) {T_cv:.3f} over {ok.size} genes > {CPM_FLOOR} CPM; "
          f"sigma median {np.median(sched[ok]):.3f}, nonzero {(sched > 0).sum()}, "
          f"{SCHED_REPS} fracpair passes", flush=True)

    epos = {t: i for i, t in enumerate(etargets)}
    targets = [t for t in sorted(set(rlab) - {CTRL}) if t in epos]
    if os.environ.get("MAXT"):                      # smoke-test path only
        targets = targets[:int(os.environ["MAXT"])]
    T = len(targets)
    print(f"{T} perturbations, arms: {[a[0] for a in ARMS]}", flush=True)

    real = {k: np.full((T, G), np.nan) for k in ("p", "q", "l")}
    out = {n: {k: np.full((T, G), np.nan) for k in ("p", "q", "l")} for n, *_ in ARMS}
    svec = {n: np.zeros((T, G), dtype=np.float32) for n, *_ in ARMS}
    npred = {n: [] for n, *_ in ARMS}
    digest = {n: [] for n, *_ in ARMS}
    keeptop_overlap = {n: [] for n, *_ in ARMS}
    base_top = None

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
        S = float(esums[k].sum())            # preserve the artifact's own library size
        h = H0[k]

        for name, law, mf, seed, rec, depth, jit in ARMS:
            if not rec:
                # 164's exact path: the profile straight off the shipped totals.
                prof_raw = (esums[k] + 0.5) / (N_EMIT * EMIT_L)
                profile = prof_raw
                totals = np.floor(N_EMIT * EMIT_L * profile).astype(np.int64)
                p_rel = profile / profile.sum()
            else:
                hq = h.copy()
                if mf:
                    hq = prod.mass_consistent_effect(hq, b, a_gate, mode="minnorm")
                lin = np.clip(np.expm1(b + hq), 0.0, None)
                p_rel = lin / lin.sum()
                profile = p_rel * (S / (N_EMIT * EMIT_L))
                totals = np.floor(N_EMIT * EMIT_L * profile).astype(np.int64)
            pcpm = 1e6 * totals / max(totals.sum(), 1)
            eff = np.abs(np.log2((pcpm + 1e-9) / (cpm + 1e-9)))
            gated = (cpm > MIN_CPM) & (b_cell > 0)
            score = np.where(gated, eff, -np.inf)
            order = np.lexsort((np.arange(G_full), -score))
            ktop = int(law[2:]) if law.startswith("hy") else KEEP_TOP
            kt = min(int(ktop), int(gated.sum()))
            sv = np.full(G_full, float(TAIL))
            sv[order[:kt]] = float(KEEP_SPREAD)
            if name == ARMS[0][0]:
                base_top = set(order[:kt].tolist())
            elif not law.startswith("hy"):
                keeptop_overlap[name].append(len(base_top & set(order[:kt].tolist())) / kt)

            rng = np.random.default_rng(seed + i)
            if law == "quasi":
                if depth == "x5":
                    # genuine 5x emitted depth: totals scale, the relative profile does not
                    libs, pp = np.full(N_EMIT, float(EMIT_L)), p_rel
                else:
                    libs, pp = np.full(N_EMIT, float(EMIT_L)), profile
                totals = np.floor(N_EMIT * EMIT_L * pp).astype(np.int64)
                vals, cols, ptr = emit_counts(
                    pp, libs, 0.1, rng, mode="quasi", spread=sv, overdispersion=float(PHI),
                )
            else:
                cidx = rng.choice(ci, N_EMIT, replace=False)
                Cc = dense(R, np.sort(cidx))       # RAW counts, own libraries, no rescaling
                dlib = Cc.sum(1)
                dlib[dlib == 0] = 1.0
                # `emit_counts` uses `L = int(min(median(lib_sizes), max_counts))`; mirror the
                # int() exactly or the totals we pin differ from the ones it repairs to
                Ld = float(int(min(np.median(dlib), 1_000_000)))
                if jit is not None:
                    sg = {"sched": sched, "sched_perm": sched_perm}.get(jit)
                    if sg is None:                 # "hkey": sigma decreasing in |H| rank
                        a = np.abs(hq)[gated_full]
                        r = np.empty(a.size)
                        r[np.argsort(a, kind="stable")] = np.arange(a.size) / max(a.size - 1, 1)
                        sg = np.zeros(G_full)
                        sg[gated_full] = 1.0 - r
                    Cc = fracpair(Cc, sg, rng, reps=SCHED_REPS)
                # pin the panel total to the artifact's own S, exactly as every other arm
                pp = p_rel * (S / (N_EMIT * Ld))
                totals = np.floor(N_EMIT * Ld * pp).astype(np.int64)
                pmask = np.zeros(G_full, dtype=bool)
                if law.startswith("hy"):
                    pmask[order[:kt]] = True
                vals, cols, ptr = emit_counts(
                    pp, dlib, 0.1, rng, mode="resample_mix", spread=sv,
                    overdispersion=float(PHI), control_counts=Cc,
                    control_profile=ctrl_rel, param_mask=pmask,
                )
                del Cc
            B = np.zeros((N_EMIT, G_full))
            for r in range(N_EMIT):
                sl = slice(ptr[r], ptr[r + 1])
                B[r, cols[sl]] = vals[sl]
            got = B.sum(0).astype(np.int64)
            assert np.array_equal(got, totals), (
                name, int(np.abs(got - totals).sum()))
            digest[name].append(hashlib.sha256(got.tobytes()).hexdigest()[:16])
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

    # ---------------------------------------------------------------- pseudobulk identity
    # 171 section 6's NON-VACUOUS form: the sha256 of the DELIVERED per-gene panel total per
    # (arm, target).  Without `emit_counts`'s exact-total repair the Poisson residual lands in
    # the totals and every arm's digest differs, so this reads differently when the thing it
    # tests is absent.  `quasi100` is expected to differ -- it emits 5x the counts on purpose.
    print("\n--- delivered per-gene panel-total digests (85 targets each) ---")
    ref_dig = digest[ARMS[1][0]] if len(ARMS) > 1 else None
    for name, *_ in ARMS:
        same = sum(int(a == bb) for a, bb in zip(digest[name], ref_dig)) if ref_dig else -1
        print(f"  {name:<15} {digest[name][0]}  matches base_mf on {same}/{len(digest[name])}")

    for name, *_ in ARMS:
        path = f"{OUT}.{name}.npz"
        np.savez_compressed(
            path, genes=genes, perts=np.array(targets),
            real_p=real["p"], real_padj=real["q"], real_lfc=real["l"],
            pred_p=out[name]["p"], pred_padj=out[name]["q"], pred_lfc=out[name]["l"],
            spread=svec[name],
        )
        ov = keeptop_overlap[name]
        print(f"wrote {path}   median n_pred {np.median(npred[name]):.0f}"
              + (f"   keep_top overlap vs incumbent {np.mean(ov):.5f}" if ov else ""),
              flush=True)
    print(f"done [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
