"""173 — realize SIGN-REPAIR arms (whitening-flipped signs set to the un-whitened consensus) on h1-test and produce their DE tables.

WHY h1-test AND NOT THE ORION PAIR
----------------------------------
157 §1: h1-test reproduces BOTH board raw `reach` points to within 0.014
(a8 0.1890 board / 0.2029 here; a15 0.1811 / 0.1828) while repaired HCT116 misses a15's by
2.6x.  It is the only surface validated on this member.  Donors are `k562_gwps` +
`orion_hek293t`; `data/donors/h1_2025` covers 85/85 of h1-test's targets and is pure oracle
leakage -- it never enters, here or in the ordering keys.

WHAT IS DIFFERENT FROM 157/164
------------------------------
Those re-emit from the built prediction's per-gene panel totals, i.e. they hold the PROFILE
fixed and move only the emission law.  Sharpening moves the profile, so the profile has to be
rebuilt.  `probe_effect.py` recovers the artifact's log-space effect `H` from its emitted
counts and proves the recovery (off-gate across-target sd 0.0015 against an on-gate 0.0656,
and it re-derives the artifact's `--min-cpm 5.0` from the data).  Here that `H` is sharpened,
optionally mass-fixed with the PRODUCTION function, re-rendered, and re-emitted through the
UNMODIFIED production `emit_counts`, keeping every other convention of `realize164.py`
byte-for-byte so 157's board calibration still applies.

The per-target panel total is preserved exactly, so each arm emits the same library sizes as
the shipped artifact.

KNOWN-ANSWER ROWS (rule 10's corollary), in increasing strength:
  * `KA_a15_raw` -- 164's exact code path, no recovery, no mass fix.  Must reproduce 157's
    median `n_pred` 6035 and `reach` 0.18283 EXACTLY (same seed, same everything).
  * `q1.00` -- the recovery round trip at the identity exponent.  Must land on `KA_a15_raw`
    inside the re-seed noise, which `q1.00_b` measures.  If the recovery were broken this row
    is where it shows, and it is the only row whose answer is known in advance.

No cell-eval2, no build, no submission.
"""
from __future__ import annotations

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
sys.path.insert(0, f"{ROOT}/experiments/169-sharpening-six-member")
from vcc2026.donor_transfer import emit_counts  # noqa: E402
from de_kernel import ctrl_tiesums, wilcoxon_z  # noqa: E402

REAL = f"{ROOT}/data/harness/h1-test.h5ad"
OUT = os.path.join(HERE, "de173")
CTRL, COL, MIN_CPM = "non-targeting", "target_gene", 5.0
N_EMIT, EMIT_L = 400, 100_000
TS = 50_000.0
# a15's shipped emission law: --overdispersion 4.0 --tail-spread 1.0 --keep-top 8000
#                             --keep-spread 0.25
TAIL, KEEP_TOP, KEEP_SPREAD, PHI = 1.0, 8000, 0.25, 4.0

# name, transform, mass_fix, seed, use_recovery, (unused)
ARMS = [
    ("KA_a15_raw",   "id",       False, 90_000,    False, "row"),   # 164's exact path: must read n_pred 6035 / reach 0.18283
    ("base_mf",      "id",       True,  90_000,    True,  "row"),   # recovery round trip + production mass fix (= 169's q1.00_mf)
    ("base_mf_b",    "id",       True,  7_000_000, True,  "row"),
    ("flip_mf",      "flip",     True,  90_000,    True,  "row"),   # sign -> cons where whitening disagrees, |H| kept
    ("flip_mf_b",    "flip",     True,  7_000_000, True,  "row"),
    ("zero_mf",      "zero",     True,  90_000,    True,  "row"),   # zero the disagreeing entries
    ("permflip_mf",  "permflip", True,  90_000,    True,  "row"),   # NULL: flip the same count of random gated entries
]
if os.environ.get("ONLY"):
    _only = set(os.environ["ONLY"].split(","))
    ARMS = [a for a in ARMS if a[0] in _only]


def load_builder():
    """The PRODUCTION `mass_consistent_effect`, not a copy of it (156's lesson, a19's cost)."""
    path = f"{ROOT}/scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("prod_builder_169", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_169"] = mod
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
    # ---- 172's un-whitened consensus (DE-universe axis) -> full axis, keyed by names -------
    sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
    from debench import Surface as _Surface
    _S = _Surface(f"{ROOT}/experiments/164-callset-artefact/de164.a15.npz")
    _F = np.load(f"{ROOT}/experiments/172-sign-confidence-key/work/h1test-feats.npz", allow_pickle=True)
    _gpos = {g: i for i, g in enumerate(egenes)}
    _de_gi = np.asarray([_gpos[g] for g in _S.genes])
    _sp = {p: i for i, p in enumerate(_S.perts)}
    CONS = np.zeros_like(H0)
    for k, t in enumerate(etargets):
        CONS[k, _de_gi] = _F["cons"][_sp[t]].astype(np.float64)
    FLIP = (np.sign(H0) != 0) & (np.sign(CONS) != 0) & (np.sign(H0) == -np.sign(CONS))
    _r = np.random.default_rng(173)
    PERM = np.zeros_like(FLIP)
    for k in range(FLIP.shape[0]):
        _idx = np.flatnonzero(a_gate); PERM[k, _r.choice(_idx, int(FLIP[k].sum()), replace=False)] = True
    print(f"whitening-flipped entries on the gate: {FLIP[:, a_gate].mean():.2%}  (energy share {(H0[FLIP]**2).sum()/(H0**2).sum():.2%})", flush=True)

    R = ad.read_h5ad(REAL)
    rlab = np.asarray(R.obs[COL].astype(str).values, dtype="U64")
    genes_full = np.asarray([str(g) for g in R.var_names])
    assert np.array_equal(genes_full, egenes), "harness / effect gene axes differ"
    G_full = genes_full.size
    ci = np.flatnonzero(rlab == CTRL)

    cpm = np.zeros(G_full)
    for i in range(0, ci.size, 2000):
        m = dense(R, ci[i:i + 2000])
        s = m.sum(1)
        s[s == 0] = 1.0
        cpm += (m * (1e6 / s[:, None])).sum(0)
    cpm /= ci.size
    keep = cpm > MIN_CPM                    # the SCORER's DE universe (mean-per-cell CPM > 5)
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
    b_cell = np.log1p(cpm)                  # realize164's gate convention, kept verbatim
    print(f"control block ready [{time.time()-t0:.0f}s]", flush=True)

    epos = {t: i for i, t in enumerate(etargets)}
    targets = [t for t in sorted(set(rlab) - {CTRL}) if t in epos]
    T = len(targets)
    print(f"{T} perturbations, arms: {[a[0] for a in ARMS]}", flush=True)

    real = {k: np.full((T, G), np.nan) for k in ("p", "q", "l")}
    out = {n: {k: np.full((T, G), np.nan) for k in ("p", "q", "l")} for n, *_ in ARMS}
    svec = {n: np.zeros((T, G), dtype=np.float32) for n, *_ in ARMS}
    npred = {n: [] for n, *_ in ARMS}
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

        for name, q, mf, seed, rec, rmode in ARMS:
            if not rec:
                # 164's exact path: the profile straight off the shipped totals.
                profile = (esums[k] + 0.5) / (N_EMIT * EMIT_L)
            else:
                if q == "id":
                    hq = h.copy()
                elif q == "flip":
                    hq = np.where(FLIP[k], -h, h)
                elif q == "zero":
                    hq = np.where(FLIP[k], 0.0, h)
                elif q == "permflip":
                    hq = np.where(PERM[k], -h, h)
                else:
                    raise ValueError(q)
                if mf:
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
            if name == ARMS[0][0]:
                base_top = set(order[:kt].tolist())
            else:
                keeptop_overlap[name].append(len(base_top & set(order[:kt].tolist())) / kt)

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
            spread=svec[name],
        )
        ov = keeptop_overlap[name]
        print(f"wrote {path}   median n_pred {np.median(npred[name]):.0f}"
              + (f"   keep_top overlap vs incumbent {np.mean(ov):.5f}" if ov else ""),
              flush=True)
    print(f"done [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
