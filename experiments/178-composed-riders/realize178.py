"""178 — realise the composed-rider arms on h1-test and produce their DE tables.

`realize173.py`, with two changes and nothing else:

  1. the effect is the REBUILT one (`work/h1arms178.npz`, from `h1rebuild178.py`) rather than
     a transform of the effect recovered from the 095 artifact, because gamma and the donor
     weights are not expressible as transforms of a frozen whitened matrix;
  2. the arm list is 178's.

Everything downstream — the production `mass_consistent_effect`, the artifact's own per-target
library total, the a15 emission law (`tail 1.0 / keep_top 8000 / keep_spread 0.25 / phi 4`),
the DE kernel, the seeds — is byte-for-byte 173's, so 157's board calibration still applies and
the numbers are directly comparable to `experiments/173-sign-repair/de173.*.npz`.

KNOWN-ANSWER ROWS (rule 10's corollary)
---------------------------------------
  `KA_a15_raw`  164's exact path, no recovery, no mass fix: median `n_pred` 6034 and `reach`
                0.18283, EXACTLY. If this row moves, the harness moved and nothing else here
                is readable.
  `rec_base_mf` the RECOVERED a15 effect through the same mass fix and emitter — i.e. 173's
                `base_mf`, re-run here. It must land on 173's stored value.
  `rb_a15` / `rb_a15_b`  the rebuilt a15-config effect and its re-seeded twin. This is the
                base every rider delta is taken against; its distance from `rec_base_mf` is
                the price of rebuilding rather than recovering, and it is REPORTED, never
                subtracted out.
  `rb_flip`     sign repair on the rebuild. 173 measured it on the recovery at precision
                +0.00313 and nmae -0.00417. Agreement in sign and rough size is the check
                that the rebuild is a usable substitute for the recovery.

RULE 9 IS IN FORCE: `fid`, `jac` and `reach` all carry `n_real`. `six178.py` reads precision
(no `n_real`), the overlap ratio, and `reach` with its twin, exactly as `six173.py` does.

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
WORK = os.path.join(HERE, "work")
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/092-call-set-truncation")
from vcc2026.donor_transfer import emit_counts  # noqa: E402
from de_kernel import ctrl_tiesums, wilcoxon_z  # noqa: E402

REAL = f"{ROOT}/data/harness/h1-test.h5ad"
EFF = f"{ROOT}/experiments/169-sharpening-six-member/work/h1test-effect.npz"
OUT = os.path.join(HERE, "de178")
CTRL, COL, MIN_CPM = "non-targeting", "target_gene", 5.0
N_EMIT, EMIT_L = 400, 100_000
TS = 50_000.0
TAIL, KEEP_TOP, KEEP_SPREAD, PHI = 1.0, 8000, 0.25, 4.0

# name, source ("art" = the recovered artifact effect, else a key in h1arms178.npz), mass fix,
# seed, use the recovered profile straight off the shipped totals
ARMS = [
    ("KA_a15_raw",   "art",        False, 90_000,    False),
    ("rec_base_mf",  "art",        True,  90_000,    True),
    ("rb_a15",       "a15",        True,  90_000,    True),
    ("rb_a15_b",     "a15",        True,  7_000_000, True),
    ("rb_g045",      "g045",       True,  90_000,    True),
    ("rb_g030",      "g030",       True,  90_000,    True),
    ("rb_flip",      "flip",       True,  90_000,    True),
    ("rb_w211",      "w211",       True,  90_000,    True),
    ("rb_flip_w211", "flip_w211",  True,  90_000,    True),
    ("rb_flip_w211_b", "flip_w211", True, 7_000_000, True),
    ("rb_g045_flip", "g045_flip",  True,  90_000,    True),
    ("rb_g045_w211", "g045_w211",  True,  90_000,    True),
    ("rb_all3_g045", "all3_g045",  True,  90_000,    True),
    ("rb_all3_g045_b", "all3_g045", True, 7_000_000, True),
    ("rb_all3_g030", "all3_g030",  True,  90_000,    True),
    ("rb_g045_b",    "g045",       True,  7_000_000, True),
    ("rb_flip_b",    "flip",       True,  7_000_000, True),
    ("rb_w211_b",    "w211",       True,  7_000_000, True),
]
if os.environ.get("ONLY"):
    _only = set(os.environ["ONLY"].split(","))
    ARMS = [a for a in ARMS if a[0] in _only]


def load_builder():
    """The PRODUCTION `mass_consistent_effect`, not a copy of it."""
    spec = importlib.util.spec_from_file_location(
        "prod_builder_178_de", f"{ROOT}/scripts/make_submission_pm_cd4.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    assert mod.TS == TS
    return mod


def dense(A, rows):
    m = A.X[rows]
    return np.asarray(m.todense() if hasattr(m, "todense") else m, dtype=np.float64)


def main() -> None:
    t0 = time.time()
    E = np.load(EFF, allow_pickle=True)
    H_art = E["H"].astype(np.float64)
    b = E["b"].astype(np.float64)
    a_gate = E["gate"]
    etargets = np.asarray([str(t) for t in E["targets"]])
    esums = E["sums"].astype(np.float64)
    egenes = np.asarray([str(g) for g in E["genes"]])

    Z = np.load(os.path.join(WORK, "h1arms178.npz"), allow_pickle=True)
    rb_targets = np.asarray([str(t) for t in Z["targets"]])
    assert np.array_equal(np.asarray([str(g) for g in Z["genes"]]), egenes), "gene axes"
    order = np.asarray([int(np.flatnonzero(rb_targets == t)[0]) for t in etargets])
    REB = {k[2:]: Z[k][order].astype(np.float64) for k in Z.files if k.startswith("H_")}
    assert np.array_equal(np.asarray(Z["gate"], bool), np.asarray(a_gate, bool)), "gates"
    print(f"recovered effect {H_art.shape}; rebuilt arms {sorted(REB)}; "
          f"artifact gate {int(np.asarray(a_gate).sum())} genes", flush=True)

    prod = load_builder()
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
    print(f"{T} perturbations, {len(ARMS)} arms: {[a[0] for a in ARMS]}", flush=True)

    real = {k: np.full((T, G), np.nan) for k in ("p", "q", "l")}
    out = {n: {k: np.full((T, G), np.nan) for k in ("p", "q", "l")} for n, *_ in ARMS}
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

        for name, src, mf, seed, rec in ARMS:
            if not rec:
                profile = (esums[k] + 0.5) / (N_EMIT * EMIT_L)
            else:
                hq = (H_art[k] if src == "art" else REB[src][k]).copy()
                if mf:
                    hq = prod.mass_consistent_effect(hq, b, a_gate, mode="minnorm")
                lin = np.clip(np.expm1(b + hq), 0.0, None)
                profile = lin / lin.sum() * (S / (N_EMIT * EMIT_L))
            totals = np.floor(N_EMIT * EMIT_L * profile).astype(np.int64)
            pcpm = 1e6 * totals / max(totals.sum(), 1)
            eff = np.abs(np.log2((pcpm + 1e-9) / (cpm + 1e-9)))
            gated = (cpm > MIN_CPM) & (b_cell > 0)
            score = np.where(gated, eff, -np.inf)
            order_g = np.lexsort((np.arange(G_full), -score))
            kt = min(int(KEEP_TOP), int(gated.sum()))
            sv = np.full(G_full, float(TAIL))
            sv[order_g[:kt]] = float(KEEP_SPREAD)

            rng = np.random.default_rng(seed + i)
            vals, cols, ptr = emit_counts(
                profile, np.full(N_EMIT, EMIT_L), 0.1, rng,
                mode="quasi", spread=sv, overdispersion=float(PHI))
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
            npred[name].append(int((qq < 0.05).sum()))
            del nb
        if (i + 1) % 10 == 0:
            print(f"  {i+1}/{T} [{time.time()-t0:.0f}s]", flush=True)

    for name, *_ in ARMS:
        path = f"{OUT}.{name}.npz"
        np.savez_compressed(
            path, genes=genes, perts=np.array(targets),
            real_p=real["p"], real_padj=real["q"], real_lfc=real["l"],
            pred_p=out[name]["p"], pred_padj=out[name]["q"], pred_lfc=out[name]["l"])
        print(f"wrote {path}   median n_pred {np.median(npred[name]):.0f}", flush=True)
    print(f"done [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
