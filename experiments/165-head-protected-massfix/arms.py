"""165 — price the mass-correction PLACEMENT family on h1-test, `reach` included.

a23 shipped `VCC_MASS_FIX=minnorm` and lost `-0.007933` scaled `reach` while gaining `pds`
+0.001737 and `mse` +0.001525.  156's `weighted.py` had already measured that the PROFILE gain
is insensitive to where the correction goes.  So the question is whether the `reach` toll is
PLACEMENT-dependent, and this script is the instrument.

WHAT VARIES, AND WHAT IS HELD
-----------------------------
Every arm is the same recovered intended effect `D` (see `recover.py`) plus a correction
`eps = lambda * u * w`, `w_g = exp(b_g + D_g)`, confined to the CPM gate, solving

    sum_g clip(expm1(b + D + eps), 0) == TS.

Only `u` changes.  In particular the EMISSION SPREAD VECTOR IS HELD FIXED across arms, exactly
as production does it: `make_submission_pm_cd4.py:2119` calls `emission_spread_vector` on
`effects[target_index]`, the effect BEFORE `profile_from_effect` applies the mass fix, so the
keep/tail split cannot move between mass-fix arms.  Recomputing it per arm would confound the
placement contrast with a call-set change.

Arms are applied as a DELTA to the frozen prediction's own emitted profile,
`emit = profile_a15 + (p_arm - p_a15_recon) * Z`, so the `a15` arm reproduces the frozen
artifact's per-gene totals bit-for-bit and the recovery's ~0.2% off-gate quantisation residual
cancels out of every contrast instead of leaking into it.

KNOWN-ANSWER ROWS (AGENTS.md 11.1 rule 10's corollary) — the run ABORTS without them:
  * `a15` (tier8000 / tail 1.0 / keep 0.25 / phi 4) must land at raw reach 0.18283
  * `a8`  (flat 0.25 / phi 4)                        must land at raw reach 0.20291
  both from 157 `power_h1test.log` at f = 1.0, which are themselves within 0.014 of the two
  board points 0.1811 and 0.1890.  They are checked in `price.py`, on the tables this writes.
  A third gate is asserted HERE: the generalised `u`-solver must reproduce production's
  `mass_consistent_effect` to machine precision at `u = 1` (minnorm) and `u = 1/w` (const).

CONTROLS THAT WOULD READ SOMETHING IF THE EFFECT WERE ABSENT (the trap of the last two days):
  * `eps_tail` is `eps_keep` with the support REVERSED.  If placement were irrelevant the two
    would agree; they cannot both be right.
  * `mn_half` removes only half the mass excess.  If the toll is driven by |eps| it halves;
    if it does not move, the toll is not the correction's size and the whole model is wrong.
  * every arm is run at two emission seeds, so the seed floor is measured, not assumed.

DONORS: the prediction is 095's h1-test artifact — `k562_gwps` + `orion_hek293t`.
`data/donors/h1_2025` covers 85/85 h1-test targets and is pure oracle leakage; it never enters.

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
from vcc2026.donor_transfer import emit_counts        # noqa: E402
from de_kernel import ctrl_tiesums, wilcoxon_z        # noqa: E402

REAL = f"{ROOT}/data/harness/h1-test.h5ad"
CTRL, COL, MIN_CPM = "non-targeting", "target_gene", 5.0
N_EMIT, EMIT_L = 400, 100_000
NL = float(N_EMIT * EMIT_L)
TS = 5.0e4
SEEDS = tuple(int(x) for x in os.environ.get("SEEDS", "90000,7000000").split(","))

# name, u-spec, emission law ("tier" = a15 shipped, "flat" = a8), fraction of the mass
# excess to remove
ARMS = [
    ("a15",       None,            "tier", 0.0),   # known-answer, and the head's source
    ("a8",        None,            "flat", 0.0),   # known-answer
    ("mn",        ("const_u", 1.0), "tier", 1.0),  # a23 as shipped: the calibration arm
    ("mn_half",   ("const_u", 1.0), "tier", 0.5),  # size control
    ("cst",       ("inv_w", None),  "tier", 1.0),  # 156's known-bad direction
    ("hp_p500",   ("head_p", 500),  "tier", 1.0),
    ("hp_p2000",  ("head_p", 2000), "tier", 1.0),
    ("hp_p6000",  ("head_p", 6000), "tier", 1.0),
    ("hp_p9000",  ("head_p", 9000), "tier", 1.0),
    ("hp_w2000",  ("head_w", 2000), "tier", 1.0),
    ("eps_keep",  ("keepblk", None), "tier", 1.0),
    ("eps_tail",  ("tailblk", None), "tier", 1.0),  # the REVERSED control
    ("inveff",    ("inv_eff", 1e-3), "tier", 1.0),
    ("eff",       ("eff", 1e-2),     "tier", 1.0),
    # `mn_c10` is NOT a placement: it is minnorm with the CORRECTION GATE tightened to
    # CPM > 10, a15's shipped `--min-cpm`. The 095 h1-test artifact was built at the builder
    # default 5.0, so its build gate and the scorer's 5-CPM DE universe coincide (10 783 vs
    # 10 779 genes) and there is NO band of genes that the scorer reads but the correction
    # cannot touch. Production HAS such a band: ~1 000 genes between 5 and 10 CPM which the
    # mass constraint lifts by exactly S/TS and no placement can reach. This arm puts that
    # band back, and is the one structural difference between this instrument and the board
    # that could explain a23's reach loss.
    ("mn_c10",    ("gate10", 10.0), "tier", 1.0),
]
if os.environ.get("ONLY"):
    keep = set(os.environ["ONLY"].split(","))
    ARMS = [a for a in ARMS if a[0] in keep]
NT = int(os.environ.get("NT", "0"))          # target cap, for smoke tests


# ----------------------------------------------------------------- the generalised solver
def mass_fix_u(h0, b, gate, u, frac=1.0, iters=64, rtol=1e-12):
    """`eps = lambda * u * w` on the live gate, Newton-solved to the mass constraint.

    `u = 1` is production's `minnorm`; `u = 1/w` is production's `const`.  `frac < 1` removes
    only that fraction of the mass excess and exists as the size control.
    """
    h = np.asarray(h0, dtype=np.float64).copy()
    lin0 = float(np.clip(np.expm1(b + h), 0.0, None).sum())
    target = lin0 - frac * (lin0 - TS)
    for _ in range(iters):
        live = gate & (b + h > 0.0)
        residual = float(np.clip(np.expm1(b + h), 0.0, None).sum() - target)
        if abs(residual) <= rtol * target:
            break
        w = np.where(live, np.exp(b + h), 0.0)
        d = np.where(live, u, 0.0) * w
        denom = float(d @ w)
        if not np.isfinite(denom) or denom <= 0.0:
            break
        h = h - (residual / denom) * d
    return h


def _production():
    path = os.path.join(ROOT, "scripts/make_submission_pm_cd4.py")
    spec = importlib.util.spec_from_file_location("prod_builder_165", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_165"] = mod
    spec.loader.exec_module(mod)
    assert mod.TS == TS and mod.MASS_FIX is None
    return mod


def dense(A, rows):
    m = A.X[rows]
    return np.asarray(m.todense() if hasattr(m, "todense") else m, dtype=np.float64)


def main() -> None:
    t0 = time.time()
    z = np.load(os.path.join(WORK, "effects.npz"), allow_pickle=True)
    D_all = z["effect"]
    b = z["b_ctrl"]
    gate = z["gate"]
    tot_all = z["totals"]
    ctrl_cpm = z["ctrl_cpm"]
    targets = [str(t) for t in z["targets"]]
    genes_full = np.asarray([str(g) for g in z["genes"]])
    G_full = genes_full.size

    # ---- GATE 3: the generalised solver must BE production's, at both shipped modes -----
    prod = _production()
    d0 = D_all[0]
    w0 = np.where(gate & (b + d0 > 0), np.exp(b + d0), 0.0)
    a = mass_fix_u(d0, b, gate, np.ones(G_full))
    bq = prod.mass_consistent_effect(d0, b, gate, mode="minnorm")
    c = mass_fix_u(d0, b, gate, np.where(w0 > 0, 1.0 / np.where(w0 > 0, w0, 1.0), 0.0))
    dq = prod.mass_consistent_effect(d0, b, gate, mode="const")
    print(f"solver check vs production: minnorm max|delta| {np.abs(a - bq).max():.3e}   "
          f"const max|delta| {np.abs(c - dq).max():.3e}")
    assert np.abs(a - bq).max() < 1e-12 and np.abs(c - dq).max() < 1e-12, \
        "ABORT: the generalised u-solver is not production's mass_consistent_effect"

    # ---- the harness: DE gene universe, sorted control block ---------------------------
    R = ad.read_h5ad(REAL)
    rlab = np.asarray(R.obs[COL].astype(str).values, dtype="U64")
    assert np.array_equal(np.asarray([str(g) for g in R.var_names]), genes_full)
    ci = np.flatnonzero(rlab == CTRL)
    cpm_cell = np.zeros(G_full)
    for i in range(0, ci.size, 2000):
        m = dense(R, ci[i:i + 2000])
        s = m.sum(1)
        s[s == 0] = 1.0
        cpm_cell += (m * (1e6 / s[:, None])).sum(0)
    cpm_cell /= ci.size
    keep = cpm_cell > MIN_CPM                      # the DE universe (157/164 convention)
    genes = genes_full[keep]
    G = int(keep.sum())
    kidx = np.flatnonzero(keep)
    b_de = np.log1p(cpm_cell)                      # 157/164's ordering gate, NOT the build's
    print(f"harness: {ci.size} controls, {G}/{G_full} genes past the {MIN_CPM}-CPM DE gate; "
          f"build gate holds {int(gate.sum())} [{time.time() - t0:.0f}s]", flush=True)

    cs = np.empty((G, ci.size))
    for i in range(0, ci.size, 2000):
        m = dense(R, ci[i:i + 2000])
        s = m.sum(1)
        s[s == 0] = 1.0
        cs[:, i:i + m.shape[0]] = np.log1p(m[:, keep] * (1e6 / s[:, None])).T
    cs.sort(axis=1)
    cs = np.ascontiguousarray(cs)
    tie = ctrl_tiesums(cs)
    ref_lin = cpm_cell[keep]
    print(f"control block ready [{time.time() - t0:.0f}s]", flush=True)

    # The head-protection rule is a BUILD-TIME choice: it is frozen once, from the first
    # a15 emission, and reused at every later seed. Re-deriving it per seed would let the
    # arm see the very draw it is scored on.
    HEAD = None
    if os.path.exists(os.path.join(WORK, "head.npz")):
        HEAD = np.load(os.path.join(WORK, "head.npz"))["head"]
        print(f"head order loaded from work/head.npz {HEAD.shape}", flush=True)

    if NT:
        targets = targets[:NT]
    T = len(targets)
    names = [a[0] for a in ARMS]
    print(f"{T} perturbations x {len(ARMS)} arms x {len(SEEDS)} seeds: {names}", flush=True)

    real = {k: np.full((T, G), np.nan) for k in ("p", "q", "l")}
    real_tot = np.zeros((T, G_full))
    out = {(n, s): {k: np.full((T, G), np.nan) for k in ("p", "q", "l")}
           for n in names for s in SEEDS}
    svec = np.zeros((T, G), dtype=np.float32)
    prof = {n: np.zeros((T, G_full), dtype=np.float32) for n in names}
    npred = {(n, s): [] for n in names for s in SEEDS}
    diag = {n: [] for n in names}                  # ||eps||/||D||, per target

    for i, t in enumerate(targets):
        # ---- reference side --------------------------------------------------------------
        Braw = dense(R, np.flatnonzero(rlab == t))
        real_tot[i] = Braw.sum(0)
        s = Braw.sum(1)
        s[s == 0] = 1.0
        Br = np.log1p(Braw[:, keep] * (1e6 / s[:, None]))
        del Braw
        zz = wilcoxon_z(cs, tie, np.ascontiguousarray(Br.T))
        real["p"][i] = 2 * norm.sf(np.abs(zz))
        real["q"][i] = false_discovery_control(real["p"][i], method="bh")
        real["l"][i] = np.log2((np.expm1(Br).mean(0) + 1e-9) / (ref_lin + 1e-9))
        del Br

        # ---- the frozen a15 profile, and the spread vector held fixed across arms --------
        tot = tot_all[i]
        prof_a15 = (tot + 0.5) / NL
        Z = float(prof_a15.sum())
        pcpm = 1e6 * tot / max(tot.sum(), 1)
        eff_rank = np.abs(np.log2((pcpm + 1e-9) / (cpm_cell + 1e-9)))
        gated_de = (cpm_cell > MIN_CPM) & (b_de > 0)
        score = np.where(gated_de, eff_rank, -np.inf)
        order = np.lexsort((np.arange(G_full), -score))
        sv_tier = np.full(G_full, 1.0)
        is_keep = np.zeros(G_full, dtype=bool)
        k8 = min(8000, int(gated_de.sum()))
        sv_tier[order[:k8]] = 0.25
        is_keep[order[:k8]] = True
        sv_flat = np.full(G_full, 0.25)
        svec[i] = sv_tier[keep]

        # ---- the recovered effect, and the a15 profile it renders to ---------------------
        Di = D_all[i]
        lin_recon = np.clip(np.expm1(b + Di), 0.0, None)
        p_recon = lin_recon / lin_recon.sum()
        w = np.where(gate & (b + Di > 0.0), np.exp(b + Di), 0.0)

        head_p = HEAD[i] if HEAD is not None else None

        for name, uspec, law, frac in ARMS:
            if uspec is None:
                emit_prof = prof_a15
            else:
                kind, par = uspec
                if kind == "const_u":
                    u = np.ones(G_full)
                elif kind == "inv_w":
                    u = np.where(w > 0, 1.0 / np.where(w > 0, w, 1.0), 0.0)
                elif kind == "head_p":
                    assert head_p is not None, "a15 must be the first arm"
                    u = np.ones(G_full)
                    u[head_p[:par]] = 0.0
                elif kind == "head_w":
                    u = np.ones(G_full)
                    u[np.lexsort((np.arange(G_full), -w))[:par]] = 0.0
                elif kind == "keepblk":
                    u = is_keep.astype(np.float64)
                elif kind == "tailblk":
                    u = (~is_keep).astype(np.float64)
                elif kind == "inv_eff":
                    u = 1.0 / (np.abs(Di) + par)
                elif kind == "eff":
                    u = np.abs(Di) + par
                elif kind == "gate10":
                    u = (ctrl_cpm > par).astype(np.float64)
                else:
                    raise ValueError(kind)
                h = mass_fix_u(Di, b, gate, u, frac=frac)
                lin = np.clip(np.expm1(b + h), 0.0, None)
                p_arm = lin / lin.sum()
                emit_prof = prof_a15 + (p_arm - p_recon) * Z
                if emit_prof.min() < 0:
                    emit_prof = np.clip(emit_prof, 0.0, None)
                diag[name].append(float(np.linalg.norm(h - Di)
                                        / max(np.linalg.norm(Di), 1e-30)))
            prof[name][i] = emit_prof / max(float(emit_prof.sum()), 1e-30)

            sv = sv_tier if law == "tier" else sv_flat
            for seed in SEEDS:
                rng = np.random.default_rng(seed + i)
                vals, cols, ptr = emit_counts(
                    emit_prof, np.full(N_EMIT, EMIT_L), 0.1, rng,
                    mode="quasi", spread=sv, overdispersion=4.0,
                )
                B = np.zeros((N_EMIT, G_full))
                for r in range(N_EMIT):
                    sl = slice(ptr[r], ptr[r + 1])
                    B[r, cols[sl]] = vals[sl]
                if name == "a15":
                    assert np.array_equal(B.sum(0).astype(np.int64), tot), \
                        "ABORT: the a15 arm does not reproduce the frozen totals"
                sc = B.sum(1)
                sc[sc == 0] = 1.0
                nb = np.log1p(B[:, keep] * (1e6 / sc[:, None]))
                del B
                zz = wilcoxon_z(cs, tie, np.ascontiguousarray(nb.T))
                pv = 2 * norm.sf(np.abs(zz))
                qv = false_discovery_control(pv, method="bh")
                o = out[(name, seed)]
                o["p"][i], o["q"][i] = pv, qv
                o["l"][i] = np.log2((np.expm1(nb).mean(0) + 1e-9) / (ref_lin + 1e-9))
                npred[(name, seed)].append(int((qv < 0.05).sum()))
                del nb
                if name == "a15" and seed == SEEDS[0] and HEAD is None:
                    # the confidence order `reach` actually walks, on the full gene axis
                    head_p = kidx[np.lexsort((np.arange(G), pv))]
        if (i + 1) % 5 == 0:
            print(f"  {i + 1}/{T} [{time.time() - t0:.0f}s]", flush=True)

    os.makedirs(WORK, exist_ok=True)
    for name in names:
        for seed in SEEDS:
            o = out[(name, seed)]
            path = os.path.join(WORK, f"de.{name}.s{seed}.npz")
            np.savez_compressed(
                path, genes=genes, perts=np.asarray(targets),
                real_p=real["p"], real_padj=real["q"], real_lfc=real["l"],
                pred_p=o["p"], pred_padj=o["q"], pred_lfc=o["l"], spread=svec,
            )
        d = diag[name]
        print(f"wrote de.{name}.*  median n_pred "
              f"{np.median(npred[(name, SEEDS[0])]):.0f}/"
              f"{np.median(npred[(name, SEEDS[1])]):.0f}"
              + (f"   ||eps||/||D|| median {np.median(d):.5f}" if d else ""), flush=True)
    np.savez_compressed(os.path.join(WORK, os.environ.get("PROF_OUT", "profiles.npz")),
                        genes=genes_full, perts=np.asarray(targets),
                        b_ctrl=b, gate=gate, real_tot=real_tot, keep=keep,
                        **{f"p_{n}": prof[n] for n in names})
    print(f"wrote profiles.npz\ndone [{time.time() - t0:.0f}s]")


if __name__ == "__main__":
    main()
