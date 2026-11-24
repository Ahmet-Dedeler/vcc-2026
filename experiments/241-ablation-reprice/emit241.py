"""241 — re-emit 200's ablation arms on h1-test with `VCC_MASS_FIX` UNSET.

200's `realize200.py` carried `mass_fix=True` on EVERY arm except `art_nomf` (its arm tuple's
fourth field).  So all four DE members of all ten components were measured on a23's renderer,
which only a23 ever shipped and which scored -0.00073 before being abandoned.  This script is
`realize200.py` with that one field flipped and the output redirected here; nothing else about
the physics changes, and `h1parts200` / `emit_counts` / the production
`mass_consistent_effect` are imported from their owners rather than copied.

KNOWN-ANSWER ROW (rule 10, aborts).  `art_nomf` is the ONE arm 200 already emitted without the
mass fix.  Re-emitting it here must reproduce `experiments/200-.../de200.art_nomf.npz`
BIT-IDENTICALLY on `pred_p`, `pred_padj`, `pred_lfc` and `spread` — otherwise this code path is
not 200's and no delta measured with it is comparable to 200's table.

Donors: `k562_gwps` + `orion_hek293t` ONLY.  `data/donors/h1_2025` covers 85/85 of h1-test's
targets and is pure oracle leakage.

    ONLY=art_nomf uv run python emit241.py
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

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
E200 = os.path.join(ROOT, "experiments/200-production-ablation")
E169 = os.path.join(ROOT, "experiments/169-sharpening-six-member")
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/092-call-set-truncation")
sys.path.insert(0, f"{ROOT}/experiments/162-retrodicting-validator")
sys.path.insert(0, E200)
from vcc2026.donor_transfer import emit_counts                 # noqa: E402
from de_kernel import ctrl_tiesums, wilcoxon_z                 # noqa: E402
import h1parts200                                              # noqa: E402

os.makedirs(WORK, exist_ok=True)
REAL = f"{ROOT}/data/harness/h1-test.h5ad"
OUT = os.path.join(HERE, "de241")
CTRL, COL, SCORER_MIN_CPM = "non-targeting", "target_gene", 5.0
N_EMIT, EMIT_L, TS = 400, 100_000, 50_000.0
E_A15 = dict(tail=1.0, keep_top=8000, keep_spread=0.25, phi=4.0)     # the shipped law
E_FLAT = dict(tail=None, keep_top=0, keep_spread=0.25, phi=4.0)      # a8's law: -tier
E_PHI0 = dict(tail=1.0, keep_top=8000, keep_spread=0.25, phi=0.0)    # -phi
SEEDS = (90_000, 7_000_000, 11_111, 22_222, 33_333)

G = dict(
    ship=dict(pm=0.50, gamma=0.65, lam=0.0, cd4_f=0.10, cd4_frac=None, weights=None),
    noWhiten=dict(pm=0.50, gamma=0.00, lam=0.0, cd4_f=0.10, cd4_frac=None, weights=None),
    noPM=dict(pm=0.00, gamma=0.65, lam=0.0, cd4_f=0.10, cd4_frac=None, weights=None),
    noCD4=dict(pm=0.50, gamma=0.65, lam=0.0, cd4_f=None, cd4_frac=None, weights=None),
    w41=dict(pm=0.50, gamma=0.65, lam=0.0, cd4_f=0.10, cd4_frac=None,
             weights={"k562_gwps": 4.0, "orion_hek293t": 1.0}),
)


def arms_all():
    """(name, geom | None, emission, mass_fix, min_cpm, seed) — mass_fix is False EVERYWHERE."""
    A = [("art_nomf_KA", None, E_A15, False, 5.0, SEEDS[0])]     # rule-10 bit-identity row
    for i, s in enumerate(SEEDS):
        tag = "" if i == 0 else f"_s{i+1}"
        A.append((f"art_flat_nomf{tag}", None, E_FLAT, False, 5.0, s))
        A.append((f"art_phi0_nomf{tag}", None, E_PHI0, False, 5.0, s))
    for i, s in enumerate(SEEDS[:3]):
        tag = "" if i == 0 else f"_s{i+1}"
        for g in ("ship", "noWhiten", "noPM", "w41"):
            A.append((f"reb_{g}_nomf{tag}", g, E_A15, False, 5.0, s))
        for mc in (0.0, 10.0):
            A.append((f"reb_ship_nomf_gate{int(mc)}{tag}", "ship", E_A15, False, mc, s))
    # one seed is enough to re-confirm the structural zero (CD4 coverage 4/85 -> frac 0)
    A.append(("reb_noCD4_nomf", "noCD4", E_A15, False, 5.0, SEEDS[0]))
    return A


def dense(A, rows):
    m = A.X[rows]
    return np.asarray(m.todense() if hasattr(m, "todense") else m, dtype=np.float64)


def load_builder():
    spec = importlib.util.spec_from_file_location(
        "prod_builder_241", f"{ROOT}/scripts/make_submission_pm_cd4.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_241"] = mod
    spec.loader.exec_module(mod)
    assert mod.TS == TS
    assert mod.MASS_FIX is None, "the module default must be off (the whole point of 241)"
    return mod


def main() -> None:
    t0 = time.time()
    want = set(os.environ["ONLY"].split(",")) if os.environ.get("ONLY") else None
    ARMS = [a for a in arms_all() if want is None or a[0] in want]
    if not ARMS:
        sys.exit("no arms selected")
    ARMS = [a for a in ARMS if not os.path.exists(f"{OUT}.{a[0]}.npz")]
    if not ARMS:
        print("every requested arm already on disk")
        return
    assert all(a[3] is False for a in ARMS), "241 emits on the SHIPPED renderer only"

    E = np.load(os.path.join(E169, "work/h1test-effect.npz"), allow_pickle=True)
    H0 = E["H"].astype(np.float64)
    b = E["b"].astype(np.float64)
    a_gate = np.asarray(E["gate"])
    etargets = np.asarray([str(t) for t in E["targets"]])
    esums = E["sums"].astype(np.float64)
    egenes = np.asarray([str(g) for g in E["genes"]])
    art_norm = float(np.median(np.linalg.norm(H0, axis=1)))
    print(f"artifact effect: {H0.shape}, gate {int(a_gate.sum())} genes at "
          f"--min-cpm {float(E['min_cpm'])}, median row norm {art_norm:.4f}", flush=True)
    load_builder()      # asserts MASS_FIX is None; never called below because mf is False

    need_geom = sorted({a[1] for a in ARMS if a[1] is not None})
    GEO: dict[str, np.ndarray] = {}
    if need_geom:
        p = h1parts200.parts()
        c4 = h1parts200.cd4(p["targets"])
        assert list(p["genes"]) == list(egenes), "harness / effect gene axes differ"
        op = h1parts200.H1Operator()
        print(f"[h1-test] operator maxeig {op.evals.max():.4f} n {op.evals.size}", flush=True)
        order = np.asarray([list(p["targets"]).index(t) for t in etargets], dtype=np.int64)
        for g in need_geom:
            M = h1parts200.effects_for(G[g], p, c4, op, art_norm)[order]
            GEO[g] = M
            cs = ((M * H0).sum(1)
                  / np.maximum(np.linalg.norm(M, axis=1) * np.linalg.norm(H0, axis=1), 1e-300))
            print(f"  geometry {g:<10} median row cosine to the artifact's H "
                  f"{float(np.median(cs)):+.5f}  median norm "
                  f"{float(np.median(np.linalg.norm(M, axis=1))):.4f}", flush=True)
        del op

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
    keep = cpm > SCORER_MIN_CPM
    genes = genes_full[keep]
    Gk = int(keep.sum())
    print(f"harness: {ci.size} controls, {Gk}/{G_full} genes past the "
          f"{SCORER_MIN_CPM}-CPM gate [{time.time()-t0:.0f}s]", flush=True)

    cs_block = np.empty((Gk, ci.size))
    for i in range(0, ci.size, 2000):
        m = dense(R, ci[i:i + 2000])
        s = m.sum(1)
        s[s == 0] = 1.0
        cs_block[:, i:i + m.shape[0]] = np.log1p(m[:, keep] * (1e6 / s[:, None])).T
    cs_block.sort(axis=1)
    cs_block = np.ascontiguousarray(cs_block)
    tie = ctrl_tiesums(cs_block)
    ref_lin = cpm[keep]
    b_cell = np.log1p(cpm)
    print(f"control block ready [{time.time()-t0:.0f}s]", flush=True)

    epos = {t: i for i, t in enumerate(etargets)}
    targets = [t for t in sorted(set(rlab.tolist()) - {CTRL}) if t in epos]
    T = len(targets)
    print(f"{T} perturbations, {len(ARMS)} arms", flush=True)

    real = {k: np.full((T, Gk), np.nan) for k in ("p", "q", "l")}
    out = {n: {k: np.full((T, Gk), np.nan) for k in ("p", "q", "l")} for n, *_ in ARMS}
    svec = {n: np.zeros((T, Gk), dtype=np.float32) for n, *_ in ARMS}
    npred = {n: [] for n, *_ in ARMS}

    for i, t in enumerate(targets):
        Br = dense(R, np.flatnonzero(rlab == t))
        s = Br.sum(1)
        s[s == 0] = 1.0
        Br = np.log1p(Br[:, keep] * (1e6 / s[:, None]))
        z = wilcoxon_z(cs_block, tie, np.ascontiguousarray(Br.T))
        real["p"][i] = 2 * norm.sf(np.abs(z))
        real["q"][i] = false_discovery_control(real["p"][i], method="bh")
        real["l"][i] = np.log2((np.expm1(Br).mean(0) + 1e-9) / (ref_lin + 1e-9))
        del Br

        k = epos[t]
        S = float(esums[k].sum())

        for name, geom, em, mf, min_cpm, seed in ARMS:
            h = H0[k] if geom is None else GEO[geom][k]
            g = a_gate if (geom is None and min_cpm == 5.0) else ((cpm > min_cpm) & (b > 0))
            hq = np.where(g, h, 0.0)
            # mass fix deliberately NOT applied: this is the shipped renderer
            lin = np.clip(np.expm1(b + hq), 0.0, None)
            profile = lin / lin.sum() * (S / (N_EMIT * EMIT_L))
            totals = np.floor(N_EMIT * EMIT_L * profile).astype(np.int64)
            pcpm = 1e6 * totals / max(totals.sum(), 1)
            eff = np.abs(np.log2((pcpm + 1e-9) / (cpm + 1e-9)))
            gated = (cpm > SCORER_MIN_CPM) & (b_cell > 0)
            score = np.where(gated, eff, -np.inf)
            sv = np.full(G_full, float(em["tail"] if em["tail"] is not None
                                       else em["keep_spread"]))
            if em["tail"] is not None:
                order_ = np.lexsort((np.arange(G_full), -score))
                kt = min(int(em["keep_top"]), int(gated.sum()))
                sv[order_[:kt]] = float(em["keep_spread"])

            rng = np.random.default_rng(seed + i)
            vals, cols, ptr = emit_counts(profile, np.full(N_EMIT, EMIT_L), 0.1, rng,
                                          mode="quasi", spread=sv,
                                          overdispersion=float(em["phi"]))
            B = np.zeros((N_EMIT, G_full))
            for r in range(N_EMIT):
                sl = slice(ptr[r], ptr[r + 1])
                B[r, cols[sl]] = vals[sl]
            assert np.array_equal(B.sum(0).astype(np.int64), totals), name
            ss = B.sum(1)
            ss[ss == 0] = 1.0
            nb = np.log1p(B[:, keep] * (1e6 / ss[:, None]))
            del B
            z = wilcoxon_z(cs_block, tie, np.ascontiguousarray(nb.T))
            pv = 2 * norm.sf(np.abs(z))
            qv = false_discovery_control(pv, method="bh")
            out[name]["p"][i], out[name]["q"][i] = pv, qv
            out[name]["l"][i] = np.log2((np.expm1(nb).mean(0) + 1e-9) / (ref_lin + 1e-9))
            svec[name][i] = sv[keep]
            npred[name].append(int((qv < 0.05).sum()))
            del nb
        if (i + 1) % 10 == 0:
            print(f"  {i+1}/{T} [{time.time()-t0:.0f}s]", flush=True)

    meta = {}
    for name, *_ in ARMS:
        path = f"{OUT}.{name}.npz"
        np.savez_compressed(
            path, genes=genes, perts=np.array(targets),
            real_p=real["p"], real_padj=real["q"], real_lfc=real["l"],
            pred_p=out[name]["p"], pred_padj=out[name]["q"], pred_lfc=out[name]["l"],
            spread=svec[name])
        meta[name] = float(np.median(npred[name]))
        print(f"wrote {os.path.basename(path)}   median n_pred {meta[name]:.0f}", flush=True)
    json.dump(meta, open(os.path.join(WORK, "emit241_npred.json"), "w"), indent=1)
    print(f"done [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
