"""313 -- run BOTH implementations of the per-target normrel weighting on ONE shared fixture.

  PYTHONDONTWRITEBYTECODE=1 uv run python parity313.py

Three layers, so that an ordering difference is never confused with a dtype difference:

  L1  PRODUCTION, the real functions: `make_submission_pm_cd4.build_v9_pre_whiten_effects`
      (-> `_normrel_pre_whiten_effects`) then `apply_perturbmatch_k562_replacement`, called on a real
      `MultiDonorModel` built from the fixture's synthetic donors.
  L2  308, the real functions: `lib308.normrel_weights` then `lib263.combine(..., norm_match=True)`,
      called on the fixture's prepared rows (panel columns zeroed, as `062.prepare_main_source` zeroes
      them) with PerturbMatch blended into the k562 row first, as `lib263.prepared_sources` does.
  L3  a clean-room float64 re-statement of each ORDERING on the same substrate and the same weights, so
      the ordering term can be measured to machine precision without the float32 casts both real
      implementations make.

Read-only outside `experiments/313-normrel-parity/`.  No bank, no build, no submission.
"""
from __future__ import annotations

import hashlib
import json
import resource
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "experiments/308-kd-weighted-donors"))
sys.path.insert(0, str(HERE))

import fixture313 as F  # noqa: E402

EXP = -1.5
LINES: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ---------------------------------------------------------------------------------- the shared substrate
def load_fixture() -> dict:
    z = np.load(WORK / "fixture313.npz", allow_pickle=False)
    return {k: z[k] for k in z.files}


def build_model(fx: dict):
    from vcc2026.donor_transfer import Donor
    from vcc2026.multi_donor import MultiDonorModel

    donors = {}
    for name in F.DONORS:
        donors[name] = Donor(
            delta=np.asarray(fx[f"delta__{name}"], dtype=np.float64),
            perts=np.asarray([str(p) for p in fx[f"perts__{name}"]]),
            genes=np.asarray([str(g) for g in fx[f"genes__{name}"]]),
            b_ctrl=np.asarray(fx[f"b_ctrl__{name}"], dtype=np.float64),
        )
    official = np.asarray([str(g) for g in fx["official"]])
    return MultiDonorModel(donors, official, gene_norm=False, eb=False)


def build_pm(fx: dict, official: list[str]) -> dict:
    """The mapping `EXP061.prepare_perturbmatch` returns, computed here on the fixture's PM matrix with the
    same recipe (column mean over every PM perturbation, median centred row norm as the scale)."""
    pm_genes = [str(g) for g in fx["pm_genes"]]
    pm_perts = [str(p) for p in fx["pm_perts"]]
    delta = np.asarray(fx["pm_delta"], dtype=np.float64)
    tindex = {g: i for i, g in enumerate(official)}
    src, tgt = [], []
    for j, g in enumerate(pm_genes):
        if g in tindex:
            src.append(j)
            tgt.append(tindex[g])
    src = np.asarray(src, dtype=np.int64)
    tgt = np.asarray(tgt, dtype=np.int64)
    mean = delta[:, src].mean(axis=0)
    scale = float(np.median(np.linalg.norm(delta[:, src] - mean, axis=1))) or 1.0
    return dict(genes=pm_genes, perts=pm_perts, pert_index={p: i for i, p in enumerate(pm_perts)},
                target_columns=tgt, source_columns=src, mean=mean, scale=scale, delta=delta)


def harness_rows(model, targets: list[str], panel: np.ndarray) -> tuple[list[np.ndarray], np.ndarray]:
    """Each donor's PREPARED row in the harness's units: `_global * _renorm[j] * v_j` on the official axis
    with the panel columns zeroed (`062.prepare_main_source:287`).  With uniform weights their mean is exactly
    the pre-whiten effect matrix divided by alpha, which KA-P1 checks against production itself."""
    G = len(model.target_genes)
    rows, avail = [], []
    for name in model.names:
        R = np.zeros((len(targets), G))
        a = np.zeros(len(targets), bool)
        for i, t in enumerate(targets):
            r = model.pert_index[name].get(t)
            if r is None:
                continue
            v = np.zeros(G)
            v[model.tcols[name]] = model.s[name][r][model.dcols[name]]
            R[i] = model._global * model._renorm[name] * v
            a[i] = True
        R[:, panel] = 0.0
        rows.append(R)
        avail.append(a)
    return rows, np.stack(avail)


def harness_pm_row(model, pm: dict, targets: list[str], panel: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """PerturbMatch as a prepared donor row in the same units: `effect_scale * (delta - mean)/scale`.
    (Production adds `alpha * k562_share * effect_scale * blend * (pm_std - incumbent)`, and the k562 slot's
    own contribution is `alpha * effect_scale/k_scale * v_k`, so `effect_scale * pm_std` is exactly the PM row
    on the same footing as `r_k562`.)"""
    G = len(model.target_genes)
    R = np.zeros((len(targets), G))
    a = np.zeros(len(targets), bool)
    for i, t in enumerate(targets):
        r = pm["pert_index"].get(t)
        if r is None:
            continue
        std = (np.asarray(pm["delta"][r, pm["source_columns"]], dtype=np.float64) - pm["mean"]) / pm["scale"]
        R[i, pm["target_columns"]] = float(model.effect_scale) * std
        a[i] = True
    R[:, panel] = 0.0
    return R, a


def f32(A: np.ndarray) -> np.ndarray:
    return np.asarray(A, dtype=np.float64).astype(np.float32).astype(np.float64)


# ---------------------------------------------------------------------------- clean-room float64 orderings
def uniform_weights(avail: np.ndarray) -> np.ndarray:
    J, n = avail.shape
    W = np.zeros((n, J))
    for t in range(n):
        js = np.flatnonzero(avail[:, t])
        if js.size:
            W[t, js] = 1.0 / js.size
    return W


def normrel_W(rows: list[np.ndarray], avail: np.ndarray, cols: np.ndarray, power: float) -> np.ndarray:
    """The 308 weight rule, re-stated in float64 (identical algebra to `lib308.normrel_weights`)."""
    J, n = avail.shape
    K = np.zeros((J, n))
    for j in range(J):
        nj = np.linalg.norm(rows[j][:, cols], axis=1)
        med = float(np.median(nj[avail[j]])) or 1.0
        K[j] = nj / med
    W = np.full((n, J), np.nan)
    for t in np.flatnonzero(avail.sum(0) >= 2):
        js = np.flatnonzero(avail[:, t])
        k = np.maximum(K[js, t], 1e-9) ** power
        W[t] = 0.0
        W[t, js] = k / k.sum()
    return W


def normrel_W_production_rule(model, targets: list[str], power: float) -> np.ndarray:
    """The PRODUCTION weight rule, re-stated: the norm is taken on the donor's OWN official-axis coverage --
    which INCLUDES the panel columns, i.e. the target's own knockdown -- of the raw `s` row, and a donor whose
    row norm or panel median is exactly zero falls back to the uniform weight."""
    J, n = len(model.names), len(targets)
    med = {}
    for name in model.names:
        vals = [float(np.linalg.norm(model.s[name][model.pert_index[name][t]][model.dcols[name]]))
                for t in targets if t in model.pert_index[name]]
        med[name] = float(np.median(vals)) if vals else 0.0
    W = np.full((n, J), np.nan)
    for i, t in enumerate(targets):
        have = [j for j, name in enumerate(model.names) if t in model.pert_index[name]]
        if not have:
            continue
        w = np.zeros(J)
        for j in have:
            name = model.names[j]
            nv = float(np.linalg.norm(model.s[name][model.pert_index[name][t]][model.dcols[name]]))
            w[j] = (nv / med[name]) ** power if (med[name] > 0.0 and nv > 0.0) else 1.0
        s = w.sum()
        if s <= 0:
            continue
        W[i] = w / s
    return W


def normrel_W_production_rule_rows(rows: list[np.ndarray], avail: np.ndarray, cols: np.ndarray,
                                   power: float) -> np.ndarray:
    """The production weight RULE (uniform fallback when the norm or the median is exactly zero) expressed on
    prepared rows, so it can be applied to the PM-blended rows a parity-correct arm would weight."""
    J, n = avail.shape
    med = np.zeros(J)
    for j in range(J):
        nj = np.linalg.norm(rows[j][:, cols], axis=1)
        med[j] = float(np.median(nj[avail[j]])) if avail[j].any() else 0.0
    W = np.full((n, J), np.nan)
    for t in range(n):
        js = np.flatnonzero(avail[:, t])
        if js.size == 0:
            continue
        w = np.zeros(J)
        for j in js:
            nv = float(np.linalg.norm(rows[j][t][cols]))
            w[j] = (nv / med[j]) ** power if (med[j] > 0.0 and nv > 0.0) else 1.0
        if w.sum() <= 0:
            continue
        W[t] = w / w.sum()
    return W


def ord_A(rows_blend: list[np.ndarray], avail: np.ndarray, W: np.ndarray, pre_blend: np.ndarray) -> np.ndarray:
    """308 / 263: PerturbMatch is already inside `rows_blend`, so the weighting and the norm match both see the
    blended k562 row.  `pre_blend` is the uniform mean of the same blended rows."""
    J, n = avail.shape
    out = pre_blend.copy()
    for t in range(n):
        if not np.isfinite(W[t]).all():
            continue
        js = np.flatnonzero(avail[:, t])
        if js.size == 0:
            continue
        w = W[t, js] / W[t, js].sum()
        acc = np.zeros(rows_blend[0].shape[1])
        for k, j in enumerate(js):
            acc += w[k] * rows_blend[j][t]
        na, nm = np.linalg.norm(acc), np.linalg.norm(pre_blend[t])
        if na > 0 and nm > 0:
            acc = acc * (nm / na)
        out[t] = acc
    return out


def ord_B(rows_raw: list[np.ndarray], avail: np.ndarray, W: np.ndarray, pre_raw: np.ndarray,
          pm_row: np.ndarray, pm_avail: np.ndarray, pm_cols: np.ndarray, k_index: int,
          blend: float) -> np.ndarray:
    """Production: weight and norm-match FIRST against the PM-free uniform row, then add the PerturbMatch
    correction at the UNIFORM k562 share, outside the norm match."""
    J, n = avail.shape
    out = pre_raw.copy()
    for t in range(n):
        js = np.flatnonzero(avail[:, t])
        if np.isfinite(W[t]).all() and js.size:
            w = W[t, js] / W[t, js].sum()
            acc = np.zeros(rows_raw[0].shape[1])
            for k, j in enumerate(js):
                acc += w[k] * rows_raw[j][t]
            nb, nc = np.linalg.norm(pre_raw[t]), np.linalg.norm(acc)
            if nb > 0 and nc > 0:
                acc = acc * (nb / nc)
            out[t] = acc
        if avail[k_index, t] and pm_avail[t]:
            share = 1.0 / float(js.size)          # the UNIFORM share, as production computes it
            out[t, pm_cols] += share * blend * (pm_row[t, pm_cols] - rows_raw[k_index][t, pm_cols])
    return out


# ------------------------------------------------------------------------------------------------- main
def main() -> None:
    fx = load_fixture()
    fx_sha = sha(WORK / "fixture313.npz")
    targets = [str(t) for t in fx["targets"]]
    official = [str(g) for g in fx["official"]]
    panel = np.asarray(fx["panel_positions"], dtype=np.int64)
    n, G = len(targets), len(official)

    model = build_model(fx)
    pm = build_pm(fx, official)

    import make_submission_pm_cd4 as MB  # noqa: E402  (imports in 1.2 s / 0.21 GiB, reads no data)

    ALPHA, BETA, BLEND = float(MB.ALPHA), float(MB.BETA), float(MB.PM_BLEND)
    uni = MB.donor_weights(model)
    k_index = model.names.index(MB.K562)
    P(f"[313] fixture sha256 {fx_sha} seed {F.SEED}; {n} targets x {G} genes x {len(model.names)} donors")
    P(f"[313] production constants: ALPHA {ALPHA} BETA {BETA} PM_BLEND {BLEND} donor_weights {uni}")
    P(f"[313] builder sha256 {sha(ROOT / 'scripts/make_submission_pm_cd4.py')[:16]}  "
      f"lib308 {sha(ROOT / 'experiments/308-kd-weighted-donors/lib308.py')[:16]}  "
      f"lib263 {sha(ROOT / 'experiments/263-donor-selection-arbiters/lib263.py')[:16]}")

    # ------------------------------------------------------------------ the two substrates
    #   the HARNESS substrate (308/263): the panel columns are zeroed, `062.prepare_main_source:287`;
    #   the PRODUCTION substrate: `model.s` keeps them, so the target's own knockdown is still in the row.
    rows_raw, avail = harness_rows(model, targets, panel)
    pm_row, pm_avail = harness_pm_row(model, pm, targets, panel)
    empty = np.asarray([], dtype=np.int64)
    rows_prod, _ = harness_rows(model, targets, empty)
    pm_row_prod, _ = harness_pm_row(model, pm, targets, empty)
    wcols = np.asarray([c for c in range(G) if c not in set(panel.tolist())], dtype=np.int64)
    pm_cols = np.asarray(pm["target_columns"], dtype=np.int64)
    shared_cols = np.asarray(sorted(set(pm_cols.tolist()) & set(wcols.tolist())
                                    & set(model.tcols[MB.K562].tolist())), dtype=np.int64)
    P(f"[313] wcols {wcols.size} (non-panel), PM columns {pm_cols.size} "
      f"(of which {len(set(pm_cols.tolist()) & set(panel.tolist()))} are PANEL columns), "
      f"263 `shared` k562∩PM∩non-panel {shared_cols.size}")
    P(f"[313] carrying donors per target: " + ", ".join(f"{t}:{int(avail[:, i].sum())}" for i, t in enumerate(targets)))

    # 263's blended k562 row and its uniform mean
    rows_blend = [R.copy() for R in rows_raw]
    both = avail[k_index] & pm_avail
    br = np.flatnonzero(both)
    rows_blend[k_index][np.ix_(br, shared_cols)] = (
        (1 - BLEND) * rows_raw[k_index][np.ix_(br, shared_cols)]
        + BLEND * pm_row[np.ix_(br, shared_cols)])

    def uniform_mean(rows):
        acc = np.zeros((n, G))
        cnt = avail.sum(0)
        for j in range(len(rows)):
            acc[avail[j]] += rows[j][avail[j]]
        cov = cnt > 0
        acc[cov] /= cnt[cov, None]
        return acc

    pre_blend64 = uniform_mean(rows_blend)
    pre_raw64 = uniform_mean(rows_raw)
    pre_prod64 = uniform_mean(rows_prod)               # production's norm-match reference `base` / alpha
    pre_blend = f32(pre_blend64)                       # the cached `pre` the harness norm-matches to
    P(f"[313] PerturbMatch blended into the k562 row for {int(both.sum())}/{n} targets")

    KA: list[dict] = []

    def ka(label, got, want, tol):
        d = abs(float(got) - float(want))
        ok = bool(d <= tol)
        KA.append(dict(label=label, got=float(got), want=float(want), abs_diff=d, tol=tol, ok=ok))
        P(f"  KA {'PASS' if ok else 'FAIL'} {label}: |d| {d:.3e} (tol {tol:.0e})")
        return ok

    # ------------------------------------------------------------------ KA-P1: the substrates agree
    MB.NORMREL_EXP = None
    prod_uniform = MB.build_v9_pre_whiten_effects(model, targets, alpha=ALPHA, beta=BETA, weights=uni)
    ka("KA-P1 production uniform pre-whiten == ALPHA * harness uniform mean (non-panel cols; float32 floor)",
       np.abs(prod_uniform[:, wcols] - ALPHA * pre_raw64[:, wcols]).max(), 0.0, 1e-6)
    f32_floor = float(np.abs(prod_uniform[:, wcols] - ALPHA * pre_raw64[:, wcols]).max())
    P(f"[313] measured float32 noise floor between the two substrates: {f32_floor:.3e} "
      f"(row scale {float(np.median(np.linalg.norm(prod_uniform[:, wcols], axis=1))):.4f})")

    # ------------------------------------------------------------------ KA-P2: exponent 0 == uniform, exactly
    MB.NORMREL_EXP = 0.0
    prod_exp0 = MB.build_v9_pre_whiten_effects(model, targets, alpha=ALPHA, beta=BETA, weights=uni)
    ka("KA-1a production exponent 0 == production uniform (bitwise)",
       np.abs(prod_exp0 - prod_uniform).max(), 0.0, 0.0)
    import lib308 as Z  # noqa: E402
    Y = Z.Y
    st = dict(pre=pre_blend, rows=rows_blend, avail=avail)
    W0 = Z.normrel_weights(rows_blend, avail, wcols, 0.0)
    a_exp0 = Y.combine(st, W0, norm_match=True)
    ka("KA-1b 308 exponent 0 == the cached uniform row (float32 combine floor)",
       np.abs(a_exp0 - pre_blend).max(), 0.0, 1e-6)
    ka("KA-1c 308 exponent-0 weights are uniform", np.abs(W0[np.isfinite(W0).all(1)]
       - uniform_weights(avail)[np.isfinite(W0).all(1)]).max(), 0.0, 1e-15)

    # ------------------------------------------------------------------ the two weight rules
    MB.NORMREL_EXP = EXP
    W_308 = Z.normrel_weights(rows_blend, avail, wcols, EXP)             # 308: blended rows, wcols only
    W_prod_rule = normrel_W_production_rule(model, targets, EXP)          # production: raw rows, full coverage
    W_clean = normrel_W(rows_blend, avail, wcols, EXP)
    ka("KA-2a clean-room 308 weight rule == lib308.normrel_weights",
       np.nanmax(np.abs(W_clean - W_308)) if np.isfinite(W_308).any() else 0.0, 0.0, 1e-15)

    # recover production's weights FROM ITS OWN OUTPUT (least squares on the donor basis)
    prod_pre = MB.build_v9_pre_whiten_effects(model, targets, alpha=ALPHA, beta=BETA, weights=uni)
    W_prod_recovered = np.full((n, len(model.names)), np.nan)
    lam_recovered = np.full(n, np.nan)
    for t in range(n):
        js = np.flatnonzero(avail[:, t])
        if js.size < 2:
            continue
        B = np.stack([rows_raw[j][t][wcols] for j in js])                 # (|js|, |wcols|)
        x, *_ = np.linalg.lstsq(B.T, prod_pre[t][wcols] / ALPHA, rcond=None)
        s = float(x.sum())
        lam_recovered[t] = s
        if s != 0:
            W_prod_recovered[t, js] = x / s
    live = np.isfinite(W_prod_rule).all(1) & np.isfinite(W_prod_recovered).all(1) & (avail.sum(0) >= 2)
    degen = np.zeros(n, bool)
    for t in range(n):                       # a combination whose norm collapses cannot be recovered by lstsq
        js = np.flatnonzero(avail[:, t])
        if js.size < 2 or not np.isfinite(W_prod_rule[t]).all():
            continue
        acc = sum(W_prod_rule[t, j] * rows_prod[j][t] for j in js)
        degen[t] = np.linalg.norm(acc) < 1e-8 * np.linalg.norm(pre_prod64[t])
    P(f"[313] degenerate combinations (norm collapse, lstsq recovery undefined): "
      f"{[targets[t] for t in np.flatnonzero(degen)]}")
    for t in np.flatnonzero(live & ~degen):
        d = float(np.abs(W_prod_recovered[t] - W_prod_rule[t]).max())
        if d > 1e-5:
            P(f"    [recovery mismatch] {targets[t]}: recovered {np.round(W_prod_recovered[t], 6).tolist()} "
              f"rule {np.round(W_prod_rule[t], 6).tolist()}")
    ka("KA-2b production weights recovered from its own output == the production weight rule",
       np.abs(W_prod_recovered[live & ~degen] - W_prod_rule[live & ~degen]).max(), 0.0, 1e-5)

    # ------------------------------------------------------------------ L1 vs L2, the real implementations
    prod_final, pm_audit = MB.apply_perturbmatch_k562_replacement(
        prod_pre, model, pm, targets, alpha=ALPHA, weights=uni)
    a_final = ALPHA * Y.combine(st, W_308, norm_match=True)

    # the INCUMBENT rows both arms sit on: uniform weights, PerturbMatch applied the same way each side
    prod_incumbent, _ = MB.apply_perturbmatch_k562_replacement(
        prod_uniform, model, pm, targets, alpha=ALPHA, weights=uni)
    a_incumbent = ALPHA * pre_blend
    # is each arm SHAPE-ONLY?  (263 `:nm` and AGENTS 8.1d: the shape-only gate is +0.0019, the
    # amplitude-changing gate +0.0089, so this decides which gate the arm has to clear.)
    dn_308 = np.linalg.norm(a_final, axis=1) / np.maximum(np.linalg.norm(a_incumbent, axis=1), 1e-300)
    dn_prod = np.linalg.norm(prod_final, axis=1) / np.maximum(np.linalg.norm(prod_incumbent, axis=1), 1e-300)
    dn_308_w = np.linalg.norm(a_final[:, wcols], axis=1) / np.maximum(
        np.linalg.norm(a_incumbent[:, wcols], axis=1), 1e-300)
    dn_prod_w = np.linalg.norm(prod_final[:, wcols], axis=1) / np.maximum(
        np.linalg.norm(prod_incumbent[:, wcols], axis=1), 1e-300)
    nd = ~np.isin(np.arange(n), [targets.index("T005")])          # the degenerate row excluded
    # the divergence that exists BEFORE 312: uniform weights on both sides, i.e. the incumbent itself.
    d_uniform = float(np.abs(prod_incumbent[:, wcols] - a_incumbent[:, wcols]).max())
    B_unif_shared = ALPHA * ord_B(rows_prod, avail, uniform_weights(avail), pre_prod64, pm_row_prod,
                                  pm_avail, shared_cols, k_index, BLEND)
    d_uniform_sharedpolicy = float(np.abs(B_unif_shared[:, wcols] - a_incumbent[:, wcols]).max())
    P("")
    P(f"[313] PRE-312 baseline: production incumbent (uniform + PM) vs the 308 harness's cached base, "
      f"non-panel columns: max|diff| {d_uniform:.4e}; with production's PM column policy replaced by "
      f"263's k562∩PM∩non-panel: {d_uniform_sharedpolicy:.4e} (float32 floor {f32_floor:.2e})")
    P("")
    P(f"[313] DELIVERED ROW NORM, arm / incumbent (the `nmae` axis; 1.0 == shape-only)")
    P(f"[313]   308  full axis: median {np.median(dn_308):.6f} range [{dn_308.min():.6f}, {dn_308.max():.6f}]; "
      f"non-panel: median {np.median(dn_308_w):.6f} range [{dn_308_w.min():.6f}, {dn_308_w.max():.6f}]")
    P(f"[313]   prod full axis: median {np.median(dn_prod):.6f} range [{dn_prod.min():.6f}, {dn_prod.max():.6f}]; "
      f"non-panel: median {np.median(dn_prod_w):.6f} range [{dn_prod_w.min():.6f}, {dn_prod_w.max():.6f}]")
    P(f"[313]   prod excl. the degenerate row: median {np.median(dn_prod[nd]):.6f} "
      f"range [{dn_prod[nd].min():.6f}, {dn_prod[nd].max():.6f}]; "
      f"panel-total norm ratio prod {np.linalg.norm(prod_final) / np.linalg.norm(prod_incumbent):.6f} "
      f"vs 308 {np.linalg.norm(a_final) / np.linalg.norm(a_incumbent):.6f}")

    D = prod_final[:, wcols] - a_final[:, wcols]
    max_abs = np.abs(D).max(1)
    nrm_p = np.linalg.norm(prod_final[:, wcols], axis=1)
    nrm_a = np.linalg.norm(a_final[:, wcols], axis=1)
    rel_norm = np.where(nrm_a > 0, (nrm_p - nrm_a) / np.maximum(nrm_a, 1e-300), np.nan)
    cosine = np.array([
        float(prod_final[t, wcols] @ a_final[t, wcols] /
              max(nrm_p[t] * nrm_a[t], 1e-300)) for t in range(n)])

    # ------------------------------------------------------------------ KA-3: PM off, weights matched
    #   With PerturbMatch absent AND the same weights, the two orderings are the SAME map.  Proved on the
    #   clean-room float64 layer (the real ones carry float32 casts at different points).
    zero_pm = np.zeros_like(pm_row)
    A_nopm = ord_A(rows_raw, avail, W_clean, pre_raw64)
    B_nopm = ord_B(rows_raw, avail, W_clean, pre_raw64, zero_pm, np.zeros(n, bool), pm_cols, k_index, BLEND)
    ka("KA-3 PerturbMatch ABSENT, identical weights: the two orderings coincide (clean float64)",
       np.abs(A_nopm - B_nopm).max(), 0.0, 1e-12)

    # and with PerturbMatch ON, same weights, the ordering term alone:
    A_pm = ord_A(rows_blend, avail, W_clean, pre_blend64)
    B_pm = ord_B(rows_raw, avail, W_clean, pre_raw64, pm_row, pm_avail, shared_cols, k_index, BLEND)
    ordering_only = np.abs(A_pm[:, wcols] - B_pm[:, wcols]).max(1)
    B_prod_clean = ord_B(rows_prod, avail, W_prod_rule, pre_prod64, pm_row_prod, pm_avail, pm_cols,
                         k_index, BLEND)
    ka("KA-4 clean-room ORD_B reproduces production's output (float32 floor)",
       np.abs(ALPHA * B_prod_clean[:, wcols] - prod_final[:, wcols]).max(), 0.0, 1e-4)
    ka("KA-5 clean-room ORD_A reproduces 308's output (float32 floor)",
       np.abs(ALPHA * A_pm[:, wcols] - a_final[:, wcols]).max(), 0.0, 1e-5)

    # ------------------------------------------------------------------ k562_share, demonstrated
    share_recovered = np.full(n, np.nan)
    for t in range(n):
        if not (avail[k_index, t] and pm_avail[t]):
            continue
        std = (np.asarray(pm["delta"][pm["pert_index"][targets[t]], pm["source_columns"]], dtype=np.float64)
               - pm["mean"]) / pm["scale"]
        incumbent = np.zeros(G)
        incumbent[model.tcols[MB.K562]] = (model.s[MB.K562][model.pert_index[MB.K562][targets[t]]][model.dcols[MB.K562]]
                                           / float(model._scale[MB.K562]))
        drive = ALPHA * float(model.effect_scale) * BLEND * (std - incumbent[pm_cols])
        applied = prod_final[t, pm_cols] - prod_pre[t, pm_cols]
        ok = np.abs(drive) > 1e-9
        share_recovered[t] = float(np.median(applied[ok] / drive[ok])) if ok.any() else np.nan

    # ------------------------------------------------------------------ the case table
    CASES = {
        "T003": "(a) PM ABSENT, 3 donors", "T019": "(a) PM ABSENT, 3 donors",
        "T005": "(d) near-zero donor row (hct116)", "T006": "(d) small donor row 0.05x (hek293t)",
        "T007": "(e) donor row 8x its own median (k562)",
        "T014": "(c) missing from ONE donor", "T015": "(c) missing from ONE donor",
        "T016": "(c) missing from TWO donors", "T017": "(c) missing from TWO donors",
        "T018": "(c) missing from k562 -- PM branch cannot fire",
    }
    rows_out = []
    for t, name in enumerate(targets):
        case = CASES.get(name, "(a) PM present, 3 donors")
        wa = W_308[t] if np.isfinite(W_308[t]).all() else np.full(len(model.names), np.nan)
        wb = W_prod_rule[t] if np.isfinite(W_prod_rule[t]).all() else np.full(len(model.names), np.nan)
        dw = float(np.nanmax(np.abs(wa - wb))) if np.isfinite(wa).any() and np.isfinite(wb).any() else float("nan")
        rows_out.append(dict(
            target=name, case=case, n_donors=int(avail[:, t].sum()),
            pm=bool(pm_avail[t] and avail[k_index, t]),
            max_abs_diff=float(max_abs[t]), rel_norm_diff=float(rel_norm[t]), cosine=float(cosine[t]),
            ordering_only_max_abs=float(ordering_only[t]),
            w_308=[None if not np.isfinite(x) else float(x) for x in wa],
            w_prod=[None if not np.isfinite(x) else float(x) for x in wb],
            max_abs_weight_diff=dw,
            k562_share_used_by_production=None if not np.isfinite(share_recovered[t]) else float(share_recovered[t]),
            k562_normrel_weight_308=None if not np.isfinite(wa[k_index]) else float(wa[k_index]),
            k562_normrel_weight_prod=None if not np.isfinite(wb[k_index]) else float(wb[k_index]),
        ))

    P("")
    P(f"{'target':7s} {'case':40s} {'nD':>2s} {'PM':>3s} {'max|diff|':>11s} {'relnorm':>9s} {'cos':>8s} "
      f"{'ordering':>11s} {'max|dW|':>8s} {'w_k 308':>8s} {'w_k prod':>8s} {'share':>7s}")
    for r in rows_out:
        P(f"{r['target']:7s} {r['case']:40s} {r['n_donors']:2d} {'Y' if r['pm'] else '-':>3s} "
          f"{r['max_abs_diff']:11.4e} {r['rel_norm_diff']:9.2e} {r['cosine']:8.5f} "
          f"{r['ordering_only_max_abs']:11.4e} {r['max_abs_weight_diff']:8.2e} "
          f"{(r['k562_normrel_weight_308'] if r['k562_normrel_weight_308'] is not None else float('nan')):8.4f} "
          f"{(r['k562_normrel_weight_prod'] if r['k562_normrel_weight_prod'] is not None else float('nan')):8.4f} "
          f"{(r['k562_share_used_by_production'] if r['k562_share_used_by_production'] is not None else float('nan')):7.4f}")

    # ------------------------------------------------------------------ attribution of the total divergence
    def maxdiff(A, B):
        return float(np.abs(A[:, wcols] - B[:, wcols]).max())

    # D1+D2 only: same weights, same substrate, only the PM ordering differs
    d_order = maxdiff(A_pm, B_pm)
    # D3+D5 only: same ordering (production's, on production's substrate), 308's weights vs production's
    B_w308 = ord_B(rows_prod, avail, W_clean, pre_prod64, pm_row_prod, pm_avail, pm_cols, k_index, BLEND)
    B_wprod = B_prod_clean
    d_weight = maxdiff(B_w308, B_wprod)
    # D5 alone: 308's weight rule but with the norm taken on the production column set (panel included)
    W_308_fullaxis = normrel_W(rows_prod, avail, np.arange(G), EXP)
    d_panelcol = float(np.nanmax(np.abs(W_308 - W_308_fullaxis)))
    # D3 alone: 308's weight rule on blended vs raw rows, same columns
    W_308_raw = normrel_W(rows_raw, avail, wcols, EXP)
    d_pmweight = float(np.nanmax(np.abs(W_308 - W_308_raw)))
    d_total = maxdiff(prod_final, a_final)
    d_total_nodegen = float(np.abs(prod_final[nd][:, wcols] - a_final[nd][:, wcols]).max())
    # D7: production's PM correction writes into every official column PerturbMatch observes; 263's blend is
    # restricted to k562 ∩ PM ∩ non-panel.  Same ordering, same weights, only the PM column set differs.
    B_pmcols = ord_B(rows_prod, avail, W_prod_rule, pre_prod64, pm_row_prod, pm_avail, pm_cols, k_index, BLEND)
    B_shared = ord_B(rows_prod, avail, W_prod_rule, pre_prod64, pm_row_prod, pm_avail, shared_cols, k_index, BLEND)
    d_pmcols = ALPHA * maxdiff(B_pmcols, B_shared)
    # what a PARITY-CORRECT production arm would emit: blend first, then weight, then norm-match -- on
    # PRODUCTION's substrate and column policy, so this is the size of the change the fix would make.
    rows_prod_blend = [R.copy() for R in rows_prod]
    rows_prod_blend[k_index][np.ix_(br, pm_cols)] = (
        (1 - BLEND) * rows_prod[k_index][np.ix_(br, pm_cols)] + BLEND * pm_row_prod[np.ix_(br, pm_cols)])
    pre_prod_blend64 = uniform_mean(rows_prod_blend)
    W_parity = normrel_W_production_rule_rows(rows_prod_blend, avail, np.arange(G), EXP)
    parity = ALPHA * ord_A(rows_prod_blend, avail, W_parity, pre_prod_blend64)
    d_parity_vs_prod = float(np.abs(parity[:, wcols] - prod_final[:, wcols]).max())
    cos_parity = np.array([
        float(parity[t, wcols] @ prod_final[t, wcols] /
              max(np.linalg.norm(parity[t, wcols]) * np.linalg.norm(prod_final[t, wcols]), 1e-300))
        for t in range(n)])
    dn_parity = np.linalg.norm(parity, axis=1) / np.maximum(np.linalg.norm(prod_incumbent, axis=1), 1e-300)

    P("")
    P(f"[313] TOTAL divergence production vs 308, non-panel columns, pre-whiten effect units: "
      f"max|diff| {d_total:.4e} (float32 floor {f32_floor:.2e}, ratio {d_total / max(f32_floor, 1e-300):.3g}x); "
      f"excluding the degenerate row {d_total_nodegen:.4e}; median row norm {float(np.median(nrm_p)):.4f}")
    P(f"[313]   D1+D2 PM ordering alone (identical weights, identical substrate): {ALPHA * d_order:.4e}")
    P(f"[313]   D3+D5 weight rule alone (identical ordering):                     {ALPHA * d_weight:.4e}")
    P(f"[313]   D7 PM column policy alone (broad official vs k562∩PM∩non-panel):  {d_pmcols:.4e}")
    P(f"[313]   D5 panel columns in the norm: max|dW| {d_panelcol:.4e}")
    P(f"[313]   D3 PM inside the norm:        max|dW| {d_pmweight:.4e}")
    P(f"[313]   worst per-target row cosine between the two outputs: {np.nanmin(cosine):.6f} "
      f"(median {np.nanmedian(cosine):.6f}, excl. degenerate min {np.nanmin(cosine[nd]):.6f}); "
      f"{int((cosine < 0.99).sum())}/{n} targets below 0.99")
    wl = np.isfinite(W_308).all(1) & np.isfinite(W_prod_rule).all(1)
    P(f"[313]   weights: max|dW| {np.abs(W_308[wl] - W_prod_rule[wl]).max():.4e} over {int(wl.sum())} weighted targets; "
      f"median |dW| {np.median(np.abs(W_308[wl] - W_prod_rule[wl])):.4e}")
    P("")
    P(f"[313] PARITY-CORRECT variant (blend -> weight -> norm-match, on production's substrate and column "
      f"policy) vs production as shipped: max|diff| {d_parity_vs_prod:.4e}, row cosine median "
      f"{np.median(cos_parity):.6f} min {cos_parity.min():.6f}; delivered norm median {np.median(dn_parity):.6f} "
      f"range [{dn_parity.min():.6f}, {dn_parity.max():.6f}]")

    out = dict(
        fixture=dict(path=str((WORK / "fixture313.npz").relative_to(ROOT)), sha256=fx_sha, seed=int(F.SEED),
                     n_targets=n, n_genes=G, donors=list(model.names), wcols=int(wcols.size),
                     pm_columns=int(pm_cols.size), shared_columns=int(shared_cols.size)),
        constants=dict(alpha=ALPHA, beta=BETA, pm_blend=BLEND, exponent=EXP,
                       donor_weights={k: float(v) for k, v in uni.items()},
                       effect_scale=float(model.effect_scale),
                       float32_noise_floor=f32_floor),
        known_answers=KA,
        per_target=rows_out,
        attribution=dict(total_max_abs_diff=d_total, total_max_abs_diff_excl_degenerate=d_total_nodegen,
                         pm_ordering_only=ALPHA * d_order, weight_rule_only=ALPHA * d_weight,
                         pm_column_policy_only=d_pmcols,
                         panel_columns_weight_delta=d_panelcol, pm_in_norm_weight_delta=d_pmweight,
                         worst_row_cosine=float(np.nanmin(cosine)), median_row_cosine=float(np.nanmedian(cosine)),
                         worst_row_cosine_excl_degenerate=float(np.nanmin(cosine[nd])),
                         n_targets_cosine_below_0p99=int((cosine < 0.99).sum()),
                         max_weight_diff=float(np.abs(W_308[wl] - W_prod_rule[wl]).max()),
                         median_weight_diff=float(np.median(np.abs(W_308[wl] - W_prod_rule[wl]))),
                         pre312_uniform_baseline=d_uniform,
                         pre312_uniform_baseline_263_pm_policy=d_uniform_sharedpolicy),
        delivered_norm=dict(
            arm_over_incumbent_308=dict(median=float(np.median(dn_308)), min=float(dn_308.min()),
                                        max=float(dn_308.max())),
            arm_over_incumbent_production=dict(median=float(np.median(dn_prod)), min=float(dn_prod.min()),
                                               max=float(dn_prod.max())),
            arm_over_incumbent_production_excl_degenerate=dict(
                median=float(np.median(dn_prod[nd])), min=float(dn_prod[nd].min()), max=float(dn_prod[nd].max())),
            panel_total_ratio_production=float(np.linalg.norm(prod_final) / np.linalg.norm(prod_incumbent)),
            panel_total_ratio_308=float(np.linalg.norm(a_final) / np.linalg.norm(a_incumbent)),
            per_target_308=[float(x) for x in dn_308],
            per_target_production=[float(x) for x in dn_prod]),
        parity_variant=dict(max_abs_diff_vs_shipped=d_parity_vs_prod,
                            median_row_cosine_vs_shipped=float(np.median(cos_parity)),
                            min_row_cosine_vs_shipped=float(cos_parity.min()),
                            delivered_norm_median=float(np.median(dn_parity)),
                            delivered_norm_min=float(dn_parity.min()),
                            delivered_norm_max=float(dn_parity.max())),
        pm_audit={k: v for k, v in pm_audit.items() if not isinstance(v, list)},
        peak_gib=peak_gib(),
    )
    (WORK / "parity313.json").write_text(json.dumps(out, indent=2))
    (WORK / "parity313.log").write_text("\n".join(LINES) + "\n")
    P("")
    P(f"[313] KA rows {sum(1 for k in KA if k['ok'])}/{len(KA)} pass; peak {peak_gib():.2f} GiB")
    if not all(k["ok"] for k in KA):
        raise SystemExit("KNOWN-ANSWER FAILURE")


if __name__ == "__main__":
    main()
