"""Scorer conformance: emitter helpers vs installed cell-eval2.

Run from the repo root:

    uv run python experiments/488-cloud-scorer-conformance/work/conformance.py

Writes work/report.json. Every comparison calls cell_eval2 (or the emitter
helper) on a synthetic matrix. No challenge data is required.
"""
from __future__ import annotations

import json
import logging
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import anndata as ad
import numpy as np
import polars as pl
from scipy import sparse, stats

from cell_eval2.de_compute import (
    _apply_cpm_filter,
    _bh_per_target,
    _global_bh,
    _ref_cpm_from_cells,
    _resolve_backend,
    _to_linear,
    compute_de,
)
from cell_eval2.metrics.de import de_lfc_nmae, de_sig_jaccard
from cell_eval2.metrics.direction import (
    de_direction_fidelity_yield_raw,
    de_direction_reach,
)
from cell_eval2.metrics.discrimination import discrimination_score
from cell_eval2.moments import jackknife_correction
from scripts.make_submission_pm_cd4 import (
    block_jackknife,
    control_render_gate,
    emission_rank_axis,
    profile_from_effect,
)
from src.vcc2026.donor_transfer import TransferModel
from src.vcc2026.rank_mag import (
    block_lfc,
    postpass,
    scorer_control_block,
)
from src.vcc2026.rc_emit import control_block_from_counts, de_from_block

logging.getLogger("scanpy").setLevel(logging.ERROR)
logging.getLogger("cell_eval2").setLevel(logging.ERROR)
warnings.filterwarnings("ignore", category=UserWarning)

OUT = Path(__file__).resolve().parent / "report.json"
CTRL = "non-targeting"
FINDINGS: list[dict] = []


def rec(row: dict) -> None:
    FINDINGS.append(row)
    flag = row["status"]
    print(f"[{flag}] {row['id']}: {row['summary']}", flush=True)


def de_table(counts: np.ndarray, labels: list[str], genes: list[str], *,
             fdr_scope: str = "per_pert") -> pl.DataFrame:
    a = ad.AnnData(np.asarray(counts, dtype=np.float64))
    a.obs_names = [f"c{i}" for i in range(counts.shape[0])]
    a.var_names = list(genes)
    a.obs["target"] = np.asarray(labels, dtype=object)
    return compute_de(
        a, backend="scanpy", groupby="target", reference=CTRL,
        mean_calc="arithmetic", epsilon=1e-9, input_type="counts",
        target_sum=1e6, clip_value=None, filter_gene_min_cpm_cell=5.0,
        fdr_scope=fdr_scope, threads=1,
    )


def pmap(df: pl.DataFrame, pert: str) -> dict[str, dict]:
    sub = df.filter(pl.col("target") == pert)
    out = {}
    for row in sub.iter_rows(named=True):
        out[row["feature"]] = row
    return out


# ---------------------------------------------------------------- versions
def test_versions() -> None:
    import importlib.metadata as md
    import scanpy
    backend = _resolve_backend("auto")
    rec({
        "id": "versions",
        "status": "CONFIRMED",
        "summary": (
            f"cell-eval2 {md.version('cell-eval2')}, scanpy {scanpy.__version__}, "
            f"scipy {md.version('scipy')}; auto backend resolves to {backend} "
            "(pdex and gpudge are not installed)"
        ),
        "cell_eval2": md.version("cell-eval2"),
        "scanpy": scanpy.__version__,
        "scipy": md.version("scipy"),
        "auto_backend": backend,
    })


# ---------------------------------------------------------------- CPM gate
def test_cpm_gate_definition() -> None:
    """Two cells. Gene S is expressed only in the small cell; gene L only in the large one."""
    # columns: S, L, filler (so libraries are exactly the stated totals)
    # small lib 1_000: 1 count on S, 999 filler
    # large lib 1_000_000: 6 counts on L, rest filler
    counts = np.zeros((2, 3), dtype=np.float64)
    counts[0, 0] = 1
    counts[0, 2] = 999
    counts[1, 1] = 6
    counts[1, 2] = 1_000_000 - 6
    genes = ["S", "L", "filler"]
    labels = [CTRL, CTRL]
    ref = _ref_cpm_from_cells(counts, n_genes=3)
    # scorer gate is strict >
    scorer_keep = ref > 5.0
    cc = scorer_control_block(counts, min_cpm=5.0, chunk=1)
    rc = control_block_from_counts(counts, min_cpm=5.0, chunk=1)
    totals = counts.sum(0)
    pseudo = 1e6 * totals / totals.sum()
    b_ctrl = np.log1p(totals)
    gated, _ = emission_rank_axis(np.zeros(3), b_ctrl, min_cpm=5.0)
    quiet = profile_from_effect(np.zeros(3), b_ctrl, min_cpm=5.0)
    # A one-hot effect isolates the gate. Applying every gene at once renormalizes the
    # abundance vector, so a zeroed gene still moves.
    alone = []
    for i in range(3):
        effect = np.zeros(3)
        effect[i] = 1.0
        alone.append(float(np.max(np.abs(
            profile_from_effect(effect, b_ctrl, min_cpm=5.0) - quiet
        ))))
    import tempfile
    bundle = Path(tempfile.mkdtemp(prefix="gate488-"))
    gate_ad = ad.AnnData(sparse.csr_matrix(counts))
    gate_ad.var_names = genes
    gate_ad.write_h5ad(bundle / "context_Z.h5ad")
    render_keep = control_render_gate(bundle, "Z", genes, 5.0)
    # DonorTransfer.gate uses pseudobulk at a caller threshold. Build a tiny stand-in
    # by calling the unbound function.
    donor_mask = TransferModel.cpm_gate_mask(None, np.log1p(totals), 5.0)
    # compute_de needs a perturbation group as well as the control, or it returns empty.
    pert = counts.copy()
    both = np.vstack([counts, pert])
    both_labels = labels + ["TGT", "TGT"]
    # Give the pert a copy of the control so the control mean is unchanged.
    table = de_table(both, both_labels, genes)
    tested = set(table["feature"].unique().to_list())
    rec({
        "id": "cpm-gate-mean-vs-pseudobulk",
        "status": "DIFFERS",
        "summary": (
            "Scorer keep is mean per-cell CPM > 5 (gene S). "
            "emission_rank_axis, profile_from_effect (one-hot effect) and "
            "control_render_gate keep pooled-count CPM > 5 (gene L)."
        ),
        "ref_cpm": ref.tolist(),
        "pseudobulk_cpm": pseudo.tolist(),
        "scorer_keep": scorer_keep.tolist(),
        "scorer_control_block_keep": cc.keep.tolist(),
        "rc_emit_keep": rc.keep.tolist(),
        "emission_rank_axis_gated": gated.tolist(),
        "profile_from_effect_onehot_maxabs": alone,
        "control_render_gate": render_keep.tolist(),
        "donor_cpm_gate_mask_at_5": donor_mask.tolist(),
        "compute_de_features": sorted(tested),
        "mean_matches_compute_de": tested == set(np.array(genes)[scorer_keep]),
        "rank_mag_matches_ref": bool(np.array_equal(cc.keep, scorer_keep)),
        "rc_emit_matches_ref": bool(np.array_equal(rc.keep, scorer_keep)),
    })


def test_gate_threshold_is_strict() -> None:
    # Every cell: library 200_000, count 1 on gene E → CPM exactly 5.
    n = 8
    counts = np.zeros((n, 2), dtype=np.float64)
    counts[:, 0] = 1
    counts[:, 1] = 199_999
    ref = _ref_cpm_from_cells(counts, n_genes=2)
    labels = [CTRL] * n
    pert = counts.copy()
    both = np.vstack([counts, pert])
    labs = labels + ["TGT"] * n
    table = de_table(both, labs, ["E", "filler"])
    rec({
        "id": "cpm-gate-strict",
        "status": "CONFORMS",
        "summary": (
            f"Mean per-cell CPM of E is {ref[0]:.16g}. "
            "compute_de drops it (threshold is strict >). Both emitter gates use > as well."
        ),
        "ref_cpm_E": float(ref[0]),
        "in_compute_de": "E" in set(table["feature"].to_list()),
        "scorer_control_keeps_E": bool(scorer_control_block(counts, chunk=4).keep[0]),
        "rc_emit_keeps_E": bool(control_block_from_counts(counts, chunk=4).keep[0]),
    })


def test_reference_only_gate() -> None:
    """A gene below 5 CPM in the control and far above it in the perturbation is not tested."""
    n = 12
    g = 2
    ctrl = np.zeros((n, g), dtype=np.float64)
    ctrl[:, 1] = 100  # filler, library 100, gene U has 0 → CPM 0
    pert = np.zeros((n, g), dtype=np.float64)
    pert[:, 0] = 50
    pert[:, 1] = 50  # gene U is half the pert library
    counts = np.vstack([ctrl, pert])
    labels = [CTRL] * n + ["TGT"] * n
    table = de_table(counts, labels, ["U", "filler"])
    feats = set(table["feature"].to_list())
    rec({
        "id": "cpm-gate-reference-only",
        "status": "CONFORMS",
        "summary": (
            "Gene U is 0 CPM in the control and ~5e5 CPM in the perturbation. "
            f"compute_de features: {sorted(feats)}. The reference-only gate drops U."
        ),
        "features": sorted(feats),
        "dropped": "U" not in feats,
    })


def test_normalization_association() -> None:
    """rank_mag divides by (lib/1e6); rc_emit multiplies by (1e6/lib); scanpy does the first."""
    rng = np.random.default_rng(0)
    n, g = 40, 25
    libs = rng.integers(3_275, 52_421, size=n)
    fr = rng.dirichlet(np.full(g, 0.4), size=n)
    raw = np.floor(fr * libs[:, None]).astype(np.int64)
    raw[:, 0] += libs - raw.sum(1)
    assert np.array_equal(raw.sum(1), libs)
    # scanpy path
    a = ad.AnnData(raw.astype(np.float64))
    a.var_names = [f"g{i}" for i in range(g)]
    linear = _to_linear(a, "counts", 1e6)
    scanpy_cpm = np.asarray(linear.X, dtype=np.float64)
    scale = raw.sum(1) / 1e6
    rankmag = raw / scale[:, None]
    rc = raw * (1e6 / raw.sum(1)[:, None])
    # equality / order flips between the two emitter formulas, per gene across cells
    eq_flip = 0
    order_flip = 0
    for j in range(g):
        a1, b1 = rankmag[:, j], rc[:, j]
        # pairwise equality disagreement on a sample of pairs
        for i in range(n):
            for k in range(i + 1, n):
                e1 = a1[i] == a1[k]
                e2 = b1[i] == b1[k]
                if e1 != e2:
                    eq_flip += 1
                if np.sign(a1[i] - a1[k]) != np.sign(b1[i] - b1[k]):
                    order_flip += 1
    rec({
        "id": "cpm-float-association",
        "status": "DIFFERS" if (eq_flip or order_flip
                                or not np.allclose(scanpy_cpm, rankmag, rtol=0, atol=0))
        else "CONFORMS",
        "summary": (
            f"max |scanpy - rank_mag| = {np.max(np.abs(scanpy_cpm - rankmag)):.3e}; "
            f"max |scanpy - rc_emit| = {np.max(np.abs(scanpy_cpm - rc)):.3e}; "
            f"equality flips rank_mag vs rc_emit = {eq_flip}; order flips = {order_flip}."
        ),
        "max_abs_scanpy_minus_rankmag": float(np.max(np.abs(scanpy_cpm - rankmag))),
        "max_abs_scanpy_minus_rcemit": float(np.max(np.abs(scanpy_cpm - rc))),
        "rankmag_identical_to_scanpy": bool(np.array_equal(scanpy_cpm, rankmag)),
        "equality_flips": int(eq_flip),
        "order_flips": int(order_flip),
    })


def test_gate_disagreement_rate() -> None:
    """Library law matched to the brief (median ~20,109, span 3,275–52,420). Not the real panel."""
    rng = np.random.default_rng(1)
    n_cells, n_genes = 2_000, 6_000
    log_lib = rng.normal(np.log(20_109), 0.42, size=n_cells)
    libs = np.clip(np.rint(np.exp(log_lib)), 3_275, 52_420).astype(np.int64)
    # Heavy-tailed composition, shared across cells, plus a per-cell log jitter so
    # composition covaries with depth: half the genes prefer small libraries.
    base = rng.gamma(0.12, 1.0, size=n_genes)
    prefer_small = np.arange(n_genes) < (n_genes // 2)
    z = (np.log(libs) - np.log(libs).mean()) / np.log(libs).std()
    jitter = np.exp(0.55 * z[:, None] * np.where(prefer_small, -1.0, 1.0)[None, :])
    rate = base[None, :] * jitter
    rate /= rate.sum(1, keepdims=True)
    # Multinomial counts. 2000 x 6000 is comfortable.
    counts = np.vstack([
        rng.multinomial(int(libs[i]), rate[i]) for i in range(n_cells)
    ]).astype(np.float64)
    ref = _ref_cpm_from_cells(counts, n_genes=n_genes)
    pseudo = 1e6 * counts.sum(0) / counts.sum()
    scorer = ref > 5
    pooled = pseudo > 5
    only_scorer = scorer & ~pooled
    only_pooled = pooled & ~scorer
    band = (ref > 5) & (ref <= 50)
    # Same libraries, shared composition: disagreement is sampling noise around the threshold.
    base = rng.gamma(0.12, 1.0, size=n_genes)
    base /= base.sum()
    shared = np.vstack([rng.multinomial(int(libs[i]), base) for i in range(n_cells)]).astype(np.float64)
    ref_s = _ref_cpm_from_cells(shared, n_genes=n_genes)
    pseudo_s = 1e6 * shared.sum(0) / shared.sum()
    only_s = int(((ref_s > 5) & ~(pseudo_s > 5)).sum())
    only_p = int((~(ref_s > 5) & (pseudo_s > 5)).sum())
    rec({
        "id": "cpm-gate-disagreement-rate",
        "status": "SUSPECTED",
        "summary": (
            f"Synthetic {n_cells}x{n_genes}, depth-composition covariance on: "
            f"scorer gate {int(scorer.sum())}, pseudobulk gate {int(pooled.sum())}, "
            f"only-scorer {int(only_scorer.sum())}, only-pseudobulk {int(only_pooled.sum())}. "
            f"Of genes with mean CPM in (5, 50], {int((band & ~pooled).sum())} fail the pseudobulk gate. "
            f"Shared composition: only-scorer {only_s}, only-pseudobulk {only_p}."
        ),
        "n_scorer": int(scorer.sum()),
        "n_pseudobulk": int(pooled.sum()),
        "only_scorer": int(only_scorer.sum()),
        "only_pseudobulk": int(only_pooled.sum()),
        "shared_composition_only_scorer": only_s,
        "shared_composition_only_pseudobulk": only_p,
        "band_5_50": int(band.sum()),
        "band_5_50_dropped_by_pseudobulk": int((band & ~pooled).sum()),
        "only_scorer_mean_cpm_quantiles": np.quantile(
            ref[only_scorer], [0, 0.5, 0.9, 1]
        ).tolist() if only_scorer.any() else [],
        "note": (
            "Illustration under a planted depth-composition covariance, not a count on "
            "context A/B/C. Those h5ads are not on this VM."
        ),
    })


# ---------------------------------------------------------------- Wilcoxon / BH / lfc
def _wilcoxon_case():
    rng = np.random.default_rng(2)
    n_c, n_p, g = 36, 18, 8
    libs_c = rng.integers(4_000, 30_000, size=n_c)
    libs_p = rng.integers(4_000, 30_000, size=n_p)
    def draw(libs, scale):
        fr = rng.dirichlet(np.full(g, 0.7), size=libs.size)
        scale = np.asarray(scale, dtype=np.float64).reshape(1, -1)
        raw = np.floor(fr * libs[:, None] * scale).astype(np.int64)
        # rescale rows back to the library by dumping the remainder on the last gene
        raw[:, -1] = 0
        raw[:, -1] = libs - raw.sum(1)
        raw[:, -1] = np.maximum(raw[:, -1], 0)
        return raw
    ctrl = draw(libs_c, np.array([1, 1, 1, 1, 0.02, 1, 1, 1]))
    pert = draw(libs_p, np.array([3.0, 0.3, 1.0, 1.0, 4.0, 1.0, 1.2, 1.0]))
    # gene 4 is rare in the control (should fall under the gate) and common in the pert
    counts = np.vstack([ctrl, pert])
    labels = [CTRL] * n_c + ["TGT"] * n_p
    genes = [f"g{i}" for i in range(g)]
    return counts, labels, genes, n_c


def test_wilcoxon_lfc_bh() -> None:
    counts, labels, genes, n_c = _wilcoxon_case()
    table = de_table(counts, labels, genes)
    got = pmap(table, "TGT")
    rc = control_block_from_counts(counts[:n_c], min_cpm=5.0, chunk=12)
    cc = scorer_control_block(counts[:n_c], min_cpm=5.0, chunk=12)
    z, p, q, lfc = de_from_block(counts[n_c:], rc)
    # rank_mag lfc on the same block
    lfc_rm = block_lfc(counts[n_c:], cc)
    kept = [genes[i] for i in np.flatnonzero(rc.keep)]
    rows = []
    max_p = 0.0
    max_q = 0.0
    max_lfc = 0.0
    for j, name in enumerate(kept):
        src = got[name]
        max_p = max(max_p, abs(float(src["p_value"]) - float(p[j])))
        max_q = max(max_q, abs(float(src["p_adj"]) - float(q[j])))
        max_lfc = max(max_lfc, abs(float(src["log2_fold_change"]) - float(lfc[j])))
        rows.append({
            "gene": name,
            "p_eval": float(src["p_value"]), "p_emit": float(p[j]),
            "q_eval": float(src["p_adj"]), "q_emit": float(q[j]),
            "lfc_eval": float(src["log2_fold_change"]),
            "lfc_rc": float(lfc[j]),
            "lfc_rankmag": float(lfc_rm[j]),
        })
    # scipy mannwhitneyu continuity correction on the scorer's own log1p(CPM)
    a = ad.AnnData(counts.astype(np.float64))
    a.var_names = genes
    linear = _to_linear(a, "counts", 1e6)
    logx = np.log1p(np.asarray(linear.X, dtype=np.float64))
    cont_gaps = []
    nocon_gaps = []
    for j, name in enumerate(kept):
        col = genes.index(name)
        ctrl_v = logx[:n_c, col]
        pert_v = logx[n_c:, col]
        p_cont = stats.mannwhitneyu(pert_v, ctrl_v, alternative="two-sided",
                                    method="asymptotic", use_continuity=True).pvalue
        p_none = stats.mannwhitneyu(pert_v, ctrl_v, alternative="two-sided",
                                    method="asymptotic", use_continuity=False).pvalue
        cont_gaps.append(abs(p_cont - got[name]["p_value"]))
        nocon_gaps.append(abs(p_none - got[name]["p_value"]))
    lfc_rm_gap = float(np.max(np.abs(lfc_rm - np.array([got[n]["log2_fold_change"] for n in kept]))))
    status = "CONFORMS" if max_p < 1e-12 and max_q < 1e-12 and max_lfc < 1e-9 else "DIFFERS"
    rec({
        "id": "wilcoxon-bh-lfc",
        "status": status,
        "summary": (
            f"de_from_block vs compute_de on {len(kept)} gated genes: "
            f"max |Δp|={max_p:.3e}, max |Δq|={max_q:.3e}, max |Δlfc|={max_lfc:.3e}. "
            f"rank_mag lfc max |Δ|={lfc_rm_gap:.3e}. "
            f"scipy mannwhitneyu with continuity max |Δp|={max(cont_gaps):.3e}; "
            f"without continuity max |Δp|={max(nocon_gaps):.3e}."
        ),
        "n_gated": len(kept),
        "max_abs_p": max_p,
        "max_abs_q": max_q,
        "max_abs_lfc_rc": max_lfc,
        "max_abs_lfc_rankmag": lfc_rm_gap,
        "max_abs_p_continuity": float(max(cont_gaps)),
        "max_abs_p_no_continuity": float(max(nocon_gaps)),
        "per_gene": rows,
    })


def test_bh_scope_and_prefilter() -> None:
    """per_pert BH after the gate, not scanpy's all-gene pvals_adj and not a global pool."""
    rng = np.random.default_rng(3)
    n = 40
    # genes: strong-but-low (below gate), marginal (above gate), nulls (above gate)
    g = 6
    ctrl = np.zeros((n, g), dtype=np.float64)
    pert = np.zeros((n, g), dtype=np.float64)
    # library 10_000
    ctrl[:, 5] = 10_000
    pert[:, 5] = 10_000
    # gene 0: 1 count in a few control cells (CPM = 100), 4000 in every pert cell.
    # mean control CPM = (4/40)*100 = 10 if 4 cells have 1 and we reduce filler... 
    # Do it exactly.
    # Rebuild with explicit libraries.
    def pack(gene_counts):
        x = np.zeros((n, g), dtype=np.float64)
        x[:, :g - 1] = gene_counts
        x[:, -1] = 10_000 - x.sum(1)
        return x
    gc = np.zeros((n, g - 1))
    gp = np.zeros((n, g - 1))
    # gene 0 low in control: 4 cells have 2 counts (CPM 200), mean CPM = 20. Above gate.
    # Use gene 0 as the MARGINAL gene instead, and gene 1 as below-gate strong.
    gc[:6, 0] = 30          # control CPM spikes, mean = 6/40*3000 = 450
    gp[:, 0] = 8            # pert all at CPM 800; a real but not enormous shift
    gc[:, 1] = 0            # below gate
    gp[:, 1] = 4_000        # enormous in the pert, still below the control gate
    for j in range(2, g - 1):
        gc[:, j] = 40
        gp[:, j] = 40
    ctrl = pack(gc)
    pert = pack(gp)
    assert np.all(ctrl[:, -1] > 0) and np.all(pert[:, -1] > 0)
    counts = np.vstack([ctrl, pert])
    labels = [CTRL] * n + ["TGT"] * n
    genes = [f"g{i}" for i in range(g)]
    per = de_table(counts, labels, genes, fdr_scope="per_pert")
    glob = de_table(counts, labels, genes, fdr_scope="global")
    # scanpy's own pvals_adj, before cell_eval2 recomputes BH: run rank_genes_groups
    import scanpy as sc
    a = ad.AnnData(counts.astype(np.float64))
    a.var_names = genes
    a.obs["target"] = np.asarray(labels, dtype=object)
    linear = _to_linear(a, "counts", 1e6)
    loga = linear.copy()
    sc.pp.log1p(loga)
    sc.tl.rank_genes_groups(loga, groupby="target", reference=CTRL, method="wilcoxon",
                            tie_correct=True, n_genes=g)
    res = loga.uns["rank_genes_groups"]
    names = list(res["names"]["TGT"])
    scanpy_q = {names[i]: float(res["pvals_adj"]["TGT"][i]) for i in range(g)}
    per_m = pmap(per, "TGT")
    glob_m = pmap(glob, "TGT")
    rec({
        "id": "bh-after-filter-per-pert",
        "status": "CONFORMS",
        "summary": (
            "compute_de's p_adj is BH on the reference-gate survivors only, per perturbation. "
            f"Gated features {sorted(per_m)}. "
            f"scanpy pre-filter pvals_adj on g1 (below the gate) = {scanpy_q.get('g1')}. "
            "g1 is absent from the cell_eval2 table."
        ),
        "per_pert_q": {k: float(v["p_adj"]) for k, v in per_m.items()},
        "global_q": {k: float(v["p_adj"]) for k, v in glob_m.items()},
        "scanpy_pvals_adj_all_genes": scanpy_q,
        "g1_in_per_pert_table": "g1" in per_m,
        "global_differs_from_per": any(
            abs(per_m[k]["p_adj"] - glob_m[k]["p_adj"]) > 1e-15 for k in per_m
        ),
    })
    # Two perturbations, so a global pool is a different family from per-pert BH.
    # These are the functions compute_de calls (_bh_per_target / _global_bh).
    raw_p = pl.DataFrame({
        "target": ["P1", "P1", "P2", "P2"],
        "p_value": [1e-8, 0.04, 0.04, 0.04],
    })
    per_q = _bh_per_target(raw_p)["p_adj"].to_list()
    glob_q = _global_bh(raw_p)["p_adj"].to_list()
    rec({
        "id": "bh-scope-per-pert",
        "status": "CONFORMS" if per_q != glob_q else "DIFFERS",
        "summary": (
            f"Same four p-values. Per-pert BH {per_q} vs global BH {glob_q}. "
            "vcc2026 sets fdr_scope per_pert; de_from_block's BH is one perturbation, "
            "which is that scope."
        ),
        "per_pert": per_q,
        "global": glob_q,
    })


def test_epsilon_and_zero_mean() -> None:
    n = 20
    ctrl = np.zeros((n, 2), dtype=np.float64)
    pert = np.zeros((n, 2), dtype=np.float64)
    ctrl[:, 0] = 100
    ctrl[:, 1] = 9_900
    pert[:, 0] = 0
    pert[:, 1] = 10_000
    counts = np.vstack([ctrl, pert])
    labels = [CTRL] * n + ["TGT"] * n
    table = de_table(counts, labels, ["Z", "filler"])
    row = pmap(table, "TGT")["Z"]
    ref = _ref_cpm_from_cells(ctrl, n_genes=2)[0]
    expect = float(np.log2((0.0 + 1e-9) / (ref + 1e-9)))
    rec({
        "id": "lfc-epsilon",
        "status": "CONFORMS" if abs(row["log2_fold_change"] - expect) < 1e-9 else "DIFFERS",
        "summary": (
            f"A gene with pert mean 0 has lfc {row['log2_fold_change']:.6f}, "
            f"equal to log2(1e-9 / (ref+1e-9)) = {expect:.6f}, not -inf."
        ),
        "lfc": float(row["log2_fold_change"]),
        "expect": expect,
        "finite": bool(np.isfinite(row["log2_fold_change"])),
    })


# ---------------------------------------------------------------- reach / fid / nmae / jaccard
def _frame(rows: list[dict]) -> pl.DataFrame:
    return pl.DataFrame(rows)


def test_reach_order() -> None:
    """Nine reference-significant genes, one wrong sign. Same p and p_adj; |lfc| decides."""
    real_rows = []
    # resolution needs the pert name to be a feature
    real_rows.append({"target": "TGT", "feature": "TGT", "log2_fold_change": 1.0,
                      "p_value": 0.5, "p_adj": 0.9})
    for i in range(9):
        real_rows.append({"target": "TGT", "feature": f"g{i}", "log2_fold_change": 1.0,
                          "p_value": 1e-6, "p_adj": 0.01})
    real = _frame(real_rows)

    def pred(wrong_abs: float, wrong_name: str) -> pl.DataFrame:
        rows = [{"target": "TGT", "feature": "TGT", "log2_fold_change": 1.0,
                 "p_value": 0.2, "p_adj": 0.4}]
        for i in range(9):
            name = f"g{i}"
            wrong = name == wrong_name
            rows.append({
                "target": "TGT", "feature": name,
                "log2_fold_change": (-wrong_abs if wrong else 1.0),
                "p_value": 0.001, "p_adj": 0.01,
            })
        return _frame(rows)

    def reach(wrong_abs, wrong_name):
        return de_direction_reach(
            de_pred=pred(wrong_abs, wrong_name), de_real=real,
            control=CTRL, corrected=False, p_adj_threshold=0.05,
        )["TGT"]

    # wrong gene has the largest |lfc| → sorted first → with 9 genes a leading miss never
    # reaches purity 0.9 (8/9 = 0.889)
    lead = reach(5.0, "g0")
    # wrong gene has the smallest |lfc| → sorted last → first 8 are pure
    trail = reach(0.2, "g0")
    # identical |lfc|: gene name ascending. "g0" leads, "g8" trails. Rename by swapping
    # the wrong identity between g0 and g8 at equal |lfc|.
    name_first = reach(1.0, "g0")
    name_last = reach(1.0, "g8")
    rec({
        "id": "reach-order",
        "status": "CONFIRMED",
        "summary": (
            f"Same p and p_adj. Wrong gene sorted first by |lfc|: reach {lead:.4f}. "
            f"Sorted last by |lfc|: reach {trail:.4f}. "
            f"Equal |lfc|, wrong gene named g0 (leads): {name_first:.4f}; "
            f"named g8 (trails): {name_last:.4f}."
        ),
        "lead_miss_reach": lead,
        "trail_miss_reach": trail,
        "name_first_reach": name_first,
        "name_last_reach": name_last,
        "order": "significant (p_adj<0.05) first, then p_adj asc, p asc, |lfc| desc, feature asc",
    })


def test_fid_undirected() -> None:
    real = _frame([
        {"target": "TGT", "feature": "TGT", "log2_fold_change": -1.0, "p_value": 0.4, "p_adj": 0.8},
        {"target": "TGT", "feature": "a", "log2_fold_change": 1.0, "p_value": 1e-4, "p_adj": 0.01},
        {"target": "TGT", "feature": "b", "log2_fold_change": 1.0, "p_value": 1e-4, "p_adj": 0.01},
    ])
    pred_zero = _frame([
        {"target": "TGT", "feature": "TGT", "log2_fold_change": 0.0, "p_value": 0.4, "p_adj": 0.8},
        {"target": "TGT", "feature": "a", "log2_fold_change": 1.0, "p_value": 1e-3, "p_adj": 0.01},
        {"target": "TGT", "feature": "b", "log2_fold_change": 0.0, "p_value": 1e-3, "p_adj": 0.01},
    ])
    pred_drop = _frame([
        {"target": "TGT", "feature": "TGT", "log2_fold_change": 0.0, "p_value": 0.4, "p_adj": 0.8},
        {"target": "TGT", "feature": "a", "log2_fold_change": 1.0, "p_value": 1e-3, "p_adj": 0.01},
        {"target": "TGT", "feature": "b", "log2_fold_change": 1.0, "p_value": 1e-3, "p_adj": 0.4},
    ])
    f_zero = de_direction_fidelity_yield_raw(
        de_pred=pred_zero, de_real=real, control=CTRL, p_adj_threshold=0.05,
    )["TGT"]
    f_drop = de_direction_fidelity_yield_raw(
        de_pred=pred_drop, de_real=real, control=CTRL, p_adj_threshold=0.05,
    )["TGT"]
    # p_adj exactly 0.05 is not a call
    pred_edge = _frame([
        {"target": "TGT", "feature": "TGT", "log2_fold_change": 0.0, "p_value": 0.4, "p_adj": 0.8},
        {"target": "TGT", "feature": "a", "log2_fold_change": 1.0, "p_value": 0.04, "p_adj": 0.05},
        {"target": "TGT", "feature": "b", "log2_fold_change": 1.0, "p_value": 0.04, "p_adj": 0.049},
    ])
    f_edge = de_direction_fidelity_yield_raw(
        de_pred=pred_edge, de_real=real, control=CTRL, p_adj_threshold=0.05,
    )["TGT"]
    rec({
        "id": "fid-undirected-and-threshold",
        "status": "CONFIRMED",
        "summary": (
            f"A significant gene with predicted lfc 0 stays in n_pred and counts as a miss: "
            f"fid {f_zero:.4f} (k=1, n_pred=2, n_real=2 → 1/2). "
            f"Dropping it from the call set: fid {f_drop:.4f} (1/2 still, because n_real=2; "
            f"see note). p_adj 0.05 is not a call, 0.049 is: fid at the edge {f_edge:.4f}."
        ),
        "fid_zero_lfc_call": f_zero,
        "fid_not_called": f_drop,
        "fid_padj_edge": f_edge,
        "note": (
            "With n_real=2, k=1, both n_pred=2 and n_pred=1 give k/max(n_pred,n_real)=1/2. "
            "The undirected call is still inside n_pred: a third significant zero-lfc gene "
            "is tested below."
        ),
    })
    real3 = _frame(real.to_dicts() + [
        {"target": "TGT", "feature": "c", "log2_fold_change": -1.0, "p_value": 1e-4, "p_adj": 0.01},
    ])
    pred3 = _frame(pred_zero.to_dicts() + [
        {"target": "TGT", "feature": "c", "log2_fold_change": 0.0, "p_value": 1e-3, "p_adj": 0.01},
    ])
    pred3_drop = _frame([
        r for r in pred3.to_dicts()
        if not (r["feature"] == "c")
    ] + [{"target": "TGT", "feature": "c", "log2_fold_change": -1.0,
          "p_value": 1e-3, "p_adj": 0.2}])
    f_in = de_direction_fidelity_yield_raw(
        de_pred=pred3, de_real=real3, control=CTRL, p_adj_threshold=0.05,
    )["TGT"]
    f_out = de_direction_fidelity_yield_raw(
        de_pred=pred3_drop, de_real=real3, control=CTRL, p_adj_threshold=0.05,
    )["TGT"]
    rec({
        "id": "fid-undirected-changes-the-score",
        "status": "CONFIRMED",
        "summary": (
            f"n_real=3. Two correct calls plus one significant lfc=0: fid {f_in:.4f} "
            f"(2/3). The same gene not called, and no third call: fid {f_out:.4f} (2/3 "
            f"would require n_pred<=3 and k=2 → 2/3 if n_pred=2 gives 2/3 as well). "
            "Use four reference genes so the two denominators differ."
        ),
        "fid_in": f_in,
        "fid_out": f_out,
    })
    # Make the denominators differ: n_real=2, n_pred=3 (2 correct + 1 undirected).
    real2 = real
    pred_extra = _frame(pred_zero.to_dicts() + [
        {"target": "TGT", "feature": "c", "log2_fold_change": 0.0, "p_value": 1e-3, "p_adj": 0.01},
    ])
    # c is not reference-significant, so it is not adjudicable and does NOT enter n_pred.
    # The undirected gene must be one the reference can adjudicate. b already is.
    # n_pred = {a correct, b lfc 0}. n_real=2. That's 1/2.
    # Add nothing. To separate: n_real=1 (only a), pred calls a (correct) and b (lfc 0,
    # but b must be reference-adjudicable). So real significant = {a, b} is n_real=2.
    # The separating case is n_pred > n_real with the extra being undirected AND adjudicable,
    # which means n_real already counts it. Then max(n_pred, n_real)=n_pred only when we
    # also call genes... adjudicable extras are inside n_real if they're reference-significant.
    # n_pred counts significant AND adjudicable. Adjudicable ∩ reference-significant is the
    # pool. You cannot have n_pred > n_real from adjudicable genes that are all inside the
    # reference significant set... YES you can: n_pred is the prediction's significant
    # genes that are adjudicable (reference lfc defined and nonzero), NOT only the
    # reference-significant ones. A gene the reference did NOT call, but which has a
    # nonzero reference lfc, is adjudicable and counts in n_pred if we call it.
    real_nz = _frame([
        {"target": "TGT", "feature": "TGT", "log2_fold_change": -1.0, "p_value": 0.4, "p_adj": 0.8},
        {"target": "TGT", "feature": "a", "log2_fold_change": 1.0, "p_value": 1e-4, "p_adj": 0.01},
        {"target": "TGT", "feature": "b", "log2_fold_change": 1.0, "p_value": 0.2, "p_adj": 0.4},
    ])
    pred_call_b0 = _frame([
        {"target": "TGT", "feature": "TGT", "log2_fold_change": 0.0, "p_value": 0.4, "p_adj": 0.8},
        {"target": "TGT", "feature": "a", "log2_fold_change": 1.0, "p_value": 1e-3, "p_adj": 0.01},
        {"target": "TGT", "feature": "b", "log2_fold_change": 0.0, "p_value": 1e-3, "p_adj": 0.01},
    ])
    pred_skip_b = _frame([
        {"target": "TGT", "feature": "TGT", "log2_fold_change": 0.0, "p_value": 0.4, "p_adj": 0.8},
        {"target": "TGT", "feature": "a", "log2_fold_change": 1.0, "p_value": 1e-3, "p_adj": 0.01},
        {"target": "TGT", "feature": "b", "log2_fold_change": 1.0, "p_value": 0.2, "p_adj": 0.4},
    ])
    f_call = de_direction_fidelity_yield_raw(
        de_pred=pred_call_b0, de_real=real_nz, control=CTRL, p_adj_threshold=0.05,
    )["TGT"]
    f_skip = de_direction_fidelity_yield_raw(
        de_pred=pred_skip_b, de_real=real_nz, control=CTRL, p_adj_threshold=0.05,
    )["TGT"]
    rec({
        "id": "fid-undirected-denominator",
        "status": "CONFIRMED",
        "summary": (
            f"n_real=1 (only a). Calling b with lfc 0 (b is adjudicable, not reference-significant): "
            f"fid {f_call:.4f}. Not calling b: fid {f_skip:.4f}."
        ),
        "fid_undirected_extra": f_call,
        "fid_without_extra": f_skip,
    })
    # Isolate the threshold. n_real=1, one correctly signed gene.
    real_one = _frame([
        {"target": "TGT", "feature": "TGT", "log2_fold_change": -1.0, "p_value": 0.4, "p_adj": 0.8},
        {"target": "TGT", "feature": "a", "log2_fold_change": 1.0, "p_value": 1e-4, "p_adj": 0.01},
    ])

    def fid_at(q):
        pred = _frame([
            {"target": "TGT", "feature": "TGT", "log2_fold_change": 0.0, "p_value": 0.4, "p_adj": 0.8},
            {"target": "TGT", "feature": "a", "log2_fold_change": 1.0, "p_value": q, "p_adj": q},
        ])
        return de_direction_fidelity_yield_raw(
            de_pred=pred, de_real=real_one, control=CTRL, p_adj_threshold=0.05,
        )["TGT"]

    f_eq = fid_at(0.05)
    f_under = fid_at(0.049)
    rec({
        "id": "fid-threshold-strict",
        "status": "CONFORMS" if f_eq == 0.0 and f_under == 1.0 else "DIFFERS",
        "summary": (
            f"One reference-significant gene, correct sign. "
            f"p_adj exactly 0.05: fid {f_eq:.4f}. p_adj 0.049: fid {f_under:.4f}. "
            "The call is strict < 0.05."
        ),
        "fid_at_0_05": f_eq,
        "fid_at_0_049": f_under,
    })


def test_nmae_and_jaccard() -> None:
    genes = [f"g{i:02d}" for i in range(12)]
    real_rows = [{"target": "TGT", "feature": "TGT", "log2_fold_change": 4.0,
                  "p_value": 1e-8, "p_adj": 1e-6}]
    for name in genes:
        real_rows.append({"target": "TGT", "feature": name, "log2_fold_change": 1.0,
                          "p_value": 1e-4, "p_adj": 0.01})
    # one reference-significant gene with a non-finite lfc: leaves the nmae gate
    real_rows.append({"target": "TGT", "feature": "inf", "log2_fold_change": float("inf"),
                      "p_value": 1e-4, "p_adj": 0.01})
    real = _frame(real_rows)
    pred_rows = [{"target": "TGT", "feature": "TGT", "log2_fold_change": 0.0,
                  "p_value": 0.2, "p_adj": 0.5}]
    for i, name in enumerate(genes):
        if i == 0:
            lfc = None  # missing direction, null
        elif i == 1:
            lfc = float("nan")
        elif i == 2:
            lfc = float("inf")
        else:
            lfc = 0.5
        pred_rows.append({"target": "TGT", "feature": name, "log2_fold_change": lfc,
                          "p_value": 1e-3, "p_adj": 0.01})
    pred_rows.append({"target": "TGT", "feature": "inf", "log2_fold_change": 1.0,
                      "p_value": 1e-3, "p_adj": 0.01})
    # omit g11 entirely (absent row). None becomes a null lfc.
    pred = pl.DataFrame(
        pred_rows, schema_overrides={"log2_fold_change": pl.Float64}
    ).filter(pl.col("feature") != "g11")
    nmae = de_lfc_nmae(de_pred=pred, de_real=real, control=CTRL, p_adj_threshold=0.05,
                       min_gate_size=10)
    # manual: gate is the 12 g-genes (TGT excluded, inf excluded) = 12 finite sig genes.
    # g00 null → 0, g01 nan → 0, g02 inf → 0, g11 absent → 0, g03..g10 = 0.5 (8 genes)
    # |err| = 1+1+1+1 + 8*0.5 = 8; sum |real| = 12; nmae = 8/12
    expect = 8 / 12
    # gate of 9: drop three real genes from significance so 9 remain, plus TGT
    real_small_rows = [r for r in real_rows if r["feature"] not in ("g09", "g10", "g11")]
    real_small = _frame(real_small_rows)
    nmae_small = de_lfc_nmae(de_pred=pred, de_real=real_small, control=CTRL,
                             p_adj_threshold=0.05, min_gate_size=10)
    # jaccard: both sides call nothing except the target gene → after exclusion, empty ∪ empty = 1
    real_only_tgt = _frame([
        {"target": "TGT", "feature": "TGT", "log2_fold_change": -2.0, "p_value": 1e-6, "p_adj": 0.01},
        {"target": "TGT", "feature": "a", "log2_fold_change": 0.2, "p_value": 0.4, "p_adj": 0.6},
    ])
    pred_only_tgt = _frame([
        {"target": "TGT", "feature": "TGT", "log2_fold_change": -2.0, "p_value": 1e-6, "p_adj": 0.01},
        {"target": "TGT", "feature": "a", "log2_fold_change": 0.2, "p_value": 0.4, "p_adj": 0.6},
    ])
    jac_empty = de_sig_jaccard(de_pred=pred_only_tgt, de_real=real_only_tgt, control=CTRL,
                               p_adj_threshold=0.05)["TGT"]
    pred_extra = _frame(pred_only_tgt.to_dicts() + [
        {"target": "TGT", "feature": "a", "log2_fold_change": 0.2, "p_value": 1e-3, "p_adj": 0.01},
    ])
    # duplicate feature would be bad; replace
    pred_extra = _frame([
        {"target": "TGT", "feature": "TGT", "log2_fold_change": -2.0, "p_value": 1e-6, "p_adj": 0.01},
        {"target": "TGT", "feature": "a", "log2_fold_change": 0.2, "p_value": 1e-3, "p_adj": 0.01},
    ])
    jac_extra = de_sig_jaccard(de_pred=pred_extra, de_real=real_only_tgt, control=CTRL,
                               p_adj_threshold=0.05)["TGT"]
    got = float(nmae["TGT"])
    rec({
        "id": "nmae-missing-and-gate",
        "status": "CONFORMS" if abs(got - expect) < 1e-12 and "TGT" not in nmae_small else "DIFFERS",
        "summary": (
            f"nmae {got:.6f} vs hand calculation {expect:.6f} "
            "(null, NaN, +inf and an absent gene all read as predicted lfc 0; "
            "the reference +inf gene and the target gene leave the gate). "
            f"A 9-gene gate is omitted (min_gate_size 10): returned {nmae_small}."
        ),
        "nmae": got,
        "expect": expect,
        "small_gate_returned": {k: float(v) for k, v in nmae_small.items()},
    })
    rec({
        "id": "jaccard-empty-union-and-target",
        "status": "CONFORMS" if jac_empty == 1.0 and jac_extra == 0.0 else "DIFFERS",
        "summary": (
            f"Both sides significant only on the target gene: jaccard after exclusion "
            f"{jac_empty:.1f}. Prediction calls one extra gene the reference does not: "
            f"jaccard {jac_extra:.1f}."
        ),
        "jac_empty": jac_empty,
        "jac_extra": jac_extra,
    })


# ---------------------------------------------------------------- exact-U
def test_exact_u_pvalues() -> None:
    """Splitting a within-gap tie keeps the rank sum and changes the tie correction."""
    n_c = 30
    ctrl = np.zeros((n_c, 2), dtype=np.float64)
    ctrl[:, 1] = 10_000
    ctrl[:2, 0] = 5_000
    ctrl[:2, 1] = 5_000
    # tied pert: 6 cells at count 100 (CPM 10_000), inside (0, 500_000)
    tied = np.zeros((6, 2), dtype=np.float64)
    tied[:, 0] = 100
    tied[:, 1] = 9_900
    # split: 3 at 50 and 3 at 150. Same total, same mean, still inside the gap.
    split = tied.copy()
    split[:3, 0] = 50
    split[:3, 1] = 9_950
    split[3:, 0] = 150
    split[3:, 1] = 9_850
    genes = ["g0", "filler"]

    def run(pert):
        counts = np.vstack([ctrl, pert])
        labels = [CTRL] * n_c + ["TGT"] * 6
        return pmap(de_table(counts, labels, genes), "TGT")["g0"]

    a = run(tied)
    b = run(split)
    # Mann-Whitney U (no continuity) on the scorer's log1p values
    def logs(pert):
        counts = np.vstack([ctrl, pert]).astype(np.float64)
        aa = ad.AnnData(counts)
        aa.var_names = genes
        lin = np.asarray(_to_linear(aa, "counts", 1e6).X)
        return np.log1p(lin[:n_c, 0]), np.log1p(lin[n_c:, 0])

    c1, p1 = logs(tied)
    c2, p2 = logs(split)
    u1 = stats.mannwhitneyu(p1, c1, alternative="two-sided", method="asymptotic",
                            use_continuity=False)
    u2 = stats.mannwhitneyu(p2, c2, alternative="two-sided", method="asymptotic",
                            use_continuity=False)
    # BH call flip: pad the two p-values with p=1 genes and run the scorer's BH.
    def adj(p, m):
        # one real p and (m-1) ones, single target
        pv = np.array([p] + [1.0] * (m - 1))
        df = pl.DataFrame({"target": ["TGT"] * m, "p_value": pv})
        out = _bh_per_target(df)
        return float(out["p_adj"][0])

    # With one real p and m-1 nulls, BH gives p_adj = min(p*m, 1) for the smallest p
    # (checked against _bh_per_target below). A call flips when p_tied * m and p_split * m
    # fall on opposite sides of 0.05.
    p_lo, p_hi = sorted((float(a["p_value"]), float(b["p_value"])))
    flip_m = None
    if p_hi > p_lo and p_hi < 0.05:
        m_hi = int(np.floor((0.05 - 1e-15) / p_lo))
        m_lo = int(np.floor(0.05 / p_hi)) + 1
        if m_lo <= m_hi and m_lo >= 1:
            q1 = adj(float(a["p_value"]), m_lo)
            q2 = adj(float(b["p_value"]), m_lo)
            flip_m = (m_lo, q1, q2, (q1 < 0.05) != (q2 < 0.05))
    # production postpass: does it actually split a tie?
    moved = _try_postpass()
    rec({
        "id": "exact-u-pvalues",
        "status": "DIFFERS",
        "summary": (
            f"Tied vs split inside one control gap: compute_de p {a['p_value']:.6e} → "
            f"{b['p_value']:.6e} (ratio {b['p_value'] / a['p_value']:.4f}), "
            f"lfc {a['log2_fold_change']:.6f} → {b['log2_fold_change']:.6f}. "
            f"scipy U {u1.statistic:.4f} → {u2.statistic:.4f}. "
            f"BH call flip search: {flip_m}. "
            f"postpass moved={moved['moved']} error={moved.get('error')}."
        ),
        "p_tied": float(a["p_value"]),
        "p_split": float(b["p_value"]),
        "p_ratio": float(b["p_value"] / a["p_value"]),
        "lfc_tied": float(a["log2_fold_change"]),
        "lfc_split": float(b["log2_fold_change"]),
        "U_tied": float(u1.statistic),
        "U_split": float(u2.statistic),
        "bh_flip": {"m": flip_m[0], "q_tied": flip_m[1], "q_split": flip_m[2]} if flip_m else None,
        "postpass": moved,
    })


def _try_postpass() -> dict:
    """A two-carrier block whose down gene starts as a within-block tie and can untie.

    Columns: g0 (the gene the pass should move), sub (control CPM 0, the payment bank),
    filler (in the gate, block lfc >= 0 so a down-only pass leaves it alone).
    """
    n_c = 16
    lib_c = 1_000_000
    ctrl = np.zeros((n_c, 3), dtype=np.int64)
    ctrl[:, 1] = 0
    ctrl[:4, 0] = 200_000
    ctrl[:, 2] = lib_c - ctrl[:, 0] - ctrl[:, 1]
    # carriers tied at CPM 1500 on g0, inside the control gap (0, 2e5).
    # filler is ~99% of each row so its block lfc stays non-negative.
    def row(lib, g0):
        sub = max(lib // 100, 20)  # payment bank, about 1%
        fill = lib - g0 - sub
        return np.array([g0, sub, fill], dtype=np.int64)
    block = np.vstack([
        row(2_000, 3),
        row(8_000, 12),
        row(2_000, 0),
        row(8_000, 0),
    ]).astype(np.int64)
    cc = scorer_control_block(ctrl.astype(np.float64), min_cpm=5.0, chunk=8)
    base = block_lfc(block.astype(np.float64), cc)
    try:
        new, diag = postpass(block, 2, cc, base, {"sel": "dn", "k": 2.5, "k_up": None, "mode": "zs"})
        pass_err = None
    except Exception as exc:  # noqa: BLE001 — the report records a failed realisation
        new, diag, pass_err = block, {"n_moved": 0}, f"{type(exc).__name__}: {exc}"
    changed = pass_err is None and not np.array_equal(new, block)
    # p of gene 0 through compute_de, before and after, embedded with the same controls
    genes = ["g0", "sub", "filler"]

    def p_of(blk):
        counts = np.vstack([ctrl, blk]).astype(np.float64)
        labels = [CTRL] * n_c + ["TGT"] * blk.shape[0]
        return pmap(de_table(counts, labels, genes), "TGT").get("g0")

    before = p_of(block)
    after = p_of(new) if changed else None
    # jackknife of the block alone (legal: every row has a positive library)
    jk_before = float(jackknife_correction(block.astype(np.float64), np.zeros(4, dtype=np.intp), 1, 5e4)[0])
    jk_after = float(jackknife_correction(new.astype(np.float64), np.zeros(4, dtype=np.intp), 1, 5e4)[0]) if changed else None
    bj_before = block_jackknife(block.astype(np.float64))
    bj_after = block_jackknife(new.astype(np.float64)) if changed else None
    flat = {k: (float(v) if isinstance(v, (int, float, np.floating, np.integer)) else v)
            for k, v in diag.items() if not str(k).startswith("_")}
    return {
        "moved": bool(changed),
        "error": pass_err,
        "diag": flat,
        "p_before": None if before is None else float(before["p_value"]),
        "p_after": None if after is None else float(after["p_value"]),
        "lfc_before": None if before is None else float(before["log2_fold_change"]),
        "lfc_after": None if after is None else float(after["log2_fold_change"]),
        "lfc_base": [float(v) for v in base],
        "keep": cc.keep.tolist(),
        "gene0_before": block[:, 0].tolist(),
        "gene0_after": new[:, 0].tolist(),
        "jk_before": jk_before,
        "jk_after": jk_after,
        "builder_jk_before": bj_before,
        "builder_jk_after": bj_after,
    }


# ---------------------------------------------------------------- bulk
def test_jackknife_and_pds() -> None:
    rng = np.random.default_rng(4)
    n, g = 25, 12
    libs = rng.integers(2_000, 8_000, size=n)
    fr = rng.dirichlet(np.full(g, 0.8), size=n)
    raw = np.floor(fr * (libs[:, None] * 0.9)).astype(np.float64)
    raw[:, 0] += libs - raw.sum(1)
    c_eval = float(jackknife_correction(raw, np.zeros(n, dtype=np.intp), 1, 50_000.0)[0])
    c_bldr = block_jackknife(raw)
    # same column sums, different allocation
    perm = raw.copy()
    perm[:-1] = raw[1:]
    perm[-1] = raw[0]
    # that rotates rows, column sums unchanged only if we swap a pair of cells' gene
    # compositions while keeping each row sum. Swap two rows' off-diagonal pattern:
    alt = raw.copy()
    alt[0, 1], alt[1, 1] = raw[1, 1], raw[0, 1]
    alt[0, 2], alt[1, 2] = raw[1, 2], raw[0, 2]
    # restore row sums by adjusting filler gene
    alt[:, -1] = 0
    alt[:, -1] = raw.sum(1) - alt.sum(1)
    c_alt = float(jackknife_correction(alt, np.zeros(n, dtype=np.intp), 1, 50_000.0)[0])
    # pds: three perts, constant prediction vs a scaled one. Panel exclusion.
    # Small counts.
    genes = [f"g{i}" for i in range(6)]
    # real: each pert knocks its own gene and lifts a private downstream gene
    def block(which):
        # Knock the pert's own gene and lift a private downstream gene (g3/g4/g5),
        # which panel exclusion does not remove.
        x = np.full((8, 6), 10.0)
        x[:, which] = 0
        x[:, 3 + which] = 80
        return x
    real = np.vstack([np.full((10, 6), 10.0), block(0), block(1), block(2)])
    labels_real = [CTRL] * 10 + ["g0"] * 8 + ["g1"] * 8 + ["g2"] * 8
    # prediction: same direction, and a count-scaled copy (composition unchanged)
    pred = real.copy()
    pred2 = real * 3
    a_real = ad.AnnData(real)
    a_real.var_names = genes
    a_real.obs["target"] = np.asarray(labels_real, dtype=object)
    a_pred = ad.AnnData(pred)
    a_pred.var_names = genes
    a_pred.obs["target"] = np.asarray(labels_real, dtype=object)
    a_pred2 = ad.AnnData(pred2)
    a_pred2.var_names = genes
    a_pred2.obs["target"] = np.asarray(labels_real, dtype=object)
    # constant prediction: paste control
    paste = np.vstack([np.full((10, 6), 10.0)] + [np.full((8, 6), 10.0) for _ in range(3)])
    a_paste = ad.AnnData(paste)
    a_paste.var_names = genes
    a_paste.obs["target"] = np.asarray(labels_real, dtype=object)
    kw = dict(pert_col="target", control=CTRL, distance="cosine",
              rank_denominator="n-1", tie_policy="midrank",
              exclude_target_gene=True, exclusion_scope="panel",
              control_source="real")
    s1 = discrimination_score(a_pred, a_real, **kw)
    s2 = discrimination_score(a_pred2, a_real, **kw)
    s0 = discrimination_score(a_paste, a_real, **kw)
    # row scope would still see other targets' genes; panel scope must not.
    kw_row = dict(kw)
    kw_row["exclusion_scope"] = "row"
    # A prediction whose only effect is on the other panel genes. Panel exclusion
    # removes every such coordinate; row exclusion leaves the other targets' genes in.
    cheat = np.vstack([np.full((10, 6), 10.0), block(0), block(1), block(2)])
    # wipe the downstream lift, keep the on-target knockdown, and knock every panel gene
    cheat[:, 3:] = 10
    cheat[10:, :3] = 0
    a_cheat = ad.AnnData(cheat)
    a_cheat.var_names = genes
    a_cheat.obs["target"] = np.asarray(labels_real, dtype=object)
    s_cheat = discrimination_score(a_cheat, a_real, **kw)
    s_cheat_row = discrimination_score(a_cheat, a_real, **kw_row)
    s_row = discrimination_score(a_paste, a_real, **kw_row)
    rec({
        "id": "jackknife-matches",
        "status": "CONFORMS" if abs(c_eval - c_bldr) / max(c_eval, 1e-30) < 1e-8 else "DIFFERS",
        "summary": (
            f"block_jackknife {c_bldr:.8e} vs jackknife_correction {c_eval:.8e} "
            f"(rel {abs(c_eval - c_bldr) / c_eval:.3e}). "
            f"Same column sums, two cells swapped: C {c_alt:.8e} "
            f"(Δ {c_alt - c_eval:.3e})."
        ),
        "builder": c_bldr,
        "cell_eval2": c_eval,
        "rel": abs(c_eval - c_bldr) / c_eval,
        "same_totals_other_allocation": c_alt,
    })
    rec({
        "id": "pds-scale-ties-panel",
        "status": "CONFORMS",
        "summary": (
            f"Count-scaling the prediction leaves PDS unchanged "
            f"({ {k: float(v) for k, v in s1.items()} } vs "
            f"{ {k: float(v) for k, v in s2.items()} }). "
            f"Control paste, panel exclusion, mean "
            f"{float(np.mean(list(s0.values()))):.4f}. "
            f"Same paste under row exclusion, mean "
            f"{float(np.mean(list(s_row.values()))):.4f}. "
            f"Panel-only cheat under panel exclusion, mean "
            f"{float(np.mean(list(s_cheat.values()))):.4f}; "
            f"under row exclusion, mean "
            f"{float(np.mean(list(s_cheat_row.values()))):.4f}."
        ),
        "pds_base": {k: float(v) for k, v in s1.items()},
        "pds_scaled": {k: float(v) for k, v in s2.items()},
        "pds_paste_panel_mean": float(np.mean(list(s0.values()))),
        "pds_paste_row_mean": float(np.mean(list(s_row.values()))),
        "pds_cheat_panel_mean": float(np.mean(list(s_cheat.values()))),
        "pds_cheat_row_mean": float(np.mean(list(s_cheat_row.values()))),
    })


def test_ballast_scaled_reallocation() -> None:
    """The 4-cell postpass move (gene0 3,12 → 1,14), padded to 400 cells.

    Column sums stay fixed. This is the jackknife the mse correction actually sees
    when ballast dominates; the 4-cell tripling is the same formula at n=4.
    """
    n = 400
    g = 8
    libs = np.full(n, 20_000, dtype=np.float64)
    libs[0] = 2_000
    libs[1] = 8_000
    before = np.zeros((n, g), dtype=np.float64)
    before[0, 0] = 3
    before[1, 0] = 12
    before[:, 1] = libs - before[:, 0]
    after = before.copy()
    after[0, 0] = 1
    after[1, 0] = 14
    after[:, 1] = libs - after[:, 0]
    assert before.sum() == after.sum()
    assert np.allclose(before.sum(0), after.sum(0))
    c0 = float(jackknife_correction(before, np.zeros(n, dtype=np.intp), 1, 50_000.0)[0])
    c1 = float(jackknife_correction(after, np.zeros(n, dtype=np.intp), 1, 50_000.0)[0])
    b0 = block_jackknife(before)
    b1 = block_jackknife(after)
    rec({
        "id": "exact-u-jackknife-at-400",
        "status": "DIFFERS" if abs(c1 - c0) > 1e-12 else "CONFORMS",
        "summary": (
            f"Same integer move as the postpass toy, n=400, column sums fixed. "
            f"jackknife_correction {c0:.8e} → {c1:.8e} (rel {(c1 - c0) / c0:.3e}). "
            f"block_jackknife {b0:.8e} → {b1:.8e}. "
            "Scaled mse on a real panel is not measured here."
        ),
        "c_before": c0,
        "c_after": c1,
        "rel": (c1 - c0) / c0,
        "builder_before": b0,
        "builder_after": b1,
    })


def _jackknife_proposed(block, ts: float = 50_000.0) -> float:
    """The formula in patches/block-jackknife-zero-remainder.patch. Not production code."""
    x = np.asarray(block, dtype=np.float64)
    p = x.sum(0)
    s = float(p.sum())
    lib = x.sum(1)
    live = p > 0
    r = s - lib
    safe = np.where(r > 0.0, r, 1.0)
    v = np.log1p(ts * (p[live][None, :] - x[:, live]) / safe[:, None])
    v = np.where((r > 0.0)[:, None], v, 0.0)
    b = np.log1p(ts * p[live] / (s if s > 0.0 else 1.0))
    dev = v - b
    n = x.shape[0]
    return float((n - 1) / n * ((dev - dev.mean(0)) ** 2).sum())


def test_jackknife_empty_remainder() -> None:
    """One cell holds every count; the other is empty. Scorer sets that row's v to 0."""
    raw = np.array([[5.0, 0.0], [0.0, 0.0]], dtype=np.float64)
    c_eval = float(jackknife_correction(raw, np.zeros(2, dtype=np.intp), 1, 50_000.0)[0])
    try:
        c_bldr = block_jackknife(raw)
        err = None
    except Exception as exc:  # noqa: BLE001
        c_bldr = None
        err = repr(exc)
    finite = c_bldr is not None and np.isfinite(c_bldr)
    rec({
        "id": "jackknife-zero-remainder",
        "status": "DIFFERS" if (not finite or abs(c_bldr - c_eval) / max(c_eval, 1e-30) > 1e-6) else "CONFORMS",
        "summary": (
            f"A cell that holds every count: cell_eval2 C={c_eval:.6e}, "
            f"block_jackknife={c_bldr}, err={err}. "
            "Legal 400-cell blocks do not hit this; the formulas differ on it."
        ),
        "cell_eval2": c_eval,
        "builder": c_bldr,
        "error": err,
    })
    proposed = _jackknife_proposed(raw)
    rng = np.random.default_rng(4)
    n, g = 25, 12
    libs = rng.integers(2_000, 8_000, size=n)
    fr = rng.dirichlet(np.full(g, 0.8), size=n)
    legal = np.floor(fr * (libs[:, None] * 0.9)).astype(np.float64)
    legal[:, 0] += libs - legal.sum(1)
    legal_eval = float(jackknife_correction(legal, np.zeros(n, dtype=np.intp), 1, 50_000.0)[0])
    legal_prop = _jackknife_proposed(legal)
    rec({
        "id": "jackknife-proposed-formula",
        "status": "CONFORMS" if abs(proposed - c_eval) < 1e-9 and abs(legal_prop - legal_eval) / legal_eval < 1e-12 else "DIFFERS",
        "summary": (
            f"Centered r==0 formula: edge {proposed:.6e} vs cell_eval2 {c_eval:.6e}; "
            f"legal 25x12 rel {abs(legal_prop - legal_eval) / legal_eval:.3e}."
        ),
        "edge_proposed": proposed,
        "edge_cell_eval2": c_eval,
        "legal_rel": abs(legal_prop - legal_eval) / legal_eval,
    })


def main() -> None:
    test_versions()
    test_cpm_gate_definition()
    test_gate_threshold_is_strict()
    test_reference_only_gate()
    test_normalization_association()
    test_gate_disagreement_rate()
    test_wilcoxon_lfc_bh()
    test_bh_scope_and_prefilter()
    test_epsilon_and_zero_mean()
    test_reach_order()
    test_fid_undirected()
    test_nmae_and_jaccard()
    test_exact_u_pvalues()
    test_ballast_scaled_reallocation()
    test_jackknife_and_pds()
    test_jackknife_empty_remainder()
    OUT.write_text(json.dumps(FINDINGS, indent=2, default=str) + "\n")
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
