"""222 mandate falsifiers 2 and 3 — the delivered profile, split by carrying row.

`experiments/210-cd4-build/profile210.py` verbatim, with ONE section added at the end (marked
`222 ADDITION`): the per-block digest, cosine and norm ratio split by whether the H1 union
carries that target.  210's own falsifier numbering is left untouched in the body below --
its "FALSIFIER 3 (digest)" is this mandate's falsifier 2 and its "FALSIFIER 4 (non-vacuity)"
is this mandate's falsifier 3.

    uv run python profile222.py <arm.h5ad> <210-cd4f034.h5ad>

Original 210 header follows.

210 falsifiers 3 and 4 — the delivered profile, and the NON-VACUITY norm.

This is the mirror image of 191/202's `digest191.py`.  Those two arms claimed the delivered
per-gene panel totals were IDENTICAL to the incumbent's on 900/900 blocks and used the
nonzero-structure digest to prove the identity was not vacuous.  **This arm claims the
opposite**: `VCC_CD4_TARGET_F` 0.10 -> 0.34 moves the post-whitening effect, so the delivered
profile MUST move, and a totals digest that MATCHES on any block would mean the residual did
not apply and the arm is void (mandate falsifier 3).

The delivered effect is reconstructed the way the scorer and the builder both read it
(`probe191.py`, `write_submission_from_effects`):

    p_hat  = per-gene panel total / sum of panel totals          (the emitted profile)
    E      = log1p(TS * p_hat) - b_ctrl,   b_ctrl = log1p(TS * pooled_ctrl / sum)

`E` is the delivered effect row -- `alpha * gated_effect` plus the renderer's mass shift --
and it is the object 207 sections 6.1/6.2 measured when it predicted

    delivered median row norm   1.0009x   (six207) / 0.99987x (amp207's median instrument)
    delivered A = mean|row|^2   0.990x    (amp207)
    ||delivered(0.34) - delivered(0)|| / ||delivered(0)||   0.3312, against 0.1023 at f = 0.10

so the arm-vs-incumbent displacement predicted here is 0.3312 - 0.1023 = 0.229 of ||d(0)||,
i.e. ~0.228 relative to the f = 0.10 artifact if the residual is collinear across doses (it
is: `fraction` enters `add_post_whiten_residual` linearly, only the preserve scalar is not).

NON-VACUITY (mandate falsifier 4).  200 and 207 both found the CD4 residual can be identically
zero while every other check still passes, because the sizing rule takes a MEDIAN over rows.
So the headline number is the relative displacement, and -- more decisive than the mean -- the
MINIMUM of the per-block displacement over all 900 blocks.  A single silently-zero block shows
up there and nowhere else.

    uv run python profile210.py <arm.h5ad> <incumbent.h5ad>
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = Path(__file__).resolve().parent
BUNDLE = ROOT / "data/challenge"
N_GENES = 18533
SCORER_MIN_CPM = 5.0
BUILDER_MIN_CPM = 10.0
sys.path.insert(0, str(ROOT / "src"))
from vcc2026.donor_transfer import TS  # noqa: E402  the builder's own pseudobulk scale


def read_obs_blocks(path: Path):
    """(context, target_gene) per row block, straight off the categorical encoding."""
    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        out = {}
        for name in ("context", "target_gene"):
            group = obs[name]
            cats = np.asarray(group["categories"][:]).astype(str)
            codes = np.asarray(group["codes"][:])
            out[name] = cats[codes]
        indptr = np.asarray(handle["X"]["indptr"][:])
        shape = tuple(handle["X"].attrs["shape"])
    labels = list(zip(out["context"], out["target_gene"]))
    blocks, start = [], 0
    for i in range(1, len(labels) + 1):
        if i == len(labels) or labels[i] != labels[start]:
            blocks.append((labels[start], start, i))
            start = i
    return blocks, indptr, shape


def scan(path: Path, blocks, indptr, tag: str):
    """One streaming pass: per-block totals, digests, per-cell CPM mean, library stats."""
    n = len(blocks)
    totals = np.zeros((n, N_GENES), dtype=np.float64)
    cpm_mean = np.zeros((n, N_GENES), dtype=np.float64)
    totals_sha, struct_sha, totals_sum, cell_stats = [], [], [], []
    t0 = time.time()
    with h5py.File(path, "r") as handle:
        data = handle["X"]["data"]
        indices = handle["X"]["indices"]
        for k, (_key, lo_row, hi_row) in enumerate(blocks):
            lo, hi = int(indptr[lo_row]), int(indptr[hi_row])
            cols = indices[lo:hi].astype(np.int64)
            vals = data[lo:hi].astype(np.float64)
            row_ptr = indptr[lo_row:hi_row + 1] - lo
            libraries = np.add.reduceat(vals, row_ptr[:-1])
            libraries[row_ptr[1:] == row_ptr[:-1]] = 0.0
            block_totals = np.bincount(cols, weights=vals, minlength=N_GENES)
            totals[k] = block_totals
            totals_sha.append(
                hashlib.sha256(block_totals.astype(np.int64).tobytes()).hexdigest()[:16]
            )
            totals_sum.append(int(block_totals.sum()))
            struct_sha.append(hashlib.sha256(
                np.bincount(cols, minlength=N_GENES).astype(np.int64).tobytes()
            ).hexdigest()[:16])
            safe = np.where(libraries > 0, libraries, 1.0)
            weights = vals * (1e6 / np.repeat(safe, np.diff(row_ptr)))
            cpm_mean[k] = np.bincount(cols, weights=weights, minlength=N_GENES) / (hi_row - lo_row)
            cell_stats.append((float(libraries.min()), float(libraries.max()),
                               float(libraries.std())))
            if (k + 1) % 150 == 0:
                print(f"  [{tag}] {k + 1}/{n} [{time.time() - t0:.0f}s]", flush=True)
    return totals, cpm_mean, totals_sha, struct_sha, totals_sum, np.asarray(cell_stats)


def control_frames():
    """`b_ctrl` and the two gates per context, exactly as `probe191.py` builds them."""
    frames = {}
    for context in ("A", "B", "C"):
        controls = ad.read_h5ad(BUNDLE / f"context_{context}.h5ad")
        matrix = controls.X
        pooled = np.asarray(matrix.sum(0), dtype=np.float64).ravel()
        n_cells = matrix.shape[0]
        cpm = np.zeros(N_GENES, dtype=np.float64)
        block = matrix.tocsr()
        for lo in range(0, n_cells, 2000):
            dense = np.asarray(block[lo:min(lo + 2000, n_cells)].todense(), dtype=np.float64)
            sums = dense.sum(1)
            sums[sums == 0] = 1.0
            cpm += (dense * (1e6 / sums[:, None])).sum(0)
        cpm /= n_cells
        del controls, matrix, block
        b_ctrl = np.log1p(TS * pooled / pooled.sum())
        expression = np.expm1(b_ctrl)
        gate_cpm = 1e6 * expression / max(float(expression.sum()), 1e-12)
        frames[context] = dict(
            b_ctrl=b_ctrl,
            scorer_gate=cpm > SCORER_MIN_CPM,
            builder_gate=(gate_cpm > BUILDER_MIN_CPM) & (b_ctrl > 0),
        )
    return frames


def delivered(totals: np.ndarray, keys, frames) -> np.ndarray:
    """E = log1p(TS * p_hat) - b_ctrl, one row per (context, target) block."""
    out = np.empty_like(totals)
    for k, (context, _target) in enumerate(keys):
        p = totals[k] / totals[k].sum()
        out[k] = np.log1p(TS * p) - frames[context]["b_ctrl"]
    return out


def main() -> None:
    arm_path, base_path = Path(sys.argv[1]), Path(sys.argv[2])
    a_blocks, a_indptr, a_shape = read_obs_blocks(arm_path)
    b_blocks, b_indptr, b_shape = read_obs_blocks(base_path)
    print(f"arm        {arm_path}  shape {a_shape}  {len(a_blocks)} blocks")
    print(f"incumbent  {base_path}  shape {b_shape}  {len(b_blocks)} blocks")
    keys = [k for k, _, _ in a_blocks]
    if keys != [k for k, _, _ in b_blocks]:
        raise SystemExit("block (context, target) ORDER differs -- comparison meaningless")
    sizes = {hi - lo for _, lo, hi in a_blocks} | {hi - lo for _, lo, hi in b_blocks}
    print(f"block sizes {sizes} (must be {{400}})")

    a_tot, a_cpm, a_tsha, a_ssha, a_tsum, a_cell = scan(arm_path, a_blocks, a_indptr, "arm")
    b_tot, b_cpm, b_tsha, b_ssha, b_tsum, b_cell = scan(base_path, b_blocks, b_indptr, "inc")

    same_tot = int(np.sum([x == y for x, y in zip(a_tsha, b_tsha)]))
    same_str = int(np.sum([x == y for x, y in zip(a_ssha, b_ssha)]))
    same_sum = int(np.sum(np.asarray(a_tsum) == np.asarray(b_tsum)))
    n = len(keys)
    print()
    print(f"per-gene PANEL TOTAL digest identical on   {same_tot}/{n}"
          "   (MUST be 0 -- this arm moves the profile)")
    print(f"panel total SUM identical on               {same_sum}/{n}")
    print(f"nonzero-STRUCTURE digest identical on      {same_str}/{n}")
    print("\nsample rows:")
    for i in list(range(3)) + [n - 2, n - 1]:
        print(f"  {keys[i]}  totals {a_tsha[i]} vs {b_tsha[i]}   "
              f"struct {a_ssha[i]} vs {b_ssha[i]}")
    f3_digest = same_tot == 0
    print(f"\nFALSIFIER 3 (digest) {'PASS' if f3_digest else 'FAIL'}: "
          f"{n - same_tot}/{n} delivered pseudobulk digests DIFFER")

    print("\nPER-CELL LIBRARY TOTALS (the weld) and the 1,000,000 format cap")
    print(f"  {'arm':<12}{'min':>10}{'max':>10}{'sd median':>12}{'sd max':>10}")
    for tag, stats in (("incumbent", b_cell), ("arm", a_cell)):
        print(f"  {tag:<12}{stats[:, 0].min():>10.0f}{stats[:, 1].max():>10.0f}"
              f"{np.median(stats[:, 2]):>12.3f}{stats[:, 2].max():>10.3f}")
    cap_ok = bool(a_cell[:, 1].max() <= 1_000_000)
    print(f"  max cell total {a_cell[:, 1].max():.0f} <= 1,000,000 cap: "
          f"{'PASS' if cap_ok else 'FAIL'}")

    print("\nDELIVERED EFFECT (falsifiers 3 and 4) -- E = log1p(TS * p_hat) - b_ctrl, "
          f"TS = {TS:g}")
    frames = control_frames()
    for context in ("A", "B", "C"):
        print(f"  ctx {context}: {int(frames[context]['scorer_gate'].sum())} genes past the "
              f"{SCORER_MIN_CPM}-CPM scorer gate, "
              f"{int(frames[context]['builder_gate'].sum())} past the builder's "
              f"{BUILDER_MIN_CPM}-CPM gate")
    a_eff = delivered(a_tot, keys, frames)
    b_eff = delivered(b_tot, keys, frames)

    a_rows = np.linalg.norm(a_eff, axis=1)
    b_rows = np.linalg.norm(b_eff, axis=1)
    d_rows = np.linalg.norm(a_eff - b_eff, axis=1)
    rel_global = float(np.linalg.norm(a_eff - b_eff) / np.linalg.norm(b_eff))
    rel_block = d_rows / b_rows
    A_arm = float(np.mean(a_rows ** 2))
    A_base = float(np.mean(b_rows ** 2))

    print(f"\n  {'quantity':<44}{'incumbent f=0.10':>20}{'arm f=0.34':>16}{'ratio':>12}")
    print(f"  {'median delivered |row|':<44}{np.median(b_rows):>20.6f}"
          f"{np.median(a_rows):>16.6f}{np.median(a_rows) / np.median(b_rows):>12.5f}")
    print(f"  {'mean delivered |row|':<44}{b_rows.mean():>20.6f}"
          f"{a_rows.mean():>16.6f}{a_rows.mean() / b_rows.mean():>12.5f}")
    print(f"  {'A = mean |row|^2':<44}{A_base:>20.6f}{A_arm:>16.6f}"
          f"{A_arm / A_base:>12.5f}")
    print(f"\n  NON-VACUITY  ||E_arm - E_inc||_F / ||E_inc||_F = {rel_global:.6f}")
    print(f"    per-block relative displacement:  min {rel_block.min():.6f}   "
          f"p01 {np.percentile(rel_block, 1):.6f}   median {np.median(rel_block):.6f}   "
          f"max {rel_block.max():.6f}")
    print(f"    blocks with displacement exactly 0: "
          f"{int(np.sum(d_rows == 0.0))}/{n}   "
          f"blocks below 1e-6 relative: {int(np.sum(rel_block < 1e-6))}/{n}")
    for context in ("A", "B", "C"):
        sel = np.array([k[0] == context for k in keys])
        print(f"    ctx {context}: ||dE||/||E|| = "
              f"{np.linalg.norm(a_eff[sel] - b_eff[sel]) / np.linalg.norm(b_eff[sel]):.6f}"
              f"   median row-norm ratio "
              f"{np.median(a_rows[sel]) / np.median(b_rows[sel]):.5f}")
    # The same three quantities restricted to the builder's own CPM gate.  Outside the gate
    # the effect is zero by construction, so E there is nothing but the emitter's integer
    # quantisation -- common-mode between the arms in energy but pure noise in the difference.
    # The gated row is therefore the sharper instrument for a 0.09% median-norm claim, and
    # both are printed rather than one being chosen after the fact.
    gate_rows = np.stack([frames[k[0]]["builder_gate"] for k in keys])
    ag = np.where(gate_rows, a_eff, 0.0)
    bg = np.where(gate_rows, b_eff, 0.0)
    ag_rows, bg_rows = np.linalg.norm(ag, axis=1), np.linalg.norm(bg, axis=1)
    dg_rows = np.linalg.norm(ag - bg, axis=1)
    print("\n  the same three, restricted to the builder's 10-CPM gate:")
    print(f"  {'median delivered |row|':<44}{np.median(bg_rows):>20.6f}"
          f"{np.median(ag_rows):>16.6f}{np.median(ag_rows) / np.median(bg_rows):>12.5f}")
    print(f"  {'A = mean |row|^2':<44}{np.mean(bg_rows ** 2):>20.6f}"
          f"{np.mean(ag_rows ** 2):>16.6f}"
          f"{np.mean(ag_rows ** 2) / np.mean(bg_rows ** 2):>12.5f}")
    print(f"  {'||dE||/||E||':<44}{'':>20}"
          f"{np.linalg.norm(ag - bg) / np.linalg.norm(bg):>16.6f}")
    print(f"    per-block, gated: min {(dg_rows / bg_rows).min():.6f}   "
          f"median {np.median(dg_rows / bg_rows):.6f}   "
          f"max {(dg_rows / bg_rows).max():.6f}")

    f4 = bool(rel_block.min() > 0.01)
    print(f"\nFALSIFIER 4 (non-vacuity) {'PASS' if f4 else 'FAIL'}: every one of the {n} "
          f"blocks moved (min relative displacement {rel_block.min():.6f})")

    print("\nCALL-SET EXPOSURE — how far the |H| head of `rc_callset`'s union can have moved")
    # `rc_callset` prescribes `top 400 by |H|` UNION `top 2,200 by control CPM`.  The CPM half
    # is a function of the CONTROL cells alone and cannot move with the dose, so the whole
    # call-set exposure of this arm is the head.  `|delivered effect|` inside the builder's
    # own CPM gate is the same ordering key up to alpha and the renderer's mass shift, so the
    # head overlap measured here bounds the call-set jaccard from below:
    #     |symmetric difference| <= 2 * (400 - overlap)
    head_overlap = np.empty(n, dtype=np.int64)
    for k, (context, _t) in enumerate(keys):
        gate = frames[context]["builder_gate"]
        idx = np.flatnonzero(gate)
        a_head = idx[np.argsort(-np.abs(a_eff[k][gate]), kind="stable")[:400]]
        b_head = idx[np.argsort(-np.abs(b_eff[k][gate]), kind="stable")[:400]]
        head_overlap[k] = np.intersect1d(a_head, b_head, assume_unique=True).size
    moved_head = 400 - head_overlap
    jac_lower = (2600.0 - moved_head) / (2600.0 + moved_head)
    print(f"  top-400 |delivered effect| head overlap with the incumbent: "
          f"min {head_overlap.min()}   median {np.median(head_overlap):.0f}   "
          f"max {head_overlap.max()}   of 400")
    print(f"  implied call-set jaccard lower bound: min {jac_lower.min():.4f}   "
          f"median {np.median(jac_lower):.4f}")

    print("\nDELIVERED log2 FOLD CHANGE (the object `nmae` reads), reference-free")
    per_block = []
    for k, (context, _t) in enumerate(keys):
        gate = frames[context]["scorer_gate"]
        delta = np.log2((a_cpm[k][gate] + 1e-9) / (b_cpm[k][gate] + 1e-9))
        per_block.append((float(np.abs(delta).max()), float(np.abs(delta).mean()),
                          float(delta.mean())))
    per_block = np.asarray(per_block)
    print(f"  |d lfc| over the DE universe: max {per_block[:, 0].max():.6f}   "
          f"mean {per_block[:, 1].mean():.6f}   median-of-block-means "
          f"{np.median(per_block[:, 1]):.6f}")
    print(f"  signed mean d lfc per block: median {np.median(per_block[:, 2]):+.6f}   "
          f"range [{per_block[:, 2].min():+.6f}, {per_block[:, 2].max():+.6f}]")

    report = dict(
        arm=str(arm_path), incumbent=str(base_path), blocks=n,
        totals_identical=same_tot, sums_identical=same_sum, structure_identical=same_str,
        cap_ok=cap_ok, max_cell_total=float(a_cell[:, 1].max()),
        cell_totals_arm=[float(a_cell[:, 0].min()), float(a_cell[:, 1].max()),
                         float(np.median(a_cell[:, 2]))],
        cell_totals_incumbent=[float(b_cell[:, 0].min()), float(b_cell[:, 1].max()),
                               float(np.median(b_cell[:, 2]))],
        median_row_norm_incumbent=float(np.median(b_rows)),
        median_row_norm_arm=float(np.median(a_rows)),
        median_row_norm_ratio=float(np.median(a_rows) / np.median(b_rows)),
        mean_row_norm_ratio=float(a_rows.mean() / b_rows.mean()),
        A_incumbent=A_base, A_arm=A_arm, A_ratio=float(A_arm / A_base),
        gated_median_row_norm_ratio=float(np.median(ag_rows) / np.median(bg_rows)),
        gated_A_ratio=float(np.mean(ag_rows ** 2) / np.mean(bg_rows ** 2)),
        gated_nonvacuity_relative_norm=float(
            np.linalg.norm(ag - bg) / np.linalg.norm(bg)),
        gated_nonvacuity_block_min=float((dg_rows / bg_rows).min()),
        gated_nonvacuity_block_median=float(np.median(dg_rows / bg_rows)),
        nonvacuity_relative_norm=rel_global,
        nonvacuity_block_min=float(rel_block.min()),
        nonvacuity_block_p01=float(np.percentile(rel_block, 1)),
        nonvacuity_block_median=float(np.median(rel_block)),
        nonvacuity_block_max=float(rel_block.max()),
        blocks_exactly_zero=int(np.sum(d_rows == 0.0)),
        head_overlap_min=int(head_overlap.min()),
        head_overlap_median=float(np.median(head_overlap)),
        head_overlap_max=int(head_overlap.max()),
        callset_jaccard_lower_min=float(jac_lower.min()),
        callset_jaccard_lower_median=float(np.median(jac_lower)),
        dlfc_max=float(per_block[:, 0].max()),
        dlfc_mean=float(per_block[:, 1].mean()),
        dlfc_signed_block_mean_median=float(np.median(per_block[:, 2])),
        falsifier3_digest=bool(f3_digest), falsifier4_nonvacuity=f4,
    )

    # -------------------------------------------------------------- 222 ADDITION
    # MANDATE FALSIFIER 2: "the profile moved on the CARRYING rows and nowhere else."
    #
    # The literal form of that -- a digest identical off the carrying rows -- is not
    # achievable for THIS arm and it is worth being exact about why.  218's construction
    # merges the two residuals into ONE add followed by production's ONE global preserve
    # rescale, so the preserve scalar is a function of both and every row is multiplied by a
    # slightly different constant than the incumbent's.  Holding CD4's own preserve fixed
    # instead would leave the delivered amplitude free to move, which is the one thing an
    # alpha-calibrated build must not do (187: alpha 0.982 sits on context A's rho cliff).
    #
    # So the sharp statement is the SCALE-CORRECTED one: off the carrying rows the arm's
    # delivered effect must equal the incumbent's times a SINGLE shared scalar -- i.e. the
    # per-block cosine must be ~1 and the per-block norm ratio must be the same constant for
    # every non-carrying block -- while the carrying blocks must move in DIRECTION.
    union = np.load(ROOT / "experiments/222-h1-union-build/work/union-on028.npz")
    panel = [line.split(",")[0] for line in
             (BUNDLE / "pert_counts.csv").read_text().splitlines()[1:] if line.strip()]
    carrying_targets = {panel[i] for i in
                        np.flatnonzero(np.asarray(union["target_condition_count"]) > 0)}
    carry = np.asarray([t in carrying_targets for _c, t in keys])
    cos_block = ((a_eff * b_eff).sum(1)
                 / (np.linalg.norm(a_eff, axis=1) * np.linalg.norm(b_eff, axis=1)))
    ratio_block = a_rows / b_rows
    print("\n" + "=" * 100)
    print(f"MANDATE FALSIFIER 2 — the union carries {len(carrying_targets)}/300 targets, "
          f"so {int(carry.sum())}/{n} of the 900 (context, target) blocks")
    print("=" * 100)
    print(f"  {'':16s}{'n':>5s}{'digest same':>13s}{'cos min':>12s}{'cos median':>13s}"
          f"{'norm-ratio min':>16s}{'norm-ratio max':>16s}")
    for tag, sel in (("carrying", carry), ("NOT carrying", ~carry)):
        same = int(np.sum([x == y for x, y, s in zip(a_tsha, b_tsha, sel) if s]))
        print(f"  {tag:16s}{int(sel.sum()):5d}{same:13d}{cos_block[sel].min():12.8f}"
              f"{np.median(cos_block[sel]):13.8f}{ratio_block[sel].min():16.8f}"
              f"{ratio_block[sel].max():16.8f}")
    # The shared scalar, read off the non-carrying blocks themselves, then divided out.
    shared = float(np.median(ratio_block[~carry]))
    corrected = np.linalg.norm(a_eff - shared * b_eff, axis=1) / b_rows
    print(f"\n  shared delivered-scale factor (median norm ratio off the carrying rows): "
          f"{shared:.8f}")
    print(f"  {'':16s}{'n':>5s}{'||dE - s*E||/||E||  min':>26s}{'median':>12s}{'max':>12s}")
    for tag, sel in (("carrying", carry), ("NOT carrying", ~carry)):
        print(f"  {tag:16s}{int(sel.sum()):5d}{corrected[sel].min():26.8f}"
              f"{np.median(corrected[sel]):12.8f}{corrected[sel].max():12.8f}")
    separated = bool(corrected[carry].min() > corrected[~carry].max())
    print(f"\n  carrying-vs-non-carrying separation on the scale-corrected displacement: "
          f"{'CLEAN' if separated else 'OVERLAPPING'} "
          f"(min carrying {corrected[carry].min():.6f} vs max non-carrying "
          f"{corrected[~carry].max():.6f}, ratio "
          f"{corrected[carry].min() / max(corrected[~carry].max(), 1e-12):.1f}x)")
    report.update(
        f2_union_carrying_targets=len(carrying_targets),
        f2_carrying_blocks=int(carry.sum()),
        f2_digest_same_carrying=int(np.sum(
            [x == y for x, y, s in zip(a_tsha, b_tsha, carry) if s])),
        f2_digest_same_non_carrying=int(np.sum(
            [x == y for x, y, s in zip(a_tsha, b_tsha, ~carry) if s])),
        f2_cos_carrying=[float(cos_block[carry].min()),
                         float(np.median(cos_block[carry])),
                         float(cos_block[carry].max())],
        f2_cos_non_carrying=[float(cos_block[~carry].min()),
                             float(np.median(cos_block[~carry])),
                             float(cos_block[~carry].max())],
        f2_shared_scale=shared,
        f2_corrected_carrying=[float(corrected[carry].min()),
                               float(np.median(corrected[carry])),
                               float(corrected[carry].max())],
        f2_corrected_non_carrying=[float(corrected[~carry].min()),
                                   float(np.median(corrected[~carry])),
                                   float(corrected[~carry].max())],
        f2_separated=separated,
    )
    (HERE / "work" / "profile222.json").write_text(json.dumps(report, indent=1))
    print(f"\nwrote {HERE / 'work' / 'profile222.json'}")


if __name__ == "__main__":
    main()
