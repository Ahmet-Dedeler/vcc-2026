"""445 — DLD-1 (223's one live source) in the BULK union, re-priced on the current bulk under ballast.

229 killed every expressible DLD-1 construction on nmae (the welded DE channel paid for the bulk change). Since 410b the
bulk is unwelded from the DE channel (ballast), so only pds/mse read a bulk composition change; 418/420 re-priced the
composition knobs (H1 dose, CD4 dose) but never an extra union source. This prices the builder-expressible construction:
`VCC_UNION_SOURCES=gse295214_h1,h1_2025,dld1_lowmoi` (one union, per-gene mean over carrying sources, carrying-row sizing),
inside the stack (CD4 x1.5) at H1-union dose g, then 410's key K3 and sharpening q 1.5, 414's ballast realisation.

KNOWN ANSWER FIRST: the H1-pair union at g 0.56 (= comp420 `stack|K3+sh1.5`) must reproduce with |d| = 0.
CONTROL: DLD-1's residual rows permuted over targets (same row set size and norms, alignment destroyed) -- must lose.

    uv run python comp445.py <fold>
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/420-bulk-stack"))

from lib420 import (ALPHA, COMPS, E062, FOLDS, GAMMA, KMAXES, LAM, PROD_OPS, UNION, WORK as W420, Donors, Operator,  # noqa: E402
                    betas_for, build_fold, cd4_frozen_and_partial, compose_arm, de_profiles, deflate, design_and_realised,
                    key_scale, pre_normrel, shared_energy_key, source_residual, union_rows, weights)
from comp420 import cmp, score_rows  # noqa: E402
import sources218  # noqa: E402

sources218.SOURCES["dld1_lowmoi"] = dict(path="experiments/223-silent-zero-audit/work/tables/dld1_lowmoi",
                                         platform="DLD-1 low-MOI CRISPRi", role="223's live source")
TOL = 1e-12
OUT = HERE / "work"


def permute_rows(src, seed):
    """Row-permuted copy of a residual (targets shuffled among the CARRYING rows only)."""
    counts = np.asarray(src.target_condition_count)
    carry = np.flatnonzero(counts > 0)
    perm = carry[np.random.default_rng(seed).permutation(carry.size)]
    eff = np.asarray(src.mean_effect).copy()
    votes = np.asarray(src.gene_votes).copy()
    eff[carry], votes[carry] = eff[perm], votes[perm]
    return E062.CD4Prepared(mode=src.mode + ":perm", mean_effect=eff, gene_votes=votes,
                            target_condition_count=counts.copy(), condition_details=[])


def compose_sep(whitened, cd4, partial, h1, extra, wcols, panel_positions, *, m_cd4, g_h1, g_extra, g_part=0.5):
    """lib418.compose_arm plus ONE more post-whiten residual (`extra`, e.g. DLD-1) sized on ITS carrying rows at its own
    dose -- the construction 223 priced live and 229 found inexpressible in a single-union builder."""
    from lib418 import CD4_F, solve_fraction, solve_fraction_carrying
    f_cd4, _ = solve_fraction(whitened, cd4, CD4_F * m_cd4)
    merged = np.where(np.asarray(cd4.gene_votes) > 0, np.asarray(cd4.mean_effect), 0.0) * f_cd4
    support = np.asarray(cd4.gene_votes) > 0
    fr = {"cd4": float(f_cd4)}
    for tag, res_, g in (("h1_union", h1, g_h1), ("cd4_partial", partial, g_part), ("extra", extra, g_extra)):
        if not g:
            continue
        f, _ = solve_fraction_carrying(whitened, res_, g)
        v = np.asarray(res_.gene_votes) > 0
        merged += np.where(v, np.asarray(res_.mean_effect), 0.0) * f
        support |= v
        fr[tag] = float(f)
    added = E062.CD4Prepared(mode="merged", mean_effect=merged, gene_votes=support.astype(np.uint8),
                             target_condition_count=np.asarray(cd4.target_condition_count), condition_details=[])
    final, info = E062.add_post_whiten_residual(whitened, added, 1.0, wcols)
    final[:, panel_positions] = 0.0
    return final, dict(fractions=fr)


def run(held: str) -> None:
    t0 = time.time()
    log = lambda *a: print(f"[{held}]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    fold = build_fold(held)
    fold["scored"] = np.asarray(fold["scored"])
    scored = fold["scored"]
    pp_pos = fold["panel_positions"]
    wcols = np.asarray(fold["wcols"])
    truth_full = np.asarray(fold["truth"], dtype=np.float64)
    truth = truth_full[:, scored]
    pre_cached = np.asarray(fold["pre"], dtype=np.float64)
    targets = list(fold["targets"])
    sources = list(FOLDS[held]["sources"])
    op = Operator(PROD_OPS / "whitener_union_pruned.npz" if FOLDS[held]["op_kind"] == "production" else FOLDS[held]["op"])
    nonpanel = np.ones(truth_full.shape[1], bool)
    nonpanel[pp_pos] = False
    betas = betas_for(held, nonpanel)
    pair = [source_residual(n, targets, "gene_norm") for n in UNION]
    dld1 = source_residual("dld1_lowmoi", targets, "gene_norm")
    n_dld1 = int((np.asarray(dld1.target_condition_count) > 0).sum())
    n_pair = int((np.asarray(union_rows(pair).target_condition_count) > 0).sum())
    dn = Donors(targets, sources)
    frozen, partial = cd4_frozen_and_partial(targets)
    cd4, _ = deflate(frozen, 5)
    pre228 = dn.pre(sources, False)
    pre228[:, pp_pos] = 0.0
    nocov = sum(dn.A[s].astype(int) for s in sources) == 0
    pre228[nocov] = pre_cached[nocov]
    pre = pre_normrel(dn, sources, "cap", -1.5, fallback_rows=pre228.copy())
    whit = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    finalPP, _ = compose_arm(whit, cd4, partial, union_rows(pair), wcols, pp_pos)
    key410, carried410, _ = shared_energy_key(targets, sources, op, wcols, cols=None)
    wK3 = weights(key410, carried410, 0.3, 3.0)[0]
    P_DE = de_profiles(finalPP)
    ref = json.loads((W420 / f"comp420-{held}.json").read_text())["arms"]
    res = dict(fold=held, n_targets=len(targets), n_carry_pair=n_pair, n_carry_dld1=n_dld1, ka={}, arms={})
    log(f"targets {len(targets)}: H1-pair carries {n_pair}, DLD-1 carries {n_dld1}, "
        f"union carries {int((np.asarray(union_rows(pair + [dld1]).target_condition_count) > 0).sum())}")
    m_cd4, g_stack = COMPS["stack"]

    def rec(name, h1res, g, check_ref=None, extra=None, g_extra=0.0):
        if extra is None:
            final, parts = compose_arm(whit, cd4, partial, h1res, wcols, pp_pos, m_cd4=m_cd4, g_h1=g)
        else:
            final, parts = compose_sep(whit, cd4, partial, h1res, extra, wcols, pp_pos, m_cd4=m_cd4, g_h1=g, g_extra=g_extra)
        bulk = key_scale(final, wK3)[0]
        dr = design_and_realised(bulk, P_DE, q=1.5, kmaxes=KMAXES)
        sd, _ = score_rows({c: dr[c]["design"] for c in "ABC"}, scored, truth, betas, nonpanel)
        entry = dict(g=g, fractions=parts["fractions"], design=sd, realised={})
        for k in KMAXES:
            sr, _ = score_rows({c: dr[c]["realised"][k] for c in "ABC"}, scored, truth, betas, nonpanel)
            entry["realised"][str(k)] = sr
        if check_ref is not None:
            r = ref[check_ref]
            d = max([cmp(sd, r["design"])] + [cmp(entry["realised"][k], r["realised"][k]) for k in entry["realised"]])
            res["ka"][name] = d
            log(f"KA {name} vs comp420 {check_ref}: worst |d| {d:.1e}")
            assert d < TOL, f"KNOWN ANSWER FAILED {held} {name}: {d}"
        res["arms"][name] = entry
        rr = entry["realised"]["16"]
        log(f"{name:26s} fr {parts['fractions'].get('h1_union', 0):.4f} design cos {sd['cos_agg']:.6f} pds@.5 {sd['pds_lam0.5']:.6f}"
            f" | R16 cos {rr['cos_agg']:.6f} pds@.5 {rr['pds_lam0.5']:.6f} A {np.mean([rr['AB'][c]['A'] for c in 'ABC']):.3f}")
        (OUT / f"comp445-{held}.json").write_text(json.dumps(res, indent=1, default=float))

    rec("pair@0.56", union_rows(pair), g_stack, check_ref="stack|K3+sh1.5")
    log("KNOWN ANSWER PASS")
    for g in (0.42, 0.56, 0.70):
        rec(f"pair+dld1@{g:g}", union_rows(pair + [dld1]), g)
    rec("CTRL:pair+dld1perm@0.56", union_rows(pair + [permute_rows(dld1, 445)]), g_stack)
    # 223's construction: DLD-1 as its OWN residual (own carrying-row sizing, own dose) on top of the unchanged stack
    for gd in (0.2, 0.35, 0.5):
        rec(f"stack+sepDLD1@{gd:g}", union_rows(pair), g_stack, extra=dld1, g_extra=gd)
    rec("CTRL:stack+sepDLD1perm@0.35", union_rows(pair), g_stack, extra=permute_rows(dld1, 445), g_extra=0.35)
    res["peak_rss_gib"] = E062.peak_rss_gib()
    (OUT / f"comp445-{held}.json").write_text(json.dumps(res, indent=1, default=float))
    log(f"done; peak RSS {res['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    run(sys.argv[1])
