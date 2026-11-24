"""213 — production-faithful staging UNDER A TRANSPORT LAW.

This is `experiments/170-gamma-mse/stage.py`, whose helpers are imported rather than copied,
with exactly one thing changed: after `load_donor_lean` reads a donor's `delta`, the law is
applied to it,

    g = law.to_law(delta, p_ctrl_donor),     p_ctrl_donor = expm1(b_ctrl)/TS

on the donor's own column axis, in row chunks, in place.  Everything after that is
production's own code -- `MultiDonorModel(gene_norm=False, eb=False)`, `donor_weights`,
`build_v9_pre_whiten_effects`, `apply_perturbmatch_k562_replacement` -- so the law is the only
difference between arms.

WHY THIS AND NOT THE PLOLO CHAIN.  PLOLO rebuilds the production effect matrix from 062's
`prepare_main_source`; measured here, its `pre` sits at per-row cosine 0.98 against
production's own `candidate_pre` (median 0.9809, min 0.892) and its operator-column median row
norm is 4.1184 against a15's recorded 4.1229.  That is fine for a DELTA instrument -- 162
validated it against six board points -- but the `mse` decomposition (096) and the h1-test DE
members read absolute geometry, so those go through the production chain.

KNOWN-ANSWER GATE (rule 10).  At `logc20` the law is the identity on `delta`, so this must
reproduce a15's shipped hashes BIT-FOR-BIT:

    baseline_pre  sha256  7678bab9ea9d1b8d951b0b24edebde4939b245e967f481d91d44ae24aac7c98d
    candidate_pre sha256  0601336fdf1ec5c533dcad69f6cb6c09d6368182bc8c69d10d516a1ecb19bf65

If either fails, the law injection is not neutral at its own incumbent and nothing downstream
may be read.  170's own trap is inherited and re-stated: `candidate_pre` reproduces only at
`--chunk-rows 128` for PerturbMatch, because its centring mean accumulates per chunk.

    uv run python stage213.py --laws logc20 logc40
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "work"
E170 = REPO / "experiments/170-gamma-mse"
sys.path.insert(0, str(E170))
sys.path.insert(0, str(HERE))

os.environ.setdefault("VCC_ALPHA", "0.982")
os.environ.setdefault("VCC_CD4_TARGET_F", "0.10")
os.environ.setdefault("VCC_MASS_FIX", "minnorm")
os.environ.setdefault("VCC_WHITEN_GAMMA", "0.65")
os.environ.setdefault("VCC_WHITEN_LAM", "0.0")

import stage as ST  # noqa: E402  (170's staging helpers; its `main` is never called)
import law213 as L  # noqa: E402

A15 = ST.A15
WHITENER = ST.WHITENER


def load_donor_law(B, path, target_genes, law, chunk=2048):
    """170's `load_donor_lean`, then the law, in row chunks, in float32 in place.

    The donor table is float32 on disk and production's own `MultiDonorModel` works in
    float32, so the transform is done at production's precision.  At the incumbent
    `to_law` returns its argument unchanged, so the round trip is the identity and the
    a15 sha256 gate below is what proves it.
    """
    d = ST.load_donor_lean(B, path, target_genes, chunk=chunk)
    if law.is_incumbent:
        return d
    p = (np.expm1(np.asarray(d.b_ctrl, dtype=np.float64)) / L.TS)
    n = d.delta.shape[0]
    for lo in range(0, n, 512):
        hi = min(lo + 512, n)
        d.delta[lo:hi] = law.to_law(
            np.asarray(d.delta[lo:hi], dtype=np.float64), p).astype(np.float32)
    return d


def stage_one(B, law, spec, genes, perts, h1_t, h1_genes, pm, donor_paths):
    with ST.lean_norm():
        donors = {n: load_donor_law(B, p, genes, law) for n, p in donor_paths.items()}
        model = B.MultiDonorModel(donors, spec.genes, gene_norm=False, eb=False)
    del donors
    gc.collect()
    weights = B.donor_weights(model)
    cov_h1 = model.coverage(list(h1_t))
    cov_chal = model.coverage(list(perts))
    ST.trim_model_rows(model, list(perts) + list(h1_t))
    gc.collect()
    operator_columns = B.production_operator_columns(model, perts)

    baseline_pre = B.build_v9_pre_whiten_effects(
        model, perts, alpha=B.ALPHA, beta=B.BETA, weights=weights)
    base_sha = B.array_sha256(baseline_pre)
    cand_chal, _ = B.apply_perturbmatch_k562_replacement(
        baseline_pre, model, pm, perts, alpha=B.ALPHA, weights=weights)
    cand_sha = B.array_sha256(cand_chal)

    base_h1 = B.build_v9_pre_whiten_effects(
        model, tuple(map(str, h1_t)), alpha=B.ALPHA, beta=B.BETA, weights=weights)
    cand_h1, pm_h1 = B.apply_perturbmatch_k562_replacement(
        base_h1, model, pm, tuple(map(str, h1_t)), alpha=B.ALPHA, weights=weights)
    del model, baseline_pre, base_h1
    gc.collect()
    audit = dict(
        law=law.key, weights=dict(weights) if isinstance(weights, dict) else list(weights),
        baseline_pre_sha256=base_sha, candidate_pre_sha256=cand_sha,
        before_all_median_norm=float(np.median(np.linalg.norm(cand_chal, axis=1))),
        before_operator_median_norm=float(
            np.median(np.linalg.norm(cand_chal[:, operator_columns], axis=1))),
        h1_donor_coverage=cov_h1, challenge_donor_coverage=cov_chal,
        h1_pm_modified=int(pm_h1["both_supported_modified_count"]),
        operator_columns=int(operator_columns.size),
        peak_rss_gib=float(B.peak_rss_gib()),
    )
    return cand_chal, cand_h1, operator_columns, audit


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--laws", nargs="*", default=["logc20", "logc40"])
    args = ap.parse_args()
    laws = L.all_laws()

    B = ST.builder()
    B.assert_frozen_contract()
    donor_paths = {
        "k562_gwps": REPO / "data/donors/k562_gwps",
        "orion_hct116": REPO / "data/donors/orion_hct116",
        "orion_hek293t": REPO / "data/donors/orion_hek293t",
    }
    hashes = B.collect_input_hashes(
        REPO / "data/challenge", ("A", "B", "C"), donor_paths,
        REPO / "data/external/perturbmatch",
        REPO / "data/external/cd4t/cloud-compact/cd4_compact",
        WHITENER, ST.WHITENER_META, whitener_id="union_pruned")
    assert hashes["whitener/whitener_union.npz"] == A15["operator_sha256"], "operator drift"

    spec = B.ChallengeSpec.load(REPO / "data/challenge")
    genes = np.asarray([str(g) for g in spec.genes])
    perts = tuple(map(str, spec.perts))
    h1_t, h1_genes = ST.h1_targets_and_genes()
    gpos = {g: i for i, g in enumerate(genes)}
    panel_positions = np.asarray([gpos[t] for t in perts], dtype=np.int64)
    h1_pos = np.asarray([gpos[t] for t in h1_t], dtype=np.int64)

    pm = B.EXP061.prepare_perturbmatch(
        REPO / "data/external/perturbmatch", tuple(genes), chunk_rows=128, max_rss_gib=9.5)
    print(f"[stage213] PerturbMatch prepared; RSS {B.peak_rss_gib():.2f} GiB", flush=True)

    if not (WORK / "cd4.npz").exists():
        _, cd4_chal = B.prepare_frozen_cd4(
            REPO / "data/external/cd4t/cloud-compact/cd4_compact",
            perts, tuple(genes), perts, chunk_rows=128, max_rss_gib=9.5)
        _, cd4_h1 = B.prepare_frozen_cd4(
            REPO / "data/external/cd4t/cloud-compact/cd4_compact",
            tuple(map(str, h1_t)), tuple(genes), tuple(map(str, h1_t)),
            chunk_rows=128, max_rss_gib=9.5)
        np.savez_compressed(
            WORK / "cd4.npz", chal_mean=cd4_chal.mean_effect, chal_votes=cd4_chal.gene_votes,
            chal_counts=cd4_chal.target_condition_count, h1_mean=cd4_h1.mean_effect,
            h1_votes=cd4_h1.gene_votes, h1_counts=cd4_h1.target_condition_count,
            mode=np.asarray([str(cd4_chal.mode)]))
        del cd4_chal, cd4_h1
        gc.collect()
    print(f"[stage213] CD4 ready; RSS {B.peak_rss_gib():.2f} GiB", flush=True)

    audits = json.loads((WORK / "stage213.json").read_text()) \
        if (WORK / "stage213.json").exists() else {}
    for key in args.laws:
        out = WORK / f"pre_{key}.npz"
        if out.exists():
            print(f"[stage213] {key} cached", flush=True)
            continue
        law = laws[key]
        cand_chal, cand_h1, opcols, audit = stage_one(
            B, law, spec, genes, perts, h1_t, h1_genes, pm, donor_paths)
        if law.is_incumbent:
            ok_b = audit["baseline_pre_sha256"] == A15["baseline_pre_sha256"]
            ok_c = audit["candidate_pre_sha256"] == A15["candidate_pre_sha256"]
            print(f"[stage213] KA baseline_pre  {'MATCH' if ok_b else 'MISMATCH'}\n"
                  f"[stage213] KA candidate_pre {'MATCH' if ok_c else 'MISMATCH'}  "
                  f"{audit['candidate_pre_sha256']}", flush=True)
            if not (ok_b and ok_c):
                raise SystemExit("a15 STAGING GATE FAILED — the law injection is not "
                                 "neutral at its own incumbent. STOP.")
        np.savez_compressed(
            out, genes=genes, perts=np.asarray(perts), cand_chal=cand_chal,
            h1_targets=h1_t, h1_genes=h1_genes, cand_h1=cand_h1,
            panel_positions=panel_positions, h1_positions=h1_pos,
            operator_columns=opcols)
        audits[key] = audit
        (WORK / "stage213.json").write_text(json.dumps(audits, indent=1))
        print(f"[stage213] {key}: before_operator_median_norm "
              f"{audit['before_operator_median_norm']!r} "
              f"(a15 4.122892101713494), RSS {B.peak_rss_gib():.2f} GiB", flush=True)
        del cand_chal, cand_h1
        gc.collect()


if __name__ == "__main__":
    main()
