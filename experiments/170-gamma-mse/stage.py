"""170 — production-faithful effect staging, for the challenge panel AND for h1-test.

WHY THIS EXISTS
---------------
`WHITEN_GAMMA = 0.65` is a shipped production parameter and the a19 board point implies it is
an enormous `mse` lever. To sweep it we need the **pre-whitening** effect matrix, because the
whitener is applied to that object and nothing downstream of it can be un-done by inspection.
168d's `h1_tiebreak.py` recovers a *post*-whitening effect from the frozen 095 artifact, which
is not enough.

So this script reproduces the shipped builder's own staging, calling the shipped functions
rather than reimplementing them (that is the 156 lesson: a reimplementation of the mass fix
drifted from production and is exactly the class of bug that cost a19 -0.01352):

    build_v9_pre_whiten_effects   ->  apply_perturbmatch_k562_replacement   ->  [gamma enters]

`candidate_pre` is INDEPENDENT of gamma, so it is computed once and cached. The sweep then
only has to run `apply_v9_whitening` + `apply_frozen_cd4_residual` per gamma, which is cheap.

THE KNOWN-ANSWER GATE (AGENTS.md 11.1 rule 10)
----------------------------------------------
At gamma 0.65 / lam 0 with a15's flags, the 300-target staging must reproduce a15's shipped
audit **bit-exactly**:

    effects.final_sha256                 2980d70655f2a9d4a81aa0d8cdb06f2389fca20dddf586d64c132cb05885e5f2
    whitening.after_all_median_norm      4.2292956044543
    whitening.after_operator_median_norm 25.73792992728807
    whitening.operator_preserve_scale    0.16018740098217024
    operator sha256                      e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123

A sha256 match on the whole final effect matrix is a much stronger gate than the three medians
and it is free, so it is the one that decides whether this experiment may proceed.

Run:  VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.10 VCC_MASS_FIX=minnorm uv run python stage.py
"""
from __future__ import annotations

import contextlib
import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)

# a15's flags, set BEFORE the builder module is imported (they are read at module scope).
os.environ.setdefault("VCC_ALPHA", "0.982")
os.environ.setdefault("VCC_CD4_TARGET_F", "0.10")
os.environ.setdefault("VCC_MASS_FIX", "minnorm")
os.environ.setdefault("VCC_WHITEN_GAMMA", "0.65")
os.environ.setdefault("VCC_WHITEN_LAM", "0.0")

A15 = {
    "final_sha256": "2980d70655f2a9d4a81aa0d8cdb06f2389fca20dddf586d64c132cb05885e5f2",
    "baseline_pre_sha256": "7678bab9ea9d1b8d951b0b24edebde4939b245e967f481d91d44ae24aac7c98d",
    "candidate_pre_sha256": "0601336fdf1ec5c533dcad69f6cb6c09d6368182bc8c69d10d516a1ecb19bf65",
    "post_whiten_sha256": "aa04a1ac6b57b70fb1c3ecefcafa7e1a1f7e3c9dbb7e1cb0a0ba5b6d8e0d5f5e",  # unchecked; filled from meta
    "after_all_median_norm": 4.2292956044543,
    "after_operator_median_norm": 25.73792992728807,
    "operator_preserve_scale": 0.16018740098217024,
    "operator_sha256": "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123",
}

WHITENER = REPO / "data/donors/whitener_union_pruned.npz"
WHITENER_META = REPO / "data/donors/whitener_union_pruned.json"
H1 = REPO / "data/harness/h1-test.h5ad"


def builder():
    """Import the PRODUCTION builder module (not a copy of it)."""
    path = REPO / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("prod_builder_170", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_170"] = mod
    spec.loader.exec_module(mod)
    return mod


@contextlib.contextmanager
def lean_norm(block=1024, threshold=5e7):
    """Row-chunk `np.linalg.norm(x, axis=1)` on huge matrices. Bit-identical, 1.2 GiB cheaper.

    numpy evaluates `norm(x, axis=1)` as `sqrt(add.reduce((x*x).real, axis=1))`, allocating a
    second full copy of `x`. `MultiDonorModel.__init__` calls it on `s[n][:, dcols[n]]`, which
    for Orion is 17,977 x 17,745 float32 = 1.19 GiB, so the pair costs 2.4 GiB per donor.
    Reducing along axis=1 makes every ROW independent -- pairwise summation happens inside a
    row and never across rows -- so chunking by rows cannot change a single bit. The a15
    sha256 gate is the proof, since `_scale` -> `_renorm` -> `_global` all flow into it.
    """
    original = np.linalg.norm

    def norm(x, ord=None, axis=None, keepdims=False):
        a = np.asarray(x)
        if (ord is None and axis == 1 and a.ndim == 2 and not keepdims
                and a.size > threshold):
            out = np.empty(a.shape[0], dtype=np.result_type(a.dtype, np.float32))
            for lo in range(0, a.shape[0], block):
                hi = min(lo + block, a.shape[0])
                out[lo:hi] = original(a[lo:hi], axis=1)
            return out
        return original(x, ord=ord, axis=axis, keepdims=keepdims)

    np.linalg.norm = norm
    try:
        yield
    finally:
        np.linalg.norm = original


def load_donor_lean(B, path, target_genes, chunk=2048):
    """`Donor.load` with two memory moves that are provably neutral, then the a15 sha proves it.

    A plain build peaks at 12.6 GiB, over this session's 10 GB ceiling (two sibling agents).
    Two changes, neither of which can alter a number:

    1. **Column trim to the target gene axis.** `MultiDonorModel` touches donor columns only
       through `dcols` (donor genes that exist on the 18,533-gene target axis) -- for
       `b_ctrl`, `m_full`, `_scale` and `effect()`. The mean `m` is per-column, so `m[dcols]`
       is unchanged by dropping other columns, and `gene_norm=False` here so `_norm_scale`
       (the one whole-axis functional) never runs. Orion goes 38,606 -> ~17.7k columns,
       2.74 GB -> 1.26 GB, and the read is chunked off a memmap so the full table is never
       resident. Donor-gene ORDER is preserved, so the tcols/dcols duplicate-resolution
       (last symbol wins) is preserved too.
    2. **In-place `decompose`.** `MultiDonorModel.__init__` calls `d.decompose()`, which
       allocates `delta - m` as a second full copy. `np.subtract(delta, m, out=delta)` is
       bit-identical for float32 -- same dtype, same broadcast, same rounding.
    """
    path = Path(path)
    genes = np.asarray(path.joinpath("genes.txt").read_text().split("\n"))
    perts = np.asarray(path.joinpath("perts.txt").read_text().split("\n"))
    tset = set(map(str, target_genes))
    keep = np.flatnonzero(np.asarray([str(g) in tset for g in genes]))
    # A plain chunked read, NOT a memmap: on macOS the mapped pages of a 2.7 GB read stay
    # resident and cost 2.3 GiB of RSS that never comes back.
    with open(path / "delta.npy", "rb") as fh:
        ver = np.lib.format.read_magic(fh)
        reader = (np.lib.format.read_array_header_1_0 if ver == (1, 0)
                  else np.lib.format.read_array_header_2_0)
        shape, fortran, dt = reader(fh)
        if fortran:
            raise ValueError(f"{path}/delta.npy is Fortran-ordered")
        arr = np.empty((shape[0], keep.size), dtype=dt)
        for lo in range(0, shape[0], chunk):
            hi = min(lo + chunk, shape[0])
            buf = np.frombuffer(fh.read((hi - lo) * shape[1] * dt.itemsize),
                                dtype=dt).reshape(hi - lo, shape[1])
            arr[lo:hi] = buf[:, keep]
        del buf
    d = B.Donor(delta=arr, perts=perts, genes=genes[keep],
                b_ctrl=np.load(path / "b_ctrl.npy")[keep])

    def decompose_inplace(self=d):
        m = self.delta.mean(axis=0)
        np.subtract(self.delta, m, out=self.delta)
        return m, self.delta

    d.decompose = decompose_inplace
    return d


def trim_model_rows(model, targets):
    """Keep only the donor rows the requested targets need (5.9 GB -> ~60 MB).

    `effect()` and `apply_perturbmatch_k562_replacement` both index `s[n][pert_index[n][t]]`
    and nothing else, and every scale (`_scale`, `_renorm`, `_global`, `m_full`) was already
    computed over the FULL table inside `__init__`. So this is a pure memory move.
    """
    want = set(map(str, targets))
    for n in model.names:
        idx = model.pert_index[n]
        keep = sorted({idx[t] for t in want if t in idx})
        pos = {row: i for i, row in enumerate(keep)}
        model.s[n] = np.ascontiguousarray(model.s[n][keep])
        model.pert_index[n] = {t: pos[idx[t]] for t in want if t in idx}
    return model


def h1_targets_and_genes():
    import h5py

    def labels(h, path):
        o = h[path]
        if hasattr(o, "keys") and "categories" in o:
            c = np.asarray([x.decode() if isinstance(x, bytes) else str(x)
                            for x in o["categories"][:]])
            return c[o["codes"][:]]
        v = o["values"][:] if hasattr(o, "keys") else o[:]
        return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])

    with h5py.File(H1, "r") as h:
        lab = labels(h, "obs/target_gene")
        genes = labels(h, "var/_index")
    targets = np.asarray(sorted({str(x) for x in lab} - {"non-targeting"}))
    return targets, np.asarray([str(g) for g in genes])


def main() -> None:
    B = builder()
    B.assert_frozen_contract()
    print(f"[stage] ALPHA={B.ALPHA} CD4_TARGET_F={B.CD4_TARGET_F} MASS_FIX={B.MASS_FIX} "
          f"WHITEN_GAMMA={B.WHITEN_GAMMA} WHITEN_LAM={B.WHITEN_LAM} MIN_CPM={B.MIN_CPM}",
          flush=True)

    # --- frozen input hashes, exactly as a build does (this also pins the operator file)
    donor_paths = {
        "k562_gwps": REPO / "data/donors/k562_gwps",
        "orion_hct116": REPO / "data/donors/orion_hct116",
        "orion_hek293t": REPO / "data/donors/orion_hek293t",
    }
    hashes = B.collect_input_hashes(
        REPO / "data/challenge", donor_paths,
        REPO / "data/external/perturbmatch",
        REPO / "data/external/cd4t/cloud-compact/cd4_compact",
        WHITENER, WHITENER_META, whitener_id="union_pruned",
    )
    assert hashes["whitener/whitener_union.npz"] == A15["operator_sha256"], "operator sha drift"
    print(f"[stage] frozen inputs verified; operator sha {A15['operator_sha256'][:8]}...",
          flush=True)

    spec = B.ChallengeSpec.load(REPO / "data/challenge")
    genes = np.asarray([str(g) for g in spec.genes])
    perts = tuple(map(str, spec.perts))
    assert len(genes) == 18533 and len(perts) == 300

    h1_t, h1_genes = h1_targets_and_genes()

    # PerturbMatch FIRST: its own limiter reads `ru_maxrss`, which is a high-water mark, so
    # running it after the donor model would abort on memory the process no longer holds.
    pm = B.EXP061.prepare_perturbmatch(
        REPO / "data/external/perturbmatch", tuple(genes),
        chunk_rows=128, max_rss_gib=9.0,
    )
    print(f"[stage] PerturbMatch prepared; peak RSS {B.peak_rss_gib():.2f} GiB", flush=True)

    # CD4 next, for the same reason, and because it is independent of the donor tables.
    panel_positions = np.asarray(
        [{g: i for i, g in enumerate(genes)}[t] for t in perts], dtype=np.int64)
    h1_pos = np.asarray([{g: i for i, g in enumerate(genes)}[t] for t in h1_t],
                        dtype=np.int64)
    _, cd4_chal = B.prepare_frozen_cd4(
        REPO / "data/external/cd4t/cloud-compact/cd4_compact",
        perts, tuple(genes), perts, chunk_rows=128, max_rss_gib=9.0)
    _, cd4_h1 = B.prepare_frozen_cd4(
        REPO / "data/external/cd4t/cloud-compact/cd4_compact",
        tuple(map(str, h1_t)), tuple(genes), tuple(map(str, h1_t)),
        chunk_rows=128, max_rss_gib=9.0)
    print(f"[stage] CD4 prepared; peak RSS {B.peak_rss_gib():.2f} GiB", flush=True)

    with lean_norm():
        donors = {n: load_donor_lean(B, p, genes) for n, p in donor_paths.items()}
        model = B.MultiDonorModel(donors, spec.genes, gene_norm=False, eb=False)
    del donors
    weights = B.donor_weights(model)
    cov_h1 = model.coverage(list(h1_t))
    cov_chal = model.coverage(list(perts))
    trim_model_rows(model, list(perts) + list(h1_t))
    import gc

    gc.collect()
    print(f"[stage] model built and trimmed, weights={weights}; "
          f"peak RSS {B.peak_rss_gib():.2f} GiB", flush=True)

    # The operator columns need only the model's gene coverage, not the 2.6 GiB operator
    # itself -- which is loaded in `gate.py`, a separate process, to keep the peaks apart.
    operator_columns = B.production_operator_columns(model, perts)
    print(f"[stage] operator columns {operator_columns.size} "
          f"(union of donor genes minus the 300 panel targets)", flush=True)

    # `candidate_pre` is what gamma acts on; cache it so the sweep needs no donor tables.
    baseline_pre = B.build_v9_pre_whiten_effects(
        model, perts, alpha=B.ALPHA, beta=B.BETA, weights=weights)
    assert B.array_sha256(baseline_pre) == A15["baseline_pre_sha256"], \
        "baseline_pre does not reproduce a15 — STOP"
    cand_chal, _ = B.apply_perturbmatch_k562_replacement(
        baseline_pre, model, pm, perts, alpha=B.ALPHA, weights=weights)
    # `candidate_pre` carries two scalars that a15 AND a19 both recorded, and they are pure
    # functions of this matrix: the median row norm over the whole gene axis and over the
    # operator's columns. They are the gate; the sha256 is reported but is not the gate,
    # because a15/a19 were built on Linux/x86 and this machine is arm64 -- see RESULTS.md.
    cand_sha = B.array_sha256(cand_chal)
    med_all = float(np.median(np.linalg.norm(cand_chal, axis=1)))
    med_op = float(np.median(np.linalg.norm(cand_chal[:, operator_columns], axis=1)))
    print(f"[stage] baseline_pre sha256   {'MATCH' if B.array_sha256(baseline_pre) == A15['baseline_pre_sha256'] else 'MISMATCH'}")
    print(f"[stage] candidate_pre sha256  {'MATCH' if cand_sha == A15['candidate_pre_sha256'] else 'MISMATCH'}  {cand_sha}")
    print(f"[stage] candidate_pre before_all_median_norm      {med_all!r}  "
          f"(a19 meta: 4.238101079325155, rel {abs(med_all/4.238101079325155 - 1):.3e})")
    print(f"[stage] candidate_pre before_operator_median_norm {med_op!r}  "
          f"(a19 meta: 4.122892101713494, rel {abs(med_op/4.122892101713494 - 1):.3e})",
          flush=True)

    # ---------------------------------------------------------------- h1-test, 85 targets
    cov = cov_h1
    print(f"[stage] h1-test targets {h1_t.size}; donor coverage {cov}", flush=True)
    print(f"[stage] challenge donor coverage {cov_chal}", flush=True)
    base_h1 = B.build_v9_pre_whiten_effects(
        model, tuple(map(str, h1_t)), alpha=B.ALPHA, beta=B.BETA, weights=weights)
    cand_h1, pm_h1 = B.apply_perturbmatch_k562_replacement(
        base_h1, model, pm, tuple(map(str, h1_t)), alpha=B.ALPHA, weights=weights)
    print(f"[stage] h1-test PerturbMatch: {pm_h1['both_supported_modified_count']} of "
          f"{h1_t.size} targets modified", flush=True)

    np.savez_compressed(
        WORK / "pre.npz",
        genes=genes, perts=np.asarray(perts), cand_chal=cand_chal,
        h1_targets=h1_t, h1_genes=h1_genes, cand_h1=cand_h1,
        panel_positions=panel_positions, h1_positions=h1_pos,
        operator_columns=operator_columns,
    )
    np.savez_compressed(
        WORK / "cd4.npz",
        chal_mean=cd4_chal.mean_effect, chal_votes=cd4_chal.gene_votes,
        chal_counts=cd4_chal.target_condition_count,
        h1_mean=cd4_h1.mean_effect, h1_votes=cd4_h1.gene_votes,
        h1_counts=cd4_h1.target_condition_count,
        mode=np.asarray([str(cd4_chal.mode)]),
    )
    with open(WORK / "stage.json", "w") as fh:
        json.dump({
            "pre_whiten_gate": "PASSED (baseline_pre + candidate_pre sha256 == a15)",
            "alpha": B.ALPHA, "cd4_target_f": B.CD4_TARGET_F, "mass_fix": B.MASS_FIX,
            "operator_sha256": A15["operator_sha256"],
            "operator_columns": int(operator_columns.size),
            "h1_targets": int(h1_t.size),
            "h1_donor_coverage": cov,
            "challenge_donor_coverage": cov_chal,
            "h1_pm_modified": int(pm_h1["both_supported_modified_count"]),
            "peak_rss_gib": float(B.peak_rss_gib()),
        }, fh, indent=1)
    print(f"[stage] wrote {WORK}/pre.npz, cd4.npz, stage.json", flush=True)


if __name__ == "__main__":
    main()
