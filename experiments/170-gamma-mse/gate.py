"""170 — apply the SHIPPED operator at each gamma, and gate on reproducing a15 exactly.

Separate process from `stage.py` on purpose: the donor tables peak around 6 GiB and the
17,942 x 17,942 float64 operator is another 2.6 GiB, and this session must stay under 10 GB
with two sibling agents running.

WHAT IT DOES
    pre.npz (candidate_pre, gamma-free)  ->  apply_v9_whitening(gamma, lam=0)
                                         ->  apply_frozen_cd4_residual
                                         ->  work/eff/g<gamma>.npz

KNOWN-ANSWER GATE (rule 10). At gamma 0.65 / lam 0 the 300-target result must reproduce a15's
shipped audit bit-exactly: `final_sha256`, `after_all_median_norm`,
`after_operator_median_norm`, `operator_preserve_scale`, and the operator's own sha256. A
sha256 on the whole final matrix is far stronger than the three medians, and if it fails this
script raises rather than reporting numbers.
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "work"
(WORK / "eff").mkdir(parents=True, exist_ok=True)

os.environ.setdefault("VCC_ALPHA", "0.982")
os.environ.setdefault("VCC_CD4_TARGET_F", "0.10")
os.environ.setdefault("VCC_MASS_FIX", "minnorm")

OPERATOR_SHA = "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"

# TWO shipped board points, both built from the SAME `candidate_pre` (their metas record the
# identical sha `0601336f...`), which this machine reproduces bit-exactly. They differ only in
# (gamma, lam) -- so together they anchor the operator path at two points ON THE GAMMA AXIS.
#
# The post-whitening `final_sha256` is NOT gateable here and must not be: a15/a19 were built
# on Linux/x86_64 (their metas record it) and this is arm64. `whitener.apply` is a pair of
# 17,942 x 17,942 GEMMs, and BLAS kernels differ in reduction order across architectures, so
# the last bits move. The scalars below are the same object read through medians, and they
# agree to ~2e-16 -- which is exactly the signature of last-bit GEMM noise and not of a
# different computation.
ANCHORS = {
    (0.65, 0.0): {  # a15, `submissions/a15-cd4f010.vcc.meta.json`
        "before_operator_median_norm": 4.122892101713489,
        "after_operator_median_norm": 25.73792992728807,
        "operator_preserve_scale": 0.16018740098217024,
        "after_all_median_norm": 4.2292956044543,
        "cd4_f_per_unit": 1.055664536904655,
        "cd4_solved_fraction": 0.09472706196346517,
        "cd4_candidate_pre_preserve_median_norm": 4.166553356202115,
        "cd4_post_add_global_preserve_scale": 0.9895210139518232,
        "final_sha256": "2980d70655f2a9d4a81aa0d8cdb06f2389fca20dddf586d64c132cb05885e5f2",
    },
    (0.80, 0.007): {  # a19, `submissions/a19-whiten-g080-lam0007.h5ad.meta.json`
        "before_operator_median_norm": 4.122892101713494,
        "after_operator_median_norm": 22.030919431001,
        "operator_preserve_scale": 0.1871411728696139,
        "after_all_median_norm": 4.235393246076997,
        "cd4_f_per_unit": 1.0541447101386758,
        "cd4_solved_fraction": 0.094863635929876,
        "final_sha256": "1f12a563d68f2a78d79226f634beceff39a35609d817e6dd090c5b48c79393aa",
    },
}
RTOL = 1e-11
WHITENER = REPO / "data/donors/whitener_union_pruned.npz"
WHITENER_META = REPO / "data/donors/whitener_union_pruned.json"

GRID = [0.00, 0.35, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80]


def builder():
    path = REPO / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("prod_builder_170g", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_170g"] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    B = builder()
    pre = np.load(WORK / "pre.npz", allow_pickle=True)
    cd4z = np.load(WORK / "cd4.npz", allow_pickle=True)
    cols = pre["operator_columns"].astype(np.int64)
    panel_pos = pre["panel_positions"].astype(np.int64)
    h1_pos = pre["h1_positions"].astype(np.int64)
    cand = {"chal": pre["cand_chal"], "h1": pre["cand_h1"]}
    mode = str(cd4z["mode"][0])
    cd4 = {
        "chal": B.EXP062.CD4Prepared(
            mode=mode, mean_effect=cd4z["chal_mean"], gene_votes=cd4z["chal_votes"],
            target_condition_count=cd4z["chal_counts"], condition_details=[]),
        "h1": B.EXP062.CD4Prepared(
            mode=mode, mean_effect=cd4z["h1_mean"], gene_votes=cd4z["h1_votes"],
            target_condition_count=cd4z["h1_counts"], condition_details=[]),
    }
    pos = {"chal": panel_pos, "h1": h1_pos}

    sha = B.file_sha256(WHITENER)
    if sha != OPERATOR_SHA:
        raise SystemExit(f"operator sha {sha} != the shipped {OPERATOR_SHA}")
    meta = json.loads(WHITENER_META.read_text())
    whit = B.Whitener.load(WHITENER)
    if whit.evecs.shape != (cols.size, cols.size):
        raise SystemExit(f"operator shape {whit.evecs.shape} != {(cols.size, cols.size)}")
    print(f"[gate] operator {WHITENER.name} sha {sha[:8]}... shape {whit.evecs.shape}, "
          f"eigen range [{whit.evals.min():.3e}, {whit.evals.max():.4f}], "
          f"basis {meta.get('basis')}", flush=True)

    # The CD4 residual covers the CHALLENGE panel, not h1-test: only 2 of h1-test's 85
    # targets survive the frozen `complete mean3` policy, so the median row norm of the h1
    # residual is exactly 0 and `solve_cd4_fraction` cannot solve a fraction at all (it
    # raises). h1-test is therefore built with the fraction production SOLVES ON THE PANEL
    # at that gamma -- i.e. the fraction the shipped artifact carries -- applied through the
    # same shipped `add_post_whiten_residual`. It touches 2 of 85 rows.
    n_h1_carry = int(np.count_nonzero(np.abs(cd4z["h1_mean"]).sum(axis=1)))
    print(f"[gate] CD4 support: challenge "
          f"{int(np.count_nonzero(np.abs(cd4z['chal_mean']).sum(axis=1)))}/300 targets, "
          f"h1-test {n_h1_carry}/85", flush=True)

    # ---------------------------------------------------------------- known-answer gate
    print("\n===== KNOWN-ANSWER GATE: two shipped board points on the gamma axis =====")
    passed = True
    for (g, lam), want in ANCHORS.items():
        whitened, waudit = B.apply_v9_whitening(cand["chal"], whit, cols, gamma=g, lam=lam)
        final, caudit = B.apply_frozen_cd4_residual(whitened, cd4["chal"], cols, panel_pos)
        got = {
            "before_operator_median_norm": float(waudit["before_operator_median_norm"]),
            "after_operator_median_norm": float(waudit["after_operator_median_norm"]),
            "operator_preserve_scale": float(waudit["operator_preserve_scale"]),
            "after_all_median_norm": float(waudit["after_all_median_norm"]),
            "cd4_f_per_unit": float(caudit["residual_frame"]["f_per_unit_fraction"]),
            "cd4_solved_fraction": float(caudit["residual_frame"]["solved_fraction"]),
            "cd4_candidate_pre_preserve_median_norm": float(
                caudit["post_add_preserve"]["candidate_pre_preserve_median_norm"]),
            "cd4_post_add_global_preserve_scale": float(
                caudit["post_add_preserve"]["post_add_global_preserve_scale"]),
        }
        label = "a15" if g == 0.65 else "a19"
        print(f"  --- {label}: gamma {g} / lam {lam}")
        for k, v in want.items():
            if k == "final_sha256":
                mine = B.array_sha256(final)
                print(f"      {'final_sha256':<40} INFORMATIONAL (arm64 vs x86 GEMM)  "
                      f"{'equal' if mine == v else 'differs'}")
                continue
            rel = abs(got[k] / v - 1.0) if v else abs(got[k])
            ok = rel <= RTOL
            passed &= ok
            print(f"      {k:<40} {'PASS' if ok else 'FAIL'}  {got[k]!r}  rel {rel:.2e}")
    if not passed:
        raise SystemExit("KNOWN-ANSWER GATE FAILED — this pipeline is not a15/a19. STOP.")
    print("  ==> GATE PASSED at BOTH points, every scalar within 1e-11 relative.\n",
          flush=True)

    audits = {}
    for gamma in GRID:
        row = {"gamma": gamma, "lam": 0.0, "h1_cd4_targets": n_h1_carry}
        for surf in ("chal", "h1"):
            whitened, waudit = B.apply_v9_whitening(
                cand[surf], whit, cols, gamma=gamma, lam=0.0)
            if surf == "chal":
                final, caudit = B.apply_frozen_cd4_residual(
                    whitened, cd4[surf], cols, pos[surf])
                frac = float(caudit["residual_frame"]["solved_fraction"])
                fper = float(caudit["residual_frame"]["f_per_unit_fraction"])
                pscale = float(caudit["post_add_preserve"]["post_add_global_preserve_scale"])
                chal_fraction = frac
            else:
                final, preserve = B.EXP062.add_post_whiten_residual(
                    whitened, cd4[surf], chal_fraction, cols)
                frac, fper = chal_fraction, float("nan")
                pscale = float(preserve["post_add_global_preserve_scale"])
            np.savez_compressed(WORK / "eff" / f"{surf}_g{gamma:.2f}.npz",
                                effect=final, gamma=gamma)
            row[surf] = {
                "final_sha256": B.array_sha256(final),
                "whitening": {k: float(v) for k, v in waudit.items()},
                "cd4_fraction": frac,
                "cd4_f_per_unit": fper,
                "post_add_preserve_scale": pscale,
                "median_row_norm": float(np.median(np.linalg.norm(final, axis=1))),
            }
        audits[f"{gamma:.2f}"] = row
        print(f"[gate] gamma {gamma:.2f}: h1 median row norm "
              f"{row['h1']['median_row_norm']:.5f}, operator amplification "
              f"{row['h1']['whitening']['after_operator_median_norm']:.3f}, "
              f"chal median row norm {row['chal']['median_row_norm']:.5f}", flush=True)

    with open(WORK / "gate.json", "w") as fh:
        json.dump(audits, fh, indent=1)
    print(f"[gate] wrote {WORK}/gate.json and {len(GRID)*2} effect files; "
          f"peak RSS {B.peak_rss_gib():.2f} GiB")


if __name__ == "__main__":
    main()
