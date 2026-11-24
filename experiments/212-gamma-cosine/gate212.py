"""212 Phase 0 — 170's operator gate, re-run unchanged, on an EXTENDED gamma grid.

This is `experiments/170-gamma-mse/gate.py` with exactly three differences, all of them
plumbing and all of them recorded here:

  1. `WORK` points at THIS experiment's work dir (nothing may be written to 170's);
  2. `pre.npz` / `cd4.npz` are READ from 170's cache -- and `candidate_pre`'s sha256 is
     re-asserted against a15's shipped audit before anything else runs, which is the whole
     point of gate 1 and is what makes the cached object usable;
  3. `GRID` gains 0.30, 0.40, 0.525, 0.575, 0.625 so the sweep can resolve the region the
     PLOLO cosine says matters. The nine gammas 170 already ran are recomputed anyway and
     gated to be **bit-identical** to 170's cached effect files -- a fourth known-answer row
     that 170 could not have.

Rule 11: the operator is `data/donors/whitener_union_pruned.npz`,
sha256 `e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123`, 17,942 columns,
basis 29,394 x 17,942. Nine near-identically-named operators exist on disk; `optable212.py`
tabulates all of them with paths, sha256 and max eigenvalues. This script asserts the sha
before loading and refuses to run if it drifts.

    VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.10 uv run python gate212.py
"""
from __future__ import annotations

import hashlib
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
E170 = REPO / "experiments/170-gamma-mse/work"

os.environ.setdefault("VCC_ALPHA", "0.982")
os.environ.setdefault("VCC_CD4_TARGET_F", "0.10")
os.environ.setdefault("VCC_MASS_FIX", "minnorm")   # 170's setting; the effect stage is
                                                   # mass-fix independent (it enters only in
                                                   # `profile_from_effect`), and the a15/a19
                                                   # anchors below prove it.

OPERATOR_SHA = "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"
CAND_CHAL_SHA = "0601336fdf1ec5c533dcad69f6cb6c09d6368182bc8c69d10d516a1ecb19bf65"

ANCHORS = {
    (0.65, 0.0): {   # a15, `submissions/a15-cd4f010.vcc.meta.json`
        "before_operator_median_norm": 4.122892101713489,
        "after_operator_median_norm": 25.73792992728807,
        "operator_preserve_scale": 0.16018740098217024,
        "after_all_median_norm": 4.2292956044543,
        "cd4_f_per_unit": 1.055664536904655,
        "cd4_solved_fraction": 0.09472706196346517,
        "cd4_candidate_pre_preserve_median_norm": 4.166553356202115,
        "cd4_post_add_global_preserve_scale": 0.9895210139518232,
    },
    (0.80, 0.007): {  # a19, `submissions/a19-whiten-g080-lam0007.h5ad.meta.json`
        "before_operator_median_norm": 4.122892101713494,
        "after_operator_median_norm": 22.030919431001,
        "operator_preserve_scale": 0.1871411728696139,
        "after_all_median_norm": 4.235393246076997,
        "cd4_f_per_unit": 1.0541447101386758,
        "cd4_solved_fraction": 0.094863635929876,
    },
}
RTOL = 1e-11
WHITENER = REPO / "data/donors/whitener_union_pruned.npz"
WHITENER_META = REPO / "data/donors/whitener_union_pruned.json"

OLD = [0.00, 0.35, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80]
NEW = [0.30, 0.40, 0.525, 0.575, 0.625]
GRID = sorted(OLD + NEW)


def builder():
    path = REPO / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("prod_builder_212g", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_212g"] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    B = builder()
    pre = np.load(E170 / "pre.npz", allow_pickle=True)
    cd4z = np.load(E170 / "cd4.npz", allow_pickle=True)
    cols = pre["operator_columns"].astype(np.int64)
    panel_pos = pre["panel_positions"].astype(np.int64)
    h1_pos = pre["h1_positions"].astype(np.int64)
    cand = {"chal": pre["cand_chal"], "h1": pre["cand_h1"]}

    # ---- GATE 1 (rule 10): the cached pre-whitening object IS a15's, to the bit.
    got = B.array_sha256(cand["chal"])
    print(f"[gate] candidate_pre sha256 {got}")
    print(f"[gate]   a15's shipped audit {CAND_CHAL_SHA}   "
          f"{'MATCH' if got == CAND_CHAL_SHA else 'MISMATCH'}")
    if got != CAND_CHAL_SHA:
        raise SystemExit("candidate_pre is not a15's — STOP. (170: this only reproduces at "
                         "--chunk-rows 128, the builder's default; PerturbMatch's centring "
                         "mean is accumulated per chunk.)")

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

    n_h1_carry = int(np.count_nonzero(np.abs(cd4z["h1_mean"]).sum(axis=1)))
    print(f"[gate] CD4 support: challenge "
          f"{int(np.count_nonzero(np.abs(cd4z['chal_mean']).sum(axis=1)))}/300 targets, "
          f"h1-test {n_h1_carry}/85", flush=True)

    # ------------------------------------------------- GATE 2: two shipped board points
    print("\n===== KNOWN-ANSWER GATE 2: a15 (g 0.65/lam 0) and a19 (g 0.80/lam 0.007) =====")
    passed = True
    for (g, lam), want in ANCHORS.items():
        whitened, waudit = B.apply_v9_whitening(cand["chal"], whit, cols, gamma=g, lam=lam)
        _, caudit = B.apply_frozen_cd4_residual(whitened, cd4["chal"], cols, panel_pos)
        got_v = {
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
        print(f"  --- {'a15' if g == 0.65 else 'a19'}: gamma {g} / lam {lam}")
        for k, v in want.items():
            rel = abs(got_v[k] / v - 1.0) if v else abs(got_v[k])
            ok = rel <= RTOL
            passed &= ok
            print(f"      {k:<44} {'PASS' if ok else 'FAIL'}  {got_v[k]!r}  rel {rel:.2e}")
    if not passed:
        raise SystemExit("KNOWN-ANSWER GATE 2 FAILED — this is not a15/a19. STOP.")
    print("  ==> GATE 2 PASSED at BOTH points.\n", flush=True)

    # ------------------------------------------------------------------- the sweep
    audits = {}
    ident = {}
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
            out = WORK / "eff" / f"{surf}_g{gamma:.3f}.npz"
            np.savez_compressed(out, effect=final, gamma=gamma)
            row[surf] = {
                "final_sha256": B.array_sha256(final),
                "whitening": {k: float(v) for k, v in waudit.items()},
                "cd4_fraction": frac, "cd4_f_per_unit": fper,
                "post_add_preserve_scale": pscale,
                "median_row_norm": float(np.median(np.linalg.norm(final, axis=1))),
            }
            # GATE 3: for the nine gammas 170 already ran, the effect array must be
            # BIT-IDENTICAL to its cached file. 170 could not run this check.
            old = E170 / "eff" / f"{surf}_g{gamma:.2f}.npz"
            if gamma in OLD and old.exists():
                ref = np.load(old)["effect"]
                same = bool(np.array_equal(ref, final))
                ident[f"{surf}_g{gamma:.2f}"] = same
                if not same:
                    ident[f"{surf}_g{gamma:.2f}_maxabs"] = float(np.max(np.abs(ref - final)))
        audits[f"{gamma:.3f}"] = row
        print(f"[gate] gamma {gamma:.3f}: h1 median row norm "
              f"{row['h1']['median_row_norm']:.5f}, operator amplification "
              f"{row['h1']['whitening']['after_operator_median_norm']:.3f}, "
              f"chal median row norm {row['chal']['median_row_norm']:.5f}", flush=True)

    n_ok = sum(1 for k, v in ident.items() if k.count("_") == 1 and v is True)
    n_tot = sum(1 for k in ident if k.count("_") == 1)
    print(f"\n[gate] GATE 3 (bit-identity vs 170's cached effects): {n_ok}/{n_tot} identical")
    if n_ok != n_tot:
        print("       NOT identical on: "
              + ", ".join(k for k, v in ident.items() if v is False))
    with open(WORK / "gate212.json", "w") as fh:
        json.dump({"audits": audits, "identity_vs_170": ident,
                   "candidate_pre_sha256": got, "operator_sha256": sha,
                   "operator_basis": meta.get("basis"),
                   "peak_rss_gib": float(B.peak_rss_gib())}, fh, indent=1)
    print(f"[gate] wrote {WORK}/gate212.json and {len(GRID)*2} effect files; "
          f"peak RSS {B.peak_rss_gib():.2f} GiB")


if __name__ == "__main__":
    main()
