"""180 gate 2 — the per-context knee alpha, validated on A/B/C where the answer is known.

Runs the PRODUCTION knee code (`make_submission_pm_cd4.load_control_frame`,
`knee_context_terms`, `solve_knee_alphas_with_columns`) against the effect matrix
`effects180.py` proved bit-identical to a23's, so this is not a reimplementation of the
builder's arithmetic -- it is the builder's arithmetic, driven offline.

Three known-answer rows, each of which reads differently if the thing under test is absent
(AGENTS.md 11.1 rule 10):

  1. the effect matrix's own sha256 must equal a23's `effects.final_sha256`;
  2. the `C_p` proxy computed from the bundle's NTC constructs must reproduce 096's cached
     16-construct vector on its own 16 constructs, to 0.0;
  3. `A`, `B`, `claim` and `rho` at the shipped alpha must reproduce 170's `chal_minnorm.json`
     row at gamma 0.65, per context.

Only then does it compute the knee alphas, re-measure at them, and price `mse`.

    VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.10 VCC_MASS_FIX=minnorm uv run python knee180.py
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "work"

os.environ.setdefault("VCC_ALPHA", "0.982")
os.environ.setdefault("VCC_CD4_TARGET_F", "0.10")
os.environ.setdefault("VCC_MASS_FIX", "minnorm")

CTX = ("A", "B", "C")
A23_FINAL_SHA = "22c09f22c749bda82d12a0862bacf25b55e709b7de697513da3e3d11288a9b97"

# 170's own numbers at gamma 0.65 (work/chal_minnorm.json), the row this must reproduce.
EXPECT_170 = {
    "A": {"A": 0.2348299009419585, "B": 0.23359310250430518, "claim": 0.23543407548420198,
          "rho": 0.9921805160272124},
    "B": {"A": 0.22477784721536484, "B": 0.22362398630982533, "claim": 0.23784695947162313,
          "rho": 0.940201156267063},
    "C": {"A": 0.23494104987562134, "B": 0.23372755463537687, "claim": 0.2650257970493955,
          "rho": 0.8819049210964724},
}
# 096's published NTC C_p, and the D/E/F transition plan's expected knee alphas.
EXPECT_CP16 = {"A": 14.543547622910769, "B": 14.771622114915147, "C": 16.371526903064655}
EXPECT_KNEE_ALPHA = {"A": 0.987, "B": 1.013, "C": 1.047}

# 096/170's board-anchored `mse` constants.
SUMD = 0.229141
MSE_MAP = (0.987862879, 0.950162566)          # score = clip((b - raw)/span, 0, 1)
A23_SCORE_MSE = 0.121581

# a23's shipped render + emission policy.
MIN_CPM = 10.0
SPREAD_KWARGS = dict(tail_spread=1.0, keep_top=8000, keep_top_frac=None, keep_spread=0.25,
                     spread_map="block", spread_lo=0.0, spread_hi=0.60)
OVERDISPERSION = 4.0
CELLS = 400


def load_builder():
    spec = importlib.util.spec_from_file_location(
        "prod_builder_180_knee", REPO / "scripts/make_submission_pm_cd4.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_180_knee"] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    os.chdir(REPO)
    M = load_builder()
    if (M.ALPHA, M.MASS_FIX, M.CD4_TARGET_F) != (0.982, "minnorm", 0.10):
        raise RuntimeError("a23's env was not applied")

    effects = np.load(WORK / "effects-a23.npy")
    digest = M.array_sha256(effects)
    print(f"[gate] effect matrix sha256 {digest}")
    if digest != A23_FINAL_SHA:
        raise SystemExit(f"KNOWN-ANSWER FAIL: {digest} != a23's {A23_FINAL_SHA}")
    print("[gate] KNOWN-ANSWER 1 PASS: the matrix is a23's own final effect matrix")

    spec = M.ChallengeSpec.load(REPO / "data/challenge")
    cols = M.panel_gene_columns(spec)
    ctrl = np.load(REPO / "experiments/096-mse-floor/work/ctrl-stats.npz", allow_pickle=True)
    cached16 = dict(
        np.load(REPO / "experiments/096-mse-floor/work/cp-proxy-ntc.npz").items()
    )
    Kc = {c: 300 * float(ctrl[f"jk_{c}"]) / (len(spec.genes) - 1) for c in CTX}

    frames, cp_audit = {}, {}
    for context in CTX:
        t0 = time.time()
        frame = M.load_control_frame(
            REPO / "data/challenge", context, spec.genes, cells_per_pert=CELLS
        )
        frames[context] = frame
        # KNOWN-ANSWER 2: 096 used the FIRST 16 ntc constructs in file order. The builder
        # uses all 46 because there is no principled "first 16" for an unseen bundle, so the
        # gate is that the first 16 of ours reproduce the cached vector exactly.
        first16 = np.sort(frame.cp_values[:16])
        target = np.sort(np.asarray(cached16[context], dtype=np.float64))
        gap = float(np.abs(first16 - target).max())
        cp_audit[context] = {
            "cp_all": frame.cp, "n_constructs": int(frame.cp_values.size),
            "cp_first16": float(frame.cp_values[:16].mean()),
            "cp_096_cached16": float(np.mean(cached16[context])),
            "max_abs_diff_first16_vs_096": gap,
            "sd_all": float(frame.cp_values.std()),
        }
        print(f"[cp] {context}: all-{frame.cp_values.size} mean {frame.cp:.6f} "
              f"(sd {frame.cp_values.std():.4f}); first-16 {frame.cp_values[:16].mean():.6f} "
              f"vs 096's {np.mean(cached16[context]):.6f}; max|diff| {gap:.3e} "
              f"[{time.time() - t0:.0f}s]", flush=True)
        if gap > 1e-9:
            raise SystemExit(f"KNOWN-ANSWER FAIL: C_p proxy does not reproduce 096 on {context}")
    print("[gate] KNOWN-ANSWER 2 PASS: the NTC C_p proxy reproduces 096 exactly")

    terms_kwargs = dict(
        min_cpm=MIN_CPM, effect_scale=1.0, shrink_scales=None, shrink_center=None,
        rank_scores=None, spread_kwargs=SPREAD_KWARGS, overdispersion=OVERDISPERSION,
        cells_per_pert=CELLS,
    )

    # --- the shipped point, against 170 ------------------------------------------------
    #
    # 170 (and 096 before it) used the FIRST 16 ntc constructs. The builder uses all 46,
    # because "first 16" is not a rule that survives contact with a bundle nobody has seen,
    # and 46 has the smaller standard error. That moves `C_p` by -0.7% / +0.2% / +0.4% and
    # therefore moves `claim`, so the gate is decomposed: `A` and `B` are OUR side and must
    # reproduce 170 outright; `claim` is `C_p` x the row-weight sum whenever the #247 cap
    # binds, so it is gated at 096's OWN `C_p` -- which is the same comparison 170 made,
    # with the estimator held fixed.
    before, gate3 = {}, {}
    for context in CTX:
        t0 = time.time()
        before[context] = M.knee_context_terms(effects, frames[context], cols, **terms_kwargs)
        got, want = before[context], EXPECT_170[context]
        weight_sum = got["weight_sum"]
        chats = np.asarray(got["C_hat_values"], dtype=np.float64)
        claim096 = float(np.mean(np.minimum(chats, EXPECT_CP16[context])) * weight_sum)
        rho096 = got["B"] / claim096
        gate3[context] = {
            "weight_sum": weight_sum, "claim_at_096_cp": claim096, "rho_at_096_cp": rho096,
            "dA": abs(got["A"] - want["A"]), "dB": abs(got["B"] - want["B"]),
            "dclaim": abs(claim096 - want["claim"]), "drho": abs(rho096 - want["rho"]),
        }
        print(f"[shipped] {context}: A {got['A']:.6f} B {got['B']:.6f} "
              f"claim(all-46) {got['claim']:.6f} rho {got['rho_unclamped']:.6f} "
              f"cap binds {got['cap_saturated_rows']}/{got['cap_rows_measured']} | "
              f"at 096's C_p claim {claim096:.6f} rho {rho096:.6f} | vs 170 "
              + " ".join(f"d{k}={gate3[context]['d' + k]:.2e}"
                         for k in ("A", "B", "claim", "rho"))
              + f"  [{time.time() - t0:.0f}s]", flush=True)
        for key, tol in (("dA", 2e-4), ("dB", 2e-4), ("dclaim", 2e-4), ("drho", 2e-3)):
            if gate3[context][key] > tol:
                raise SystemExit(
                    f"KNOWN-ANSWER FAIL: {context}.{key} = {gate3[context][key]:.3e} > {tol}"
                )
    print("[gate] KNOWN-ANSWER 3 PASS: A / B reproduce 170 outright, and claim / rho "
          "reproduce it at 096's own C_p")

    # --- the knee ----------------------------------------------------------------------
    solved = M.solve_knee_alphas_with_columns(
        spec, effects, frames, cols, alpha_ref=0.982, **terms_kwargs
    )
    after, sensitivity = {}, {}
    for context in CTX:
        scale = solved[context]["knee_scale"]
        t0 = time.time()
        after[context] = M.knee_context_terms(
            effects * scale, frames[context], cols, **terms_kwargs
        )
        # The same knee at 096's `C_p` -- the estimator the plan's 0.987 / 1.013 / 1.047
        # were computed under -- plus the standard error the proxy itself carries. alpha
        # goes as sqrt(claim), so d(alpha)/alpha = 0.5 * d(C_p)/C_p.
        values = frames[context].cp_values
        se = float(values.std(ddof=1) / np.sqrt(values.size))
        alpha096 = 0.982 * float(np.sqrt(
            gate3[context]["claim_at_096_cp"] / before[context]["B"]
        ))
        sensitivity[context] = {
            "alpha_all46": solved[context]["alpha"],
            "alpha_at_096_cp": alpha096,
            "C_p_all46": frames[context].cp,
            "C_p_standard_error": se,
            "C_p_relative_se": se / frames[context].cp,
            "alpha_relative_se": 0.5 * se / frames[context].cp,
            "alpha_1se": 0.5 * se / frames[context].cp * solved[context]["alpha"],
        }
        print(f"[knee] {context}: alpha {solved[context]['alpha']:.6f} (all-46 C_p) / "
              f"{alpha096:.6f} (096's C_p; plan expected ~{EXPECT_KNEE_ALPHA[context]}) "
              f"+-{sensitivity[context]['alpha_1se']:.4f} (1 s.e. of C_p) -> rho_unclamped "
              f"{after[context]['rho_unclamped']:.6f} (was "
              f"{before[context]['rho_unclamped']:.6f})  [{time.time() - t0:.0f}s]",
              flush=True)

    # --- the delivered geometry: `X`'s multiplier, and `pds`'s invariance ---------------
    #
    # `X = sum_p w_p <dh_p, D_p>` is LINEAR in the delivered effect, so it does not survive
    # an amplitude change unchanged -- holding it fixed while `A` and `B` move by s^2 is the
    # single easiest way to get this pricing wrong. What multiplies it is the per-row ratio
    # `kappa_p = <dh_p(s), dh_p(1)> / ||dh_p(1)||^2`, which is a constant times the identity
    # only to the extent the renderer is homogeneous (a23's minimum-norm mass fix is not).
    #
    # The SAME two numbers settle `pds`. The scored member is `pds_cosine` with
    # `exclusion_scope="panel"`: every distance it ranks is a cosine between our delivered
    # delta row and a reference delta row, in one fixed feature space with all 300 panel
    # target genes dropped. A cosine is invariant to a positive per-row scale, so if
    # `cos(dh_p(s), dh_p(1)) == 1` on the panel-excluded axis for every row, then EVERY
    # entry of the n x n distance matrix is unchanged, every rank is unchanged, and
    # `d pds = 0` exactly -- for any reference, without needing one. This is a proof, and
    # `1 - cos` is the size of the hole in it.
    panel_mask = np.ones(len(spec.genes), dtype=bool)
    panel_mask[np.asarray(sorted(set(cols.values())), dtype=np.int64)] = False
    print(f"\n[geom] panel-excluded feature space: {int(panel_mask.sum())} of "
          f"{len(spec.genes)} genes ({len(spec.genes) - int(panel_mask.sum())} panel "
          f"targets dropped, as pds_cosine exclusion_scope='panel' does)")

    def delivered_delta(matrix, frame):
        rows = np.empty_like(matrix)
        for index in range(matrix.shape[0]):
            rows[index] = np.log1p(M.TS * M.profile_from_effect(
                matrix[index], frame.b_ctrl, min_cpm=MIN_CPM
            )) - frame.b_ctrl
        return rows

    geometry = {}
    for context in CTX:
        t0 = time.time()
        base = delivered_delta(effects, frames[context])
        moved = delivered_delta(effects * solved[context]["knee_scale"], frames[context])
        num = np.einsum("ij,ij->i", moved, base)
        den = np.einsum("ij,ij->i", base, base)
        kappa = num / den
        cos_full = num / np.sqrt(den * np.einsum("ij,ij->i", moved, moved))
        bp, mp = base[:, panel_mask], moved[:, panel_mask]
        cos_panel = (
            np.einsum("ij,ij->i", mp, bp)
            / np.sqrt(np.einsum("ij,ij->i", bp, bp) * np.einsum("ij,ij->i", mp, mp))
        )
        geometry[context] = {
            "knee_scale": solved[context]["knee_scale"],
            "kappa_mean": float(kappa.mean()), "kappa_min": float(kappa.min()),
            "kappa_max": float(kappa.max()), "kappa_sd": float(kappa.std()),
            "kappa_over_scale_max_dev": float(
                np.abs(kappa / solved[context]["knee_scale"] - 1.0).max()
            ),
            "cos_full_min": float(cos_full.min()),
            "one_minus_cos_full_max": float((1.0 - cos_full).max()),
            "cos_panel_excluded_min": float(cos_panel.min()),
            "one_minus_cos_panel_excluded_max": float((1.0 - cos_panel).max()),
            "rows": int(kappa.size),
        }
        del base, moved, bp, mp
        print(f"[geom] {context}: kappa {kappa.mean():.6f} +-{kappa.std():.2e} "
              f"(scale {solved[context]['knee_scale']:.6f}, max dev "
              f"{geometry[context]['kappa_over_scale_max_dev']:.2e}); "
              f"max(1-cos) full {geometry[context]['one_minus_cos_full_max']:.3e}, "
              f"panel-excluded "
              f"{geometry[context]['one_minus_cos_panel_excluded_max']:.3e} "
              f"[{time.time() - t0:.0f}s]", flush=True)

    # --- `mse`, through 096's claim route ----------------------------------------------
    target_raw = MSE_MAP[0] - A23_SCORE_MSE * MSE_MAP[1]
    num_before = float(np.mean([
        before[c]["A"] + Kc[c] - min(before[c]["claim"], before[c]["B"]) for c in CTX
    ]))
    X = (num_before - SUMD * (target_raw - 1.0)) / 2.0
    num_after = float(np.mean([
        after[c]["A"] + Kc[c] - min(after[c]["claim"], after[c]["B"]) for c in CTX
    ]))
    x_after = float(np.mean([geometry[c]["kappa_mean"] * X for c in CTX]))
    raw_a23 = 1.0 + (num_before - 2.0 * X) / SUMD
    raw_knee = 1.0 + (num_after - 2.0 * x_after) / SUMD
    # The same arithmetic with `X` held fixed, so the size of that error is on the record.
    raw_knee_frozen_x = 1.0 + (num_after - 2.0 * X) / SUMD
    d_score = (raw_a23 - raw_knee) / MSE_MAP[1]
    print(f"\n[mse] board anchor a23 score_mse {A23_SCORE_MSE} -> raw {target_raw:.6f}; "
          f"solved X {X:.6f} (170 solved {0.017931763405396428:.6f})")
    print(f"[mse] X at the knee {x_after:.6f} (mean kappa "
          f"{np.mean([geometry[c]['kappa_mean'] for c in CTX]):.6f})")
    print(f"[mse] raw a23 {raw_a23:.6f} -> knee {raw_knee:.6f}   "
          f"d score_mse {d_score:+.6f}   d overall {d_score / 6:+.6f}")
    print(f"[mse] with X WRONGLY held fixed: raw {raw_knee_frozen_x:.6f}, d overall "
          f"{(raw_a23 - raw_knee_frozen_x) / MSE_MAP[1] / 6:+.6f}  <- not the answer")
    for context in CTX:
        n0 = before[context]["A"] + Kc[context] - min(
            before[context]["claim"], before[context]["B"])
        n1 = after[context]["A"] + Kc[context] - min(
            after[context]["claim"], after[context]["B"])
        d_ctx = ((n1 - 2 * geometry[context]["kappa_mean"] * X) - (n0 - 2 * X)) / SUMD
        print(f"[mse]   {context}: num {n0:+.6f} -> {n1:+.6f}, d raw mse {d_ctx:+.6f}, "
              f"d score_mse {-d_ctx / MSE_MAP[1]:+.6f}")

    out = {
        "effects_sha256": digest,
        "cp": cp_audit,
        "Kc": Kc,
        "gate3_vs_170": gate3,
        "cp_sensitivity": sensitivity,
        "per_context": {
            c: {
                "before": before[c],
                "knee_scale": solved[c]["knee_scale"],
                "alpha": solved[c]["alpha"],
                "after": after[c],
                "expected_alpha_plan": EXPECT_KNEE_ALPHA[c],
            }
            for c in CTX
        },
        "geometry": geometry,
        "mse": {
            "sumD": SUMD, "X": X, "X_at_knee": x_after,
            "raw_a23": raw_a23, "raw_knee": raw_knee,
            "raw_knee_with_X_frozen": raw_knee_frozen_x,
            "d_score_mse": d_score, "d_overall": d_score / 6,
            "board_anchor_raw": target_raw,
        },
    }
    (WORK / "knee180.json").write_text(json.dumps(out, indent=1))

    fails = []
    for context in CTX:
        if abs(after[context]["rho_unclamped"] - 1.0) > 0.01:
            fails.append(f"{context} rho_unclamped {after[context]['rho_unclamped']:.4f}")
        # The plan's 0.987 / 1.013 / 1.047 were computed under 096's own `C_p`, so that is
        # the number the plan's tolerance applies to. The shipped alpha uses all 46
        # constructs and is reported with its own C_p standard error alongside.
        if abs(sensitivity[context]["alpha_at_096_cp"] - EXPECT_KNEE_ALPHA[context]) > 0.002:
            fails.append(
                f"{context} alpha(096 C_p) {sensitivity[context]['alpha_at_096_cp']:.4f} "
                f"vs plan {EXPECT_KNEE_ALPHA[context]}"
            )
    print("\n[gate2] " + ("PASS" if not fails else "FAIL: " + "; ".join(fails)))
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
