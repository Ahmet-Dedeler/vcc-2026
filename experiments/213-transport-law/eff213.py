"""213 — apply the SHIPPED operator + CD4 residual to a law's staged `candidate_pre`.

This is `experiments/170-gamma-mse/gate.py` with the law, not gamma, as the free variable:
gamma stays at production's 0.65 and lam at 0.0.  Everything is the production builder's own
`apply_v9_whitening` / `apply_frozen_cd4_residual`.

KNOWN-ANSWER GATE.  At `logc20` the four scalars a15 recorded in its own meta must reproduce:

    before_operator_median_norm   4.122892101713489
    after_operator_median_norm   25.73792992728807
    operator_preserve_scale       0.16018740098217024
    after_all_median_norm         4.2292956044543
    cd4_solved_fraction           0.09472706196346517

`final_sha256` is NOT gateable (a15 was built on x86_64, this is arm64, and `whitener.apply`
is a pair of 17,942^2 GEMMs whose reduction order differs) -- 170 established that and the
scalars agree to ~2e-16, which is last-bit GEMM noise and not a different computation.

    uv run python eff213.py --laws logc20 logc40
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
(WORK / "eff").mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REPO / "experiments/170-gamma-mse"))
sys.path.insert(0, str(HERE))

os.environ.setdefault("VCC_ALPHA", "0.982")
os.environ.setdefault("VCC_CD4_TARGET_F", "0.10")
os.environ.setdefault("VCC_MASS_FIX", "minnorm")

import gate as G  # noqa: E402  (170's gate module; its `main` is never called)

GAMMA, LAM = 0.65, 0.0
ANCHOR = G.ANCHORS[(0.65, 0.0)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--laws", nargs="*", default=["logc20", "logc40"])
    args = ap.parse_args()

    B = G.builder()
    sha = B.file_sha256(G.WHITENER)
    if sha != G.OPERATOR_SHA:
        raise SystemExit(f"operator sha {sha} != the shipped {G.OPERATOR_SHA}")
    whit = B.Whitener.load(G.WHITENER)
    cd4z = np.load(WORK / "cd4.npz", allow_pickle=True)
    mode = str(cd4z["mode"][0])
    audits = json.loads((WORK / "eff213.json").read_text()) \
        if (WORK / "eff213.json").exists() else {}

    for key in args.laws:
        src = WORK / f"pre_{key}.npz"
        if not src.exists():
            print(f"[eff213] {key}: no staged pre_{key}.npz — run stage213.py first")
            continue
        pre = np.load(src, allow_pickle=True)
        cols = pre["operator_columns"].astype(np.int64)
        panel_pos = pre["panel_positions"].astype(np.int64)
        h1_pos = pre["h1_positions"].astype(np.int64)
        cd4 = {
            "chal": B.EXP062.CD4Prepared(
                mode=mode, mean_effect=cd4z["chal_mean"], gene_votes=cd4z["chal_votes"],
                target_condition_count=cd4z["chal_counts"], condition_details=[]),
            "h1": B.EXP062.CD4Prepared(
                mode=mode, mean_effect=cd4z["h1_mean"], gene_votes=cd4z["h1_votes"],
                target_condition_count=cd4z["h1_counts"], condition_details=[]),
        }
        rec = {}
        for frame, cand, pos in (("chal", pre["cand_chal"], panel_pos),
                                 ("h1", pre["cand_h1"], h1_pos)):
            out = WORK / "eff" / f"{frame}_{key}.npz"
            whitened, waudit = B.apply_v9_whitening(cand, whit, cols, gamma=GAMMA, lam=LAM)
            if frame == "chal":
                final, caudit = B.apply_frozen_cd4_residual(whitened, cd4["chal"], cols,
                                                            panel_pos)
                frac = float(caudit["residual_frame"]["solved_fraction"])
                rec.update({k: float(waudit[k]) for k in (
                    "before_operator_median_norm", "after_operator_median_norm",
                    "operator_preserve_scale", "after_all_median_norm")})
                rec["cd4_solved_fraction"] = frac
                rec["cd4_f_per_unit"] = float(
                    caudit["residual_frame"]["f_per_unit_fraction"])
            else:
                # 170: the CD4 residual covers the challenge panel, not h1-test (2/85
                # targets survive the frozen `complete mean3` policy and the median row
                # norm is 0, so `solve_cd4_fraction` cannot solve at all).  h1-test is
                # built with the fraction production SOLVES ON THE PANEL, through the same
                # shipped `add_post_whiten_residual`.
                # 170's own h1 branch does NOT zero the panel columns here
                # (`apply_frozen_cd4_residual` does that for the challenge frame only),
                # so this reproduces its effect files exactly.
                final, _ = B.EXP062.add_post_whiten_residual(
                    whitened, cd4["h1"], rec["cd4_solved_fraction"], cols)
            np.savez_compressed(out, effect=final)
            rec[f"{frame}_final_sha256"] = B.array_sha256(final)
            rec[f"{frame}_median_row_norm"] = float(np.median(np.linalg.norm(final, axis=1)))
            del whitened, final
            gc.collect()
        audits[key] = rec
        (WORK / "eff213.json").write_text(json.dumps(audits, indent=1))
        if key in ("logc20", "bc0"):
            print(f"\n[eff213] KNOWN-ANSWER GATE for {key} (a15's own meta):")
            ok = True
            for k, want in ANCHOR.items():
                if k == "final_sha256":
                    continue
                if k not in rec:
                    continue
                rel = abs(rec[k] / want - 1.0)
                good = rel <= 1e-11
                ok &= good
                print(f"   {'PASS' if good else 'FAIL'}  {k:<36} {rec[k]!r}  "
                      f"want {want!r}  rel {rel:.2e}")
            print(f"   final_sha256 {rec['chal_final_sha256']} "
                  f"(a15 {ANCHOR['final_sha256']}; arm64 vs x86_64, NOT gateable)")
            if not ok:
                raise SystemExit("a15 OPERATOR GATE FAILED. STOP.")
        print(f"[eff213] {key}: after_all_median_norm {rec['after_all_median_norm']:.9f}  "
              f"cd4_f {rec['cd4_solved_fraction']:.8f}  "
              f"chal row norm {rec['chal_median_row_norm']:.6f}", flush=True)
        del pre
        gc.collect()


if __name__ == "__main__":
    main()
