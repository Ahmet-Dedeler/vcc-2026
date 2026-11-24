"""170 — the six-member gamma sweep on h1-test, end to end, nothing frozen.

At every gamma the production effect matrix (`work/eff/h1_g*.npz`, built by `gate.py` from
the bit-exact `candidate_pre`) is rendered through the shipped `profile_from_effect`
(CPM gate 10, MASS_FIX=minnorm), 400 cells are emitted per target through the shipped
`emit_counts` under a15's tier-8000 policy, the delete-1 jackknife is recomputed on those
cells, and all six members are scored:

  mse    096's `claim` route -- `cell_eval2.metrics.delta.mse_unbiased_capped` /
         `distance_unbiased` on the emitted group sums, with the full A / X / B / claim /
         rho / sumD decomposition. NOT a cosine proxy.
  pds    `cell_eval2.metrics.discrimination_score`, vcc2026 settings.
  nmae   153's `debench.Surface.nmae` -- exact, deterministic in the mean profile.
  reach  }  153's `debench.Surface` on the realised DE tables. These carry `n_real`
  fid    }  (rule 9), so ONLY the deltas are reported and only as ordering evidence.
  jac    }

KNOWN-ANSWER GATE (rule 10): before any arm is scored, `msekit.score` is run on the FROZEN
095 h1-test artifact and must return 2.0551018299736388, the value the full
`031-somatic -> cell-eval2 run` pipeline wrote to
`095/out-h1test/.../run/agg_results.csv`. That anchors this exact mse code path to a real
scorer run. It aborts if it does not reproduce.

    uv run python h1sweep.py [gamma,gamma,...]
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from h1kit import (  # noqa: E402
    H1, NCELL, TS, DEKit, Surface170, builder, deliver, emit_and_score, pds_exact,
    to_h1_axis)

REPO = HERE.parents[1]
WORK = HERE / "work"

MSE_ANCHOR = 2.0551018299736388
FROZEN_PRED = str(REPO / "experiments/095-clean-harness/out-h1test/runs/"
                         "g0.65_q1.0_a1.08_spread0.25/pred.h5ad")


def known_answer_gate():
    from msekit import group_stats, score

    rp, rs, rc, rj, rg = group_stats(H1)
    pp, ps, pc, pj, pg = group_stats(FROZEN_PRED)
    if not np.array_equal(rg, pg):
        raise SystemExit("frozen 095 prediction and h1-test gene axes differ")
    m, _, _ = score(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    rel = abs(m / MSE_ANCHOR - 1.0)
    print(f"[gate] msekit on the frozen 095 h1-test artifact -> {m!r}\n"
          f"       cell-eval2's own agg_results.csv           -> {MSE_ANCHOR!r}   "
          f"rel {rel:.2e}", flush=True)
    if rel > 1e-12:
        raise SystemExit("MSE KNOWN-ANSWER GATE FAILED — the mse route is not the scorer's.")
    print("  ==> mse route GATE PASSED.\n", flush=True)


def main(grid, seed=0, out_name="h1sweep.json") -> None:
    from debench import Surface
    from decompose import terms
    from msekit import score

    known_answer_gate()
    B = builder()
    S = Surface170()
    pre = np.load(WORK / "pre.npz", allow_pickle=True)
    prod_genes = pre["genes"]
    h1_genes = pre["h1_genes"]
    print(f"[h1] {S.targets.size} targets, {int(S.gate.sum())} gated genes "
          f"(CPM > 10 on h1-test's own control), emitted library median "
          f"{np.median(S.libs):.0f}", flush=True)

    de = DEKit()
    real = de.real_table(S.targets)
    de.release_reference()
    nr = np.asarray([(real["real_padj"][i] < 0.05).sum() for i in range(S.targets.size)])
    print(f"[h1] DE universe {de.genes.size} genes (control CPM > 5); reference n_real "
          f"median {np.median(nr):.0f} mean {nr.mean():.1f}  -- RULE 9: fid/jac/reach are "
          f"reported as DELTAS only", flush=True)

    rows = []
    for gamma in grid:
        t0 = time.time()
        eff_prod = np.load(WORK / "eff" / f"h1_g{gamma:.2f}.npz")["effect"]
        eff, have = to_h1_axis(eff_prod, prod_genes, h1_genes)
        # order the rows the same way the reference groups them
        pos = {str(t): i for i, t in enumerate(pre["h1_targets"])}
        order = np.asarray([pos[str(t)] for t in S.targets])
        eff = eff[order]

        profiles, dh = deliver(B, eff, S.b_ctrl)
        sums, jk, (pp_, pq_, pl_) = emit_and_score(
            B, profiles, eff, S.b_ctrl, S.libs, de, seed=seed)

        pred_perts = S.targets
        mse, _, _ = score(pred_perts, sums, np.full(pred_perts.size, float(NCELL)), jk,
                          S.rp, S.rs, S.rc, S.rj, S.genes)
        t = terms(pred_perts, sums, np.full(pred_perts.size, float(NCELL)), jk,
                  S.rp, S.rs, S.rc, S.rj, S.genes)
        pds = pds_exact(pred_perts, sums, S)

        np.savez_compressed(WORK / f"de_g{gamma:.2f}_s{seed}.npz",
                            genes=de.genes, perts=S.targets, **real,
                            pred_p=pp_, pred_padj=pq_, pred_lfc=pl_)
        surf = Surface(str(WORK / f"de_g{gamma:.2f}_s{seed}.npz"))
        nr, npd = surf.counts()

        row = dict(
            gamma=gamma, mse=float(mse), pds=float(pds),
            nmae=float(surf.nmae()), reach=float(surf.reach()),
            fid=float(surf.fid()), jac=float(surf.jac()),
            A=float(t["A"]), X=float(t["X"]), B=float(t["budget"]),
            claim=float(t["claim"]), rho=float(t["rho"]), sumD=float(t["sumD"]),
            R=float(t["R"]), K=float(t["K"]), Kc=float(t["Kc"]),
            cos_agg=float(t["X"] / np.sqrt(t["A"] * t["R"])),
            cos_deconf=float(t["X"] / np.sqrt(t["A"] * t["sumD"])),
            n_pred=float(np.mean(npd)), n_real=float(np.mean(nr)),
            sign_acc=float(surf.sign_accuracy()),
            row_norm=float(np.median(np.linalg.norm(dh, axis=1))),
            rr=float(_rowrow(dh)),
            secs=round(time.time() - t0, 1),
        )
        rows.append(row)
        print(f"[h1] gamma {gamma:.2f}  mse {row['mse']:.5f}  pds {row['pds']:.5f}  "
              f"nmae {row['nmae']:.5f}  reach {row['reach']:.5f}  fid {row['fid']:.5f}  "
              f"jac {row['jac']:.5f}  | A {row['A']:.5f} X {row['X']:.6f} "
              f"B {row['B']:.5f} claim {row['claim']:.5f} rho {row['rho']:.4f} "
              f"[{row['secs']:.0f}s]", flush=True)
        with open(WORK / out_name, "w") as fh:
            json.dump(rows, fh, indent=1)
    print(f"[h1] wrote {WORK}/{out_name}; peak RSS {B.peak_rss_gib():.2f} GiB")


def _rowrow(D):
    U = D / np.maximum(np.linalg.norm(D, axis=1, keepdims=True), 1e-300)
    G = U @ U.T
    n = G.shape[0]
    return G[~np.eye(n, dtype=bool)].mean()


if __name__ == "__main__":
    g = ([float(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1
         else [0.35, 0.45, 0.55, 0.65, 0.75, 0.80])
    sd = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    nm = sys.argv[3] if len(sys.argv) > 3 else ("h1sweep.json" if sd == 0
                                                else f"h1sweep_s{sd}.json")
    main(g, sd, nm)
