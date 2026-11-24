"""213 Phase 3 — the six members on h1-test, per transport law.

`experiments/170-gamma-mse/h1sweep.py`'s `main`, transcribed with the law replacing gamma and
the RENDER replaced by the law's own; every scoring primitive (`msekit.score`,
`decompose.terms`, `h1kit.pds_exact`, `h1kit.emit_and_score`, `debench.Surface`) is imported
from 170/153/096 and unchanged.

Donors on h1-test are `k562_gwps` + `orion_hek293t` + `orion_hct116` -- production's own three
-- and `data/donors/h1_2025`, which covers 85/85 of the targets, is NOT among them (6.4's
oracle-leakage trap).

`nmae` is `debench.Surface.nmae`, whose pool is the reference's own BH < 0.05 significant set
with at least 10 genes per row: the significance-selected pool 188 identified as the
challenge's construction.  A top-N-by-magnitude pool INVERTS its sign.

Rule 9: h1-test's `n_real` is 765 mean / ~90 median against the challenge's ~89 median, so
`fid`, `jac` and `reach`'s normalisation are read as DELTAS ONLY and only as ordering evidence.

KNOWN-ANSWER GATE: `msekit.score` on the frozen 095 h1-test artifact must return
2.0551018299736388 (170's own gate, `cell-eval2`'s `agg_results.csv` value), and `logc20` must
reproduce 170's own gamma-0.65 row.

    uv run python h1sweep213.py --laws logc20 logc40
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "work"
E170 = REPO / "experiments/170-gamma-mse"
sys.path.insert(0, str(E170))
sys.path.insert(0, str(HERE))

import h1kit  # noqa: E402
for _n in ("de_real_h1.npz",):
    if not (WORK / _n).exists():
        shutil.copyfile(E170 / "work" / _n, WORK / _n)
h1kit.WORK = WORK

import h1sweep as H  # noqa: E402
H.WORK = WORK

import law213 as L  # noqa: E402

# 170's own gamma-0.65 row: the known-answer target for `logc20`.
KA170 = dict(mse=1.0159689545093225, pds=0.8753501400560225, nmae=0.959881122307442,
             reach=0.17751503650931924, fid=0.5081565598953068, jac=0.06424535489807914,
             A=0.08349223403695392, X=0.007733929547645832, B=0.08204237493686707,
             claim=0.06700840536519416, cos_deconf=0.06841128111448559)


def deliver_law(law, effect_h1, b_ctrl, gate, mass_fix):
    """`h1kit.deliver`, with the law's render and the law's mass fix."""
    p_ctrl = np.expm1(np.asarray(b_ctrl, dtype=np.float64)) / h1kit.TS
    P = np.empty_like(effect_h1)
    for i in range(effect_h1.shape[0]):
        g = np.where(gate, effect_h1[i], 0.0)
        if mass_fix is not None:
            g = L.mass_consistent_law(g, p_ctrl, b_ctrl, gate, law, mode=mass_fix)
        p = np.clip(law.render(g, p_ctrl, b_ctrl), 0.0, None)
        P[i] = p / float(p.sum())
    return P, np.log1p(h1kit.TS * P) - b_ctrl[None, :]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--laws", nargs="*", default=["logc20", "logc40"])
    ap.add_argument("--seed", type=int, default=0)
    # BETA == 0 in the production builder, so the staged effect is EXACTLY linear in
    # VCC_ALPHA (the whitener renormalises to the pre-whitening median row norm, which is
    # scale-equivariant, and the CD4 fraction is a ratio of norms, which is scale-
    # invariant).  `--escale s` therefore reproduces alpha = 0.982*s exactly, without
    # re-staging.
    ap.add_argument("--escale", type=float, default=1.0)
    args = ap.parse_args()

    from debench import Surface
    from decompose import terms
    from msekit import score

    H.known_answer_gate()
    B = h1kit.builder()                      # forces VCC_MASS_FIX=minnorm (a23 geometry)
    S = h1kit.Surface170()
    laws = L.all_laws()

    de = h1kit.DEKit()
    real = de.real_table(S.targets)
    de.release_reference()
    nr = np.asarray([(real["real_padj"][i] < 0.05).sum() for i in range(S.targets.size)])
    print(f"[h1] {S.targets.size} targets, {int(S.gate.sum())} gated genes; reference "
          f"n_real median {np.median(nr):.0f} mean {nr.mean():.1f}", flush=True)

    out_path = WORK / "h1sweep213.json"
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["law"], r.get("escale", 1.0), r.get("seed", 0)) for r in rows}

    for key in args.laws:
        if (key, args.escale, args.seed) in done:
            continue
        law = laws[key]
        t0 = time.time()
        z = np.load(WORK / f"pre_{key}.npz", allow_pickle=True)
        eff_prod = np.load(WORK / "eff" / f"h1_{key}.npz")["effect"]
        eff, have = h1kit.to_h1_axis(eff_prod, z["genes"], z["h1_genes"])
        pos = {str(t): i for i, t in enumerate(z["h1_targets"])}
        eff = eff[np.asarray([pos[str(t)] for t in S.targets])] * args.escale

        profiles, dh = deliver_law(law, eff, S.b_ctrl, S.gate, B.MASS_FIX)
        sums, jk, (pp_, pq_, pl_) = h1kit.emit_and_score(
            B, profiles, eff, S.b_ctrl, S.libs, de, seed=args.seed)
        mse, _, _ = score(S.targets, sums, np.full(S.targets.size, float(h1kit.NCELL)), jk,
                          S.rp, S.rs, S.rc, S.rj, S.genes)
        t = terms(S.targets, sums, np.full(S.targets.size, float(h1kit.NCELL)), jk,
                  S.rp, S.rs, S.rc, S.rj, S.genes)
        pds = h1kit.pds_exact(S.targets, sums, S)
        np.savez_compressed(WORK / f"de_{key}_e{args.escale:g}_s{args.seed}.npz", genes=de.genes,
                            perts=S.targets, **real, pred_p=pp_, pred_padj=pq_, pred_lfc=pl_)
        surf = Surface(str(WORK / f"de_{key}_e{args.escale:g}_s{args.seed}.npz"))
        nrl, npd = surf.counts()
        row = dict(law=key, seed=args.seed, escale=args.escale,
                   alpha=0.982 * args.escale, mse=float(mse), pds=float(pds),
                   nmae=float(surf.nmae()), reach=float(surf.reach()),
                   fid=float(surf.fid()), jac=float(surf.jac()),
                   A=float(t["A"]), X=float(t["X"]), B=float(t["budget"]),
                   claim=float(t["claim"]), rho=float(t["rho"]), sumD=float(t["sumD"]),
                   R=float(t["R"]), K=float(t["K"]), Kc=float(t["Kc"]),
                   cos_deconf=float(t["X"] / np.sqrt(t["A"] * t["sumD"])),
                   n_pred=float(np.mean(npd)), n_real=float(np.mean(nrl)),
                   sign_acc=float(surf.sign_accuracy()),
                   row_norm=float(np.median(np.linalg.norm(dh, axis=1))),
                   secs=round(time.time() - t0, 1))
        rows.append(row)
        out_path.write_text(json.dumps(rows, indent=1))
        print(f"[h1] {key}  mse {row['mse']:.5f}  pds {row['pds']:.5f}  nmae "
              f"{row['nmae']:.5f}  reach {row['reach']:.5f}  fid {row['fid']:.5f}  jac "
              f"{row['jac']:.5f} | A {row['A']:.5f} X {row['X']:.6f} cos "
              f"{row['cos_deconf']:.6f} n_pred {row['n_pred']:.0f} [{row['secs']:.0f}s]",
              flush=True)
        if key in ("logc20", "bc0") and args.escale == 1.0:
            print("\n[h1] KNOWN-ANSWER: `logc20` == 170's own gamma-0.65 row")
            ok = True
            for k, want in KA170.items():
                rel = abs(row[k] / want - 1.0)
                good = rel <= 1e-9
                ok &= good
                print(f"   {'PASS' if good else 'FAIL'}  {k:<12} {row[k]!r}  want {want!r}"
                      f"  rel {rel:.2e}")
            if not ok:
                print("   NOTE: any deviation here is the h1-test instrument, not the law; "
                      "it must be explained before the law rows are read.")
        del profiles, dh, sums, eff, eff_prod


if __name__ == "__main__":
    main()
