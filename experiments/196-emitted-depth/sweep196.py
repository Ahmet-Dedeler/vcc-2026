"""196 Phase 1 — the EMITTED DEPTH sweep on h1-test, through the production emitters.

One variable moves: `L`, the per-cell library the emitter writes.  Everything else is a15's
shipped emission policy, held bit-identical across depths:

  * the spread vector's rank axis (`|log2(pcpm/cpm)|`) is computed ONCE at m = 1 and reused,
    so `keep_top 8000 / keep_spread 0.25 / tail_spread 1.0` selects the SAME 8,000 genes at
    every depth -- a depth arm, not a depth-and-policy arm;
  * `overdispersion` 4.0, `DISPERSION` 0.1, and the same per-target seed;
  * the intended profile is 196's reconstruction (`recon196.Recon`), whose `m = 1` row
    re-delivers a15's own per-gene panel totals BIT-IDENTICALLY -- asserted per target, and
    the run aborts if it ever fails.

KNOWN-ANSWER ROWS (rule 10, checked by `six196.py`): the `m = 1` quasi arm must reproduce
185/191's a15 row on h1-test -- reach 0.18283, median n_pred 6035 (raw BH) / 6034 (debench),
jac 0.06824, nmae 0.96659 -- and the `m = 1` rc_hybrid arm must reproduce reach 0.21681,
n_pred 6043, jac 0.06827, nmae 0.96670.

    DEPTHS=1,2,4,8 ARMS=quasi,rc_hybrid MAXT=8 uv run python sweep196.py
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import false_discovery_control, norm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/185-rank-controlled-emission"))
sys.path.insert(0, str(ROOT / "src"))
import harness185 as H                                     # noqa: E402
import recon196                                            # noqa: E402
import vcc2026.rc_emit as RC                                # noqa: E402
from vcc2026.donor_transfer import emit_counts              # noqa: E402

WORK = HERE / "work"
n, PHI, DISPERSION = 400, 4.0, 0.1
DEPTHS = tuple(float(x) for x in os.environ.get("DEPTHS", "1,2,4,8").split(","))
ARMS = tuple(os.environ.get("ARMS", "quasi,rc_hybrid").split(","))
MAXT = int(os.environ.get("MAXT", "0"))
TAG = os.environ.get("OUTTAG", "de196")


def to_block(vals, cols, ptr, G):
    B = np.zeros((n, G), np.int64)
    for r in range(n):
        sl = slice(ptr[r], ptr[r + 1])
        B[r, cols[sl]] = vals[sl]
    return B


def main() -> None:
    t0 = time.time()
    C = H.load_controls()
    cs, tie, keep, cpm, ref_lin = C["cs"], C["tie"], C["keep"], C["cpm"], C["ref_lin"]
    genes_full = C["genes"]
    G_full, G = genes_full.size, cs.shape[0]
    gated_scorer = keep                                   # the SCORER's 5-CPM gate
    control = RC.RcControlBlock(keep=keep, cs=cs, tie=tie, ref_lin=ref_lin,
                                n_controls=cs.shape[1], min_cpm=H.MIN_CPM)

    D = recon196.load_h1()
    assert np.array_equal(D["genes"], genes_full)
    targets = list(D["targets"])
    sums = D["sums"]
    Hd = H.load_H()
    hpos = {t: i for i, t in enumerate(Hd["targets"])}
    if MAXT:
        targets, sums = targets[:MAXT], sums[:MAXT]
    T_n = len(targets)
    keys = [f"{a}_m{m:g}" for a in ARMS for m in DEPTHS]
    print(f"h1-test: {cs.shape[1]} controls, {G}/{G_full} scorer-gated, "
          f"{int(D['gate'].sum())} builder-gated, {T_n} perturbations; arms {keys} "
          f"[{time.time()-t0:.0f}s]", flush=True)

    real = {k: np.full((T_n, G), np.nan) for k in ("p", "q", "l")}
    out = {k: {j: np.full((T_n, G), np.nan) for j in ("p", "q", "l", "z")} for k in keys}
    absH_all = np.zeros((T_n, G), np.float32)
    diag = {k: [] for k in keys}
    npred = {k: [] for k in keys}
    digest = {k: [] for k in keys}
    failed = {}

    import anndata as ad
    R = ad.read_h5ad(H.REAL, backed="r")
    rlab = np.asarray(R.obs[H.COL].astype(str).values, dtype="U64")
    R.file.close()
    sys.path.insert(0, str(ROOT / "experiments/185-rank-controlled-emission"))
    from realize185 import RowReader                        # noqa: E402
    rd = RowReader(H.REAL)

    for i, t in enumerate(targets):
        Br = rd.rows(np.flatnonzero(rlab == t), G_full)
        s = Br.sum(1)
        s[s == 0] = 1.0
        Br = np.log1p(Br[:, keep] * (1e6 / s[:, None]))
        z = H.wilcoxon_z(cs, tie, np.ascontiguousarray(Br.T))
        real["p"][i] = 2 * norm.sf(np.abs(z))
        real["q"][i] = false_discovery_control(real["p"][i], method="bh")
        real["l"][i] = np.log2((np.expm1(Br).mean(0) + 1e-9) / (ref_lin + 1e-9))
        del Br

        Rec = recon196.Recon(sums[i], D["r"], D["gate"], n=n, seed=1000 + i)
        assert np.array_equal(Rec.totals(1.0), sums[i]), f"{t}: m=1 identity broken"
        # a15's spread rank axis, frozen at m = 1 so only the depth moves
        M1 = sums[i].astype(np.int64)
        S1 = int(M1.sum())
        pcpm = 1e6 * M1 / max(S1, 1)
        eff = np.abs(np.log2((pcpm + 1e-9) / (cpm + 1e-9)))
        order = np.lexsort((np.arange(G_full), -np.where(gated_scorer, eff, -np.inf)))
        sv = np.full(G_full, 1.0)
        sv[order[:min(8000, int(gated_scorer.sum()))]] = 0.25
        assert float(sv[~keep].max()) == 1.0
        keyH = Hd["H"][hpos[t]]
        absH_all[i] = np.abs(keyH)[keep]

        for m in DEPTHS:
            Lm = Rec.emit_L(m)
            prof = Rec.profile(m)
            want = Rec.totals(m)
            libs = np.full(n, Lm)
            for arm in ARMS:
                k = f"{arm}_m{m:g}"
                rng = np.random.default_rng(90_000 + i)
                if arm == "quasi":
                    vals, cols, ptr = emit_counts(prof, libs, DISPERSION, rng, mode="quasi",
                                                  spread=sv, overdispersion=PHI)
                    B = to_block(vals, cols, ptr, G_full)
                    dg = {}
                else:
                    # `rc_hybrid`'s greedy placement can fail to weld a row at depth -- that
                    # is itself a result, so it is recorded rather than allowed to abort the
                    # run.  The perturbation is then dropped from EVERY arm's cohort by
                    # `six196.py`, so the comparison stays paired.
                    try:
                        B, dg = RC.rc_hybrid_block(
                            prof, libs, DISPERSION, rng,
                            pass2_rng=np.random.default_rng(90_000 + i), spread=sv,
                            overdispersion=PHI, control=control, key=keyH,
                            jitter_seed=1000 + i, emit_counts=emit_counts)
                    except AssertionError as exc:
                        failed.setdefault(k, []).append((t, str(exc)))
                        digest[k].append("FAILED")
                        npred[k].append(-1)
                        diag[k].append(dict(nnz=0, max_cell=0, min_cell=0, cell_total_sd=0.0,
                                            emit_L=Lm, panel_S=int(want.sum()), failed=1))
                        print(f"    !! {k} {t}: {exc}", flush=True)
                        continue
                delivered = B.sum(0).astype(np.int64)
                assert np.array_equal(delivered, want), (k, t)
                digest[k].append(hashlib.sha256(delivered.tobytes()).hexdigest()[:16])
                rows = B.sum(1)
                dg.update(nnz=int((B > 0).sum()), max_cell=int(rows.max()),
                          min_cell=int(rows.min()), cell_total_sd=float(rows.std()),
                          emit_L=Lm, panel_S=int(delivered.sum()), failed=0)
                zz, p, q, lfc = H.de_from_block(B, keep, cs, tie, ref_lin)
                out[k]["z"][i], out[k]["p"][i] = zz, p
                out[k]["q"][i], out[k]["l"][i] = q, lfc
                npred[k].append(int((q < 0.05).sum()))
                diag[k].append(dg)
                del B
        if (i + 1) % 5 == 0 or i == 0:
            print(f"  {i+1}/{T_n} [{time.time()-t0:.0f}s]  " +
                  " ".join(f"{k}={npred[k][-1]}" for k in keys), flush=True)

    if failed:
        print("\nEMISSION FAILURES (rc_hybrid could not weld a row):")
        for k, v in failed.items():
            print(f"  {k:<16} {len(v)}/{T_n} targets: {[x[0] for x in v][:8]}")
    print("\nPSEUDOBULK: delivered totals equal floor(m * x) on every target, by assertion.")
    ref = digest[keys[0]]
    for k in keys:
        same = sum(a == b for a, b in zip(digest[k], ref))
        print(f"  {k:<16} digest identical to {keys[0]} on {same}/{T_n}")
    for k in keys:
        np.savez_compressed(WORK / f"{TAG}.{k}.npz", genes=genes_full[keep],
                            perts=np.array(targets), real_p=real["p"],
                            real_padj=real["q"], real_lfc=real["l"], pred_p=out[k]["p"],
                            pred_padj=out[k]["q"], pred_lfc=out[k]["l"], pred_z=out[k]["z"],
                            absH=absH_all, digest=np.array(digest[k]))
        d = diag[k]
        summ = {j: float(np.median([x[j] for x in d])) for j in d[0]}
        print(f"wrote {TAG}.{k}.npz  median n_pred {np.median(npred[k]):.0f}  {summ}",
              flush=True)
    json.dump(dict(npred=npred, diag=diag, depths=list(DEPTHS), arms=list(ARMS),
                   failed={k: [x[0] for x in v] for k, v in failed.items()}),
              open(WORK / f"{TAG}_diag.json", "w"), default=float)
    print(f"done [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
