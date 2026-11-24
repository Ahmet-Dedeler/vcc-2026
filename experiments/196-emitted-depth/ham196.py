"""196 — the ROUNDING RULE at fixed depth: `floor` against largest-remainder (Hamilton).

The depth sweep found that raising `L` buys `nmae` (the delivered profile stops being
systematically short) and pays for it on `mse` (the jackknife `Chat` falls out of `claim`'s
saturation and raw `mse` goes 0.882 -> 1.008).  This script asks whether the `nmae` half is
obtainable WITHOUT the depth, i.e. whether the defect is the depth or the ROUNDING.

`emit_counts` delivers `totals = floor(n L p)`.  Every gene is therefore short by its own
fractional part -- mean 0.5 counts -- and 0.5 counts is 0.9% of a 56-count gene and 0.02% of
a 2,500-count one.  The shortfall is a fixed number of COUNTS and the metrics read RATIOS, so
it is a size-correlated TILT, not a common mode: `work/round196.md` measures it at -0.035
log2 on h1-test's lowest gated decile against -0.0004 on its highest.

Largest-remainder allocation removes exactly that: floor everywhere, then hand the
`S - sum(floor)` leftover counts to the largest fractional remainders.  The per-gene error is
then centred on zero and uncorrelated with the gene's size.  It costs nothing -- the same `L`,
the same nnz to within one count per gene, the same emitted cell-to-cell spread, so the same
`Chat` and the same saturated `claim`.

ONE VARIABLE.  The spread vector, the seeds, `overdispersion`, `DISPERSION`, the reference
tables and the reconstruction are all the depth sweep's, bit for bit; `total_rule` is the only
difference.  Its `floor` arms must reproduce `work/de196.quasi_m1.npz` and
`work/rc196.rc_hybrid_m1.npz` exactly, which is the known-answer row (rule 10) -- and is why
this script re-emits them instead of reading them off disk.

    ARMS=quasi,rc_hybrid uv run python ham196.py
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
ARMS = tuple(os.environ.get("ARMS", "quasi,rc_hybrid").split(","))
RULES = tuple(os.environ.get("RULES", "floor,ham").split(","))
MAXT = int(os.environ.get("MAXT", "0"))
TAG = os.environ.get("OUTTAG", "ham196")


def hamilton(x, S=None):
    """Largest-remainder allocation of `S` counts over the intended per-gene means `x`.

    `S` defaults to `round(sum(x))` = `round(n L)`, the library the emitter is asked for.
    Ties on the remainder go to the lower gene index, so the rule is deterministic.
    """
    T = np.floor(x).astype(np.int64)
    S = int(round(float(x.sum()))) if S is None else int(S)
    k = S - int(T.sum())
    if k > 0:
        rem = x - T
        idx = np.lexsort((np.arange(x.size), -rem))[:k]
        T[idx] += 1
    elif k < 0:
        rem = x - np.floor(x)
        idx = np.lexsort((np.arange(x.size), rem))
        idx = idx[T[idx] > 0][:-k]
        T[idx] -= 1
    return T


def to_block(vals, cols, ptr, G):
    B = np.zeros((n, G), np.int64)
    for r in range(n):
        sl = slice(ptr[r], ptr[r + 1])
        B[r, cols[sl]] = vals[sl]
    return B


def tilt_stats(T, x, gate_scorer, cpm):
    """Delivered-vs-intended error over the SCORER-gated pool, and its size tilt."""
    g = gate_scorer & (x > 0)
    e = (T[g] - x[g]) / x[g]                       # signed relative error
    lo = cpm[g] <= np.quantile(cpm[g], 0.10)
    hi = cpm[g] >= np.quantile(cpm[g], 0.90)
    l2 = np.log2((T[g] + 1e-9) / x[g])
    return dict(mean_rel=float(e.mean()), mean_abs_rel=float(np.abs(e).mean()),
                common_log2=float(np.log2(T[g].sum() / x[g].sum())),
                d1_log2=float(l2[lo].mean()), d10_log2=float(l2[hi].mean()),
                tilt_log2=float(l2[lo].mean() - l2[hi].mean()),
                zero=int((T[g] == 0).sum()))


def main() -> None:
    t0 = time.time()
    C = H.load_controls()
    cs, tie, keep, cpm, ref_lin = C["cs"], C["tie"], C["keep"], C["cpm"], C["ref_lin"]
    genes_full = C["genes"]
    G_full, G = genes_full.size, cs.shape[0]
    gated_scorer = keep
    control = RC.RcControlBlock(keep=keep, cs=cs, tie=tie, ref_lin=ref_lin,
                                n_controls=cs.shape[1], min_cpm=H.MIN_CPM)
    assert cpm.size == G_full, (cpm.size, G_full)     # `load_controls` gives the full axis
    cpm_full = cpm

    D = recon196.load_h1()
    assert np.array_equal(D["genes"], genes_full)
    targets, sums = list(D["targets"]), D["sums"]
    Hd = H.load_H()
    hpos = {t: i for i, t in enumerate(Hd["targets"])}
    if MAXT:
        targets, sums = targets[:MAXT], sums[:MAXT]
    T_n = len(targets)
    keys = [f"{a}_{r}" for a in ARMS for r in RULES]
    print(f"h1-test: {cs.shape[1]} controls, {G}/{G_full} scorer-gated, {T_n} perturbations; "
          f"arms {keys} [{time.time()-t0:.0f}s]", flush=True)

    real = {k: np.full((T_n, G), np.nan) for k in ("p", "q", "l")}
    out = {k: {j: np.full((T_n, G), np.nan) for j in ("p", "q", "l", "z")} for k in keys}
    absH_all = np.zeros((T_n, G), np.float32)
    diag = {k: [] for k in keys}
    npred = {k: [] for k in keys}
    digest = {k: [] for k in keys}
    tilt = {r: [] for r in RULES}
    failed = {}

    import anndata as ad
    R = ad.read_h5ad(H.REAL, backed="r")
    rlab = np.asarray(R.obs[H.COL].astype(str).values, dtype="U64")
    R.file.close()
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
        Lp = Rec.emit_L(1.0)
        libs = np.full(n, Lp)
        # a15's spread rank axis, frozen on the SHIPPED (floor) totals so only the rule moves
        M1 = sums[i].astype(np.int64)
        pcpm = 1e6 * M1 / max(int(M1.sum()), 1)
        eff = np.abs(np.log2((pcpm + 1e-9) / (cpm_full + 1e-9)))
        order = np.lexsort((np.arange(G_full), -np.where(gated_scorer, eff, -np.inf)))
        sv = np.full(G_full, 1.0)
        sv[order[:min(8000, int(gated_scorer.sum()))]] = 0.25
        assert float(sv[~keep].max()) == 1.0
        keyH = Hd["H"][hpos[t]]
        absH_all[i] = np.abs(keyH)[keep]

        want = {"floor": Rec.totals(1.0), "ham": hamilton(Rec.x)}
        for r in RULES:
            tilt[r].append(tilt_stats(want[r], Rec.x, gated_scorer, cpm_full))
        for r in RULES:
            Tw = want[r]
            prof = (Tw + 0.5) / (n * Lp)
            for arm in ARMS:
                k = f"{arm}_{r}"
                rng = np.random.default_rng(90_000 + i)
                if arm == "quasi":
                    vals, cols, ptr = emit_counts(prof, libs, DISPERSION, rng, mode="quasi",
                                                  spread=sv, overdispersion=PHI)
                    B = to_block(vals, cols, ptr, G_full)
                    dg = {}
                else:
                    try:
                        B, dg = RC.rc_hybrid_block(
                            prof, libs, DISPERSION, rng,
                            pass2_rng=np.random.default_rng(90_000 + i), spread=sv,
                            overdispersion=PHI, control=control, key=keyH,
                            jitter_seed=1000 + i, emit_counts=emit_counts)
                    except AssertionError as exc:
                        failed.setdefault(k, []).append(t)
                        digest[k].append("FAILED")
                        npred[k].append(-1)
                        diag[k].append(dict(nnz=0, max_cell=0, cell_total_sd=0.0, failed=1))
                        print(f"    !! {k} {t}: {exc}", flush=True)
                        continue
                delivered = B.sum(0).astype(np.int64)
                assert np.array_equal(delivered, Tw), (k, t)
                digest[k].append(hashlib.sha256(delivered.tobytes()).hexdigest()[:16])
                rows = B.sum(1)
                dg.update(nnz=int((B > 0).sum()), max_cell=int(rows.max()),
                          cell_total_sd=float(rows.std()), panel_S=int(delivered.sum()),
                          failed=0)
                zz, p, q, lfc = H.de_from_block(B, keep, cs, tie, ref_lin)
                out[k]["z"][i], out[k]["p"][i] = zz, p
                out[k]["q"][i], out[k]["l"][i] = q, lfc
                npred[k].append(int((q < 0.05).sum()))
                diag[k].append(dg)
                del B
        if (i + 1) % 5 == 0 or i == 0:
            print(f"  {i+1}/{T_n} [{time.time()-t0:.0f}s]  " +
                  " ".join(f"{k}={npred[k][-1]}" for k in keys), flush=True)

    print("\nDELIVERED-vs-INTENDED over the scorer-gated pool (median over perturbations)")
    print(f"{'rule':<8}{'mean rel':>11}{'mean |rel|':>12}{'common log2':>13}"
          f"{'dec1 log2':>11}{'dec10 log2':>12}{'TILT log2':>11}{'zeroed':>8}")
    tsum = {}
    for r in RULES:
        med = {j: float(np.median([x[j] for x in tilt[r]])) for j in tilt[r][0]}
        tsum[r] = med
        print(f"{r:<8}{med['mean_rel']:>+11.5f}{med['mean_abs_rel']:>12.5f}"
              f"{med['common_log2']:>+13.6f}{med['d1_log2']:>+11.5f}{med['d10_log2']:>+12.5f}"
              f"{med['tilt_log2']:>+11.5f}{med['zero']:>8.0f}")

    if failed:
        print("\nEMISSION FAILURES:")
        for k, v in failed.items():
            print(f"  {k:<18} {len(v)}/{T_n}: {v[:8]}")

    print("\nKNOWN-ANSWER ROW (rule 10): the `floor` arms must be BIT-IDENTICAL to the depth "
          "sweep's own m=1 arms.")
    ok = True
    for k, ref in (("quasi_floor", WORK / "de196.quasi_m1.npz"),
                   ("rc_hybrid_floor", WORK / "rc196.rc_hybrid_m1.npz")):
        if k not in keys or not ref.exists():
            continue
        z = np.load(ref, allow_pickle=True)
        same_d = bool(np.array_equal(np.asarray(z["digest"]).astype(str)[:T_n],
                                     np.array(digest[k], dtype=str)))
        same_p = bool(np.array_equal(np.nan_to_num(z["pred_p"][:T_n], nan=-1.0),
                                     np.nan_to_num(out[k]["p"], nan=-1.0)))
        same_r = bool(np.array_equal(np.nan_to_num(z["real_p"][:T_n], nan=-1.0),
                                     np.nan_to_num(real["p"], nan=-1.0)))
        print(f"  {k:<18} totals digest {'MATCH' if same_d else 'DIFFER'}   pred p "
              f"{'MATCH' if same_p else 'DIFFER'}   reference p "
              f"{'MATCH' if same_r else 'DIFFER'}")
        ok &= same_d and same_p and same_r
    if not ok and not os.environ.get("NOABORT"):
        raise SystemExit("ABORT (rule 10): the floor arms did not reproduce the depth sweep")
    print("  known-answer rows PASS" if ok else "")

    for k in keys:
        np.savez_compressed(WORK / f"{TAG}.{k}.npz", genes=genes_full[keep],
                            perts=np.array(targets), real_p=real["p"],
                            real_padj=real["q"], real_lfc=real["l"], pred_p=out[k]["p"],
                            pred_padj=out[k]["q"], pred_lfc=out[k]["l"], pred_z=out[k]["z"],
                            absH=absH_all, digest=np.array(digest[k]))
        d = [x for x in diag[k] if not x.get("failed")]
        summ = {j: float(np.median([x[j] for x in d])) for j in ("nnz", "max_cell",
                                                                 "cell_total_sd", "panel_S")}
        print(f"wrote {TAG}.{k}.npz  median n_pred {np.median(npred[k]):.0f}  {summ}",
              flush=True)
    json.dump(dict(npred=npred, diag=diag, arms=list(ARMS), rules=list(RULES),
                   tilt=tsum, failed=failed), open(WORK / f"{TAG}_diag.json", "w"),
              default=float)
    print(f"done [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
