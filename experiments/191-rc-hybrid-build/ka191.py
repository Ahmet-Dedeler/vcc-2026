"""191 gate 0 — the KNOWN-ANSWER row: drive the PRODUCTION emitter over h1-test and
reproduce 185's `rc_hybrid` row before it is trusted on the challenge.

    a15        reach 0.18283   n_pred 6034 (debench convention)   nmae 0.96659
    rc_hybrid  reach 0.21681   n_pred 6043                        nmae 0.96670
               jac 0.06827     P(1) 0.8780                        twin spread ~7e-5

Everything on the h1-test side -- the control block, a15's per-gene totals, the `|H|` key,
the seeds, the spread policy -- is 185's, read from 185's own cached work (read-only).  The
EMITTER is `vcc2026.rc_emit`, i.e. the object the challenge build uses.  If the row does not
reproduce, the port is a different emitter and the preregistration does not apply to it.

Two seedings are run:
  `rc_hybrid`      185's exact bookkeeping (pass 1 and pass 2 each open default_rng(90000+i))
  `rc_hybrid_seq`  production's (one stream, pass 2 continues where pass 1 stopped)
The second is the arm the challenge build emits; the pair measures the seeding difference
against 185's own twin spread.

Donors: 095's h1-test artifact (k562_gwps + orion_hek293t).  `h1_2025` never enters.
No cell-eval2, no build, no submission.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time

import numpy as np
from scipy.stats import false_discovery_control, norm

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/185-rank-controlled-emission")
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")

import harness185 as H  # noqa: E402
from vcc2026 import rc_emit as RC  # noqa: E402
from vcc2026.donor_transfer import emit_counts  # noqa: E402

n, EMIT_L, PHI = 400, H.EMIT_L, 4.0
OUT = os.path.join(WORK, "de191.%s.npz")
MAXT = int(os.environ.get("MAXT", "0"))


class RowReader:
    """Dense rows out of a CSR h5ad without loading it (realize185.py's)."""

    def __init__(self, path):
        import h5py

        self.f = h5py.File(path, "r")
        self.indptr = self.f["X"]["indptr"][:]
        self.data = self.f["X"]["data"]
        self.indices = self.f["X"]["indices"]

    def rows(self, rows, G):
        rows = np.sort(np.asarray(rows))
        out = np.zeros((rows.size, G), np.float64)
        for k, r in enumerate(rows):
            a, b = int(self.indptr[r]), int(self.indptr[r + 1])
            out[k, self.indices[a:b]] = self.data[a:b]
        return out


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
    b_cell = np.log1p(cpm)
    gated = (cpm > H.MIN_CPM) & (b_cell > 0)
    assert gated.sum() == G and np.array_equal(np.flatnonzero(gated), np.flatnonzero(keep))
    # The production emitter's own view of the same control block.  Built from 185's cache
    # rather than re-read, so this script cannot disagree with 185 about the surface.
    control = RC.RcControlBlock(
        keep=keep, cs=cs, tie=tie, ref_lin=ref_lin, n_controls=cs.shape[1],
        min_cpm=H.MIN_CPM,
    )

    TOT = H.load_totals()
    assert np.array_equal(TOT["genes"], genes_full)
    tmask = TOT["targets"] != "non-targeting"
    targets = list(TOT["targets"][tmask])
    totals_all = TOT["totals"][tmask]
    if MAXT:
        targets, totals_all = targets[:MAXT], totals_all[:MAXT]
    Hd = H.load_H()
    assert np.array_equal(Hd["genes"], genes_full)
    hpos = {t: i for i, t in enumerate(Hd["targets"])}
    T_n = len(targets)

    import anndata as ad

    R = ad.read_h5ad(H.REAL, backed="r")
    rlab = np.asarray(R.obs[H.COL].astype(str).values, dtype="U64")
    assert np.array_equal(np.asarray([str(g) for g in R.var_names]), genes_full)
    R.file.close()
    rd = RowReader(H.REAL)

    ARMS = ["a15", "rc_hybrid", "rc_hybrid_seq"]
    real = {k: np.full((T_n, G), np.nan) for k in ("p", "q", "l")}
    out = {a: {k: np.full((T_n, G), np.nan) for k in ("p", "q", "l", "z")} for a in ARMS}
    npred = {a: [] for a in ARMS}
    digest = {a: [] for a in ARMS}
    diag = {a: [] for a in ARMS}
    absH_all = np.zeros((T_n, G), np.float32)
    print(f"h1-test: {cs.shape[1]} controls, {G}/{G_full} gated, {T_n} perturbations "
          f"[{time.time() - t0:.0f}s]", flush=True)

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

        Mfull = totals_all[i].astype(np.int64)
        S = int(Mfull.sum())
        profile = (Mfull + 0.5) / (n * EMIT_L)      # floor(n L profile) == Mfull exactly
        pcpm = 1e6 * Mfull / max(S, 1)
        eff = np.abs(np.log2((pcpm + 1e-9) / (cpm + 1e-9)))
        score = np.where(gated, eff, -np.inf)
        order = np.lexsort((np.arange(G_full), -score))
        sv = np.full(G_full, 1.0)                   # a15: tail 1.0, keep block 0.25
        sv[order[:min(8000, int(gated.sum()))]] = 0.25
        assert float(sv[~keep].max()) == 1.0, "h1-test sub-gate spread is not a15's tail law"
        key = Hd["H"][hpos[t]]
        absH_all[i] = np.abs(key)[keep]
        libs = np.full(n, EMIT_L)

        for name in ARMS:
            if name == "a15":
                rng = np.random.default_rng(90_000 + i)
                vals, cols, ptr = emit_counts(profile, libs, 0.1, rng, mode="quasi",
                                              spread=sv, overdispersion=PHI)
                B = to_block(vals, cols, ptr, G_full)
                dg = {}
            else:
                if name == "rc_hybrid":             # 185's bookkeeping: two fresh streams
                    rng = np.random.default_rng(90_000 + i)
                    rng2 = np.random.default_rng(90_000 + i)
                else:                               # production's: one stream
                    rng = rng2 = np.random.default_rng(90_000 + i)
                B, dg = RC.rc_hybrid_block(
                    profile, libs, 0.1, rng, pass2_rng=rng2, spread=sv,
                    overdispersion=PHI, control=control, key=key, jitter_seed=1000 + i,
                    emit_counts=emit_counts,
                )
            delivered = B.sum(0).astype(np.int64)
            assert np.array_equal(delivered, Mfull), name
            digest[name].append(hashlib.sha256(delivered.tobytes()).hexdigest()[:16])
            zz, p, q, lfc = H.de_from_block(B, keep, cs, tie, ref_lin)
            out[name]["z"][i], out[name]["p"][i] = zz, p
            out[name]["q"][i], out[name]["l"][i] = q, lfc
            npred[name].append(int((q < 0.05).sum()))
            dg["cell_total_sd"] = float(B.sum(1).std())
            diag[name].append(dg)
            del B
        if (i + 1) % 10 == 0 or i == 0:
            print(f"  {i + 1}/{T_n} [{time.time() - t0:.0f}s]  n_pred " +
                  " ".join(f"{a}={npred[a][-1]}" for a in ARMS), flush=True)

    print("\nPSEUDOBULK BIT-IDENTITY (sha256 of the delivered per-gene panel total)")
    ref = digest["a15"]
    for a in ARMS:
        same = sum(x == y for x, y in zip(digest[a], ref))
        print(f"  {a:<16} identical to a15 on {same}/{T_n}   first digest {digest[a][0]}")

    for a in ARMS:
        np.savez_compressed(
            OUT % a, genes=genes_full[keep], perts=np.array(targets),
            real_p=real["p"], real_padj=real["q"], real_lfc=real["l"],
            pred_p=out[a]["p"], pred_padj=out[a]["q"], pred_lfc=out[a]["l"],
            pred_z=out[a]["z"], absH=absH_all, digest=np.array(digest[a]),
        )
        d = diag[a]
        summ = {k: float(np.median([x[k] for x in d])) for k in d[0]} if d and d[0] else {}
        print(f"wrote {OUT % a}   median n_pred {np.median(npred[a]):.0f}   diag {summ}",
              flush=True)
    with open(os.path.join(WORK, "ka191_diag.json"), "w") as fh:
        json.dump(dict(npred=npred, digest=digest, diag=diag), fh)
    print(f"done [{time.time() - t0:.0f}s]")


if __name__ == "__main__":
    main()
