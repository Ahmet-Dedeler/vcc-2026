"""202 gate 0b — drive the PATCHED emitter over h1-test, EMITTED, not prescribed.

192 priced the arm by *prescribing* a call set on 185's frozen DE tables; `ka202_prescribed`
reproduced that grid row for row.  This script asks the different question the build has to
answer: does the emitter actually DELIVER the prescribed set, in the real regime, at the
identities the arm rests on?

Two arms, both on production's seeding (one stream, which is what the challenge writer does):

  `rc_hybrid`   MUST reproduce 191's `work/de191.rc_hybrid_seq.npz` EXACTLY -- same seeds,
                same control block, same totals.  This is the known-answer row that says the
                202 patch did not move the arm the board already scored at +0.0093325.
  `rc_callset`  the arm.

Everything on the h1-test side is 185's cache, read-only, exactly as `ka191.py` read it, and
the reference DE table is 191's own (`de191.a15.npz`), so this script cannot disagree with
191 about the surface.  The emitter is loaded from `work/rc_emit_candidate.py` by path, so
`src/vcc2026/rc_emit.py` is untouched while another experiment's build is running.

No cell-eval2, no build, no submission.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import sys
import time

import numpy as np
from scipy.stats import false_discovery_control, norm  # noqa: F401

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/185-rank-controlled-emission")

import harness185 as H  # noqa: E402
from vcc2026.donor_transfer import emit_counts  # noqa: E402

CAND = os.environ.get("RC_EMIT_CANDIDATE", os.path.join(WORK, "rc_emit_candidate.py"))
REF191 = f"{ROOT}/experiments/191-rc-hybrid-build/work"
n, EMIT_L, PHI = 400, H.EMIT_L, 4.0
N_HEAD, N_CPM, Z_CAP = 400, 2200, 1.90
MAXT = int(os.environ.get("MAXT", "0"))


def load_candidate():
    spec = importlib.util.spec_from_file_location("rc_emit_candidate", CAND)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class RowReader:
    """Dense rows out of a CSR h5ad without loading it (191's, which is 185's)."""

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


def main() -> None:
    RC = load_candidate()
    t0 = time.time()
    C = H.load_controls()
    cs, tie, keep, cpm, ref_lin = C["cs"], C["tie"], C["keep"], C["cpm"], C["ref_lin"]
    genes_full = C["genes"]
    G_full, G = genes_full.size, cs.shape[0]
    b_cell = np.log1p(cpm)
    gated = (cpm > H.MIN_CPM) & (b_cell > 0)
    assert gated.sum() == G and np.array_equal(np.flatnonzero(gated), np.flatnonzero(keep))
    control = RC.RcControlBlock(keep=keep, cs=cs, tie=tie, ref_lin=ref_lin,
                                n_controls=cs.shape[1], min_cpm=H.MIN_CPM)
    call_set = RC.CallSet(N_HEAD, N_CPM, Z_CAP)

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

    # 191's own tables: the reference side, and the row this run must reproduce.
    ref = np.load(os.path.join(REF191, "de191.a15.npz"), allow_pickle=True)
    seq = np.load(os.path.join(REF191, "de191.rc_hybrid_seq.npz"), allow_pickle=True)
    assert list(map(str, ref["perts"]))[:T_n] == targets
    assert np.array_equal(ref["genes"], genes_full[keep])

    ARMS = ["rc_hybrid", "rc_callset"]
    out = {a: {k: np.full((T_n, G), np.nan) for k in ("p", "q", "l", "z")} for a in ARMS}
    npred = {a: [] for a in ARMS}
    digest = {a: [] for a in ARMS}
    diag = {a: [] for a in ARMS}
    extra = []
    absH_all = np.zeros((T_n, G), np.float32)
    print(f"h1-test: {cs.shape[1]} controls, {G}/{G_full} gated, {T_n} perturbations; "
          f"call set = top {N_HEAD} |H| U top {N_CPM} cpm = {call_set.size}, "
          f"z_cap {Z_CAP} [{time.time() - t0:.0f}s]", flush=True)

    for i, t in enumerate(targets):
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
            rng = np.random.default_rng(90_000 + i)   # production's single stream
            B, dg = RC.rc_hybrid_block(
                profile, libs, 0.1, rng, pass2_rng=rng, spread=sv, overdispersion=PHI,
                control=control, key=key, jitter_seed=1000 + i, emit_counts=emit_counts,
                call_set=(call_set if name == "rc_callset" else None),
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

        # the composition of the set that was actually delivered, against the one prescribed
        key_abs = np.abs(np.asarray(key, dtype=np.float64)[keep])
        head = np.zeros(G, bool)
        head[np.argsort(-key_abs, kind="stable")[:N_HEAD]] = True
        want = head.copy()
        want[np.argsort(-np.where(head, -np.inf, ref_lin), kind="stable")[:N_CPM]] = True
        raw_cpm = np.zeros(G, bool)
        raw_cpm[np.argsort(-ref_lin, kind="stable")[:N_CPM]] = True
        got = out["rc_callset"]["q"][i] < 0.05
        inc = out["rc_hybrid"]["q"][i] < 0.05
        extra.append(dict(
            target=t,
            prescribed=int(want.sum()), delivered=int(got.sum()),
            jaccard=float((got & want).sum() / max((got | want).sum(), 1)),
            extra_called=int((got & ~want).sum()), missed=int((~got & want).sum()),
            halves_shared=int((head & raw_cpm).sum()),
            overlap_rc_hybrid=int((want & inc).sum()),
            jac_vs_rc_hybrid=float((want & inc).sum() / max((want | inc).sum(), 1)),
            n_pred_rc_hybrid=int(inc.sum()),
        ))
        if (i + 1) % 10 == 0 or i == 0:
            e = extra[-1]
            print(f"  {i + 1}/{T_n} [{time.time() - t0:.0f}s]  n_pred " +
                  " ".join(f"{a}={npred[a][-1]}" for a in ARMS) +
                  f"   delivered/prescribed {e['delivered']}/{e['prescribed']} "
                  f"jaccard {e['jaccard']:.4f}", flush=True)

    print("\nKNOWN-ANSWER ROW 0b — `rc_hybrid` through the PATCHED emitter vs 191's table")
    fields = [("p", "pred_p"), ("q", "pred_padj"), ("l", "pred_lfc"), ("z", "pred_z")]
    ok = True
    for mine, theirs in fields:
        same = np.array_equal(out["rc_hybrid"][mine], seq[theirs][:T_n])
        ok &= same
        worst = float(np.nanmax(np.abs(out["rc_hybrid"][mine] - seq[theirs][:T_n])))
        print(f"  {theirs:<10} bit-identical on all {T_n} x {G}: "
              f"{'YES' if same else 'NO'}   max |diff| {worst:.3g}")
    same_digest = sum(x == y for x, y in
                      zip(digest["rc_hybrid"], list(map(str, seq["digest"][:T_n]))))
    print(f"  per-gene panel-total digest identical to 191's on {same_digest}/{T_n}")
    print(f"  KNOWN-ANSWER ROW 0b: {'PASS' if ok and same_digest == T_n else 'FAIL'}")

    print("\nPSEUDOBULK BIT-IDENTITY (sha256 of the delivered per-gene panel total vs a15)")
    a15 = list(map(str, ref["digest"][:T_n]))
    for a in ARMS:
        same = sum(x == y for x, y in zip(digest[a], a15))
        print(f"  {a:<12} identical to a15 on {same}/{T_n}   first digest {digest[a][0]}")

    d = np.array([e["delivered"] for e in extra])
    j = np.array([e["jaccard"] for e in extra])
    u = np.array([e["extra_called"] for e in extra])
    print(f"\nDELIVERY: median n_pred {np.median(d):.0f} against a prescribed "
          f"{call_set.size} ({100 * (np.median(d) / call_set.size - 1):+.2f}%), "
          f"mean {d.mean():.0f} ({100 * (d.mean() / call_set.size - 1):+.2f}%), "
          f"range [{d.min()}, {d.max()}]")
    print(f"  jaccard(delivered, prescribed) median {np.median(j):.4f} "
          f"mean {j.mean():.4f} min {j.min():.4f}")
    print(f"  UNSUPPRESSIBLE (the |z| floor the delivered profile sets): median "
          f"{np.median(u):.0f}  mean {u.mean():.0f}  p90 {np.quantile(u, 0.9):.0f}  "
          f"max {u.max()}  ({int((u > 0.25 * call_set.size).sum())} of {len(u)} "
          f"perturbations above 25% of the prescribed size)")
    print(f"  missed (prescribed but not called) median "
          f"{np.median([e['missed'] for e in extra]):.0f}  "
          f"max {max(e['missed'] for e in extra)}")
    print(f"  halves shared (|H| head n top-{N_CPM} cpm) median "
          f"{np.median([e['halves_shared'] for e in extra]):.0f}")
    print(f"  overlap with `rc_hybrid`'s own call set: median "
          f"{np.median([e['overlap_rc_hybrid'] for e in extra]):.0f} genes, jaccard "
          f"{np.median([e['jac_vs_rc_hybrid'] for e in extra]):.4f} "
          f"(rc_hybrid n_pred median {np.median([e['n_pred_rc_hybrid'] for e in extra]):.0f})")

    for a in ARMS:
        path = os.path.join(WORK, f"de202.{a}.npz")
        np.savez_compressed(
            path, genes=genes_full[keep], perts=np.array(targets),
            real_p=ref["real_p"][:T_n], real_padj=ref["real_padj"][:T_n],
            real_lfc=ref["real_lfc"][:T_n],
            pred_p=out[a]["p"], pred_padj=out[a]["q"], pred_lfc=out[a]["l"],
            pred_z=out[a]["z"], absH=absH_all, digest=np.array(digest[a]),
        )
        dd = diag[a]
        summ = {k: float(np.median([x[k] for x in dd])) for k in dd[0]}
        print(f"\nwrote {path}   median n_pred {np.median(npred[a]):.0f}")
        print(f"  diag {json.dumps(summ, sort_keys=True)}")
    with open(os.path.join(WORK, "ka202_emit.json"), "w") as fh:
        json.dump(dict(npred=npred, digest=digest, diag=diag, composition=extra,
                       known_answer_0b=bool(ok and same_digest == T_n)), fh, default=float)
    print(f"\ndone [{time.time() - t0:.0f}s]")
    if not (ok and same_digest == T_n):
        raise SystemExit("known-answer row 0b FAILED — the patch moved `rc_hybrid`")


if __name__ == "__main__":
    main()
