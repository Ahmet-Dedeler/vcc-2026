"""241 — the PerturbMatch row on the EMITTER we actually ship (`rc_callset`), not `quasi`.

200 (and therefore this experiment's main table) emits with the `quasi` law.  The shipped
package emits with `rc_callset` (202; board 0.1918638, `jac` +0.0276 / `reach` +0.0226 scaled
over `rc_hybrid`).  214 §1.1 established that `rc_hybrid`/`rc_callset`/a23 all share ONE
effect matrix, so the emitter cannot move any component's `pds` or `mse` — but it re-draws the
CALL SET, so it can and does move `fid`, `jac` and `reach`.  PerturbMatch's verdict lives
almost entirely in those three plus `nmae`, so it is the row worth re-emitting.

FALSIFIER (rule 10, aborts).  `rc_callset` is pseudobulk-preserving by construction, and that
is the identity every "emission arms leave `pds`/`mse` pinned" claim rests on.  Each block
asserts the delivered per-gene panel totals equal `floor(n L p)` EXACTLY — the same totals
`de241.reb_*_nomf.npz` delivered — so the `quasi` -> `rc_callset` comparison is single-variable
by construction rather than by argument.

    ARMS=rc_ship,rc_noPM SEEDS=3 uv run python emitrc241.py
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
E200 = os.path.join(ROOT, "experiments/200-production-ablation")
E169 = os.path.join(ROOT, "experiments/169-sharpening-six-member")
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/185-rank-controlled-emission")
sys.path.insert(0, E200)
import harness185 as H                                          # noqa: E402
import vcc2026.rc_emit as RC                                    # noqa: E402
from vcc2026.donor_transfer import emit_counts                  # noqa: E402
import h1parts200                                               # noqa: E402

os.makedirs(WORK, exist_ok=True)
OUT = os.path.join(HERE, "de241")
N_EMIT, EMIT_L = 400, 100_000
PHI, DISPERSION = 4.0, 0.1
KEEP_TOP, KEEP_SPREAD, TAIL = 8000, 0.25, 1.0
N_HEAD, N_CPM, Z_CAP = 400, 2200, 1.90       # 202's shipped call set
SEEDS = (90_000, 7_000_000, 11_111)
G = dict(
    ship=dict(pm=0.50, gamma=0.65, lam=0.0, cd4_f=0.10, cd4_frac=None, weights=None),
    noPM=dict(pm=0.00, gamma=0.65, lam=0.0, cd4_f=0.10, cd4_frac=None, weights=None),
)


def main() -> None:
    t0 = time.time()
    nseed = int(os.environ.get("SEEDS", "3"))
    geoms = os.environ.get("GEOMS", "ship,noPM").split(",")
    ARMS = []
    for g in geoms:
        for i, s in enumerate(SEEDS[:nseed]):
            tag = "" if i == 0 else f"_s{i+1}"
            ARMS.append((f"rc_{g}_nomf{tag}", g, s))
    ARMS = [a for a in ARMS if not os.path.exists(f"{OUT}.{a[0]}.npz")]
    if not ARMS:
        print("every requested arm already on disk")
        return

    C = H.load_controls()
    cs, tie, keep, cpm, ref_lin = C["cs"], C["tie"], C["keep"], C["cpm"], C["ref_lin"]
    genes_full = np.asarray([str(g) for g in C["genes"]])
    G_full, Gk = genes_full.size, int(keep.sum())
    control = RC.RcControlBlock(keep=keep, cs=cs, tie=tie, ref_lin=ref_lin,
                                n_controls=cs.shape[1], min_cpm=H.MIN_CPM)
    call_set = RC.CallSet(N_HEAD, N_CPM, Z_CAP)
    print(f"h1-test: {cs.shape[1]} controls, {Gk}/{G_full} scorer-gated; call set "
          f"{call_set.size}, z_cap {Z_CAP}; arms {[a[0] for a in ARMS]} "
          f"[{time.time()-t0:.0f}s]", flush=True)

    E = np.load(os.path.join(E169, "work/h1test-effect.npz"), allow_pickle=True)
    H0 = E["H"].astype(np.float64)
    b = E["b"].astype(np.float64)
    a_gate = np.asarray(E["gate"])
    etargets = np.asarray([str(t) for t in E["targets"]])
    esums = E["sums"].astype(np.float64)
    assert np.array_equal(np.asarray([str(g) for g in E["genes"]]), genes_full)
    art_norm = float(np.median(np.linalg.norm(H0, axis=1)))

    p = h1parts200.parts()
    c4 = h1parts200.cd4(p["targets"])
    op = h1parts200.H1Operator()
    order = np.asarray([list(p["targets"]).index(t) for t in etargets], dtype=np.int64)
    GEO = {}
    for g in sorted({a[1] for a in ARMS}):
        GEO[g] = h1parts200.effects_for(G[g], p, c4, op, art_norm)[order]
        print(f"  geometry {g:<8} median norm "
              f"{float(np.median(np.linalg.norm(GEO[g], axis=1))):.4f}", flush=True)
    del op

    # the reference side is identical for every arm -- take it from the quasi tables
    ref = np.load(os.path.join(HERE, "de241.reb_ship_nomf.npz"), allow_pickle=True)
    targets = [str(t) for t in ref["perts"]]
    assert np.array_equal(np.asarray([str(g) for g in ref["genes"]]), genes_full[keep])
    epos = {t: i for i, t in enumerate(etargets)}
    T = len(targets)
    b_cell = np.log1p(cpm)

    out = {n: {k: np.full((T, Gk), np.nan) for k in ("p", "q", "l")} for n, *_ in ARMS}
    npred = {n: [] for n, *_ in ARMS}
    digest = {n: [] for n, *_ in ARMS}
    quasi_digest = []

    for i, t in enumerate(targets):
        k = epos[t]
        S = float(esums[k].sum())
        for name, geom, seed in ARMS:
            h = GEO[geom][k]
            g = a_gate
            hq = np.where(g, h, 0.0)
            lin = np.clip(np.expm1(b + hq), 0.0, None)
            profile = lin / lin.sum() * (S / (N_EMIT * EMIT_L))
            totals = np.floor(N_EMIT * EMIT_L * profile).astype(np.int64)
            pcpm = 1e6 * totals / max(totals.sum(), 1)
            eff = np.abs(np.log2((pcpm + 1e-9) / (cpm + 1e-9)))
            gated = (cpm > H.MIN_CPM) & (b_cell > 0)
            score = np.where(gated, eff, -np.inf)
            sv = np.full(G_full, TAIL)
            o = np.lexsort((np.arange(G_full), -score))
            sv[o[:min(KEEP_TOP, int(gated.sum()))]] = KEEP_SPREAD

            rng = np.random.default_rng(seed + i)
            B, _dg = RC.rc_hybrid_block(
                profile, np.full(N_EMIT, EMIT_L), DISPERSION, rng, pass2_rng=rng,
                spread=sv, overdispersion=PHI, control=control, key=hq,
                jitter_seed=1000 + i, emit_counts=emit_counts, call_set=call_set)
            delivered = B.sum(0).astype(np.int64)
            # FALSIFIER: pseudobulk-preserving, asserted per block
            assert np.array_equal(delivered, totals), (name, t)
            digest[name].append(hashlib.sha256(delivered.tobytes()).hexdigest()[:16])
            _z, pv, qv, lfc = H.de_from_block(B, keep, cs, tie, ref_lin)
            out[name]["p"][i], out[name]["q"][i], out[name]["l"][i] = pv, qv, lfc
            npred[name].append(int((qv < 0.05).sum()))
            del B
        if (i + 1) % 10 == 0:
            print(f"  {i+1}/{T} [{time.time()-t0:.0f}s]  "
                  + " ".join(f"{n}={npred[n][-1]}" for n, *_ in ARMS[:2]), flush=True)

    meta = {}
    for name, *_ in ARMS:
        np.savez_compressed(
            f"{OUT}.{name}.npz", genes=genes_full[keep], perts=np.array(targets),
            real_p=ref["real_p"], real_padj=ref["real_padj"], real_lfc=ref["real_lfc"],
            pred_p=out[name]["p"], pred_padj=out[name]["q"], pred_lfc=out[name]["l"],
            spread=np.zeros((T, Gk), dtype=np.float32))
        meta[name] = dict(median_npred=float(np.median(npred[name])),
                          totals_digest=hashlib.sha256(
                              "".join(digest[name]).encode()).hexdigest()[:16])
        print(f"wrote de241.{name}.npz   median n_pred {meta[name]['median_npred']:.0f}  "
              f"totals digest {meta[name]['totals_digest']}", flush=True)
    json.dump(meta, open(os.path.join(WORK, "emitrc241.json"), "w"), indent=1)
    print(f"done [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
