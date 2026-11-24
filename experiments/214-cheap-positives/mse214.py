"""214 — `mse` for arms A (drop PerturbMatch) and B (Hamilton), marginal to `cd4f034`.

TWO CHANNELS, PRICED SEPARATELY BECAUSE THEY ARE DIFFERENT ARM CLASSES (rule 13, 200's
mse200 docstring):

  A  is a GEOMETRY arm.  It moves `A`, `B`, `claim` and `X` in 096's identity

        mse - 1 = ( s^2 A - 2 s X + Kc - min(s^2 B, claim) ) / sumD

     `X` needs the withheld truth, so 200 held it at the null (rX = 1) and transferred
     `A`/`B`/`claim` as RATIOS against its own SHIP row.  That is reproduced here, and its
     PerturbMatch row is a known-answer gate.  Because `rX = 1` deliberately refuses to
     credit any direction change, the claim route is a PARTIAL price for a geometry arm; the
     cosine channel is priced separately from PLOLO through 197's identity in `price214.py`.

  B  is an EMISSION-LAW arm that changes the delivered per-gene totals.  It is priced the
     way 208 priced it: `claim` on the challenge panel with the INCUMBENT's emitter
     (`rc_callset`) over `cd4f034`'s OWN delivered profile, plus the rounding-energy term
     converted through the PLOLO cross ratio measured at the incumbent's configuration.

KNOWN-ANSWER GATES (each aborts): 164 §1's four `claim` levels, the A/B round trip at s = 1,
a15 raw `mse` 0.88228 / score 0.11112, and 200's own PerturbMatch geometry row.

    MAXT=8 uv run python mse214.py
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
P200 = ROOT / "experiments/200-production-ablation"
P196 = ROOT / "experiments/196-emitted-depth"
P162 = ROOT / "experiments/162-retrodicting-validator"
for p in (str(P162), str(P200), str(P196), str(ROOT / "src"),
          str(ROOT / "experiments/096-mse-floor"), str(HERE)):
    sys.path.insert(0, p)

import anndata as ad                                                     # noqa: E402
import scipy.sparse as sp                                                # noqa: E402
import mech_alpha_chal as M                                              # noqa: E402
import mse200                                                            # noqa: E402
import recon196                                                          # noqa: E402
import vcc2026.rc_emit as RC                                             # noqa: E402
from ham196 import hamilton                                              # noqa: E402
from mse196 import DISPERSION, KNEE, PHI, S15, TS, n, raw_mse, score, sD  # noqa: E402
from mse196 import K as CAP_K                                            # noqa: E402
from vcc2026.donor_transfer import emit_counts                           # noqa: E402
from pds214 import BASE, arm, deliver_rounded, effects_for               # noqa: E402

MAXT = int(os.environ.get("MAXT", "8"))
N_HEAD, N_CPM, Z_CAP = 400, 2200, 1.90
KEEP_TOP, KEEP_SPREAD = 8000, 0.25
TAIL = 1.0
STATS = str(ROOT / "experiments/096-mse-floor/work/v11-stats.npz")

# The four geometry configurations.  `inc`/`A` are the live pair; `mf_*` reproduce 200's own
# PerturbMatch row, which was measured on a23's renderer at f = 0.10.
GEOM = {
    "inc_f034_pm50": arm(cd4_f=0.34, pm=0.50, mass=None),
    "A_f034_pm00": arm(cd4_f=0.34, pm=0.00, mass=None),
    "ka200_SHIP": arm(cd4_f=0.10, pm=0.50, mass="minnorm"),
    "ka200_noPM": arm(cd4_f=0.10, pm=0.00, mass="minnorm"),
    "f010_pm50": arm(cd4_f=0.10, pm=0.50, mass=None),
    "f010_pm00": arm(cd4_f=0.10, pm=0.00, mass=None),
}


def cd4f034_chal(ctx):
    """`recon196.load_chal`, with the INCUMBENT's own delivered totals swapped in."""
    D = recon196.load_chal(ctx)
    z = np.load(WORK / "totals-cd4f034.npz", allow_pickle=True)
    blocks = np.asarray(z["blocks"]).astype(str)
    sel = np.flatnonzero(np.array([b.split("|")[0] for b in blocks]) == ctx)
    assert np.array_equal(np.asarray(z["genes"]).astype(str), D["genes"])
    assert [blocks[i].split("|")[1] for i in sel] == list(D["targets"])
    D["sums"] = z["totals"][sel]
    return D


def plolo_cross(cfg):
    """The cross term at a given configuration: how much of the rounding error's own energy
    lands in squared error.  196 §4b measured 1.4-2.5x for `floor`, 0.88-1.25x for `ham`."""
    from folds200 import frames, parts
    from plolo import FOLDS, Operator, cd4_for, deliver
    out = {}
    for f in FOLDS:
        p = parts(f)
        p["cd4"] = cd4_for(p["targets"])
        op = Operator(p["op_path"])
        E, _ = effects_for(p, cfg, op)
        del op
        fname, b, gate = frames(p, cfg["min_cpm"])[0]        # chalA, 196's convention
        inf = deliver(E, b, gate, cfg["mass"])
        T, c = np.asarray(p["truth"], dtype=np.float64), p["scored"]
        u_inf = float(np.mean(((inf[:, c] - T[:, c]) ** 2).sum(1) / c.size))
        for rule in ("floor", "ham"):
            D, _rel = deliver_rounded(E, b, gate, cfg["mass"], 400 * 20_109, rule)
            eps = D - inf
            u = float(np.mean(((D[:, c] - T[:, c]) ** 2).sum(1) / c.size))
            e2 = float(np.mean((eps[:, c] ** 2).sum(1) / c.size))
            out[(f, rule)] = dict(dU=u - u_inf, eps2=e2, ratio=(u - u_inf) / e2 if e2 else 0.0)
        print("  [%s] " % f + "  ".join(
            f"{r}: dU {out[(f,r)]['dU']:.3e} eps2 {out[(f,r)]['eps2']:.3e} "
            f"ratio {out[(f,r)]['ratio']:.3f}" for r in ("floor", "ham")), flush=True)
        del p
    return out


def main() -> None:
    from cell_eval2.metrics.delta import PRED_TRACE_CAP_K
    from cell_eval2.moments import jackknife_correction
    t0 = time.time()
    out: dict = {}

    # ------------------------------------------------------ Phase 0: 200's own claim gates
    D = M.load(STATS)
    cpn = dict(np.load(ROOT / "experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
    bs = np.load(ROOT / "experiments/096-mse-floor/work/board_solve.npz")
    sDb, A0, B0 = float(bs["sumD"]), float(bs["v11_A"]), float(bs["v11_B"])
    Kc, X0, Kcap = float(bs["v11_Kc"]), float(bs["X"]), float(bs["K"])
    raws = {c: D[c]["dh"] for c in M.CTX}
    print(f"board solve: sumD {sDb:.6f}  X {X0:.6f}  Kc {Kc:.6f}  cap K {Kcap:.5f}")
    print(f"v11 A {A0:.6f}  B {B0:.6f} | a15 amplitude s = {S15:.5f}, "
          f"s^2 B = {S15**2 * B0:.5f}\n")
    print("--- PHASE 0: 096/164's `claim` known-answer gates ---")
    bad = []
    GATES = [("flat 0.25/phi0", None, 0, 0.25, 0.0, 0.0, 0.11860),
             ("flat 1.00/phi0", None, 0, 1.00, 0.0, 0.0, 0.19216),
             ("flat 0.25/phi4 (a8)", None, 0, 0.25, 4.0, 4.0, 0.16478),
             ("a15 as shipped", TAIL, KEEP_TOP, KEEP_SPREAD, PHI, PHI, 0.24581)]
    for name, tail, kt, ks, pk, pt, want in GATES:
        got, ch, cp = mse200.claim_of(D, raws, tail, kt, ks, pk, pt, cpn)
        ok = abs(got - want) / want < 0.03
        print(f"  {name:<24} claim {got:.5f}  expect {want:.5f}  "
              f"({100*(got-want)/want:+.2f}%)  {'PASS' if ok else 'FAIL'}")
        if not ok:
            bad.append(name)
    A_chk, B_chk = mse200.AB(D, raws, 1.0)
    print(f"  A/B round trip at s=1: A {A_chk:.6f} ({100*(A_chk/A0-1):+.3f}%)  "
          f"B {B_chk:.6f} ({100*(B_chk/B0-1):+.3f}%)")
    a15_claim, _, _ = mse200.claim_of(D, raws, TAIL, KEEP_TOP, KEEP_SPREAD, PHI, PHI, cpn)
    a15_raw = 1 + (S15**2 * A0 - 2 * S15 * X0 + Kc - min(S15**2 * B0, a15_claim)) / sDb
    print(f"  a15 raw mse {a15_raw:.5f} (expect 0.88228)  score {mse200.SCORE(a15_raw):.5f} "
          f"(expect 0.11112)")
    if abs(a15_raw - 0.88228) > 3e-4 or abs(mse200.SCORE(a15_raw) - 0.11112) > 3e-4:
        bad.append("a15 raw/score")
    if bad and not os.environ.get("NOABORT"):
        raise SystemExit(f"ABORT (rule 10): {bad}")
    print("  ALL CLAIM GATES PASS\n")
    out["anchor"] = dict(sumD=sDb, A0=A0, B0=B0, Kc=Kc, X0=X0, K=Kcap, S15=S15,
                         a15_claim=a15_claim, a15_raw=a15_raw,
                         a15_score=mse200.SCORE(a15_raw))

    # ------------------------------------------------------ ARM A: the geometry route
    print("--- ARM A: geometry, A/B/claim as RATIOS (rX = 1 null), 200's route ---")
    import chalparts200
    from plolo import PROD_OPS, Operator
    p = chalparts200.parts()
    c4 = chalparts200.cd4()
    tgt_index = {t: i for i, t in enumerate(p["targets"])}
    order = {c: np.asarray([tgt_index[t] for t in D[c]["tgt"]], dtype=np.int64)
             for c in M.CTX}
    op = Operator(PROD_OPS / "whitener_union_pruned.npz")
    grec = {}
    for name, cfg in GEOM.items():
        E = mse200.geometry_effects(cfg, p, c4, op)
        eff = {}
        for c in M.CTX:
            b = D[c]["bctrl"]
            g = mse200.gate_of(b, cfg["min_cpm"])
            h = np.where(g[None, :], E[order[c]], 0.0)
            if cfg["mass"] is not None:
                h = np.vstack([mse200.mass_fix_row(h[i], b, g, cfg["mass"])
                               for i in range(h.shape[0])])
            eff[c] = h
        Aq, Bq = mse200.AB(D, eff, 1.0)
        cl, ch, cp = mse200.claim_of(D, eff, TAIL, KEEP_TOP, KEEP_SPREAD, PHI, PHI, cpn,
                                     min_cpm=cfg["min_cpm"])
        grec[name] = dict(A=Aq, B=Bq, claim=cl, chat=ch)
        print(f"    {name:<16} A {Aq:.5f}  B {Bq:.5f}  claim {cl:.5f}  C_hat {ch:.2f}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
    del op

    def geo_mse(rec, base):
        rA, rB, rC = rec["A"]/base["A"], rec["B"]/base["B"], rec["claim"]/base["claim"]
        m = 1 + (S15**2 * A0 * rA - 2 * S15 * X0 + Kc
                 - min(S15**2 * B0 * rB, a15_claim * rC)) / sDb
        return m, rA, rB, rC

    print(f"\n  {'pair':<34}{'rA':>9}{'rB':>9}{'rClaim':>9}{'raw mse':>10}{'score':>9}"
          f"{'d overall':>11}")
    pairs = [("200's own PerturbMatch row (a23 renderer, f 0.10)", "ka200_noPM",
              "ka200_SHIP"),
             ("PerturbMatch, incumbent renderer, f 0.10", "f010_pm00", "f010_pm50"),
             ("ARM A: drop PerturbMatch, f 0.34", "A_f034_pm00", "inc_f034_pm50")]
    prec = {}
    for label, a, b in pairs:
        mb, *_ = geo_mse(grec[b], grec[b])
        ma, rA, rB, rC = geo_mse(grec[a], grec[b])
        d = (mse200.SCORE(ma) - mse200.SCORE(mb)) / 6
        prec[label] = dict(raw=ma, base_raw=mb, rA=rA, rB=rB, rClaim=rC, d_overall=d)
        print(f"  {label:<34}{rA:>9.5f}{rB:>9.5f}{rC:>9.5f}{ma:>10.5f}"
              f"{mse200.SCORE(ma):>9.5f}{d:>+11.6f}")
    ka = prec["200's own PerturbMatch row (a23 renderer, f 0.10)"]["d_overall"]
    ok = abs(ka - 0.00025551338264093865) < 2e-5
    print(f"\n  KNOWN-ANSWER (200 rank200.json, PerturbMatch `mse` = -0.000256, sign "
          f"convention flipped): {ka:+.6f}  {'PASS' if ok else 'FAIL'}")
    if not ok and not os.environ.get("NOABORT"):
        raise SystemExit("ABORT (rule 10): 200's PerturbMatch geometry row")
    out["geometry"] = grec
    out["geometry_pairs"] = prec

    # ------------------------------------------------------ ARM B: the emission route
    print(f"\n--- ARM B: Hamilton, `claim` with the INCUMBENT emitter over `cd4f034`'s own "
          f"delivered profile ---")
    print("  PLOLO cross ratio at the incumbent's configuration (f 0.34, mass None):")
    cross = plolo_cross(GEOM["inc_f034_pm50"])
    xr = {r: float(np.mean([v["ratio"] for k, v in cross.items() if k[1] == r]))
          for r in ("floor", "ham")}
    print(f"    mean cross ratio: floor {xr['floor']:.3f}   ham {xr['ham']:.3f}")

    effects = np.load(ROOT / "experiments/178-composed-riders/work/eff.a15.npy")
    perts300 = np.asarray([ln.split(",")[0] for ln in
                           (ROOT / "data/challenge/pert_counts.csv")
                           .read_text().splitlines()[1:]])
    pos300 = {str(q): i for i, q in enumerate(perts300)}
    call_set = RC.CallSet(N_HEAD, N_CPM, Z_CAP)
    chat, cp, nnz, maxcell, prof = {}, {}, {}, {}, {}
    for ctx in "ABC":
        Dc = cd4f034_chal(ctx)
        G = Dc["genes"].size
        assert effects.shape[1] == G
        Ad = ad.read_h5ad(ROOT / f"data/challenge/context_{ctx}.h5ad")
        X = Ad.X.tocsr() if not sp.isspmatrix_csr(Ad.X) else Ad.X
        ids = np.asarray(Ad.obs["ntc_id"].astype(str).values)
        cps = [float(jackknife_correction(sp.csr_matrix(X[np.flatnonzero(ids == nm)]),
                                          np.zeros(int((ids == nm).sum()), np.intp), 1, TS,
                                          chunk=64)[0])
               for nm in dict.fromkeys(ids.tolist())]
        cp[ctx] = float(np.mean(cps))
        control = RC.control_block_from_counts(X, min_cpm=5.0)
        del Ad, X
        print(f"  ctx {ctx}: C_p {cp[ctx]:.4f}, scorer gate {int(control.keep.sum())}/{G} "
              f"[{time.time()-t0:.0f}s]", flush=True)
        sel = np.linspace(0, Dc["sums"].shape[0] - 1, MAXT).astype(int)
        acc = {r: [] for r in ("floor", "ham")}
        for i in sel:
            R = recon196.Recon(Dc["sums"][i], Dc["r"], Dc["gate"], n=n, seed=1000 + int(i))
            Lp = R.emit_L(1.0)
            M1 = Dc["sums"][i].astype(np.int64)
            pcpm = 1e6 * M1 / max(int(M1.sum()), 1)
            e = np.abs(np.log2((pcpm + 1e-9) / (Dc["cpm"] + 1e-9)))
            rank = np.lexsort((np.arange(G), -np.where(Dc["gate"], e, -np.inf)))
            sv = np.full(G, 1.0)
            sv[rank[:min(KEEP_TOP, int(Dc["gate"].sum()))]] = KEEP_SPREAD
            key = np.where(Dc["gate"], effects[pos300[str(Dc["targets"][i])]], 0.0)
            si = np.flatnonzero(Dc["genes"] == Dc["targets"][i])
            for rule in ("floor", "ham"):
                Tw = R.totals(1.0) if rule == "floor" else hamilton(R.x)
                eps = np.log1p(TS * Tw / Tw.sum()) - np.log1p(TS * R.x / R.x.sum())
                if si.size:
                    eps[si[0]] = 0.0
                acc[rule].append(float((eps ** 2).sum() / (G - 1)))
                pr = (Tw + 0.5) / (n * Lp)
                rng = np.random.default_rng(90_000 + int(i))
                B, _dg = RC.rc_hybrid_block(
                    pr, np.full(n, Lp), DISPERSION, rng, pass2_rng=rng, spread=sv,
                    overdispersion=PHI, control=control, key=key,
                    jitter_seed=1000 + int(i), emit_counts=emit_counts, call_set=call_set)
                blk = sp.csr_matrix(B.astype(np.float64))
                del B
                assert np.array_equal(np.asarray(blk.sum(0)).ravel().astype(np.int64),
                                      Tw.astype(np.int64)), (ctx, rule, i)
                k = (ctx, rule)
                chat.setdefault(k, []).append(float(jackknife_correction(
                    blk, np.zeros(n, np.intp), 1, TS, chunk=64)[0]))
                nnz.setdefault(k, []).append(int(blk.nnz))
                maxcell.setdefault(k, []).append(int(np.asarray(blk.sum(1)).ravel().max()))
                del blk
        for rule in ("floor", "ham"):
            prof[(ctx, rule)] = float(np.mean(acc[rule]))
        print(f"    " + "  ".join(
            f"{r}: Chat {np.mean(chat[(ctx,r)]):.3f} (/C_p "
            f"{np.mean(chat[(ctx,r)])/cp[ctx]:.3f})  round E {prof[(ctx,r)]:.3e}"
            for r in ("floor", "ham")) + f"  [{time.time()-t0:.0f}s]", flush=True)
        del control

    print(f"\n  {'rule':>6}{'Chat/C_p':>10}{'rows<cap':>10}{'claim':>10}{'branch':>11}"
          f"{'d claim':>10}{'round E':>11}{'d prof':>10}{'raw mse':>10}{'score':>9}"
          f"{'d overall':>11}{'nnz/blk':>12}{'max cell':>10}")
    base_raw = raw_mse(S15, 0.24581)
    e_floor = float(np.mean([prof[(c, "floor")] for c in "ABC"]))
    erec = {}
    for rule in ("floor", "ham"):
        ratios, claims, below = [], [], 0
        for ctx in "ABC":
            ch = np.asarray(chat[(ctx, rule)])
            cap = PRED_TRACE_CAP_K * cp[ctx]
            ratios.append(float(ch.mean() / cp[ctx]))
            claims.append(float(np.mean(np.minimum(ch, cap)) / cp[ctx]))
            below += int((ch < cap).sum())
        ratio = float(np.mean(ratios))
        claim = float(np.mean(claims)) * CAP_K
        branch = "SATURATED" if claim >= KNEE else "linear"
        d_claim = (min(KNEE, claim) - min(KNEE, 0.24581)) / sD * -1.0
        energy = float(np.mean([prof[(c, rule)] for c in "ABC"]))
        d_prof = (energy - e_floor) * xr[rule] / sD
        rawv = base_raw + d_claim + d_prof
        mn = float(np.mean([np.mean(nnz[(c, rule)]) for c in "ABC"]))
        mx = float(np.max([np.max(maxcell[(c, rule)]) for c in "ABC"]))
        erec[rule] = dict(ratio=ratio, claim=claim, branch=branch, d_claim=d_claim,
                          energy=energy, cross=xr[rule], d_prof=d_prof, raw=rawv,
                          score=score(rawv), rows_below_cap=below, nnz=mn, max_cell=mx)
        print(f"  {rule:>6}{ratio:>10.3f}{below:>10d}{claim:>10.5f}{branch:>11}"
              f"{d_claim:>+10.5f}{energy:>11.3e}{d_prof:>+10.5f}{rawv:>10.5f}"
              f"{score(rawv):>9.5f}{(score(rawv)-score(base_raw))/6:>+11.6f}"
              f"{mn:>12,.0f}{mx:>10,.0f}")
    dB = (erec["ham"]["score"] - erec["floor"]["score"]) / 6
    print(f"\n  ARM B `mse` marginal to cd4f034: {dB:+.6f} overall  "
          f"(nnz {erec['ham']['nnz']/erec['floor']['nnz']:.4f}x)")
    out["emission"] = erec
    out["armB_d_overall"] = dB
    out["cross"] = {f"{k[0]}|{k[1]}": v for k, v in cross.items()}
    out["prof"] = {f"{k[0]}|{k[1]}": v for k, v in prof.items()}
    out["C_p"] = cp
    json.dump(out, open(WORK / "mse214.json", "w"), indent=1, default=float)
    print(f"\nwrote {WORK/'mse214.json'} [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
