"""208 — the two free riders, emitted on h1-test on top of `rc_callset` (the incumbent).

FOUR ARMS, two binary variables, factorial so the composition can be MEASURED rather than
inferred (178: "gains compose sub-additively, costs do not"):

  rule   `floor`  `totals = floor(n L p)`, what `emit_counts` ships           (196)
         `ham`    largest-remainder over the same `n L p`                     (rider A)
  spread `tier`   a15's tail 1.0 / top-8000 block 0.25                        (incumbent)
         `flat`   flat 1.00 everywhere                                        (rider B)

Everything else is `rc_callset`'s, bit for bit: the emitter is
`rc_emit.rc_hybrid_block(..., call_set=CallSet(400, 2200, 1.90))`, the seeds are
production's single stream (`rng = default_rng(90000 + i)`, `pass2_rng = rng`), the jitter
seed is `1000 + i`, `phi = 4.0`, `DISPERSION = 0.1`, the reference side is 191's own
`de191.a15.npz`, and the h1-test control block is 185's cache.

KNOWN-ANSWER ROW (rule 10, aborts).  `floor_tier` IS `rc_callset`.  It must reproduce
`experiments/202-callset-rider/work/de202.rc_callset.npz` BIT-IDENTICALLY on all
85 x 10,779 entries of `pred_p`, `pred_padj`, `pred_lfc`, `pred_z` and on the per-gene
panel-total digest -- or this code path is not the one the board scored.

THE ROUNDING RULE IS NOT PSEUDOBULK-PRESERVING.  `ham` moves the delivered per-gene panel
totals themselves (by at most one count each), so `pds`, `mse` and `nmae` are NOT pinned for
any arm carrying it and are priced separately (`pds208.py`, `mse208.py`).  `flat` alone IS
pseudobulk-preserving -- `emit_counts` computes `totals` before `spread` is ever read -- and
this script asserts that identity per block rather than claiming it.

    ARMS=floor_tier,ham_tier,floor_flat,ham_flat uv run python emit208.py
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/196-emitted-depth"))
sys.path.insert(0, str(ROOT / "experiments/185-rank-controlled-emission"))
sys.path.insert(0, str(ROOT / "src"))

import harness185 as H  # noqa: E402
import recon196  # noqa: E402
import vcc2026.rc_emit as RC  # noqa: E402
from ham196 import hamilton  # noqa: E402
from vcc2026.donor_transfer import emit_counts  # noqa: E402

REF191 = ROOT / "experiments/191-rc-hybrid-build/work"
KA_TABLE = ROOT / "experiments/202-callset-rider/work/de202.rc_callset.npz"
n, PHI, DISPERSION = 400, 4.0, 0.1
N_HEAD, N_CPM, Z_CAP = 400, 2200, 1.90
KEEP_TOP, KEEP_SPREAD = 8000, 0.25
ARMS = tuple(os.environ.get(
    "ARMS", "floor_tier,ham_tier,floor_flat,ham_flat").split(","))
MAXT = int(os.environ.get("MAXT", "0"))
TAG = os.environ.get("OUTTAG", "de208")


def main() -> None:
    t0 = time.time()
    C = H.load_controls()
    cs, tie, keep, cpm, ref_lin = C["cs"], C["tie"], C["keep"], C["cpm"], C["ref_lin"]
    genes_full = C["genes"]
    G_full, G = genes_full.size, cs.shape[0]
    control = RC.RcControlBlock(keep=keep, cs=cs, tie=tie, ref_lin=ref_lin,
                                n_controls=cs.shape[1], min_cpm=H.MIN_CPM)
    call_set = RC.CallSet(N_HEAD, N_CPM, Z_CAP)

    D = recon196.load_h1()
    assert np.array_equal(D["genes"], genes_full)
    targets, sums = list(D["targets"]), D["sums"]
    TOT = H.load_totals()
    tm = TOT["targets"] != "non-targeting"
    assert list(TOT["targets"][tm]) == targets
    assert np.array_equal(TOT["totals"][tm].astype(np.int64), sums.astype(np.int64)), \
        "196's reconstruction source and 185's a15 totals disagree"
    Hd = H.load_H()
    hpos = {t: i for i, t in enumerate(Hd["targets"])}
    if MAXT:
        targets, sums = targets[:MAXT], sums[:MAXT]
    T_n = len(targets)

    ref = np.load(REF191 / "de191.a15.npz", allow_pickle=True)
    assert list(map(str, ref["perts"]))[:T_n] == targets
    assert np.array_equal(ref["genes"], genes_full[keep])

    out = {a: {k: np.full((T_n, G), np.nan) for k in ("p", "q", "l", "z")} for a in ARMS}
    npred = {a: [] for a in ARMS}
    digest = {a: [] for a in ARMS}
    diag = {a: [] for a in ARMS}
    absH_all = np.zeros((T_n, G), np.float32)
    tilt = {r: [] for r in ("floor", "ham")}
    tot_move = []
    print(f"h1-test: {cs.shape[1]} controls, {G}/{G_full} scorer-gated, {T_n} perturbations; "
          f"arms {list(ARMS)}; call set top {N_HEAD} |H| U top {N_CPM} cpm = "
          f"{call_set.size}, z_cap {Z_CAP} [{time.time()-t0:.0f}s]", flush=True)

    for i, t in enumerate(targets):
        Rec = recon196.Recon(sums[i], D["r"], D["gate"], n=n, seed=1000 + i)
        assert np.array_equal(Rec.totals(1.0), sums[i].astype(np.int64)), f"{t}: m=1 broken"
        Lp = Rec.emit_L(1.0)
        libs = np.full(n, Lp)

        # a15's spread rank axis, frozen on the SHIPPED (floor) totals, so the ROUNDING
        # rule never moves which genes the tier selects -- one variable at a time.
        M1 = sums[i].astype(np.int64)
        pcpm = 1e6 * M1 / max(int(M1.sum()), 1)
        eff = np.abs(np.log2((pcpm + 1e-9) / (cpm + 1e-9)))
        order = np.lexsort((np.arange(G_full), -np.where(keep, eff, -np.inf)))
        sv_tier = np.full(G_full, 1.0)
        sv_tier[order[:min(KEEP_TOP, int(keep.sum()))]] = KEEP_SPREAD
        assert float(sv_tier[~keep].max()) == 1.0
        sv_flat = np.full(G_full, 1.0)
        SV = {"tier": sv_tier, "flat": sv_flat}

        key = Hd["H"][hpos[t]]
        absH_all[i] = np.abs(key)[keep]
        want = {"floor": Rec.totals(1.0), "ham": hamilton(Rec.x)}
        for r in ("floor", "ham"):
            tilt[r].append(tilt_stats(want[r], Rec.x, keep, cpm))
        dT = want["ham"] - want["floor"]
        tot_move.append(dict(max_abs=int(np.abs(dT).max()),
                             moved=int((dT != 0).sum()),
                             moved_gated=int((dT[keep] != 0).sum()),
                             sum_floor=int(want["floor"].sum()),
                             sum_ham=int(want["ham"].sum())))

        for a in ARMS:
            rule, spread_name = a.split("_")
            Tw = want[rule]
            prof = (Tw + 0.5) / (n * Lp)
            rng = np.random.default_rng(90_000 + i)
            B, dg = RC.rc_hybrid_block(
                prof, libs, DISPERSION, rng, pass2_rng=rng, spread=SV[spread_name],
                overdispersion=PHI, control=control, key=key, jitter_seed=1000 + i,
                emit_counts=emit_counts, call_set=call_set)
            delivered = B.sum(0).astype(np.int64)
            assert np.array_equal(delivered, Tw), (a, t)
            digest[a].append(hashlib.sha256(delivered.tobytes()).hexdigest()[:16])
            rows = B.sum(1)
            dg.update(nnz=int((B > 0).sum()), max_cell=int(rows.max()),
                      cell_total_sd=float(rows.std()), panel_S=int(delivered.sum()),
                      emit_L=int(Lp))
            zz, p, q, lfc = H.de_from_block(B, keep, cs, tie, ref_lin)
            out[a]["z"][i], out[a]["p"][i] = zz, p
            out[a]["q"][i], out[a]["l"][i] = q, lfc
            npred[a].append(int((q < 0.05).sum()))
            diag[a].append(dg)
            del B
        if (i + 1) % 10 == 0 or i == 0:
            print(f"  {i+1}/{T_n} [{time.time()-t0:.0f}s]  " +
                  " ".join(f"{a}={npred[a][-1]}" for a in ARMS), flush=True)

    # ---- the pseudobulk-identity claims, ASSERTED per block, not argued -----------------
    print("\nPER-GENE PANEL TOTAL DIGESTS (falsifier 1, h1-test side)")
    base = digest.get("floor_tier")
    for a in ARMS:
        if base is None:
            break
        same = sum(x == y for x, y in zip(digest[a], base))
        print(f"  {a:<12} identical to floor_tier on {same}/{T_n} blocks")
    mv = {k: float(np.median([x[k] for x in tot_move])) for k in tot_move[0]}
    print(f"  ham-vs-floor totals: max |dT| per block median {mv['max_abs']:.0f}, "
          f"worst {max(x['max_abs'] for x in tot_move)}; genes moved median "
          f"{mv['moved']:.0f} (of {G_full}), gated {mv['moved_gated']:.0f} (of {G}); "
          f"panel S {mv['sum_floor']:,.0f} -> {mv['sum_ham']:,.0f}")

    print("\nDELIVERED-vs-INTENDED over the scorer-gated pool (median over perturbations)")
    print(f"{'rule':<8}{'mean rel':>11}{'mean |rel|':>12}{'common log2':>13}"
          f"{'dec1 log2':>11}{'dec10 log2':>12}{'TILT log2':>11}")
    tsum = {}
    for r in ("floor", "ham"):
        med = {j: float(np.median([x[j] for x in tilt[r]])) for j in tilt[r][0]}
        tsum[r] = med
        print(f"{r:<8}{med['mean_rel']:>+11.5f}{med['mean_abs_rel']:>12.5f}"
              f"{med['common_log2']:>+13.6f}{med['d1_log2']:>+11.5f}"
              f"{med['d10_log2']:>+12.5f}{med['tilt_log2']:>+11.5f}")

    # ---- KNOWN-ANSWER ROW (rule 10) ----------------------------------------------------
    ok = True
    if "floor_tier" in ARMS and KA_TABLE.exists():
        print("\nKNOWN-ANSWER ROW (rule 10): `floor_tier` IS `rc_callset` and must reproduce "
              "202's own emitted table bit for bit.")
        z = np.load(KA_TABLE, allow_pickle=True)
        for mine, theirs in (("p", "pred_p"), ("q", "pred_padj"), ("l", "pred_lfc"),
                             ("z", "pred_z")):
            same = bool(np.array_equal(out["floor_tier"][mine], z[theirs][:T_n]))
            worst = float(np.nanmax(np.abs(out["floor_tier"][mine] - z[theirs][:T_n])))
            print(f"  {theirs:<10} bit-identical on {T_n} x {G}: "
                  f"{'YES' if same else 'NO'}   max |diff| {worst:.3g}")
            ok &= same
        sd = sum(x == y for x, y in zip(digest["floor_tier"],
                                        list(map(str, z["digest"][:T_n]))))
        print(f"  per-gene panel-total digest identical on {sd}/{T_n}")
        ok &= sd == T_n
        print(f"  KNOWN-ANSWER ROW: {'PASS' if ok else 'FAIL'}")
    if not ok and not os.environ.get("NOABORT"):
        raise SystemExit("ABORT (rule 10): floor_tier is not rc_callset")

    for a in ARMS:
        np.savez_compressed(
            WORK / f"{TAG}.{a}.npz", genes=genes_full[keep], perts=np.array(targets),
            real_p=ref["real_p"][:T_n], real_padj=ref["real_padj"][:T_n],
            real_lfc=ref["real_lfc"][:T_n], pred_p=out[a]["p"],
            pred_padj=out[a]["q"], pred_lfc=out[a]["l"], pred_z=out[a]["z"],
            absH=absH_all, digest=np.array(digest[a]))
        d = diag[a]
        summ = {k: float(np.median([x[k] for x in d])) for k in
                ("nnz", "max_cell", "cell_total_sd", "panel_S", "emit_L")}
        print(f"wrote {TAG}.{a}.npz  median n_pred {np.median(npred[a]):.0f}  {summ}",
              flush=True)
    json.dump(dict(npred=npred, diag=diag, arms=list(ARMS), tilt=tsum,
                   totals_move=tot_move, digest=digest),
              open(WORK / f"{TAG}_diag.json", "w"), default=float)
    print(f"done [{time.time()-t0:.0f}s]")


def tilt_stats(T, x, gate_scorer, cpm):
    """Delivered-vs-intended error over the SCORER-gated pool, and its size tilt (196)."""
    g = gate_scorer & (x > 0)
    e = (T[g] - x[g]) / x[g]
    lo = cpm[g] <= np.quantile(cpm[g], 0.10)
    hi = cpm[g] >= np.quantile(cpm[g], 0.90)
    l2 = np.log2((T[g] + 1e-9) / x[g])
    return dict(mean_rel=float(e.mean()), mean_abs_rel=float(np.abs(e).mean()),
                common_log2=float(np.log2(T[g].sum() / x[g].sum())),
                d1_log2=float(l2[lo].mean()), d10_log2=float(l2[hi].mean()),
                tilt_log2=float(l2[lo].mean() - l2[hi].mean()))


if __name__ == "__main__":
    main()
