"""214 — h1-test DE tables for arms A (drop PerturbMatch) and B (Hamilton), on the
INCUMBENT's emitter.

WHY THIS IS NOT 200's TABLE.  200 priced PerturbMatch's four DE members with the `quasi`
emitter, because that was the incumbent then.  `rc_hybrid` and `rc_callset` have shipped
since, and `rc_callset` PRESCRIBES the call set (`top 400 by |H| U top 2,200 by control
CPM`).  208 measured what that does to a rider: Hamilton's `reach` half was +0.00607 raw
against `rc_hybrid` and +0.00006 against `rc_callset`.  A rider's value is not a property of
the rider (AGENTS §11.1 rule 13's corollary), so PerturbMatch must be re-priced on the
emitter that ships.

TWO FAMILIES.
  ARTIFACT (`art_*`)  the profile is a15's own delivered h1-test block, reconstructed by
                      196's `Recon`, and the ordering key is 169's recovered `H`.  These
                      carry the known-answer rows: `art_floor` must reproduce
                      `de202.rc_callset.npz` BIT-IDENTICALLY, and `art_ham` must reproduce
                      208's `de208.ham_tier`.
  REBUILD (`reb_*`)   the effect matrix is rebuilt from the donor banks (200's
                      `h1parts200`) with PerturbMatch on / off, renormalised to the
                      artifact's own median row norm so the ablation is a pure SHAPE change,
                      and each arm supplies its OWN ordering key.  This is the only way to
                      price a geometry component's DE members at all.

THE RENDERER IS THE INCUMBENT'S.  `submissions/210-cd4f034.h5ad.meta.json` records
`mass_fix: null`, so `mass_consistent_effect` is NOT applied here.  200's rebuild family used
a23's mass fix; that flag has not shipped since a23.

THE CD4 DOSE IS INVISIBLE HERE, AND THAT IS MEASURED, NOT ASSUMED.  200 found h1-test has
CD4 coverage on 4 of 85 targets, so production's median-based sizing returns fraction 0 and
`de200.reb_noCD4.npz` is byte-identical to `de200.reb_ship.npz`.  The h1-test effect matrix
is therefore the same for `rc_callset` (f = 0.10) and `cd4f034` (f = 0.34), and this script
asserts it rather than relying on 200's finding.

Donors: `k562_gwps` + `orion_hek293t` ONLY.  `data/donors/h1_2025` covers 85/85 and is pure
oracle leakage.

    uv run python emit214.py
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
for p in ("experiments/196-emitted-depth", "experiments/185-rank-controlled-emission",
          "experiments/200-production-ablation", "experiments/162-retrodicting-validator",
          "src"):
    sys.path.insert(0, str(ROOT / p))

import h1parts200                                                        # noqa: E402
import harness185 as H                                                   # noqa: E402
import recon196                                                          # noqa: E402
import vcc2026.rc_emit as RC                                             # noqa: E402
from ham196 import hamilton                                              # noqa: E402
from vcc2026.donor_transfer import emit_counts                           # noqa: E402

REF191 = ROOT / "experiments/191-rc-hybrid-build/work"
KA_TABLE = ROOT / "experiments/202-callset-rider/work/de202.rc_callset.npz"
E169 = ROOT / "experiments/169-sharpening-six-member/work/h1test-effect.npz"

n, PHI, DISPERSION = 400, 4.0, 0.1
N_HEAD, N_CPM, Z_CAP = 400, 2200, 1.90       # rc_callset's prescribed call set
KEEP_TOP, KEEP_SPREAD = 8000, 0.25           # a15's tiered spread policy (the incumbent's)
BUILDER_MIN_CPM = 10.0                       # the SHIPPED builder gate (a15..cd4f034)
SEEDS = (90_000, 7_000_000, 11_111)          # 200's matched-seed ladder
MAXT = int(os.environ.get("MAXT", "0"))
TAG = os.environ.get("OUTTAG", "de214")

# (name, geometry | "ARTIFACT", rule, seed index).  Seed index 0 reproduces 208 exactly.
GEOM = dict(
    pm50=dict(pm=0.50, gamma=0.65, lam=0.0, cd4_f=0.34, cd4_frac=None, weights=None),
    pm00=dict(pm=0.00, gamma=0.65, lam=0.0, cd4_f=0.34, cd4_frac=None, weights=None),
    pm50_f010=dict(pm=0.50, gamma=0.65, lam=0.0, cd4_f=0.10, cd4_frac=None, weights=None),
)


def arms_all():
    A = [("art_floor", "ARTIFACT", "floor", 0), ("art_ham", "ARTIFACT", "ham", 0)]
    for k in range(len(SEEDS)):
        tag = "" if k == 0 else f"_s{k+1}"
        for g in ("pm50", "pm00"):
            for rule in ("floor", "ham"):
                A.append((f"reb_{g}_{rule}{tag}", g, rule, k))
    # the CD4-dose invisibility check on h1-test, one seed
    A.append(("reb_pm50_f010_floor", "pm50_f010", "floor", 0))
    return A


def tilt_stats(T, x, gate_scorer, cpm):
    """208's `tilt_stats`, verbatim: delivered-vs-intended error and its size tilt."""
    g = gate_scorer & (x > 0)
    e = (T[g] - x[g]) / x[g]
    lo = cpm[g] <= np.quantile(cpm[g], 0.10)
    hi = cpm[g] >= np.quantile(cpm[g], 0.90)
    l2 = np.log2((T[g] + 1e-9) / x[g])
    return dict(mean_rel=float(e.mean()), mean_abs_rel=float(np.abs(e).mean()),
                common_log2=float(np.log2(T[g].sum() / x[g].sum())),
                d1_log2=float(l2[lo].mean()), d10_log2=float(l2[hi].mean()),
                tilt_log2=float(l2[lo].mean() - l2[hi].mean()))


def main() -> None:
    t0 = time.time()
    ARMS = arms_all()
    if os.environ.get("ONLY"):
        keep_arms = set(os.environ["ONLY"].split(","))
        ARMS = [a for a in ARMS if a[0] in keep_arms]
    names = [a[0] for a in ARMS]

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
    assert np.array_equal(TOT["totals"][tm].astype(np.int64), sums.astype(np.int64))
    Hd = H.load_H()
    hpos = {t: i for i, t in enumerate(Hd["targets"])}
    H0 = Hd["H"]
    art_norm = float(np.median(np.linalg.norm(H0, axis=1)))
    b_eff = np.asarray(Hd["b"], dtype=np.float64)

    # ------------------------------------------------------------------ geometry rebuilds
    need = sorted({a[1] for a in ARMS if a[1] != "ARTIFACT"})
    GEO: dict[str, np.ndarray] = {}
    if need:
        p = h1parts200.parts()
        c4 = h1parts200.cd4(p["targets"])
        assert list(p["genes"]) == list(map(str, Hd["genes"]))
        op = h1parts200.H1Operator()
        print(f"[h1-test] operator maxeig {op.evals.max():.4f} n {op.evals.size}; "
              f"artifact median row norm {art_norm:.4f}", flush=True)
        order = np.asarray([list(p["targets"]).index(t) for t in Hd["targets"]],
                           dtype=np.int64)
        for g in need:
            M = h1parts200.effects_for(GEOM[g], p, c4, op, art_norm)[order]
            GEO[g] = M
            cosv = ((M * H0).sum(1) / np.maximum(
                np.linalg.norm(M, axis=1) * np.linalg.norm(H0, axis=1), 1e-300))
            print(f"  geometry {g:<10} median row cosine to the artifact's H "
                  f"{float(np.median(cosv)):+.5f}  median norm "
                  f"{float(np.median(np.linalg.norm(M, axis=1))):.4f}", flush=True)
        del op
        # the CD4-dose invisibility assertion (200's finding, re-proved not assumed)
        if "pm50" in GEO and "pm50_f010" in GEO:
            same = bool(np.array_equal(GEO["pm50"], GEO["pm50_f010"]))
            mx = float(np.abs(GEO["pm50"] - GEO["pm50_f010"]).max())
            print(f"  CD4 dose f=0.34 vs f=0.10 on h1-test: identical {same}, "
                  f"max |diff| {mx:.3e}   (200: the median-based sizing returns 0 on a "
                  f"surface with 4/85 CD4-covered targets)", flush=True)

    if MAXT:
        targets, sums = targets[:MAXT], sums[:MAXT]
    T_n = len(targets)
    ref = np.load(REF191 / "de191.a15.npz", allow_pickle=True)
    assert list(map(str, ref["perts"]))[:T_n] == targets
    assert np.array_equal(ref["genes"], genes_full[keep])

    out = {a: {k: np.full((T_n, G), np.nan) for k in ("p", "q", "l", "z")} for a in names}
    npred = {a: [] for a in names}
    digest = {a: [] for a in names}
    diag = {a: [] for a in names}
    absH = {a: np.zeros((T_n, G), np.float32) for a in names}
    tilt = {a: [] for a in names}
    print(f"h1-test: {cs.shape[1]} controls, {G}/{G_full} scorer-gated, {T_n} "
          f"perturbations; {len(ARMS)} arms; call set top {N_HEAD} |H| U top {N_CPM} cpm "
          f"= {call_set.size}, z_cap {Z_CAP}; builder gate {BUILDER_MIN_CPM} CPM; "
          f"mass_fix None (the incumbent's renderer) [{time.time()-t0:.0f}s]", flush=True)

    b_gate = (cpm > BUILDER_MIN_CPM) & (b_eff > 0)
    for i, t in enumerate(targets):
        Rec = recon196.Recon(sums[i], D["r"], D["gate"], n=n, seed=1000 + i)
        assert np.array_equal(Rec.totals(1.0), sums[i].astype(np.int64)), f"{t}: m=1 broken"
        Lp = Rec.emit_L(1.0)
        libs = np.full(n, Lp)
        k_art = hpos[t]

        # a15's spread rank axis, frozen on the SHIPPED (floor, artifact) totals, so neither
        # the rounding rule nor the geometry moves WHICH genes the tier selects.
        M1 = sums[i].astype(np.int64)
        pcpm = 1e6 * M1 / max(int(M1.sum()), 1)
        effr = np.abs(np.log2((pcpm + 1e-9) / (cpm + 1e-9)))
        rank = np.lexsort((np.arange(G_full), -np.where(keep, effr, -np.inf)))
        sv = np.full(G_full, 1.0)
        sv[rank[:min(KEEP_TOP, int(keep.sum()))]] = KEEP_SPREAD

        # the intended x = n L p for each geometry
        X = {"ARTIFACT": Rec.x}
        KEY = {"ARTIFACT": H0[k_art]}
        for g in GEO:
            h = np.where(b_gate, GEO[g][k_art], 0.0)      # mass_fix None: the incumbent
            lin = np.clip(np.expm1(b_eff + h), 0.0, None)
            X[g] = (n * Lp) * lin / lin.sum()
            KEY[g] = GEO[g][k_art]

        for name, geom, rule, ks in ARMS:
            x = X[geom]
            Tw = (np.floor(x).astype(np.int64) if rule == "floor"
                  else hamilton(x).astype(np.int64))
            tilt[name].append(tilt_stats(Tw, x, keep, cpm))
            prof = (Tw + 0.5) / (n * Lp)
            rng = np.random.default_rng(SEEDS[ks] + i)
            B, dg = RC.rc_hybrid_block(
                prof, libs, DISPERSION, rng, pass2_rng=rng, spread=sv,
                overdispersion=PHI, control=control, key=KEY[geom],
                jitter_seed=1000 + i + 100_000 * ks,
                emit_counts=emit_counts, call_set=call_set)
            delivered = B.sum(0).astype(np.int64)
            assert np.array_equal(delivered, Tw), (name, t)
            digest[name].append(hashlib.sha256(delivered.tobytes()).hexdigest()[:16])
            rows = B.sum(1)
            dg.update(nnz=int((B > 0).sum()), max_cell=int(rows.max()),
                      cell_total_sd=float(rows.std()), panel_S=int(delivered.sum()),
                      emit_L=int(Lp))
            zz, pv, qv, lfc = H.de_from_block(B, keep, cs, tie, ref_lin)
            out[name]["z"][i], out[name]["p"][i] = zz, pv
            out[name]["q"][i], out[name]["l"][i] = qv, lfc
            npred[name].append(int((qv < 0.05).sum()))
            diag[name].append(dg)
            absH[name][i] = np.abs(KEY[geom])[keep]
            del B
        if (i + 1) % 10 == 0 or i == 0:
            print(f"  {i+1}/{T_n} [{time.time()-t0:.0f}s]  " +
                  " ".join(f"{a}={npred[a][-1]}" for a in names[:4]), flush=True)

    # ---- known-answer row (rule 10) ------------------------------------------------------
    ok = True
    if "art_floor" in names and KA_TABLE.exists():
        print("\nKNOWN-ANSWER ROW (rule 10): `art_floor` IS `rc_callset` and must reproduce "
              "202's emitted table bit for bit.")
        z = np.load(KA_TABLE, allow_pickle=True)
        for mine, theirs in (("p", "pred_p"), ("q", "pred_padj"), ("l", "pred_lfc"),
                             ("z", "pred_z")):
            same = bool(np.array_equal(out["art_floor"][mine], z[theirs][:T_n]))
            worst = float(np.nanmax(np.abs(out["art_floor"][mine] - z[theirs][:T_n])))
            print(f"  {theirs:<10} bit-identical on {T_n} x {G}: "
                  f"{'YES' if same else 'NO'}   max |diff| {worst:.3g}")
            ok &= same
        sd = sum(a == b for a, b in zip(digest["art_floor"],
                                        list(map(str, z["digest"][:T_n]))))
        print(f"  per-gene panel-total digest identical on {sd}/{T_n}")
        ok &= sd == T_n
        print(f"  KNOWN-ANSWER ROW: {'PASS' if ok else 'FAIL'}")
    if not ok and not os.environ.get("NOABORT"):
        raise SystemExit("ABORT (rule 10): art_floor is not rc_callset")

    print("\nPER-GENE PANEL TOTAL DIGESTS vs `art_floor` (rule 13's class test)")
    base = digest.get("art_floor")
    if base:
        for a in names:
            same = sum(x == y for x, y in zip(digest[a], base))
            print(f"  {a:<24} identical to art_floor on {same}/{T_n} blocks")

    print("\nDELIVERED-vs-INTENDED TILT over the scorer-gated pool (median over blocks)")
    print(f"{'arm':<24}{'mean rel':>11}{'common log2':>13}{'TILT log2':>12}")
    tsum = {}
    for a in names:
        med = {j: float(np.median([x[j] for x in tilt[a]])) for j in tilt[a][0]}
        tsum[a] = med
        print(f"{a:<24}{med['mean_rel']:>+11.5f}{med['common_log2']:>+13.6f}"
              f"{med['tilt_log2']:>+12.5f}")

    for a in names:
        np.savez_compressed(
            WORK / f"{TAG}.{a}.npz", genes=genes_full[keep], perts=np.array(targets),
            real_p=ref["real_p"][:T_n], real_padj=ref["real_padj"][:T_n],
            real_lfc=ref["real_lfc"][:T_n], pred_p=out[a]["p"], pred_padj=out[a]["q"],
            pred_lfc=out[a]["l"], pred_z=out[a]["z"], absH=absH[a],
            digest=np.array(digest[a]))
        d = diag[a]
        summ = {k: float(np.median([x[k] for x in d]))
                for k in ("nnz", "max_cell", "cell_total_sd", "panel_S", "emit_L")}
        print(f"wrote {TAG}.{a}.npz  median n_pred {np.median(npred[a]):.0f}  {summ}",
              flush=True)
    json.dump(dict(npred=npred, diag=diag, arms=names, tilt=tsum, digest=digest),
              open(WORK / f"{TAG}_diag.json", "w"), default=float)
    print(f"done [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
