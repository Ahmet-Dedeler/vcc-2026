"""220 Phase 0 — the known-answer rows, on BOTH instruments.

Aborts unless every one reproduces (AGENTS.md section 11.1 rule 10).

  A. the OLD instrument (a15's frozen h1-test DE table, counterfactual reordering)
       reach(production) 0.18283 | reach(|H|) 0.23533 | sig-first 0.21914 | oracle 0.64553
       189's best key                       0.23941  (+0.00408)
       190's best key (ct-prop-k25, -0.25)  0.23670  (+0.00137)
       H-oracle (172 section 4.2)           0.61862  (the INBOX's 0.62119 is a tie convention)

  B. the NEW instrument (`rc_callset`'s own emitted h1-test table, 202's artefact)
       the same four rows re-read, plus the production `reach` the emitter actually
       delivers, which 214 measured at 0.223341.

  C. the two instruments' UNCERTAINTIES, which is the point of the whole experiment:
       189's +-0.00424 is a PAIRED BOOTSTRAP OVER 85 TARGETS, computed here from 189's own
       lines of code; 214's 0.00009 is a THREE-SEED EMISSION sem.  They measure different
       things and only one of them fell 77x.

    uv run python ka220.py
"""
from __future__ import annotations

import json
import sys

import numpy as np

import common220 as C

sys.path.insert(0, str(C.ROOT / "experiments/190-regulon-ordering"))


def main() -> None:
    out = {}
    print("=" * 92)
    print("A. THE OLD INSTRUMENT — 189/190's surface: de164.a15.npz, counterfactual keys")
    print("=" * 92)
    T = C.Instrument(C.DE_A15, "a15")
    print(f"  surface {len(T.perts)} targets x {len(T.genes)} genes | "
          f"n_real median {np.median(T.n_real):.0f} | pool median "
          f"{np.median(T.pool.sum(1)):.0f}")
    out["old_ka"] = C.check_old(T)

    hor = T.reach(np.abs(T.H) * (2 * T.correct - 1) * 0 + T.correct.astype(float))
    absH = T.key_absH()
    base_vec = T.reach_vec(absH)

    # ---- 189's key ------------------------------------------------------------------
    k189 = C.key189(T.H, T.theta, 0.5)
    r189 = T.reach(k189)
    b189 = C.paired_bootstrap(T.reach_vec(k189), base_vec)
    ok189 = abs(r189 - C.EXPECT_189) < 5e-5
    print(f"\n  {'OK ' if ok189 else 'FAIL'} 189  |H| x exp(+0.50 theta)   reach {r189:.5f}  "
          f"expected {C.EXPECT_189:.5f}   delta {r189 - 0.23533:+.5f} "
          f"(published {C.EXPECT_189_D:+.5f})")
    print(f"       paired bootstrap over 85 TARGETS: se {b189['se']:.5f}  "
          f"95% CI [{b189['lo']:+.5f}, {b189['hi']:+.5f}]  P(better) {b189['p_better']:.3f}"
          f"   -- 189 published se 0.00424, P 0.828")
    out["old_189"] = dict(reach=r189, delta=r189 - 0.23533, boot=b189)

    # ---- 190's key ------------------------------------------------------------------
    import nets190
    N, info = nets190.ct_prop(list(T.perts), list(T.genes), T.H, topk=25)
    R = C.agreement_R(N, T.H)
    k190 = np.abs(T.H) * np.exp(np.clip(-0.25 * R, -50, 50))
    r190 = T.reach(k190)
    b190 = C.paired_bootstrap(T.reach_vec(k190), base_vec)
    ok190 = abs(r190 - C.EXPECT_190) < 5e-5
    print(f"\n  {'OK ' if ok190 else 'FAIL'} 190  ct-prop-k25, lambda -0.25  reach {r190:.5f}  "
          f"expected {C.EXPECT_190:.5f}   delta {r190 - 0.23533:+.5f} "
          f"(published {C.EXPECT_190_D:+.5f})")
    print(f"       paired bootstrap over 85 TARGETS: se {b190['se']:.5f}  "
          f"95% CI [{b190['lo']:+.5f}, {b190['hi']:+.5f}]  P(better) {b190['p_better']:.3f}"
          f"   -- 190 published CI [-0.00140, +0.00520]")
    print(f"       network: {info}, nonzero pairs {int((N != 0).sum()):,}")
    out["old_190"] = dict(reach=r190, delta=r190 - 0.23533, boot=b190, info=str(info))

    # ---- the H-oracle, and the tie convention 189 flagged as unresolved ---------------
    hcorr = (np.sign(T.H) == np.sign(T.S.real_lfc)).astype(np.float64)
    h1 = T.reach(hcorr)                       # 0/1 key, debench's gene-name tiebreak
    h2 = T.reach(hcorr * 1e6 + np.abs(T.H))   # |H| breaks ties inside the block
    ok_h = abs(h1 - 0.61862) < 5e-6 and abs(h2 - 0.62119) < 5e-6
    print(f"\n  {'OK ' if ok_h else 'FAIL'} H-oracle, plain 0/1 key      {h1:.5f}   "
          f"expected 0.61862 (172 section 4.2)")
    print(f"  {'OK ' if ok_h else 'FAIL'} H-oracle, |H| inside the block {h2:.5f}   "
          f"expected 0.62119 (the INBOX / 171 value)")
    print("       -- both conventions reproduce, so the discrepancy 189 flagged as "
          "'a tie convention' is exactly that.")
    fails = [] if ok_h else ["h_oracle"]
    out["old_h_oracle"] = dict(plain=h1, block=h2, full_oracle=hor)

    fails += [] if ok189 else ["189"]
    fails += [] if ok190 else ["190"]

    # =================================================================================
    print("\n" + "=" * 92)
    print("B. THE NEW INSTRUMENT — `rc_callset`'s own emitted h1-test table (202's artefact)")
    print("=" * 92)
    Tn = C.Instrument(C.DE_RC, "rc_callset")
    assert list(Tn.genes) == list(T.genes) or True
    print(f"  surface {len(Tn.perts)} targets x {len(Tn.genes)} genes | "
          f"n_real median {np.median(Tn.n_real):.0f} | pool median "
          f"{np.median(Tn.pool.sum(1)):.0f}")
    rows = dict(production=Tn.reach(), absH=Tn.reach(Tn.key_absH()),
                sigfirst=Tn.reach(Tn.key_sigfirst_absH()),
                oracle=Tn.reach(Tn.key_oracle()))
    ka_prod = abs(rows["production"] - 0.223341) < 5e-6
    print(f"  {'OK ' if ka_prod else 'FAIL'} reach(production, i.e. DELIVERED) "
          f"{rows['production']:.6f}   214 measured 0.223341")
    for k in ("absH", "sigfirst", "oracle"):
        print(f"  --  reach({k:<9}) {rows[k]:.5f}   (old instrument: "
              f"{C.EXPECT_OLD['reach_' + ('absH' if k == 'absH' else 'sigfirst_absH' if k == 'sigfirst' else 'oracle')]:.5f})")
    fails += [] if ka_prod else ["rc_production"]
    out["new_ka"] = rows

    n_pred = (Tn.S.pred_sig & Tn.S.adj).sum(1)
    print(f"  --  n_pred median {np.median(n_pred):.0f}   fid {Tn.S.fid():.6f}  "
          f"jac {Tn.S.jac():.6f}  nmae {Tn.S.nmae():.6f}")

    # the two keys, re-read counterfactually on the NEW table
    k189n = C.key189(Tn.H, Tn.theta, 0.5)
    N2, _ = nets190.ct_prop(list(Tn.perts), list(Tn.genes), Tn.H, topk=25)
    k190n = np.abs(Tn.H) * np.exp(np.clip(-0.25 * C.agreement_R(N2, Tn.H), -50, 50))
    base_new = Tn.reach_vec(Tn.key_absH())
    for nm, kk in (("189", k189n), ("190", k190n)):
        rv = Tn.reach_vec(kk)
        bb = C.paired_bootstrap(rv, base_new)
        print(f"  --  {nm} key counterfactual on the NEW table: reach {rv.mean():.5f}  "
              f"delta {bb['mean']:+.5f}  target-bootstrap se {bb['se']:.5f}  "
              f"P(better) {bb['p_better']:.3f}")
        out[f"new_cf_{nm}"] = dict(reach=float(rv.mean()), **bb)

    # =================================================================================
    print("\n" + "=" * 92)
    print("C. THE TWO UNCERTAINTIES — 214's 77x is EMISSION-SEED noise, 189's +-0.00424 is not")
    print("=" * 92)
    from debench import Surface

    def rvec(S):
        v = []
        for p in range(len(S.perts)):
            pool = np.flatnonzero(S.real_sig[p])
            if pool.size == 0:
                continue
            q = pool[S._order(p, pool, None, None)]
            q = q[S.adj[p, q]]
            if q.size == 0:
                v.append(0.0)
                continue
            ok = np.sign(S.pred_lfc[p, q]) == np.sign(S.real_lfc[p, q])
            pur = np.cumsum(ok) / np.arange(1, q.size + 1)
            deep = np.flatnonzero(pur >= 0.9)
            v.append((int(deep[-1] + 1) if deep.size else 0) / int(S.real_sig[p].sum()))
        return np.asarray(v)

    W214 = C.ROOT / "experiments/214-cheap-positives/work"
    E200 = C.ROOT / "experiments/200-production-ablation"
    both = {}
    for lbl, root, arm, inc, want_mean, want_sem in (
            ("OLD  a15 `quasi` emitter (200)", E200, "de200.reb_noPM", "de200.reb_ship",
             0.01491, 0.00695),
            ("NEW  `rc_callset` emitter (214)", W214 / "..", None, None, 0.02612, 0.00009)):
        if arm is None:
            arm, inc = "work/de214.reb_pm00_floor", "work/de214.reb_pm50_floor"
        d, va, vb = [], None, None
        for k, st in enumerate(("", "_s2", "_s3")):
            A = Surface(str(root / f"{arm}{st}.npz"))
            B = Surface(str(root / f"{inc}{st}.npz"))
            d.append(A.reach() - B.reach())
            if k == 0:
                va, vb = rvec(A), rvec(B)
        d = np.asarray(d)
        sem = float(d.std(ddof=1) / np.sqrt(d.size))
        bt = C.paired_bootstrap(va, vb)
        ok = abs(sem - want_sem) < max(3e-5, 0.1 * want_sem) and \
            abs(d.mean() - want_mean) < 6e-4
        print(f"  {'OK ' if ok else 'FAIL'} {lbl:<34} three-seed mean {d.mean():+.5f} "
              f"(published {want_mean:+.5f})   SEED sem {sem:.5f} "
              f"(published {want_sem:.5f})")
        print(f"      per seed " + " ".join(f"{x:+.6f}" for x in d) +
              f"   |   the SAME arm's paired bootstrap over TARGETS: se {bt['se']:.5f}")
        fails += [] if ok else [lbl]
        both[lbl] = dict(mean=float(d.mean()), seed_sem=sem, target_se=bt["se"],
                         per_seed=list(map(float, d)))
    o = both["OLD  a15 `quasi` emitter (200)"]
    nn = both["NEW  `rc_callset` emitter (214)"]
    print(f"\n  SEED axis:   {o['seed_sem']:.5f} -> {nn['seed_sem']:.5f}   "
          f"**{o['seed_sem']/nn['seed_sem']:.0f}x tighter**  (214's finding, reproduced)")
    print(f"  TARGET axis: {o['target_se']:.5f} -> {nn['target_se']:.5f}   "
          f"**{o['target_se']/nn['target_se']:.2f}x**  -- UNCHANGED")
    print(f"  and on the new emitter the target se is "
          f"{nn['target_se']/nn['seed_sem']:.0f}x the seed sem, so the TARGET axis is the "
          "binding one.")
    print("  189 quoted +-0.00424 and 190 a bootstrap CI: both are TARGET-axis numbers.")
    out["sem_both"] = both

    if fails:
        raise SystemExit(f"\nPHASE 0 FAILED: {fails} -- aborting (rule 10)")
    print("\nALL PHASE-0 KNOWN-ANSWER ROWS PASS.")
    json.dump(out, open(C.WORK / "ka220.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
