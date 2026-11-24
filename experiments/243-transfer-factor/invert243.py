"""243 Phase 2 — the board side is an INVERSION, not a measurement.  Put the arms' own
amplitude terms into 197's identity and read the cosine each arm actually carries.

197 §1, per context c (187's correction of 096):

    MSE_c - 1 = ( A_c + Kc_c - min(claim_c, B_c) - 2 X_c ) / sumD,      raw = mean_c MSE_c

with  A_c = sum_p w_p ||dh_p||^2  (the arm's delivered energy, own-target coordinate dropped),
      B_c = cell_eval2's `_across_pert_budget` on the arm's delivered pseudobulk,
      Kc_c = the reference's control jackknife term (096 (e)),
      claim_c <= cap_c = (sum_p w_p) C_p  (hard cap; 211: the rc_callset emitter sits ON it),
      X_c  = sum_p w_p <dh_p, d_p>  (alignment with the withheld truth),
      sumD = sum_p w_p ||d_p||^2   (reference property, 187 [KA-i]).

Writing r_c = sqrt(A_c / sumD) and c = X / (sqrt(sumD) * mean_c sqrt(A_c)):

    raw = 1 + mean_c[ r_c^2 - (min(claim_c, B_c) - Kc_c)/sumD ] - 2 c * mean_c sqrt(A_c)/sqrt(sumD)

i.e. 197's "raw = r^2 - 2 r c + 1" plus the budget/jackknife correction.  The one-sided FLOOR
`cos >= (1 + Kc/sumD - raw) / (2 sqrt(K/sumD))` replaces the arm's own (A, B, claim) with the
reference constants, so between two arms it attributes EVERY change in raw to c.  This script
measures A_c and B_c on the three shipped artifacts (rc_callset = a15's profile, cd4f034,
h1union) and reads the cosine with the amplitude terms in place.

KNOWN-ANSWER ROWS
  * b_ctrl from 096's ctrl-stats.npz == 162's cached bctrl-challenge-{A,B,C}.npy
  * rc_callset and cd4f034 A_c / B_c / cap_c == 211's bx211.json (1e-9)
  * rc_callset A_c == 197's a15 terms (196's independently streamed totals)
  * a14 and a15 measured cosines == 197's 0.074175 / 0.079244; rc_callset == 0.079249
  * my streamed cd4f034 totals == 214's cached totals on 900/900 blocks (totals243.py)

    uv run python invert243.py
"""
from __future__ import annotations

import json

import numpy as np
from cell_eval2.metrics.delta import _across_pert_budget, _row_weights

from lib243 import (COEFF, CP46, CTX, FLOOR0, INV197, KC, N_GENES, PDS_PER_COS_215, ROOT, SUMD,
                    TS, WORK, board_row, floor_cos, rss_gb)

TOTALS = {
    "rc_callset": ROOT / "experiments/214-cheap-positives/work/totals-rc_callset.npz",
    "cd4f034": WORK / "totals-cd4f034.npz",          # streamed here, == 214's cache (KA)
    "h1union": WORK / "totals-h1union.npz",          # streamed here (the one new measurement)
}
BX211 = json.loads((ROOT / "experiments/211-amplitude-reopened/work/bx211.json").read_text())
LINES = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LINES.append(s)


def bctrl():
    z = np.load(ROOT / "experiments/096-mse-floor/work/ctrl-stats.npz", allow_pickle=True)
    genes = np.asarray(z["genes"]).astype(str)
    out = {}
    for c in CTX:
        s = z[f"sum_{c}"].astype(np.float64)
        out[c] = np.log1p(TS * s / s.sum())
        cached = np.load(ROOT / f"experiments/162-retrodicting-validator/work/bctrl-challenge-{c}.npy")
        d = float(np.abs(out[c] - cached).max())
        P(f"  b_ctrl {c}: 096 ctrl-stats vs 162 cache max|d| {d:.2e}  {'PASS' if d < 1e-9 else 'FAIL'}")
        assert d < 1e-9
    return genes, out


def terms(tot, blocks, genes, b):
    """211 bx211.py's per-context (A, B, Kc, cap) on per-block panel totals, verbatim."""
    gpos = {g: i for i, g in enumerate(genes)}
    per = {}
    for c in CTX:
        rows = [i for i, k in enumerate(blocks) if k.split("|")[0] == c]
        T = tot[rows].astype(np.float64)
        mu = np.log1p(TS * T / T.sum(1, keepdims=True))
        dh = mu - b[c][None, :]
        tgts = [blocks[i].split("|")[1] for i in rows]
        cols = {j: gpos[t] for j, t in enumerate(tgts) if t in gpos}
        w = _row_weights(len(rows), N_GENES, cols)
        a = np.einsum("ij,ij->i", dh, dh)
        for j, col in cols.items():
            a[j] -= dh[j, col] ** 2
        per[c] = dict(A=float((a * w).sum()), B=float(_across_pert_budget(mu, cols, w)),
                      Kc=KC[c], cap=float(w.sum()) * CP46[c], n=len(rows), resolved=len(cols),
                      wsum=float(w.sum()))
    return per


def invert(per, raw, claim_key="cap"):
    """X_bar from one board raw; cosine on 216/197's functional X/(sqrt(sumD) mean_c sqrt(A_c))."""
    num = float(np.mean([per[c]["A"] + per[c]["Kc"] - min(per[c][claim_key], per[c]["B"])
                         for c in CTX]))
    X = (num - (raw - 1.0) * SUMD) / 2.0
    msa = float(np.mean([np.sqrt(per[c]["A"]) for c in CTX]))
    abar = float(np.mean([per[c]["A"] for c in CTX]))
    return dict(X=X, num=num, cos=X / (np.sqrt(SUMD) * msa), cos_abar=X / np.sqrt(abar * SUMD),
                mean_sqrtA=msa, Abar=abar,
                legs={c: ("claim" if per[c][claim_key] < per[c]["B"] else "B") for c in CTX},
                r={c: float(np.sqrt(per[c]["A"] / SUMD)) for c in CTX})


def floor_sensitivity(per, X, claim_key="cap", ds=0.005):
    """d(floor cos)/d(ln s) for a global amplitude rescale s of the delivered effect
    (A -> s^2 A, B -> s^2 B, X -> s X) at s = 1, with X_c split in proportion to sqrt(A_c)."""
    msa = float(np.mean([np.sqrt(per[c]["A"]) for c in CTX]))

    def raw_at(s):
        return 1.0 + float(np.mean([
            s * s * per[c]["A"] + per[c]["Kc"] - min(per[c][claim_key], s * s * per[c]["B"])
            - 2.0 * s * X * np.sqrt(per[c]["A"]) / msa for c in CTX])) / SUMD
    fp, fm = floor_cos(raw_at(1 + ds)), floor_cos(raw_at(1 - ds))
    return (fp - fm) / (np.log(1 + ds) - np.log(1 - ds)), raw_at(1.0)


def main() -> None:
    genes, b = bctrl()
    arms = {}
    for name, path in TOTALS.items():
        z = np.load(path, allow_pickle=True)
        assert np.array_equal(np.asarray(z["genes"]).astype(str), genes), f"{name}: gene axis"
        blocks = np.asarray(z["blocks"]).astype(str)
        arms[name] = terms(np.asarray(z["totals"]), blocks, genes, b)
        P(f"  {name}: {blocks.size} blocks from {path.relative_to(ROOT)}  rss {rss_gb():.2f} GB")

    # ---- known-answer rows -----------------------------------------------------------------
    P("\n=== KA: A_c / B_c / cap_c against 211's bx211.json and 197's a15 terms ===")
    ok = True
    for name in ("rc_callset", "cd4f034"):
        for c in CTX:
            ref = BX211["per_context"][name][c]
            for k, rk in (("A", "A"), ("B", "B"), ("cap", "claim_cap")):
                d = abs(arms[name][c][k] - ref[rk])
                ok &= d < 1e-9
                if d >= 1e-9:
                    P(f"  FAIL {name} {c} {k}: {arms[name][c][k]:.12f} vs 211 {ref[rk]:.12f}")
    for c in CTX:
        d = abs(arms["rc_callset"][c]["A"] - INV197["ours"]["a15"]["per"][c]["A"])
        ok &= d < 1e-9
    P(f"  {'PASS' if ok else 'FAIL'}  rc_callset and cd4f034 (A, B, cap) == bx211.json on 3 contexts; "
      f"rc_callset A_c == 197's a15 A_c (196's totals)")
    assert ok

    # a14 / a15 from 197's own stored terms (a14 = a12's artifact rescaled, quasi emitter,
    # claim = proxy; a15 = 196's totals + 191's claim) -> reproduce 197's measured cosines
    ka = {}
    for arm, want in (("a14", 0.07417525298473132), ("a15", 0.07924410082529984)):
        per = {c: dict(INV197["ours"][arm]["per"][c]) for c in CTX}
        r = invert(per, INV197["ours"][arm]["raw"], claim_key="claim")
        ka[arm] = r
        P(f"  {'PASS' if abs(r['cos_abar'] - want) < 1e-9 else 'FAIL'}  197 {arm}: cos(Abar form) "
          f"{r['cos_abar']:.9f} vs {want:.9f};  cos(mean-sqrtA form) {r['cos']:.9f}")
        assert abs(r["cos_abar"] - want) < 1e-9

    # ---- the three rc-emitter arms, claim at the cap (211 chat211: 0 rows under the cap) ---
    P("\n=== the arms, inverted with their OWN amplitude terms (claim at the cap) ===")
    P(f"  {'arm':11s}{'raw mse':>13s}{'A_A':>10s}{'A_B':>10s}{'A_C':>10s}{'B/A':>8s}"
      f"{'legs':>10s}{'X':>10s}{'cos meas':>10s}{'floor':>9s}{'tight':>7s}{'r mean':>8s}")
    res = {}
    for name in ("rc_callset", "cd4f034", "h1union"):
        bd = board_row(name)
        per = arms[name]
        r = invert(per, bd["raw"]["mse"])
        fl = floor_cos(bd["raw"]["mse"])
        sens, raw_chk = floor_sensitivity(per, r["X"])
        assert abs(raw_chk - bd["raw"]["mse"]) < 1e-12
        res[name] = dict(board=bd, per=per, inv=r, floor=fl, tight=fl / r["cos"],
                         floor_sens_per_ln_s=sens, r_mean=float(np.mean(list(r["r"].values()))),
                         BA=float(np.mean([per[c]["B"] / per[c]["A"] for c in CTX])))
        P(f"  {name:11s}{bd['raw']['mse']:>13.10f}{per['A']['A']:>10.6f}{per['B']['A']:>10.6f}"
          f"{per['C']['A']:>10.6f}{res[name]['BA']:>8.5f}"
          f"{''.join(r['legs'][c][0] for c in CTX):>10s}{r['X']:>10.6f}{r['cos']:>10.6f}"
          f"{fl:>9.6f}{fl/r['cos']:>7.3f}{res[name]['r_mean']:>8.4f}")
    P(f"  (a15 from 197's own terms: cos {ka['a15']['cos']:.6f}; rc_callset above must equal it "
      f"-> |d| {abs(ka['a15']['cos'] - res['rc_callset']['inv']['cos']):.2e})")
    P(f"  (a14 from 197's own terms: cos {ka['a14']['cos']:.6f}, floor "
      f"{floor_cos(INV197['ours']['a14']['raw']):.6f})")

    # ---- the decomposition of each board step --------------------------------------------
    P("\n=== each board step: floor d(cos) vs measured d(cos), and where the raw moved ===")
    steps = [("a14->a15", None, "rc_callset"), ("cd4f034", "rc_callset", "cd4f034"),
             ("h1union", "cd4f034", "h1union")]
    dec = {}
    for step, a, bname in steps:
        if a is None:
            pa = {c: dict(INV197["ours"]["a14"]["per"][c]) for c in CTX}
            ra = ka["a14"]
            raw_a = INV197["ours"]["a14"]["raw"]
            ck_a = "claim"
            pds_a = board_row("a14")["raw"]["pds"]
        else:
            pa, ra, raw_a, ck_a = arms[a], res[a]["inv"], res[a]["board"]["raw"]["mse"], "cap"
            pds_a = res[a]["board"]["raw"]["pds"]
        pb, rb, raw_b = arms[bname], res[bname]["inv"], res[bname]["board"]["raw"]["mse"]
        pds_b = res[bname]["board"]["raw"]["pds"]
        term_a = float(np.mean([pa[c]["A"] + pa[c]["Kc"] - min(pa[c][ck_a], pa[c]["B"]) for c in CTX]))
        term_b = float(np.mean([pb[c]["A"] + pb[c]["Kc"] - min(pb[c]["cap"], pb[c]["B"]) for c in CTX]))
        d_raw = raw_b - raw_a
        d_raw_term = (term_b - term_a) / SUMD             # amplitude + budget + jackknife leg
        d_raw_align = -2.0 * (rb["X"] - ra["X"]) / SUMD    # alignment leg
        assert abs(d_raw_term + d_raw_align - d_raw) < 1e-12
        fl_a, fl_b = floor_cos(raw_a), floor_cos(raw_b)
        d_cos_floor = fl_b - fl_a
        d_cos_meas = rb["cos"] - ra["cos"]
        # the floor's d(cos) split: the part that is really amplitude/budget, in floor units
        d_cos_floor_from_term = -d_raw_term / COEFF
        d_cos_floor_from_align = -d_raw_align / COEFF
        # 215's amplitude-free second equation
        d_cos_pds = (pds_b - pds_a) / PDS_PER_COS_215
        d_cos_mse_215 = -d_raw / COEFF
        dec[step] = dict(
            raw_a=raw_a, raw_b=raw_b, d_raw=d_raw, d_raw_term=d_raw_term, d_raw_align=d_raw_align,
            cos_meas_a=ra["cos"], cos_meas_b=rb["cos"], rel_meas=rb["cos"] / ra["cos"] - 1.0,
            floor_a=fl_a, floor_b=fl_b, rel_floor=fl_b / fl_a - 1.0,
            d_cos_floor=d_cos_floor, d_cos_floor_from_term=d_cos_floor_from_term,
            d_cos_floor_from_align=d_cos_floor_from_align, d_cos_meas=d_cos_meas,
            amp_ratio=float(np.mean([np.sqrt(pb[c]["A"] / pa[c]["A"]) for c in CTX])),
            mean_sqrtA_a=ra["mean_sqrtA"], mean_sqrtA_b=rb["mean_sqrtA"],
            d_pds_board=pds_b - pds_a, d_cos_pds_215=d_cos_pds, d_cos_mse_215=d_cos_mse_215,
            rho_215=d_cos_pds / d_cos_mse_215 if d_cos_mse_215 else float("nan"),
            rel_pds_route=d_cos_pds / ra["cos"])
        P(f"\n  {step}:  raw {raw_a:.7f} -> {raw_b:.7f}  (d {d_raw:+.7f})")
        P(f"    d raw from the amplitude/budget leg  {d_raw_term:+.7f}   from the alignment leg "
          f"{d_raw_align:+.7f}")
        P(f"    delivered amplitude sqrt(A_b/A_a) - 1   {100*(dec[step]['amp_ratio']-1):+.3f}%")
        P(f"    FLOOR    cos {fl_a:.6f} -> {fl_b:.6f}   d {d_cos_floor:+.6f}  rel {100*(fl_b/fl_a-1):+.3f}%"
          f"   [of which amplitude/budget {d_cos_floor_from_term:+.6f}, alignment "
          f"{d_cos_floor_from_align:+.6f}]")
        P(f"    MEASURED cos {ra['cos']:.6f} -> {rb['cos']:.6f}   d {d_cos_meas:+.6f}  rel "
          f"{100*(rb['cos']/ra['cos']-1):+.3f}%")
        P(f"    215 pds route: d raw pds {pds_b-pds_a:+.6f} -> d cos {d_cos_pds:+.6f}  rel "
          f"{100*d_cos_pds/ra['cos']:+.3f}%   rho = {dec[step]['rho_215']:+.3f}")

    # ---- the slack: what a 1% amplitude change does to the floor, per arm -----------------
    P("\n=== floor sensitivity to a global amplitude rescale, per arm ===")
    P(f"  {'arm':11s}{'d floor / d ln s':>18s}{'per +1% amp':>13s}{'as % of cos':>13s}{'legs':>6s}")
    for name in ("rc_callset", "cd4f034", "h1union"):
        s = res[name]["floor_sens_per_ln_s"]
        P(f"  {name:11s}{s:>+18.5f}{0.01*s:>+13.6f}{100*0.01*s/res[name]['inv']['cos']:>+12.2f}%"
          f"{''.join(res[name]['inv']['legs'][c][0] for c in CTX):>6s}")

    # ---- metas, for the record -------------------------------------------------------------
    m210 = json.loads((ROOT / "submissions/210-cd4f034.h5ad.meta.json").read_text())
    m222 = json.loads((ROOT / "submissions/222-h1union-g028.h5ad.meta.json").read_text())
    ps210 = m210["effects"]["cd4"]["post_add_preserve"]["post_add_global_preserve_scale"]
    ps222 = m222["effects"]["h1_union"]["merged_preserve"]["post_add_global_preserve_scale"]
    P(f"\n  metas: post_add_global_preserve_scale {ps210:.6f} -> {ps222:.6f} ({100*(ps222/ps210-1):+.3f}%);"
      f" after_all_median_norm {m210['effects']['whitening']['after_all_median_norm']:.6f} -> "
      f"{m222['effects']['whitening']['after_all_median_norm']:.6f}; alpha {m210['config']['alpha']} -> "
      f"{m222['config']['alpha']}")
    P(f"  measured delivered energy ratio sqrt(A_h1union / A_cd4f034) - 1 = "
      f"{100*(dec['h1union']['amp_ratio']-1):+.3f}%  (median-norm preserve says {100*(ps222/ps210-1):+.3f}%)")

    out = dict(arms={k: dict(per=v["per"], inv=v["inv"], floor=v["floor"], tight=v["tight"],
                             floor_sens_per_ln_s=v["floor_sens_per_ln_s"], r_mean=v["r_mean"],
                             BA=v["BA"], board=v["board"]) for k, v in res.items()},
               a14=dict(per=INV197["ours"]["a14"]["per"], inv=ka["a14"],
                        raw=INV197["ours"]["a14"]["raw"]),
               steps=dec, sumD=SUMD, floor0=FLOOR0, coeff=COEFF, pds_per_cos_215=PDS_PER_COS_215,
               preserve_scale={"cd4f034": ps210, "h1union": ps222}, peak_rss_gb=rss_gb())
    (WORK / "invert243.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "invert243.log").write_text("\n".join(LINES))
    P(f"\nwrote {WORK/'invert243.json'}  peak rss {rss_gb():.2f} GB")


if __name__ == "__main__":
    main()
