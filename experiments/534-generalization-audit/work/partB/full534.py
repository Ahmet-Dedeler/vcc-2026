"""534 Part B -- the FULL panel (300 targets, the three production donors, the panel's own 3-donor operator): 438's
full438.py restricted to the preregistered designs. Render-space energies of every design as REALISED by the ballast
(KMAX 16), the cross-panel energy anchor, and each design's re-level factor (ntclevel's rule) -> level534.json, which
comp534 --level reads for the levelled twins.

Energy route (420's r417a, 438 section 2.3): an arm's A_c/B_c on this panel = 417a's MEASURED A_c/B_c on A/B/C x (this
panel's realised energy of the arm / A/B/C's realised energy of 417a). Claims: ntclevel417's f = 1 claim from the bundle's
own NTC cells (sd 0.01). Keys on the full panel: the 410 key clip (0.3, 3); LODO-V1 at (0.5, 3); and the percentile cap
(def_lodo_key.py --hi-pct 96.64697522769474: the numpy-linear percentile of key/median over carried targets).

KNOWN ANSWERS (--ka, efp or abc; abort before writing on failure):
  KA-COMP   the PP constants == the builder's meta (efp: 432's ef_b1; abc: 410b) (<= 5e-5 rel; 438's KA-COMP)
  KA-KEY    the full-panel 410 key == the builder's artifact (efp c3e429dd, abc d7f76a22), w identical (438's KA-KEY)
  KA-CLAIMS my claims == 438's (equal)
  KA-E438   every design 438's full438.json has (h1c1|none, M0-M3, C-a): realised A/B per slot, the anchored AB, AB_level,
            rho, rho_level and the level factors == 438's full438.json / level438.json, |d| <= 1e-12 relative

    uv run python full534.py <panel> [--op PATH] [--ka]
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard534  # noqa: E402,F401
import argparse  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

import designs534 as D  # noqa: E402
import lib534 as L  # noqa: E402
from lib534 import ALPHA, GAMMA, LAM, SLOTS, ROOT, W420, E062  # noqa: E402

TOL = 1e-12
W438 = ROOT / "experiments/438-october-stack-repricing/work"
KEY_ART = {"efp": (L.R432 / "key/pert_scale_ef432_lo0.3_hi3.npz", "c3e429ddff3c39e51626980c24c4206718f43535552c956ecc247d16f29d36de"),
           "abc": (ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.3_hi3.npz",
                   "d7f76a22f62ad6991d7220ce8a65218db38b2d71dcfe3c3781ec6a8f70409e0e")}


def meta_constants(meta):
    """full438.meta_constants, verbatim."""
    e = json.loads(Path(meta).read_text())["effects"]
    return dict(before_operator_median_norm=e["whitening"]["before_operator_median_norm"],
                operator_preserve_scale=e["whitening"]["operator_preserve_scale"],
                before_all_median_norm=e["whitening"]["before_all_median_norm"],
                after_all_median_norm=e["whitening"]["after_all_median_norm"], frac_cd4=e["cd4"]["fraction"],
                frac_h1=e["h1_union"]["residual_frame"]["applied_fraction"],
                frac_partial=e["cd4_partial"]["residual_frame"]["applied_fraction"],
                s=e["cd4"]["post_add_preserve"]["post_add_global_preserve_scale"],
                defl_energy_removed=e["cd4"]["deflation"]["energy_removed_frac"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("panel")
    ap.add_argument("--op", default=None, help="3-donor operator (default: the panel's spec full_op)")
    ap.add_argument("--sd", type=float, default=None, help="re-level sd (default 0.002 on A/B/C, 0.01 elsewhere; 438)")
    ap.add_argument("--ka", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    log = lambda *x: print(f"[full534 {a.panel}]", *x, f"({time.time() - t0:.0f}s)", flush=True)
    if a.ka:
        assert a.panel in ("efp", "abc")
    ctx = L.bind(a.panel)
    targets = list(ctx.panel)
    pp_pos, nonpanel = ctx.pp_pos, ctx.nonpanel
    op_path = Path(a.op) if a.op else ctx.full_op
    assert op_path is not None, f"{a.panel}: no 3-donor operator"
    wcols = L.E064.operator_columns([L.bank(n) for n in L.MAIN], L.OFFICIAL, ctx.panel_set)
    op = L.plolo.Operator(op_path)
    op._check_axis(wcols)
    assert not (set(targets) & {L.OFFICIAL[i] for i in wcols}), "panel genes on the operator axis"
    res = dict(panel=a.panel, rules=ctx.rules, op=str(op_path.resolve().relative_to(ROOT)), op_sha256=op.sha256,
               op_axis_checked=op.axis_checked, op_axis_sha256=getattr(op, "axis_sha256", None), ka={}, designs={},
               contexts=list(ctx.contexts))
    log(f"op {op.name} {op.sha256[:12]} ({op.axis_checked}); rules {ctx.rules}; bundle {ctx.bundle}")

    # ------------------------------------------------------------------ composition (full438's lines)
    dn = L.Donors438(targets, L.MAIN)
    ncar = sum(dn.A[s].astype(int) for s in L.MAIN)
    fb = None
    if (ncar == 0).any():
        src = [L.prepared(n, targets, pp_pos) for n in L.MAIN]
        fb, _ = E062.build_baseline(src)
    res["carriers"] = {str(k): int((ncar == k).sum()) for k in np.unique(ncar)}
    pre = L.pre_normrel(dn, L.MAIN, "cap", -1.5, fallback_rows=fb)
    cd4, partial, h1, daud = L.residuals(ctx, targets)
    res["deflation"] = daud
    rows_of = lambda r: int((np.linalg.norm(np.where(np.asarray(r.gene_votes) > 0, np.asarray(r.mean_effect), 0.0), axis=1) > 0).sum())
    res["coverage"] = dict(donors={s: int(dn.A[s].sum()) for s in L.MAIN}, h1_union_rows=rows_of(h1),
                           h1_sources={n: rows_of(L.source_residual(n, targets, "gene_norm")) for n in L.UNION},
                           cd4_complete_rows=rows_of(cd4), cd4_partial_rows=rows_of(partial))
    log(f"coverage {json.dumps(res['coverage'])}; carriers {res['carriers']}")
    whit = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    FC = L.Finals(ctx, whit, cd4, partial, h1, wcols, pp_pos)
    for c in sorted({d[1] for d in D.DESIGNS}):
        FC.get(c)
    res["comp"] = FC.meta
    from lodo415 import whiten_split
    cpres, _ = whiten_split(op, pre * ALPHA, wcols)
    parts = FC.parts("h1c1")
    mine = dict(before_operator_median_norm=float(np.median(np.linalg.norm((pre * ALPHA)[:, wcols], axis=1))),
                operator_preserve_scale=cpres,
                before_all_median_norm=float(np.median(np.linalg.norm(pre * ALPHA, axis=1))),
                after_all_median_norm=float(np.median(np.linalg.norm(whit, axis=1))),
                frac_cd4=parts["fractions"]["cd4"], frac_h1=parts["fractions"]["h1_union"],
                frac_partial=parts["fractions"]["cd4_partial"], s=parts["s"], defl_energy_removed=daud["energy_removed_frac"])
    res["constants"] = mine
    if a.ka:
        ref = meta_constants(ctx.spec["meta"])
        rel = {k: mine[k] / ref[k] - 1 for k in ref}
        worst = max(abs(v) for v in rel.values())
        res["ka"]["composition_vs_meta"] = dict(meta=str(Path(ctx.spec["meta"]).relative_to(ROOT)), rel=rel, worst=worst)
        log(f"KA-COMP vs {Path(ctx.spec['meta']).name}: worst rel {worst:.2e}")
        assert worst < 5e-5, f"KA-COMP FAILED: {rel}"

    # ------------------------------------------------------------------ keys on the full panel
    key410, carried410, nd = L.shared_energy_key(targets, L.MAIN, op, wcols, cols=None)
    carried = ncar >= 2
    bases = L.lodo_bases(op, dn, L.MAIN, "cap", -1.5, pre, wcols)
    d_rows = {s: L.donor_rows(s, targets)[0] for s in L.MAIN}
    kv1 = L.lodo_key(FC.final("h1c1"), FC.parts("h1c1"), bases, d_rows, {s: dn.A[s] for s in L.MAIN}, carried, L.MAIN, pp_pos)
    del bases, d_rows
    kvv = np.nan_to_num(kv1, nan=0.0)
    hi = float(np.percentile(kvv[carried] / float(np.median(kvv[carried])), D.P_PROD))
    W = {"K3": L.weights(key410, carried410, 0.3, 3.0)[0], "V1_5_3": L.weights(kvv, carried, 0.5, 3.0)[0],
         "V1pct": L.weights(kvv, carried, 0.5, hi)[0]}
    if a.ka:
        path, sha = KEY_ART[a.panel]
        assert L.sha256(path) == sha, f"key artifact sha {path}"
        z = np.load(path)
        assert [str(t) for t in z["targets"]] == targets
        assert np.array_equal(z["w"], W["K3"]), "KA-KEY: full-panel 410 key != the builder's artifact"
        res["ka"]["key410_vs_artifact"] = f"{path.name} {sha[:12]}: w identical"
        log(f"KA-KEY: 410 key w == {path.name} ({sha[:12]})")
    res["keys"] = {k: dict(q=np.quantile(w, [0, .1, .25, .5, .75, .9, 1]).tolist(),
                           frac_floor=float(np.mean(w[carried] <= 0.5 + 1e-12)), w_max=float(np.max(w))) for k, w in W.items()}
    res["keys"]["carried"] = int(carried.sum())
    res["keys"]["pct_hi"] = hi
    np.savez(ctx.work / "keys534-full.npz", targets=np.asarray(targets), key410=np.asarray(key410), carried410=np.asarray(carried410),
             kv_lodo=kvv, carried_lodo=carried, **{f"w_{k}": w for k, w in W.items()})
    log(f"full-panel keys: carried {int(carried.sum())}; percentile cap hi {hi:.4f}")

    # ------------------------------------------------------------------ energies of every design
    P_DE = L.carrier_profiles(FC.final("h1c1"))

    def energies(name, comp, key, q, scale=None):
        bulk = FC.final(comp) if key is None else L.key_scale(FC.final(comp), W[key])[0]
        dr = L.realised_rows(bulk, P_DE, q=q, scale=scale, kmaxes=(16,), zero_cols=pp_pos)
        e = dict(design={c: L.energy(dr[c]["design"], nonpanel) for c in SLOTS},
                 realised={c: L.energy(dr[c]["realised"][16], nonpanel) for c in SLOTS})
        K = np.concatenate([dr[c]["stats"][16]["K"] for c in SLOTS])
        cosr = np.concatenate([L.lib420.rowcos(dr[c]["design"][:, nonpanel], dr[c]["realised"][16][:, nonpanel]) for c in SLOTS])
        e["ballast"] = dict(frac_capped=float(np.mean(K >= 16)), K_mean=float(np.mean(K)),
                            cos_realised_vs_design_q1=float(np.quantile(cosr, 0.01)))
        return e

    ENERGY = [d for d in D.DESIGNS if not d[0].startswith("CTRL:")]
    for name, comp, key, q, _ in ENERGY:
        res["designs"][name] = energies(name, comp, key, q)
    log(f"{len(res['designs'])} designs rendered")

    # ------------------------------------------------------------------ the anchor, rho, the level factors
    a420 = json.loads((W420 / "abc420.json").read_text())
    M0 = {c: a420["measured"]["417a"][c] for c in SLOTS}
    E0 = a420["designs"]["PP|K3+sh1.25"]["realised"]["16"]
    if a.panel.startswith("abc"):                # 438: A/B/C uses 187's claims (analyse418.Cliff) and sd 0.002
        import analyse418
        claims = {c: float(analyse418.Cliff().claim[c]) for c in SLOTS}
        a.sd = 0.002 if a.sd is None else a.sd
    else:
        nc = L.ntc_claims(ctx)
        claims = {c: float(nc[c]["claim"]) for c in SLOTS}
        a.sd = 0.01 if a.sd is None else a.sd
    res["claims"], res["anchor"] = claims, dict(M0_417a=M0, E_abc_417a_realised=E0)

    def anchored(e):
        return {c: dict(A=M0[c]["A"] * e[c]["A"] / E0[c]["A"], B=M0[c]["B"] * e[c]["B"] / E0[c]["B"]) for c in SLOTS}

    for name, e in res["designs"].items():
        e["AB"] = anchored(e["realised"])
        e["rho"] = {c: e["AB"][c]["B"] / claims[c] for c in SLOTS}
    rho_t = L.level_rho(a.sd)
    fac = {}
    for name, comp, key, q, _ in ENERGY:
        e = res["designs"][name]
        fac[name] = {c: float(np.sqrt(rho_t * claims[c] / e["AB"][c]["B"])) for c in SLOTS}
        el = energies(name, comp, key, q, scale=fac[name])
        e["AB_level"] = anchored(el["realised"])
        e["rho_level"] = {c: e["AB_level"][c]["B"] / claims[c] for c in SLOTS}
        e["ballast_level"] = el["ballast"]
        log(f"{name:28s} rho " + "/".join(f"{e['rho'][c]:.4f}" for c in SLOTS) + " -> factors "
            + "/".join(f"{fac[name][c]:.4f}" for c in SLOTS) + " -> rho after " + "/".join(f"{e['rho_level'][c]:.4f}" for c in SLOTS))
    for n in D.NAMES:                            # the controls use their unpermuted design's factors (designs534)
        fac.setdefault(n, fac[D.energy_name(n)])
    lv = dict(panel=a.panel, sd=a.sd, rho_target=rho_t, claims=claims, factors=fac, controls=D.CTRL_ENERGY,
              rule="ntclevel417: s = sqrt(rho* claim / B) with rho* = exp(2 sd Phi^-1(0.148)); B from the anchored route")
    res["level"] = dict(sd=a.sd, rho_target=rho_t)

    if a.ka:                                     # KA-CLAIMS / KA-E438 against 438's full panel (efp or abc)
        f438 = json.loads((W438 / a.panel / "full438.json").read_text())
        l438 = json.loads((W438 / a.panel / "level438.json").read_text())
        assert f438["op_sha256"] == op.sha256, "KA: not 438's efp operator"
        dcl = max(abs(claims[c] - f438["claims"][c]) for c in SLOTS)
        assert dcl == 0.0, f"KA-CLAIMS: {dcl}"
        worst = 0.0
        for name in ("h1c1|none", "h1c1|K3", "h1c1|K3+sh1.25", "h2c1|K3+sh1.25", "h2c1|K3+sh1.5", "h2c1.5|K3+sh1.75"):
            e, r = res["designs"][name], f438["designs"][name]
            for blk in ("AB", "AB_level"):
                worst = max(worst, max(abs(e[blk][c][x] / r[blk][c][x] - 1) for c in SLOTS for x in ("A", "B")))
            for blk in ("rho", "rho_level"):
                worst = max(worst, max(abs(e[blk][c] / r[blk][c] - 1) for c in SLOTS))
            worst = max(worst, max(abs(e["realised"][c][x] / r["realised"][c][x] - 1) for c in SLOTS for x in ("A", "B")))
            worst = max(worst, max(abs(fac[name][c] / l438["factors"][name][c] - 1) for c in SLOTS))
        res["ka"]["energies_levels_vs_438"] = worst
        res["ka"]["claims_vs_438"] = dcl
        log(f"KA-CLAIMS |d| {dcl:.1e}; KA-E438 (realised A/B, AB, AB_level, rho, rho_level, level factors; 6 designs) worst rel {worst:.1e}")
        assert worst <= TOL, f"KA-E438 FAILED: {worst}"

    (ctx.work / "level534.json").write_text(json.dumps(lv, indent=1))
    res["peak_rss_gib"] = E062.peak_rss_gib()
    res["seconds"] = time.time() - t0
    (ctx.work / "full534.json").write_text(json.dumps(res, indent=1, default=float))
    log(f"done; peak RSS {res['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main()
