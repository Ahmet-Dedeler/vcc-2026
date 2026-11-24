"""438 on the FULL panel (300 targets, the three production donors, the panel's 3-donor operator): render-space energies
of every grid design as DESIGNED and as REALISED by the ballast (KMAX 16), the cross-panel energy anchor, the re-level
factors each design would get on the day, and the known answers that make those numbers usable.

The ENERGY ROUTE (420's r417a, extended across panels): an arm's A_c/B_c on this panel =
    417a's MEASURED A_c/B_c on A/B/C  x  (this panel's realised energy of the arm / A/B/C's realised energy of 417a)
so rho_c = B_c / claim_c is predicted with no artifact of this panel. Claims: 187's on A/B/C; on a D/E/F bundle the
ntclevel417 claim with f = 1 (Cp46 x sum w from the bundle's own controls) -- what the day's re-level will use.

KNOWN ANSWERS (abort on failure):
  KA-PP     A/B/C: the PP constants reproduce the 410b build meta (<= 1e-5 rel; 418's KA-PP) and pre == sens415.Donors.pre
  KA-COMP   D/E/F: the PP constants (whitening medians, operator preserve, CD4/H1/partial fractions, post-add preserve,
            pinned-deflation energy) reproduce the builder's own meta on the same bundle (432's b0 / ef_b1), <= 1e-5 rel
  KA-KEY    the full-panel 410 key == the artifact the builder used on that panel (A/B/C d7f76a22, real 4f1f0209,
            efp c3e429dd): identical w; A/B/C's LODO-V1 == 415pp's key vector
  KA-E420   A/B/C design + realised energies of the 420 grid rows == abc420.json (1e-12 rel)
VALIDATION (reported, not gated): predicted rho vs 432's MEASURED rho (ntclevel) on real's b0 (417a analogue) and b1
(the stack), and efp's ef_b1 (the stack).

    uv run python full438.py <panel> [--op PATH] [--tag T] [--sd SD]      -> work/<tag>/full438.json, level438.json
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

import lib438 as L
from lib438 import ALPHA, GAMMA, LAM, SLOTS, ROOT, W420, E062, sens415

REF438 = "h1c1|K3+sh1.25"
ABC420_ROWS = {"PP|K3": "h1c1|K3", "PP|K3+sh1.25": "h1c1|K3+sh1.25", "PP|K3+sh1.5": "h1c1|K3+sh1.5",
               "stack|K3+sh1.25": "h2c1.5|K3+sh1.25", "stack|K3+sh1.5": "h2c1.5|K3+sh1.5",
               "stack|V1_5_3+sh1.5": "h2c1.5|V1_5_3+sh1.5"}
KEY_ART = {"abc": (ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.3_hi3.npz",
                   "d7f76a22f62ad6991d7220ce8a65218db38b2d71dcfe3c3781ec6a8f70409e0e"),
           "real": (L.R432 / "key/pert_scale_def432_lo0.3_hi3.npz",
                    "4f1f02090c46baaa5e8e257e42b4047db8d05e4f914d0491facec85d6aed692a"),
           "efp": (L.R432 / "key/pert_scale_ef432_lo0.3_hi3.npz",
                   "c3e429ddff3c39e51626980c24c4206718f43535552c956ecc247d16f29d36de")}


def meta_constants(meta):
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
    ap.add_argument("--op", default=None, help="3-donor operator (default: the panel's registry entry)")
    ap.add_argument("--tag", default=None)
    ap.add_argument("--sd", type=float, default=None, help="re-level sd (default 0.002 on A/B/C, 0.01 elsewhere)")
    ap.add_argument("--no-level", action="store_true")
    ap.add_argument("--resume", action="store_true", help="add only the grid designs a finished run does not have")
    a = ap.parse_args()
    t0 = time.time()
    log = lambda *x: print(f"[full438 {a.tag or a.panel}]", *x, f"({time.time() - t0:.0f}s)", flush=True)
    ctx = L.bind(a.panel, tag=a.tag)
    targets = list(ctx.panel)
    pp_pos, nonpanel = ctx.pp_pos, ctx.nonpanel
    op_path = Path(a.op) if a.op else ctx.full_op
    assert op_path is not None, f"{a.panel}: no 3-donor operator registered; pass --op (a same-axis surrogate)"
    surrogate = a.op is not None and Path(a.op).resolve() != (ctx.full_op.resolve() if ctx.full_op else None)
    wcols = L.E064.operator_columns([L.bank(n) for n in L.MAIN], L.OFFICIAL, ctx.panel_set)
    op = L.plolo.Operator(op_path)
    op._check_axis(wcols)
    assert not (set(targets) & {L.OFFICIAL[i] for i in wcols}), "panel genes on the operator axis"
    res = dict(panel=a.panel, tag=a.tag or a.panel, rules=ctx.rules, op=str(op_path.relative_to(ROOT)),
               op_sha256=op.sha256, surrogate=surrogate, ka={}, designs={}, contexts=list(ctx.contexts))
    log(f"op {op.name} {op.sha256[:12]} ({op.axis_checked}); surrogate={surrogate}; rules {ctx.rules}")

    # ------------------------------------------------------------------ composition
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
    log(f"coverage {json.dumps(res['coverage'])}")
    whit = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    FC = L.Finals(ctx, whit, cd4, partial, h1, wcols, pp_pos)
    need = sorted({L.parse(n)[0] for n, *_ in L.grid()})
    for c in need:
        FC.get(c)
    res["comp"] = FC.meta
    parts = {"h1c1": FC.parts("h1c1")}
    from lodo415 import whiten_split
    cpres, _ = whiten_split(op, pre * ALPHA, wcols)
    mine = dict(before_operator_median_norm=float(np.median(np.linalg.norm((pre * ALPHA)[:, wcols], axis=1))),
                operator_preserve_scale=cpres,
                before_all_median_norm=float(np.median(np.linalg.norm(pre * ALPHA, axis=1))),
                after_all_median_norm=float(np.median(np.linalg.norm(whit, axis=1))),
                frac_cd4=parts["h1c1"]["fractions"]["cd4"], frac_h1=parts["h1c1"]["fractions"]["h1_union"],
                frac_partial=parts["h1c1"]["fractions"]["cd4_partial"], s=parts["h1c1"]["s"],
                defl_energy_removed=daud["energy_removed_frac"])
    res["constants"] = mine
    meta = ctx.spec.get("meta")
    if meta and not surrogate:
        ref = meta_constants(meta)
        rel = {k: mine[k] / ref[k] - 1 for k in ref}
        worst = max(abs(v) for v in rel.values())
        res["ka"]["composition_vs_meta"] = dict(meta=str(Path(meta).relative_to(ROOT)), rel=rel, worst=worst)
        log(f"KA-COMP vs {Path(meta).name}: worst rel {worst:.2e} " + json.dumps({k: f"{v:+.1e}" for k, v in rel.items()}))
        # tolerance: A/B/C's own KA-PP sits at 8.87e-6 (415/418: the builder's pre differs from this reconstruction at
        # ~1e-6 and whitening amplifies it); efp reads 8.9e-6, real 2.6e-5. A construction error (wrong sizing rule,
        # basis, dose) moves these constants by >= 1e-3.
        assert worst < 5e-5, f"KA-COMP FAILED: {rel}"
    if a.panel == "abc" and not surrogate:
        ref = sens415.Donors(targets, L.MAIN)
        assert np.array_equal(pre, ref.pre(L.MAIN, True)), "KA-PP pre != sens415.Donors.pre"
        log("KA-PP: pre == sens415.Donors.pre; constants within 1e-5 of the 410b meta")

    # ------------------------------------------------------------------ keys on the full panel
    key410, carried410, nd = L.shared_energy_key(targets, L.MAIN, op, wcols, cols=None)
    carried = ncar >= 2
    bases = L.lodo_bases(op, dn, L.MAIN, "cap", -1.5, pre, wcols)
    d_rows = {s: L.donor_rows(s, targets)[0] for s in L.MAIN}
    kv1 = L.lodo_key(FC.final("h1c1"), FC.parts("h1c1"), bases, d_rows, {s: dn.A[s] for s in L.MAIN}, carried, L.MAIN, pp_pos)
    W = {"K3": L.weights(key410, carried410, 0.3, 3.0)[0], "K5_3": L.weights(key410, carried410, 0.5, 3.0)[0],
         "V1": L.weights(np.nan_to_num(kv1, nan=0.0), carried, 0.3, 3.0)[0],
         "V1_5_3": L.weights(np.nan_to_num(kv1, nan=0.0), carried, 0.5, 3.0)[0]}
    if a.panel in KEY_ART and not surrogate:
        path, sha = KEY_ART[a.panel]
        assert L.sha256(path) == sha, f"key artifact sha {path}"
        z = np.load(path)
        assert [str(t) for t in z["targets"]] == targets
        assert np.array_equal(z["w"], W["K3"]), "KA-KEY: full-panel 410 key != the builder's artifact"
        res["ka"]["key410_vs_artifact"] = f"{path.name} {sha[:12]}: w identical"
        log(f"KA-KEY: 410 key w == {path.name} ({sha[:12]})")
    if a.panel == "abc" and not surrogate:
        z = np.load(ROOT / "experiments/415-lodo-key/work/pert_scale_415pp_lo0.3_hi3.npz")
        dk = float(np.abs(np.nan_to_num(kv1, nan=0.0) - z["key"]).max())
        dw = float(np.abs(W["V1"] - z["w"]).max())
        assert dk < 1e-12 and dw < 1e-12, f"KA-KEY V1 vs 415pp: key {dk} w {dw}"
        res["ka"]["v1_vs_415pp"] = dict(key=dk, w=dw)
        log(f"KA-KEY: LODO-V1 == 415pp key (max |d| key {dk:.1e}, w {dw:.1e}; 415 used sens415's LOO path)")
    res["keys"] = {k: dict(q=np.quantile(w, [0, .1, .25, .5, .75, .9, 1]).tolist(),
                           frac_floor=float(np.mean(w[carried] <= L.KEYS[k][1] + 1e-12)),
                           frac_cap=float(np.mean(w[carried] >= 3.0 - 1e-12))) for k, w in W.items()}
    res["keys"]["carried"] = int(carried.sum())

    # ------------------------------------------------------------------ energies of every design
    P_DE = L.carrier_profiles(FC.final("h1c1"))
    bulks = {}

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

    old = None
    if a.resume and (ctx.work / "full438.json").exists():
        old = json.loads((ctx.work / "full438.json").read_text())
        assert old["op_sha256"] == op.sha256 and old["rules"] == ctx.rules
        e = energies(REF438, *L.parse(REF438))
        dd = max(abs(e["realised"][c][x] / old["designs"][REF438]["realised"][c][x] - 1) for c in SLOTS for x in ("A", "B"))
        assert dd == 0.0, f"RESUME: re-built setup does not reproduce the stored 417a-analogue energy: {dd}"
        res["designs"] = old["designs"]
        res["ka"] = old.get("ka", {})
        log(f"RESUME: {len(old['designs'])} designs loaded; the re-built setup reproduces the 417a analogue exactly")
    for name, comp, key, q in L.grid():
        if name not in res["designs"]:
            res["designs"][name] = energies(name, comp, key, q)
    log(f"{len(res['designs'])} designs rendered")
    if a.panel == "abc" and not surrogate:
        a420 = json.loads((W420 / "abc420.json").read_text())["designs"]
        worst = 0.0
        for r420, mine_n in ABC420_ROWS.items():
            for route, k in (("design", None), ("realised", "16")):
                ref = a420[r420][route] if k is None else a420[r420][route][k]
                got = res["designs"][mine_n][route]
                worst = max(worst, max(abs(got[c][x] / ref[c][x] - 1) for c in SLOTS for x in ("A", "B")))
        res["ka"]["energies_vs_abc420"] = worst
        log(f"KA-E420: design + realised energies of {len(ABC420_ROWS)} rows == abc420.json, worst rel {worst:.1e}")
        assert worst < 1e-12, "KA-E420 FAILED"

    # ------------------------------------------------------------------ the anchor, rho, validation
    a420 = json.loads((W420 / "abc420.json").read_text())
    M0 = {c: a420["measured"]["417a"][c] for c in SLOTS}
    E0 = a420["designs"]["PP|K3+sh1.25"]["realised"]["16"]
    if a.panel.startswith("abc"):
        import analyse418
        cl = analyse418.Cliff()
        claims = {c: float(cl.claim[c]) for c in SLOTS}
        sd = 0.002 if a.sd is None else a.sd
    else:
        nc = L.ntc_claims(ctx)
        claims = {c: float(nc[c]["claim"]) for c in SLOTS}
        sd = 0.01 if a.sd is None else a.sd
    res["claims"], res["anchor"] = claims, dict(M0_417a=M0, E_abc_417a_realised=E0)

    def anchored(e):
        return {c: dict(A=M0[c]["A"] * e[c]["A"] / E0[c]["A"], B=M0[c]["B"] * e[c]["B"] / E0[c]["B"]) for c in SLOTS}

    for name, e in res["designs"].items():
        e["AB"] = anchored(e["realised"])
        e["rho"] = {c: e["AB"][c]["B"] / claims[c] for c in SLOTS}
    ref = res["designs"]["h1c1|K3+sh1.25"]
    log("417a analogue rho " + "/".join(f"{ref['rho'][c]:.4f}" for c in SLOTS) + " | stack rho "
        + "/".join(f"{res['designs']['h2c1.5|K3+sh1.25']['rho'][c]:.4f}" for c in SLOTS))
    for tag, rel in (ctx.spec.get("measured") or {}).items():
        m = json.loads((L.R432 / rel).read_text())["per_context"]
        design = "h1c1|K3+sh1.25" if tag == "b0" else "h2c1.5|K3+sh1.25"
        pred = res["designs"][design]
        v = {}
        for s, cx in zip(SLOTS, ctx.contexts):
            v[cx] = dict(measured_B=m[cx]["B"], predicted_B=pred["AB"][s]["B"], measured_A=m[cx]["A"],
                         predicted_A=pred["AB"][s]["A"], measured_rho=m[cx]["rho"], predicted_rho=pred["rho"][s],
                         err_B=pred["AB"][s]["B"] / m[cx]["B"] - 1, err_A=pred["AB"][s]["A"] / m[cx]["A"] - 1,
                         claim_measured=m[cx]["claim"], claim_mine=claims[s])
        res.setdefault("validation", {})[tag] = dict(design=design, per_context=v)
        log(f"VALIDATION {tag} ({design}): rho measured " + "/".join(f"{v[cx]['measured_rho']:.4f}" for cx in ctx.contexts)
            + " predicted " + "/".join(f"{v[cx]['predicted_rho']:.4f}" for cx in ctx.contexts)
            + " | err B " + "/".join(f"{100 * v[cx]['err_B']:+.2f}%" for cx in ctx.contexts)
            + " err A " + "/".join(f"{100 * v[cx]['err_A']:+.2f}%" for cx in ctx.contexts))

    # ------------------------------------------------------------------ re-level factors (ntclevel's rule) + their energies
    if not a.no_level:
        rho_t = L.level_rho(sd)
        fac = {}
        for name, comp, key, q in L.grid():
            e = res["designs"][name]
            fac[name] = {c: float(np.sqrt(rho_t * claims[c] / e["AB"][c]["B"])) for c in SLOTS}
            if "AB_level" in e and old is not None:
                continue
            el = energies(name, comp, key, q, scale=fac[name])
            e["AB_level"] = anchored(el["realised"])
            e["rho_level"] = {c: e["AB_level"][c]["B"] / claims[c] for c in SLOTS}
            e["ballast_level"] = el["ballast"]
            if name in ("h1c1|K3+sh1.25", "h2c1.5|K3+sh1.25", "h2c1|K3+sh1.25", "h1c1.5|K3+sh1.25", "h2c1.5|K3+sh1.5"):
                log(f"level {name}: factors " + "/".join(f"{fac[name][c]:.4f}" for c in SLOTS) + " -> rho after "
                    + "/".join(f"{e['rho_level'][c]:.4f}" for c in SLOTS) + f" (target {rho_t:.4f})")
        lv = dict(panel=a.panel, tag=a.tag or a.panel, sd=sd, rho_target=rho_t, claims=claims, factors=fac,
                  rule="ntclevel417: s = sqrt(rho* claim / B) with rho* = exp(2 sd Phi^-1(0.148)); B from the anchored route")
        (ctx.work / "level438.json").write_text(json.dumps(lv, indent=1))
        res["level"] = dict(sd=sd, rho_target=rho_t)
    res["peak_rss_gib"] = E062.peak_rss_gib()
    (ctx.work / "full438.json").write_text(json.dumps(res, indent=1, default=float))
    log(f"done; peak RSS {res['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main()
