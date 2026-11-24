"""534 Part C -- the FULL panel (the panel's 3-donor operator) for abc / real / efp: each B design's render-space energy on
its own base (its own donor mean, compositions, full-panel LODO key and carriers), the anchored route, and its ntclevel
re-level factor -> level534c.json for comp534c --level. 438's full438 arithmetic, as Part B's full534.

KNOWN ANSWERS (abort before writing):
  KA-COMP   production PP constants == the builder's meta (abc 410b, real 432 b0, efp 432 ef_b1), <= 5e-5 rel (438's value)
  KA-KEY    the full-panel 410 key == the builder's artifact (abc d7f76a22, real 4f1f0209, efp c3e429dd), w identical
  KA-E438   production h1c1|none, h1c1|K3, h2c1|V1_5_3+sh1.5: realised A/B, anchored AB / AB_level, rho / rho_level and
            the level factor == 438's full438.json / level438.json (relative |d| <= 1e-12)
  KA-534B   B0 == Part B's h2c1|V1_5_3+sh1.75 (efp, abc): the same quantities and factor
  KA-UNIFORM  B2's base == the uniform mean, bit for bit

    uv run python full534c.py <panel>
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

import designs534c as D  # noqa: E402
import lib534c as L  # noqa: E402
from lib534c import ALPHA, GAMMA, LAM, SLOTS, ROOT, W420, E062  # noqa: E402

TOL = 1e-12
W438 = ROOT / "experiments/438-october-stack-repricing/work"
WPB = ROOT / "experiments/534-generalization-audit/work/partB"
KEY_ART = {"abc": (ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.3_hi3.npz",
                   "d7f76a22f62ad6991d7220ce8a65218db38b2d71dcfe3c3781ec6a8f70409e0e"),
           "real": (L.R432 / "key/pert_scale_def432_lo0.3_hi3.npz", "4f1f02090c46baaa5e8e257e42b4047db8d05e4f914d0491facec85d6aed692a"),
           "efp": (L.R432 / "key/pert_scale_ef432_lo0.3_hi3.npz", "c3e429ddff3c39e51626980c24c4206718f43535552c956ecc247d16f29d36de")}
B_NAME = f"{D.COMP}|{D.KEY}+sh{D.Q:g}"


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
    ap.add_argument("--pk", action="store_true", help="secondary: add B1pk..B3pk (toggled base, production key)")
    a = ap.parse_args()
    sfx = "pk" if a.pk else ""
    assert a.panel in ("abc", "real", "efp"), "efr has no 3-donor operator: the fold route (foldlevel534c)"
    t0 = time.time()
    log = lambda *x: print(f"[full534c {a.panel}]", *x, f"({time.time() - t0:.0f}s)", flush=True)
    L.set_toggles(**D.PROD)
    ctx = L.bind(a.panel)
    targets = list(ctx.panel)
    pp_pos, nonpanel = ctx.pp_pos, ctx.nonpanel
    wcols = L.E064.operator_columns([L.bank(n) for n in L.MAIN], L.OFFICIAL, ctx.panel_set)
    op = L.plolo.Operator(ctx.full_op)
    op._check_axis(wcols)
    assert not (set(targets) & {L.OFFICIAL[i] for i in wcols}), "panel genes on the operator axis"
    res = dict(panel=a.panel, rules=ctx.rules, op=str(Path(ctx.full_op).resolve().relative_to(ROOT)), op_sha256=op.sha256,
               op_axis_checked=op.axis_checked, ka={}, designs={}, states={})
    log(f"op {op.name} {op.sha256[:12]} ({op.axis_checked}); rules {ctx.rules}")
    dn = L.Donors438(targets, L.MAIN)
    ncar = sum(dn.A[s].astype(int) for s in L.MAIN)
    fb = None
    if (ncar == 0).any():
        src = [L.prepared(n, targets, pp_pos) for n in L.MAIN]
        fb, _ = E062.build_baseline(src)
    cd4, partial, h1, daud = L.residuals(ctx, targets)
    key410, carried410, _ = L.shared_energy_key(targets, L.MAIN, op, wcols, cols=None)
    W3 = L.weights(key410, carried410, 0.3, 3.0)[0]
    path, sha = KEY_ART[a.panel]                                              # KA-KEY
    assert L.sha256(path) == sha, f"key artifact sha {path}"
    zk = np.load(path)
    assert [str(t) for t in zk["targets"]] == targets and np.array_equal(zk["w"], W3), "KA-KEY"
    res["ka"]["key410_vs_artifact"] = f"{path.name} {sha[:12]}: w identical"
    d_rows = {s: L.donor_rows(s, targets)[0] for s in L.MAIN}
    carried = ncar >= 2
    uniform = L.pre_normrel(dn, L.MAIN, None, None, fallback_rows=fb)

    def build_state(tag):
        L.set_toggles(**D.BASES[tag])
        pre = L.pre_base(dn, L.MAIN, fallback_rows=fb)
        whit = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
        FC = L.Finals(ctx, whit, cd4, partial, h1, wcols, pp_pos, pinned=("h1c1", "h2c1"), size=1)
        FC.get("h1c1")
        FC.get("h2c1")
        bases = L.lodo_bases_base(op, dn, L.MAIN, pre, wcols)
        kv1 = L.lodo_key(FC.final("h1c1"), FC.parts("h1c1"), bases, d_rows, {s: dn.A[s] for s in L.MAIN}, carried, L.MAIN, pp_pos)
        del bases
        W = {"K3": W3, "V1_5_3": L.weights(np.nan_to_num(kv1, nan=0.0), carried, 0.5, 3.0)[0]}
        res["states"][tag] = dict(toggles=L.toggles(), comp=FC.meta, pre_equals_uniform=bool(np.array_equal(pre, uniform)))
        return dict(tag=tag, pre=pre, whit=whit, FC=FC, W=W, P_DE=L.carrier_profiles(FC.final("h1c1")))

    def energies(S, comp, key, q, scale=None):
        assert L.toggles() == D.BASES[S["tag"]]
        bulk = S["FC"].final(comp) if key is None else L.key_scale(S["FC"].final(comp), S["W"][key])[0]
        dr = L.realised_rows(bulk, S["P_DE"], q=q, scale=scale, kmaxes=(16,), zero_cols=pp_pos)
        e = dict(design={c: L.energy(dr[c]["design"], nonpanel) for c in SLOTS},
                 realised={c: L.energy(dr[c]["realised"][16], nonpanel) for c in SLOTS})
        K = np.concatenate([dr[c]["stats"][16]["K"] for c in SLOTS])
        e["ballast"] = dict(frac_capped=float(np.mean(K >= 16)), K_mean=float(np.mean(K)))
        return e

    a420 = json.loads((W420 / "abc420.json").read_text())
    M0 = {c: a420["measured"]["417a"][c] for c in SLOTS}
    E0 = a420["designs"]["PP|K3+sh1.25"]["realised"]["16"]
    if a.panel == "abc":                                    # 438: A/B/C uses 187's claims (analyse418.Cliff), sd 0.002
        import analyse418
        claims, sd = {c: float(analyse418.Cliff().claim[c]) for c in SLOTS}, 0.002
    else:
        nc = L.ntc_claims(ctx)
        claims, sd = {c: float(nc[c]["claim"]) for c in SLOTS}, 0.01
    rho_t = L.level_rho(sd)
    res["claims"], res["level"] = claims, dict(sd=sd, rho_target=rho_t)
    anchored = lambda e: {c: dict(A=M0[c]["A"] * e[c]["A"] / E0[c]["A"], B=M0[c]["B"] * e[c]["B"] / E0[c]["B"]) for c in SLOTS}
    fac = {}

    def design(S, name, comp, key, q):
        e = energies(S, comp, key, q)
        e["AB"] = anchored(e["realised"])
        e["rho"] = {c: e["AB"][c]["B"] / claims[c] for c in SLOTS}
        fac[name] = {c: float(np.sqrt(rho_t * claims[c] / e["AB"][c]["B"])) for c in SLOTS}
        el = energies(S, comp, key, q, scale=fac[name])
        e["AB_level"] = anchored(el["realised"])
        e["rho_level"] = {c: e["AB_level"][c]["B"] / claims[c] for c in SLOTS}
        res["designs"][name] = e
        log(f"{name:24s} rho " + "/".join(f"{e['rho'][c]:.4f}" for c in SLOTS) + " -> factors "
            + "/".join(f"{fac[name][c]:.4f}" for c in SLOTS) + " -> rho after " + "/".join(f"{e['rho_level'][c]:.4f}" for c in SLOTS))

    S0 = build_state("B0")                                   # production base: KA designs and B0
    parts, whit = S0["FC"].parts("h1c1"), S0["whit"]
    from lodo415 import whiten_split
    cpres, _ = whiten_split(op, S0["pre"] * ALPHA, wcols)
    mine = dict(before_operator_median_norm=float(np.median(np.linalg.norm((S0["pre"] * ALPHA)[:, wcols], axis=1))),
                operator_preserve_scale=cpres, before_all_median_norm=float(np.median(np.linalg.norm(S0["pre"] * ALPHA, axis=1))),
                after_all_median_norm=float(np.median(np.linalg.norm(whit, axis=1))),
                frac_cd4=parts["fractions"]["cd4"], frac_h1=parts["fractions"]["h1_union"],
                frac_partial=parts["fractions"]["cd4_partial"], s=parts["s"], defl_energy_removed=daud["energy_removed_frac"])
    ref = meta_constants(ctx.spec["meta"])                    # KA-COMP
    rel = {k: mine[k] / ref[k] - 1 for k in ref}
    res["ka"]["composition_vs_meta"] = dict(meta=str(Path(ctx.spec["meta"]).relative_to(ROOT)), worst=max(abs(v) for v in rel.values()))
    log(f"KA-COMP vs {Path(ctx.spec['meta']).name}: worst rel {res['ka']['composition_vs_meta']['worst']:.2e}")
    assert res["ka"]["composition_vs_meta"]["worst"] < 5e-5, f"KA-COMP FAILED: {rel}"
    for name, comp, key, q in (("h1c1|none", "h1c1", None, None), ("h1c1|K3", "h1c1", "K3", None),
                               (D.NEAREST_438, "h2c1", "V1_5_3", 1.5), ("B0", D.COMP, D.KEY, D.Q)):
        design(S0, name, comp, key, q)
    W0v = S0["W"]["V1_5_3"]
    del S0

    # ------------------------------------------------------------------ known answers on the production base, BEFORE any
    # B1-B3 energy exists (PREREGISTRATION-C section 3)
    f438 = json.loads((W438 / a.panel / "full438.json").read_text())
    l438 = json.loads((W438 / a.panel / "level438.json").read_text())
    assert f438["op_sha256"] == op.sha256 and max(abs(claims[c] - f438["claims"][c]) for c in SLOTS) == 0.0, "KA-CLAIMS"

    def cmpE(e, r, f, fr):
        w = max(abs(e[b][c][x] / r[b][c][x] - 1) for b in ("AB", "AB_level") for c in SLOTS for x in ("A", "B"))
        w = max(w, max(abs(e[b][c] / r[b][c] - 1) for b in ("rho", "rho_level") for c in SLOTS))
        w = max(w, max(abs(e["realised"][c][x] / r["realised"][c][x] - 1) for c in SLOTS for x in ("A", "B")))
        return max(w, max(abs(f[c] / fr[c] - 1) for c in SLOTS))

    worst = {n: cmpE(res["designs"][n], f438["designs"][n], fac[n], l438["factors"][n]) for n in ("h1c1|none", "h1c1|K3", D.NEAREST_438)}
    if a.panel in ("efp", "abc"):
        pbf = json.loads((WPB / a.panel / "full534.json").read_text())
        pbl = json.loads((WPB / a.panel / "level534.json").read_text())
        worst["534B B0"] = cmpE(res["designs"]["B0"], pbf["designs"][B_NAME], fac["B0"], pbl["factors"][B_NAME])
    for k, v in worst.items():
        log(f"KA-E {k:24s} worst rel {v:.1e}")
        assert v <= TOL, f"KA-E FAILED {k}: {v}"
    for tag in ("B1", "B2", "B3"):                           # the toggled bases
        S = build_state(tag)
        if tag == "B2":
            assert res["states"]["B2"]["pre_equals_uniform"], "KA-UNIFORM"
        design(S, tag, D.COMP, D.KEY, D.Q)
        if a.pk:
            design(dict(S, W={"K3": W3, "V1_5_3": W0v}), f"{tag}pk", D.COMP, D.KEY, D.Q)
        del S
    L.set_toggles(**D.PROD)
    res["ka"]["energies"] = worst
    if a.pk:                                                 # KA-PRIMARY: the primary run's energies and factors, exactly
        pf = json.loads((ctx.work / "full534c.json").read_text())["designs"]
        pl = json.loads((ctx.work / "level534c.json").read_text())["factors"]
        wp = {n: cmpE(res["designs"][n], pf[n], fac[n], pl[n]) for n in pf}
        assert max(wp.values()) == 0.0, f"KA-PRIMARY {wp}"
        res["ka"]["primary"] = wp
        log(f"KA-PRIMARY: {len(wp)} designs' energies and factors == the primary run's")
    (ctx.work / f"level534c{sfx}.json").write_text(json.dumps(dict(panel=a.panel, route="full", sd=sd, rho_target=rho_t, claims=claims,
                                                               factors=fac), indent=1))
    res["peak_rss_gib"] = E062.peak_rss_gib()
    res["seconds"] = time.time() - t0
    (ctx.work / f"full534c{sfx}.json").write_text(json.dumps(res, indent=1, default=float))
    log(f"done; peak RSS {res['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main()
