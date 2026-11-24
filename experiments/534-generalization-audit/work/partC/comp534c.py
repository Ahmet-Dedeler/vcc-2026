"""534 Part C -- one PLOLO fold of one panel: B0..B3 (designs534c), each scored as DESIGNED and as REALISED (KMAX 16), plus
each design's re-levelled twin when a level file is given. One process per fold (062's CD4 guard reads the process's
PEAK RSS, so a second fold in the same process could trip it; Part B's abc run did).

The arithmetic is comp438's through lib534c (lib438 with WORK under work/partC and three base toggles):
  production state (B0): bind -> build_fold438 (the panel's held-out operator, axis asserted) -> Donors438 ->
    residuals -> pre_base (cap, -1.5) -> op.apply -> Finals (h1c1, h2c1) -> the fold's 410 key and its LODO-V1 key
    (lodo_bases_base on this base) -> carriers from this base's h1c1 -> betas from this base's PP|K3 energy (438)
  toggled states (B1, B2, B3): the same, with the base re-derived under the design's toggles: its donor mean, whitened
    rows, compositions, LODO key and carriers follow the toggles; the fold's calibration constants (betas, the pricer's
    anchors h1c1|K3 / h1c1|none) stay the production base's.
Rows: h1c1|none, h1c1|K3, h2c1|V1_5_3+sh1.5 (the nearest design 438 stored), B0..B3 = h2c1|V1_5_3+sh1.75 per base; on
`abc` also each toggled base's h1c1|K3 design row (418's stored toggled compositions).

KNOWN ANSWERS (all abort before anything is written, |d| <= 1e-12 with L.cmp's fields):
  KA-FOLD     my fold cache (built from scratch) == 438's fold-<held>-held.npz, bit for bit
  KA-438      production h1c1|none, h1c1|K3, h2c1|V1_5_3+sh1.5 (design, realised 16, and each levelled twin when levels
              are given) == 438's comp438-<held>-held.json: the toggles at production reproduce the unpatched lib438
  KA-546      B0 == 546o's `base q1.75` (design, realised 16) on real / efp / efr (unpatched lib438, the exact B0 design)
  KA-534B     B0 and B0:L == Part B's h2c1|V1_5_3+sh1.75 rows on efp and abc (the unpatched lib534)
  KA-UNIFORM  B2's base (exponent 0) == lib418.pre_normrel's uniform mean (mode None), bit for bit
  KA-418      on abc: B1 / B2 / B3's h1c1|K3 design rows == 418's f2:plain-1.5|K3 / e:exp0|K3 / d:g0|K3
  KA-RESUME   a second pass (the fold route's levelled pass) reproduces the first pass's unlevelled rows exactly

    uv run python comp534c.py <panel> <held> [--level work/partC/<panel>/level534c.json]
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard534  # noqa: E402,F401
import argparse  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import designs534c as D  # noqa: E402
import lib534c as L  # noqa: E402
from lib534c import ALPHA, GAMMA, LAM, SLOTS, ROOT, E062  # noqa: E402

TOL = 1e-12
W438 = ROOT / "experiments/438-october-stack-repricing/work"
W546 = ROOT / "experiments/546-october-cap-pertarget/work"
W418 = ROOT / "experiments/418-bulk-composition/work"
WPB = ROOT / "experiments/534-generalization-audit/work/partB"
B_NAME = f"{D.COMP}|{D.KEY}+sh{D.Q:g}"


def run(panel, held, level_file=None, pk=False):
    t0 = time.time()
    log = lambda *a: print(f"[comp534c {panel}/{held}]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    L.set_toggles(**D.PROD)
    ctx = L.bind(panel)
    fold = L.build_fold438(ctx, held, "held")
    scored, pp_pos, wcols = fold["scored"], fold["panel_positions"], fold["wcols"]
    truth_full, pre_cached = fold["truth"], fold["pre"]
    truth = truth_full[:, scored]
    targets, sources, op = list(fold["targets"]), list(fold["sources"]), fold["op"]
    nonpanel = ctx.nonpanel
    sfx = "pk" if pk else ""
    out_json = ctx.work / f"comp534c{sfx}-{held}.json"
    prev = json.loads(out_json.read_text()) if out_json.exists() else None
    res = dict(panel=panel, fold=held, rules=ctx.rules, n_targets=len(targets), sources=sources, audit=fold["audit"],
               op=str(fold["op_path"].relative_to(ROOT)), op_sha256=op.sha256, op_axis_checked=op.axis_checked,
               ka={}, arms={}, states={}, level_file=str(level_file) if level_file else None)
    log(f"{len(targets)} targets; op {op.name} {op.sha256[:12]} ({op.axis_checked})")
    z = np.load(W438 / panel / f"fold-{held}-held.npz", allow_pickle=True)          # KA-FOLD
    assert tuple(str(t) for t in z["targets"]) == tuple(targets), "KA-FOLD targets"
    for k in ("pre", "truth"):
        assert np.array_equal(np.asarray(z[k], np.float64), fold[k]), f"KA-FOLD {k}"
    for k in ("b_ctrl", "gate", "scored", "wcols"):
        assert np.array_equal(np.asarray(z[k]), np.asarray(fold[k])), f"KA-FOLD {k}"
    res["ka"]["fold"] = f"bit-identical to 438's {panel} fold cache"
    del z
    if panel in D.FOLD_ROUTE:                            # the fold route's claims (bundle-only; cached for foldlevel534c)
        res["claims"] = {c: v["claim"] for c, v in L.ntc_claims(ctx).items()}

    dn = L.Donors438(targets, sources)
    cd4, partial, h1, daud = L.residuals(ctx, targets)   # toggle-independent; first (062's guard reads the peak RSS)
    res["deflation"] = daud
    key410, carried410, _ = L.shared_energy_key(targets, sources, op, wcols, cols=None)
    W3 = L.weights(key410, carried410, 0.3, 3.0)[0]
    d_rows = {s: L.donor_rows(s, targets)[0] for s in sources}
    A_ = {s: dn.A[s] for s in sources}
    carried_all = sum(A_[s].astype(int) for s in sources) >= 2
    lv = json.loads(open(level_file).read()) if level_file else None
    uniform = L.pre_normrel(dn, sources, None, None, fallback_rows=pre_cached)

    def build_state(tag):
        L.set_toggles(**D.BASES[tag])
        pre = L.pre_base(dn, sources, fallback_rows=pre_cached)
        whit = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
        FC = L.Finals(ctx, whit, cd4, partial, h1, wcols, pp_pos, pinned=("h1c1", "h2c1"), size=1)
        FC.get("h1c1")
        FC.get("h2c1")
        bases = L.lodo_bases_base(op, dn, sources, pre, wcols)
        kv1 = L.lodo_key(FC.final("h1c1"), FC.parts("h1c1"), bases, d_rows, A_, carried_all, sources, pp_pos)
        del bases
        kvv = np.nan_to_num(kv1, nan=0.0)
        W = {"K3": W3, "V1_5_3": L.weights(kvv, carried_all, 0.5, 3.0)[0]}
        P_DE = L.carrier_profiles(FC.final("h1c1"))
        res["states"][tag] = dict(toggles=L.toggles(), comp=FC.meta,
                                  key_q=np.quantile(W["V1_5_3"], [0, .1, .5, .9, 1]).tolist(),
                                  key_frac_floor=float(np.mean(W["V1_5_3"][carried_all] <= 0.5 + 1e-12)),
                                  pre_equals_uniform=bool(np.array_equal(pre, uniform)),
                                  pre_rel_diff_vs_production=None)
        return dict(tag=tag, pre=pre, FC=FC, W=W, P_DE=P_DE)

    arrays = {}

    def rec(name, S, comp, key, q, bet, scale=None, realised=True):
        assert L.toggles() == D.BASES[S["tag"]], (name, L.toggles())
        bulk = S["FC"].final(comp) if key is None else L.key_scale(S["FC"].final(comp), S["W"][key])[0]
        dr = L.realised_rows(bulk, S["P_DE"], q=q, scale=scale, kmaxes=(16,) if realised else ())
        sd, ad = L.row_score({c: dr[c]["design"] for c in SLOTS}, scored, truth, bet, nonpanel)
        entry = dict(state=S["tag"], comp=comp, key=key, q=q, scale=scale, design=sd, realised={})
        arrays[f"{name}@D"] = ad
        if realised:
            sr, ar = L.row_score({c: dr[c]["realised"][16] for c in SLOTS}, scored, truth, bet, nonpanel)
            K = np.concatenate([dr[c]["stats"][16]["K"] for c in SLOTS])
            sr["frac_capped"], sr["K_mean"] = float(np.mean(K >= 16)), float(np.mean(K))
            entry["realised"]["16"] = sr
            arrays[f"{name}@R16"] = ar
        res["arms"][name] = entry
        r16 = entry["realised"].get("16")
        log(f"{name:28s} cos {sd['cos_agg']:.6f} pds@.5 {sd['pds_lam0.5']:.6f}"
            + (f" | R16 cos {r16['cos_agg']:.6f} pds@.5 {r16['pds_lam0.5']:.6f}" if r16 else ""))

    fac = (lambda n: {c: float(lv["factors"][n][c]) for c in SLOTS}) if lv else None
    # ------------------------------------------------------------------ production base (B0) and the fold's constants
    S0 = build_state("B0")
    bulk_k3 = L.key_scale(S0["FC"].final("h1c1"), W3)[0]
    dk3 = L.realised_rows(bulk_k3, S0["P_DE"], q=None, kmaxes=(), zero_cols=None)
    bet = L.betas(ctx, held, nonpanel, A_pp_k3={c: L.energy(dk3[c]["design"], nonpanel)["A"] for c in SLOTS})
    del dk3, bulk_k3
    rec("h1c1|none", S0, "h1c1", None, None, bet, realised=False)
    rec("h1c1|K3", S0, "h1c1", "K3", None, bet)
    rec(D.NEAREST_438, S0, "h2c1", "V1_5_3", 1.5, bet)
    rec("B0", S0, D.COMP, D.KEY, D.Q, bet)
    if lv:
        for n, comp, key, q in (("h1c1|none", "h1c1", None, None), ("h1c1|K3", "h1c1", "K3", None),
                                (D.NEAREST_438, "h2c1", "V1_5_3", 1.5), ("B0", D.COMP, D.KEY, D.Q)):
            if n in lv["factors"]:                       # the fold route has no factor for the design-only h1c1|none
                rec(f"{n}:L", S0, comp, key, q, bet, scale=fac(n))
    prod_pre = S0["pre"]
    W0v = S0["W"]["V1_5_3"]                              # the production base's LODO key (the pk variant)
    del S0

    # ------------------------------------------------------------------ known answers on the production base, BEFORE any
    # B1-B3 row exists (PREREGISTRATION-C section 3): the toggles at production reproduce the unpatched copies
    worst = {}

    def chk(tag, mine, ref, fields=("cos_agg", "pds", "pds_lam0.5", "pds_lam1")):
        d = L.cmp(mine, ref, fields)
        worst[tag] = d
        log(f"KA {tag:52s} |d| {d:.1e}")
        assert d <= TOL, f"KA FAILED {tag}: {d}"

    c438 = json.loads((W438 / panel / f"comp438-{held}-held.json").read_text())["arms"]
    for n in ("h1c1|none", "h1c1|K3", D.NEAREST_438):
        for tw in (("", ":L") if lv else ("",)):
            if tw and (n + tw) not in c438:              # 438 levelled only what its level file priced (efr: no h1c1|none)
                assert (n + tw) not in res["arms"], n + tw
                continue
            mine, ref = res["arms"][n + tw], c438[n + tw]
            chk(f"438 {n + tw} design", mine["design"], ref["design"])
            if ref.get("realised") and mine["realised"]:
                chk(f"438 {n + tw} realised16", mine["realised"]["16"], ref["realised"]["16"])
    if panel in ("real", "efp", "efr"):
        r546 = json.loads((W546 / f"comp546o-{panel}-{held}.json").read_text())["arms"]["base q1.75"]
        chk("546o base q1.75 == B0 design", res["arms"]["B0"]["design"], r546["design"])
        chk("546o base q1.75 == B0 realised16", res["arms"]["B0"]["realised"]["16"], r546["realised"]["16"])
    if panel in ("efp", "abc"):
        pb = json.loads((WPB / panel / f"comp534-{held}.json").read_text())["arms"]
        for tw in (("", ":L") if lv else ("",)):
            chk(f"534B {B_NAME + tw} == B0{tw} design", res["arms"]["B0" + tw]["design"], pb[B_NAME + tw]["design"])
            chk(f"534B {B_NAME + tw} == B0{tw} realised16", res["arms"]["B0" + tw]["realised"]["16"],
                pb[B_NAME + tw]["realised"]["16"])
    assert not res["states"]["B0"]["pre_equals_uniform"]
    log(f"production known answers pass ({len(worst)} checks); now the toggled bases")

    # ------------------------------------------------------------------ toggled bases: each one's own known answer first
    a418 = json.loads((W418 / f"comp418-{held}.json").read_text())["arms"] if panel == "abc" else None
    for tag in ("B1", "B2", "B3"):
        S = build_state(tag)
        res["states"][tag]["pre_rel_diff_vs_production"] = float(np.linalg.norm(S["pre"] - prod_pre) / np.linalg.norm(prod_pre))
        if tag == "B2":
            assert res["states"]["B2"]["pre_equals_uniform"], "KA-UNIFORM: B2's base (exponent 0) != the uniform mean"
            worst["B2 base == uniform mean (bit for bit)"] = 0.0
        else:
            assert not res["states"][tag]["pre_equals_uniform"], tag
        if panel == "abc":                               # KA-418: the toggle reproduces 418's stored toggled composition
            rec(f"{tag}:h1c1|K3", S, "h1c1", "K3", None, bet, realised=False)
            chk(f"418 {D.KA418[tag]} == {tag}:h1c1|K3 design", res["arms"][f"{tag}:h1c1|K3"]["design"], a418[D.KA418[tag]])
        rec(tag, S, D.COMP, D.KEY, D.Q, bet)
        if lv:
            rec(f"{tag}:L", S, D.COMP, D.KEY, D.Q, bet, scale=fac(tag))
        if pk:                                           # secondary: the toggled base with the PRODUCTION base's LODO key
            Spk = dict(S, W={"K3": W3, "V1_5_3": W0v})
            rec(f"{tag}pk", Spk, D.COMP, D.KEY, D.Q, bet)
            if lv:
                rec(f"{tag}pk:L", Spk, D.COMP, D.KEY, D.Q, bet, scale=fac(f"{tag}pk"))
        del S
    L.set_toggles(**D.PROD)
    if prev and lv and not prev.get("level_file"):      # KA-RESUME (the fold route's second pass)
        for n, e in prev["arms"].items():
            d = L.cmp(res["arms"][n]["design"], e["design"])
            if e["realised"]:
                d = max(d, L.cmp(res["arms"][n]["realised"]["16"], e["realised"]["16"]))
            worst[f"resume {n}"] = d
            assert d == 0.0, f"KA-RESUME {n}: {d}"
        log(f"KA-RESUME: the levelled pass reproduces the first pass's {len(prev['arms'])} rows exactly")
    res["ka"]["rows"] = worst
    res["ka"]["worst"] = max(worst.values())
    log(f"ALL KNOWN ANSWERS PASS on {panel}/{held}: {len(worst)} checks, worst |d| {res['ka']['worst']:.1e}")

    if pk:                                               # KA-PRIMARY: every row the primary run has is reproduced exactly
        prim = json.loads((ctx.work / f"comp534c-{held}.json").read_text())["arms"]
        for n, e in res["arms"].items():
            if n in prim and bool(prim[n]["realised"]) == bool(e["realised"]):
                d = L.cmp(e["design"], prim[n]["design"])
                if e["realised"]:
                    d = max(d, L.cmp(e["realised"]["16"], prim[n]["realised"]["16"]))
                assert d == 0.0, f"KA-PRIMARY {n}: {d}"
                res["ka"].setdefault("primary", {})[n] = d
        log(f"KA-PRIMARY: {len(res['ka']['primary'])} rows == the primary run's, bit for bit")
    res["peak_rss_gib"] = E062.peak_rss_gib()
    res["seconds"] = time.time() - t0
    out_json.write_text(json.dumps(res, indent=1, default=float))
    np.savez_compressed(ctx.work / f"arrays534c{sfx}-{held}.npz", names=np.asarray(list(arrays)), targets=np.asarray(targets),
                        d=np.asarray((truth * truth).sum(1)), **{f"{n}__{k}": v for n, a in arrays.items() for k, v in a.items()})
    log(f"done; {len(res['arms'])} arms; peak RSS {res['peak_rss_gib']:.2f} GiB -> {out_json.name}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("panel")
    ap.add_argument("fold")
    ap.add_argument("--level", default=None)
    ap.add_argument("--pk", action="store_true", help="secondary: add B1pk..B3pk (toggled base, production key)")
    a = ap.parse_args()
    run(a.panel, a.fold, a.level, a.pk)
