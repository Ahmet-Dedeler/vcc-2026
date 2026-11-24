"""551 (C10b) -- comp534.py with 551's designs. Additions, nothing else changed:
  * per-target designs (S_P2 / P2g on the M5 or M6 base, and their q-permuted controls): the fold's own exponent vector
    (S_P2: tercile of the unkeyed h2c1 row's non-panel norm on the fold's targets; P2g: carried-only tercile of the base
    key's FOLD weights, floor LOW), rendered by lib551.realised_rows_pt; levelled twins at their own full551 factors;
  * KA-546o (efp): V1p546 (P 96.647) at q 1.75 / + P2g / + S_P2 == 546o's `p96.6` / `p96.6 + P2g` / `p96.6 + S_P2`
    (design and realised 16, |d| <= 1e-12; K_mean too);
  * KA-493h (efp): V1_5_3 + P2g / + S_P2 == 493h's `P2g_carried` / `S_P2` (design and realised 16);
  * KA-534 (efp): every 534 design and its levelled twin == 534's efp comp534-<fold>.json, |d| <= 1e-12;
  * KA-REPRO (an R6 panel, INSIDE the C10b access, before any new design's number is read or written): every 534 design
    and its levelled twin == C10a's stored comp534-<fold>.json on that panel, |d| <= 1e-12 (the re-fitted operators and
    rebuilt fold caches are 534's); plus the fold cache itself == 534's cached fold, bit for bit.
Outputs: work/<panel>/comp551-<fold>.json and arrays551-<fold>.npz. comp534's docstring follows.

534 Part B -- one PLOLO fold of one panel: 438's comp438.py restricted to the preregistered designs (designs534), each
scored as DESIGNED and as REALISED by the ballast (KMAX 16), plus each design's re-levelled twin (factors from full534).

The arithmetic is comp438's, line for line, through lib534 (lib438 with its caches pointed here):
  bind(panel) -> build_fold438 (targets = panel n held coverage; the fold's HELD-OUT operator, axis asserted) ->
  Donors438 -> residuals (the panel's CD4 banks, D/E/F rules; H1 union with the panel's columns masked) ->
  pre_normrel -> op.apply -> Finals -> keys re-fitted on the fold (410 key; LODO-V1 on h1c1; the percentile cap on the
  fold's own LODO key) -> key_scale -> sharpen -> render -> ballast realise (KMAX 16) -> row_score (betas at lam 0.5).

KNOWN ANSWERS (--ka: efp, or abc for the supplementary A/B/C column (438 and 493 rows only); abort before any
file is written if one fails, |d| <= 1e-12):
  KA-FOLD   my fold cache (targets, pre, truth, b_ctrl, gate, scored, wcols) == 438's efp fold-<held>-held.npz, bit for bit
  KA-438    h1c1|none (design), M0 h1c1|K3, M1, M2, M3 and C-a h2c1.5|K3+sh1.75 == comp438-<held>-held.json (design and
            realised 16: cos_agg, pds, pds lam 0.5/1, render A/B per frame, frac_capped), and their levelled twins
  KA-493    M4 h2c1|K3+sh1.75 == 493's `base q1.75` (lib438 path, K3, h2c1, q 1.75; design and realised 16)
  KA-546    M5 h2c1|V1_5_3+sh1.75 == 546o's `base q1.75` (design and realised 16)
  KA-522    M6 h2c1|V1pct+sh1.75 == 522's `pctPROD` row (cos_agg, pds, pds lam 0.5, render A/B, K_mean; hi == w_max)

    uv run python comp551.py <panel> <held>[,<held>...] [--level work/<panel>/level551.json] [--ka | --repro]
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard551  # noqa: E402,F401
import argparse  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import designs551 as D  # noqa: E402
import lib551 as L  # noqa: E402
from lib551 import ALPHA, GAMMA, LAM, SLOTS, ROOT, E062  # noqa: E402

TOL = 1e-12
W438 = ROOT / "experiments/438-october-stack-repricing/work"
W493 = ROOT / "experiments/493-pertarget-q-october/work"
W546 = ROOT / "experiments/546-october-cap-pertarget/work"
W522 = ROOT / "experiments/522-october-bulk-key/work"
W534 = ROOT / "experiments/534-generalization-audit/work/partB"


def pct_hi(kv, carried, p):
    """522's fold522.pct_hi / def_lodo_key.py --hi-pct: numpy-linear p-th percentile of key/median over carried."""
    return float(np.percentile(kv[carried] / float(np.median(kv[carried])), p))


def run(panel, held, level_file=None, ka=False, repro=False, cfg=None):
    t0 = time.time()
    log = lambda *a: print(f"[comp551 {panel}/{held}]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    if ka:
        assert panel in ("efp", "abc"), "the known answers run on efp (and abc, the supplementary A/B/C column)"
    ctx = L.bind(panel)
    fold = L.build_fold438(ctx, held, "held")
    scored, pp_pos, wcols = fold["scored"], fold["panel_positions"], fold["wcols"]
    truth_full, pre_cached = fold["truth"], fold["pre"]
    truth = truth_full[:, scored]
    targets, sources, op = list(fold["targets"]), list(fold["sources"]), fold["op"]
    nonpanel = ctx.nonpanel
    res = dict(panel=panel, fold=held, config=cfg, rules=ctx.rules, n_targets=len(targets), sources=sources, audit=fold["audit"],
               op=str(fold["op_path"].relative_to(ROOT)), op_sha256=op.sha256, op_axis_checked=op.axis_checked,
               op_axis_sha256=getattr(op, "axis_sha256", None), ka={}, arms={}, comp={}, keys={})
    log(f"{len(targets)} targets; op {op.name} {op.sha256[:12]} ({op.axis_checked}); audit {json.dumps(fold['audit'])}")
    if ka:                                                   # KA-FOLD: my fold cache from scratch == 438's
        z = np.load(W438 / panel / f"fold-{held}-held.npz", allow_pickle=True)
        assert tuple(str(t) for t in z["targets"]) == tuple(targets), "KA-FOLD targets"
        for k in ("pre", "truth"):
            assert np.array_equal(np.asarray(z[k], np.float64), fold[k]), f"KA-FOLD {k}"
        for k in ("b_ctrl", "gate", "scored", "wcols"):
            assert np.array_equal(np.asarray(z[k]), np.asarray(fold[k])), f"KA-FOLD {k}"
        res["ka"]["fold"] = f"bit-identical to 438's {panel} fold cache (targets, pre, truth, b_ctrl, gate, scored, wcols)"
        log("KA-FOLD: my fold cache == 438's, bit for bit")
    if repro:                                                # KA-REPRO (fold): the rebuilt R6 fold cache == 534's (C10a's)
        z = np.load(W534 / panel / f"fold-{held}-held.npz", allow_pickle=True)
        assert tuple(str(t) for t in z["targets"]) == tuple(targets), "KA-REPRO fold targets"
        for k in ("pre", "truth"):
            assert np.array_equal(np.asarray(z[k], np.float64), fold[k]), f"KA-REPRO fold {k}"
        for k in ("b_ctrl", "gate", "scored", "wcols"):
            assert np.array_equal(np.asarray(z[k]), np.asarray(fold[k])), f"KA-REPRO fold {k}"
        res["ka"]["repro_fold"] = "bit-identical to 534's (C10a's) fold cache on this R6 panel"
        log("KA-REPRO fold: my fold cache == C10a's, bit for bit")

    dn = L.Donors438(targets, sources)
    cd4, partial, h1, daud = L.residuals(ctx, targets)
    res["deflation"] = daud
    rows_of = lambda r: int((np.linalg.norm(np.where(np.asarray(r.gene_votes) > 0, np.asarray(r.mean_effect), 0.0), axis=1) > 0).sum())
    res["coverage"] = dict(h1_union_rows=rows_of(h1), cd4_complete_rows=rows_of(cd4), cd4_partial_rows=rows_of(partial))
    pre = L.pre_normrel(dn, sources, "cap", -1.5, fallback_rows=pre_cached)
    whit = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    FC = L.Finals(ctx, whit, cd4, partial, h1, wcols, pp_pos)
    DESIGNS = D.designs_for(cfg, ka=ka and panel == "efp") if (ka and panel == "efp") else D.designs_for(cfg)
    for c in sorted({d[1] for d in DESIGNS}):
        FC.get(c)
    res["comp"] = FC.meta
    log("compositions: " + ", ".join(f"{c} f={FC.meta[c]['fractions']['cd4']:.4f}/{FC.meta[c]['fractions']['h1_union']:.4f}"
                                      for c in sorted(FC.meta)) + f"; coverage {json.dumps(res['coverage'])}")

    # ============================================================== keys (comp438's lines + the percentile cap)
    d_rows = {s: L.donor_rows(s, targets)[0] for s in sources}
    A_ = {s: dn.A[s] for s in sources}
    carried_all = sum(A_[s].astype(int) for s in sources) >= 2
    key410, carried410, _ = L.shared_energy_key(targets, sources, op, wcols, cols=None)
    bases = L.lodo_bases(op, dn, sources, "cap", -1.5, pre, wcols)
    kv1 = L.lodo_key(FC.final("h1c1"), FC.parts("h1c1"), bases, d_rows, A_, carried_all, sources, pp_pos)
    del bases, d_rows
    kvv = np.nan_to_num(kv1, nan=0.0)
    hi = pct_hi(kvv, carried_all, D.P_PROD)
    W = {"K3": L.weights(key410, carried410, 0.3, 3.0)[0],
         "V1_5_3": L.weights(kvv, carried_all, 0.5, 3.0)[0],
         "V1pct": L.weights(kvv, carried_all, 0.5, hi)[0]}
    rng = np.random.default_rng(D.CTRL_SEED)                  # 438's CTRL construction, on the V1_5_3 weights
    W["V1_5_3shuf"] = W["V1_5_3"][rng.permutation(len(W["V1_5_3"]))]
    if ka and panel == "efp":                                 # 546o's rounded P (KA rows only)
        hi546 = pct_hi(kvv, carried_all, D.P_546)
        W["V1p546"] = L.weights(kvv, carried_all, 0.5, hi546)[0]
    res["keys"] = {k: dict(q=np.quantile(w, [0, .1, .25, .5, .75, .9, 1]).tolist(),
                           frac_floor=float(np.mean(w[carried_all] <= 0.5 + 1e-12)), w_max=float(np.max(w)))
                   for k, w in W.items()}
    res["keys"]["carried_lodo"] = int(carried_all.sum())
    res["keys"]["carried_410"] = int(np.asarray(carried410).sum())
    res["keys"]["pct_hi"] = hi
    res["keys"]["n_above_cap3_under_pct"] = int((W["V1pct"] > 3.0 + 1e-12).sum())
    log(f"keys: carried LODO {int(carried_all.sum())}, 410 {int(np.asarray(carried410).sum())}; percentile cap hi {hi:.4f} "
        f"({res['keys']['n_above_cap3_under_pct']} targets above 3)")

    P_DE = L.carrier_profiles(FC.final("h1c1"))
    bulk_k3 = L.key_scale(FC.final("h1c1"), W["K3"])[0]
    dk3 = L.realised_rows(bulk_k3, P_DE, q=None, kmaxes=(), zero_cols=None)
    A_k3 = {c: L.energy(dk3[c]["design"], nonpanel)["A"] for c in SLOTS}
    bet = L.betas(ctx, held, nonpanel, A_pp_k3=A_k3)
    arrays = {}
    strength = np.linalg.norm(FC.final("h2c1")[:, nonpanel], axis=1)       # 493h/546o: the unkeyed h2c1 row's norm
    QV = {}

    def qvec(name):
        """551: the design's per-target exponent vector on this fold (None for a scalar-q design)."""
        _, comp, key, q, pt, perm, _ = D.BY_NAME[name]
        if pt is None:
            return None
        if name not in QV:
            qv = L.q_s_p2(strength)[0] if pt == "S_P2" else L.q_p2g(W[key], carried_all)[0]
            if perm:
                qv = qv[np.random.default_rng(D.PT_CTRL_SEED).permutation(len(qv))]
            QV[name] = qv
        return QV[name]

    def rec(name, comp, key, q, scale=None, realised=True):
        base = name[:-2] if name.endswith(":L") else name
        bulk = FC.final(comp) if key is None else L.key_scale(FC.final(comp), W[key])[0]
        qv = qvec(base)
        if qv is None:
            dr = L.realised_rows(bulk, P_DE, q=q, scale=scale, kmaxes=(16,) if realised else ())
        else:
            dr = L.realised_rows_pt(bulk, P_DE, qv, scale=scale, kmaxes=(16,) if realised else ())
        sd, ad = L.row_score({c: dr[c]["design"] for c in SLOTS}, scored, truth, bet, nonpanel)
        entry = dict(comp=comp, key=key, q=q, pt=D.BY_NAME[base][4], perm=D.BY_NAME[base][5], scale=scale, design=sd,
                     realised={})
        arrays[f"{name}@D"] = ad
        if realised:
            sr, ar = L.row_score({c: dr[c]["realised"][16] for c in SLOTS}, scored, truth, bet, nonpanel)
            K = np.concatenate([dr[c]["stats"][16]["K"] for c in SLOTS])
            sr["frac_capped"], sr["K_mean"] = float(np.mean(K >= 16)), float(np.mean(K))
            entry["realised"]["16"] = sr
            arrays[f"{name}@R16"] = ar
        res["arms"][name] = entry
        r16 = entry["realised"].get("16")
        log(f"{name:32s} cos {sd['cos_agg']:.6f} pds@.5 {sd['pds_lam0.5']:.6f}"
            + (f" | R16 cos {r16['cos_agg']:.6f} pds@.5 {r16['pds_lam0.5']:.6f} capped {100 * r16['frac_capped']:.0f}%" if r16 else ""))
        return entry

    for name, comp, key, q, pt, perm, realised in DESIGNS:
        rec(name, comp, key, q, realised=realised)
    if level_file:
        lv = json.loads(open(level_file).read())
        res["level_file"] = str(level_file)
        res["level_sha256"] = L.sha256(level_file)
        for name, comp, key, q, pt, perm, realised in DESIGNS:
            fac = {c: float(lv["factors"][D.energy_name(name)][c]) for c in SLOTS}
            rec(f"{name}:L", comp, key, q, scale=fac)
    res["pt_counts"] = {n: np.unique(v, return_counts=True)[1].tolist() for n, v in QV.items()}

    # ============================================================== known answers (efp), before anything is written
    if ka:
        worst = {}

        def chk(tag, mine, ref, fields=("cos_agg", "pds", "pds_lam0.5", "pds_lam1")):
            d = L.cmp(mine, ref, fields)
            worst[tag] = d
            log(f"KA {tag:58s} |d| {d:.1e}")
            assert d <= TOL, f"KA FAILED {tag}: {d}"

        c438 = json.loads((W438 / panel / f"comp438-{held}-held.json").read_text())["arms"]
        for n in ("h1c1|none", "h1c1|K3", "h1c1|K3+sh1.25", "h2c1|K3+sh1.25", "h2c1|K3+sh1.5", "h2c1.5|K3+sh1.75"):
            for tw in ("", ":L") if level_file else ("",):
                mine, ref = res["arms"][n + tw], c438[n + tw]
                chk(f"438 {n + tw} design", mine["design"], ref["design"])
                if ref.get("realised") and not (n == "h1c1|none" and tw == ""):
                    chk(f"438 {n + tw} realised16", mine["realised"]["16"], ref["realised"]["16"])
                    for f in ("frac_capped", "K_mean"):
                        if f in ref["realised"]["16"]:
                            d = abs(mine["realised"]["16"][f] - ref["realised"]["16"][f])
                            assert d <= TOL, f"KA {n + tw} {f}: {d}"
        r493 = json.loads((W493 / f"comp493-{panel}-{held}.json").read_text())["arms"]["base q1.75"]
        chk("493 base q1.75 == M4 design", res["arms"][D.M["M4"]]["design"], r493["design"])
        chk("493 base q1.75 == M4 realised16", res["arms"][D.M["M4"]]["realised"]["16"], r493["realised"]["16"])
    if ka and panel == "efp":
        r546 = json.loads((W546 / f"comp546o-efp-{held}.json").read_text())["arms"]["base q1.75"]
        chk("546o base q1.75 == M5 design", res["arms"][D.M["M5"]]["design"], r546["design"])
        chk("546o base q1.75 == M5 realised16", res["arms"][D.M["M5"]]["realised"]["16"], r546["realised"]["16"])
        a522 = json.loads((W522 / f"arms522-efp-{held}.json").read_text())
        p = a522["arms"]["pctPROD"]
        mine = res["arms"][D.M["M6"]]["realised"]["16"]
        d = max(abs(mine["cos_agg"] - p["cos_agg"]), abs(mine["pds"] - p["pds"]), abs(mine["pds_lam0.5"] - p["pds_lam05"]),
                abs(mine["K_mean"] - p["K_mean"]), abs(hi - a522["his"]["hi_pct_prod"]), abs(float(np.max(W["V1pct"])) - p["w_max"]),
                *[abs(mine["AB"][c][x] / p["AB"][c][x] - 1) for c in SLOTS for x in ("A", "B")])
        worst["522 pctPROD == M6 realised16 (cos_agg, pds, pds lam .5, K_mean, AB, hi)"] = d
        log(f"KA 522 pctPROD == M6 realised16 (+hi, AB, K_mean)                   |d| {d:.1e}")
        assert d <= TOL, f"KA FAILED 522 pctPROD: {d}"
        c3 = a522["arms"]["cap3"]
        m5 = res["arms"][D.M["M5"]]["realised"]["16"]
        d = max(abs(m5["cos_agg"] - c3["cos_agg"]), abs(m5["pds"] - c3["pds"]))
        worst["522 cap3 == M5 realised16 (cos_agg, pds)"] = d
        assert d <= TOL, f"KA FAILED 522 cap3: {d}"
    if ka and panel == "efp":                                # 551: KA-546o, KA-493h, KA-534
        a546 = json.loads((W546 / f"comp546o-efp-{held}.json").read_text())["arms"]
        for mine_n, ref_n in (("h2c1|V1p546+sh1.75", "p96.6"), ("h2c1|V1p546+ptP2g", "p96.6 + P2g"),
                              ("h2c1|V1p546+ptS_P2", "p96.6 + S_P2")):
            chk(f"546o {ref_n} == {mine_n} design", res["arms"][mine_n]["design"], a546[ref_n]["design"])
            chk(f"546o {ref_n} == {mine_n} realised16", res["arms"][mine_n]["realised"]["16"], a546[ref_n]["realised"]["16"])
            d = abs(res["arms"][mine_n]["realised"]["16"]["K_mean"] - a546[ref_n]["realised"]["16"]["K_mean"])
            worst[f"546o {ref_n} K_mean"] = d
            assert d <= TOL, f"KA 546o {ref_n} K_mean: {d}"
        a493 = json.loads((W493 / f"comp493h-efp-{held}.json").read_text())["arms"]
        for mine_n, ref_n in (("h2c1|V1_5_3+ptP2g", "P2g_carried"), ("h2c1|V1_5_3+ptS_P2", "S_P2")):
            chk(f"493h {ref_n} == {mine_n} design", res["arms"][mine_n]["design"], a493[ref_n]["design"])
            chk(f"493h {ref_n} == {mine_n} realised16", res["arms"][mine_n]["realised"]["16"], a493[ref_n]["realised"]["16"])
        c534 = json.loads((W534 / "efp" / f"comp534-{held}.json").read_text())["arms"]
        for n, ref in c534.items():
            chk(f"534 {n} design", res["arms"][n]["design"], ref["design"])
            if ref.get("realised"):
                chk(f"534 {n} realised16", res["arms"][n]["realised"]["16"], ref["realised"]["16"])
    if repro:                                                # KA-REPRO: every 534 design and twin == C10a's stored values
        worst = {}

        def chk(tag, mine, ref, fields=("cos_agg", "pds", "pds_lam0.5", "pds_lam1")):
            d = L.cmp(mine, ref, fields)
            worst[tag] = d
            assert d <= TOL, f"KA-REPRO FAILED {tag}: {d}"
        c534 = json.loads((W534 / panel / f"comp534-{held}.json").read_text())
        assert c534["op_sha256"] == op.sha256, "KA-REPRO: the re-fitted fold operator is not C10a's"
        for n, ref in c534["arms"].items():
            chk(f"C10a {n} design", res["arms"][n]["design"], ref["design"])
            if ref.get("realised"):
                chk(f"C10a {n} realised16", res["arms"][n]["realised"]["16"], ref["realised"]["16"])
        res["ka"]["repro"] = dict(rows=len(worst), worst=max(worst.values()), op_sha256=op.sha256)
        log(f"KA-REPRO: {len(worst)} C10a rows reproduced on {panel}/{held}, worst |d| {max(worst.values()):.1e}")
    if ka:
        res["ka"]["rows"] = worst
        res["ka"]["worst"] = max(worst.values())
        log(f"ALL KNOWN ANSWERS PASS on {panel}/{held}: worst |d| {res['ka']['worst']:.1e}")

    res["peak_rss_gib"] = E062.peak_rss_gib()
    res["seconds"] = time.time() - t0
    out = ctx.work / f"comp551-{held}.json"
    out.write_text(json.dumps(res, indent=1, default=float))
    np.savez_compressed(ctx.work / f"arrays551-{held}.npz", names=np.asarray(list(arrays)), targets=np.asarray(targets),
                        d=np.asarray((truth * truth).sum(1)), **{f"{n}__{k}": v for n, a in arrays.items() for k, v in a.items()})
    log(f"done; {len(res['arms'])} arms; peak RSS {res['peak_rss_gib']:.2f} GiB -> {out.name}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("panel")
    ap.add_argument("fold", help="one fold, or several comma-separated: run in sequence in ONE process (one gate opening)")
    ap.add_argument("--level", default=None)
    ap.add_argument("--ka", action="store_true")
    ap.add_argument("--repro", action="store_true", help="R6 panel inside the C10b access: reproduce C10a's rows first")
    ap.add_argument("--config", default=None, help="F, one of designs551.CONFIGS (required except for --ka on efp)")
    a = ap.parse_args()
    assert not (a.ka and a.repro)
    if a.panel.startswith("S"):
        assert a.repro and a.config and not a.ka, "an R6 panel runs only inside the C10b access: --repro --config F"
        cfgrec = json.loads((L.WORK / "c10b-config.json").read_text())
        assert cfgrec["config"] == a.config, f"--config {a.config} != the recorded configuration {cfgrec['config']}"
    import gc
    for fold in a.fold.split(","):
        done_json = L.WORK / a.panel / f"comp551-{fold}.json"
        if done_json.exists() and "peak_rss_gib" in json.loads(done_json.read_text()) and \
                bool(json.loads(done_json.read_text())["ka"]) == a.ka:
            print(f"[comp551 {a.panel}/{fold}] already priced -> skip", flush=True)
            continue
        run(a.panel, fold, a.level, a.ka, a.repro, a.config)
        gc.collect()
