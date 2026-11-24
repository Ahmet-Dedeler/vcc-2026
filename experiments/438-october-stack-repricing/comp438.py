"""438 -- 420's comp420 on ANY panel, one PLOLO fold per process: the stack grid (H1 x CD4 multipliers x key x q), scored
as DESIGNED and as REALISED by the ballast (KMAX 16; 16/24/32 on the known-answer rows), plus each design's re-levelled
variant when a level file (full438's output) is given.

KNOWN ANSWERS on the real A/B/C panel (abort before any grid number is written):
  KA-ROWS   my panel-agnostic donor rows (cached bank statistics) == 415's full-<donor>.npz, bit for bit (4 banks)
  KA-FOLD   my fold (targets, pre, truth, b_ctrl, gate, scored, wcols) == 420's cached fold, bit for bit
  KA-DONORS Donors438 (F, A, normrel medians, PM correction) == sens415.Donors on the fold, bit for bit
  KA-PP     compose(rules abc, m = 1) == sens415.compose_pp (final rows + merged residual), bit for bit
  KA-KEYS   K3 and LODO-V1 weights == 418's saved W3 / V1__PP, bit for bit
  KA-418    PP|none (design) == comp418's PP|none
  KA-420    PP|K3, PP|K3+sh1.25, h1x2|K3+sh1.25, cd4x1.5|K3+sh1.25, stack|K3+sh1.25 == comp420-<fold>.json, design and
            realised KMAX 16/24/32 (cos_agg, pds, pds lam 0.5/1, render A/B per frame): |d| = 0      (the mandate's (b))
  KA-REAL   my row scorer on the designed rows == lib418.score on those rows (1e-12)

    uv run python comp438.py <panel> <fold> [--fold-ops held|prod] [--level work/<panel>/level438.json] [--tag T]
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np

import lib438 as L
from lib438 import ALPHA, GAMMA, LAM, SLOTS, W418, W420, ROOT, E062, lib418, sens415

ALIAS420 = {"PP|K3": "h1c1|K3", "PP|K3+sh1.25": "h1c1|K3+sh1.25", "h1x2|K3+sh1.25": "h2c1|K3+sh1.25",
            "cd4x1.5|K3+sh1.25": "h1c1.5|K3+sh1.25", "stack|K3+sh1.25": "h2c1.5|K3+sh1.25"}


def run(panel, held, fold_ops="held", level_file=None, tag=None, resume=False):
    t0 = time.time()
    log = lambda *a: print(f"[{panel}/{held}/{fold_ops}]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    ctx = L.bind(panel, tag=tag)
    fold = L.build_fold438(ctx, held, fold_ops)
    scored, pp_pos, wcols = fold["scored"], fold["panel_positions"], fold["wcols"]
    truth_full, pre_cached = fold["truth"], fold["pre"]
    truth = truth_full[:, scored]
    targets, sources, op = list(fold["targets"]), list(fold["sources"]), fold["op"]
    nonpanel = ctx.nonpanel
    res = dict(panel=panel, fold=held, fold_ops=fold_ops, rules=ctx.rules, n_targets=len(targets), sources=sources,
               audit=fold["audit"], op=str(fold["op_path"].relative_to(ROOT)), op_sha256=op.sha256, ka={}, arms={},
               comp={}, keys={})
    log(f"{len(targets)} targets; op {op.name} ({op.axis_checked}); audit {json.dumps(fold['audit'])}")
    abc = panel == "abc" and fold_ops == "held"

    # ============================================================== KA-ROWS / KA-FOLD / KA-DONORS (A/B/C only)
    if abc:
        for s in (*L.MAIN, "perturbmatch"):
            z = np.load(ROOT / f"experiments/415-lodo-key/work/full-{s}.npz")
            f, a, *_ = L.rows(s, list(ctx.panel))
            assert np.array_equal(f, z["effects"]) and np.array_equal(a, z["available"]), f"KA-ROWS {s}"
        z = np.load(W420 / f"plolo/fold-{held}.npz", allow_pickle=True)
        assert tuple(str(t) for t in z["targets"]) == tuple(targets), "KA-FOLD targets"
        for k in ("pre", "truth"):
            assert np.array_equal(np.asarray(z[k], np.float64), fold[k]), f"KA-FOLD {k}"
        for k in ("b_ctrl", "gate", "scored", "wcols"):
            assert np.array_equal(np.asarray(z[k]), np.asarray(fold[k])), f"KA-FOLD {k}"
        res["ka"]["rows_fold"] = "bit-identical (4 banks; pre, truth, b_ctrl, gate, scored, wcols, targets)"
        log("KA-ROWS (4 banks) and KA-FOLD vs 420's cached fold: bit-identical")

    dn = L.Donors438(targets, sources)
    if abc:
        ref = sens415.Donors(targets, sources)
        for s in sources:
            assert np.array_equal(dn.F[s], ref.F[s]) and np.array_equal(dn.A[s], ref.A[s]) and dn.med[s] == ref.med[s]
        assert np.array_equal(dn.pmcorr["prod"], ref.pmcorr["prod"]), "KA-DONORS pmcorr"
        log("KA-DONORS: Donors438 == sens415.Donors (rows, availability, medians, PM correction)")
    cd4, partial, h1, daud = L.residuals(ctx, targets)
    res["deflation"] = daud
    rows_of = lambda r: int((np.linalg.norm(np.where(np.asarray(r.gene_votes) > 0, np.asarray(r.mean_effect), 0.0), axis=1) > 0).sum())
    res["coverage"] = dict(h1_union_rows=rows_of(h1), cd4_complete_rows=rows_of(cd4), cd4_partial_rows=rows_of(partial))
    pre = L.pre_normrel(dn, sources, "cap", -1.5, fallback_rows=pre_cached)
    whit = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    FC = L.Finals(ctx, whit, cd4, partial, h1, wcols, pp_pos)
    need = sorted({L.parse(n)[0] for n, *_ in L.grid()})
    for c in need:
        FC.get(c)
    res["comp"] = FC.meta
    if abc:
        frozen, part_ = sens415.cd4_frozen_and_partial(targets)
        cd4_def, _ = sens415.deflate(frozen, 5)
        f415, _, _, radd415, _ = sens415.compose_pp(op, pre, cd4_def, part_, h1, wcols, pp_pos)
        assert np.array_equal(FC.final("h1c1"), f415) and np.array_equal(FC.parts("h1c1")["radd"], radd415), "KA-PP"
        log("KA-PP: compose(abc rules, m=1) == sens415.compose_pp")
    log("compositions: " + ", ".join(f"{c} f={FC.meta[c]['fractions']['cd4']:.4f}/{FC.meta[c]['fractions']['h1_union']:.4f}"
                                      for c in ("h1c1", "h2c1", "h1c1.5", "h2c1.5")))

    # ============================================================== keys
    d_rows = {s: L.donor_rows(s, targets)[0] for s in sources}
    A_ = {s: dn.A[s] for s in sources}
    carried_all = sum(A_[s].astype(int) for s in sources) >= 2
    key410, carried410, _ = L.shared_energy_key(targets, sources, op, wcols, cols=None)
    bases = L.lodo_bases(op, dn, sources, "cap", -1.5, pre, wcols)
    kv1 = L.lodo_key(FC.final("h1c1"), FC.parts("h1c1"), bases, d_rows, A_, carried_all, sources, pp_pos)
    W = {"K3": L.weights(key410, carried410, 0.3, 3.0)[0], "K5_3": L.weights(key410, carried410, 0.5, 3.0)[0],
         "V1": L.weights(np.nan_to_num(kv1, nan=0.0), carried_all, 0.3, 3.0)[0],
         "V1_5_3": L.weights(np.nan_to_num(kv1, nan=0.0), carried_all, 0.5, 3.0)[0]}
    if abc:
        z418 = np.load(W418 / f"arrays418-{held}.npz")
        assert np.array_equal(W["V1"], z418["V1__PP"]) and np.array_equal(W["K3"], z418["W3"]), "KA-KEYS"
        log("KA-KEYS: K3 and LODO-V1 weights == 418's saved W3 / V1__PP")
    res["keys"] = {k: dict(q=np.quantile(w, [0, .1, .25, .5, .75, .9, 1]).tolist(),
                           frac_floor=float(np.mean(w <= L.KEYS[k][1] + 1e-12)), frac_cap=float(np.mean(w >= 3.0 - 1e-12)))
                   for k, w in W.items()}
    res["keys"]["carried"] = int(carried_all.sum())
    res["keys"]["spearman_K3_V1"] = float(np.corrcoef(np.argsort(np.argsort(W["K3"])), np.argsort(np.argsort(W["V1"])))[0, 1])

    P_DE = L.carrier_profiles(FC.final("h1c1"))
    # the board-emulated beta needs the fold's PP|K3 design energy (on A/B/C it is lib418.betas_for verbatim)
    bulk_k3 = L.key_scale(FC.final("h1c1"), W["K3"])[0]
    dk3 = L.realised_rows(bulk_k3, P_DE, q=None, kmaxes=(), zero_cols=None)
    A_k3 = {c: L.energy(dk3[c]["design"], nonpanel)["A"] for c in SLOTS}
    bet = L.betas(ctx, held, nonpanel, A_pp_k3=A_k3)
    arrays = {}
    old_json = ctx.work / f"comp438-{held}-{fold_ops}.json"
    if resume and old_json.exists():                 # add only the designs a finished run does not have (grid extension)
        old = json.loads(old_json.read_text())
        assert old["fold"] == held and old["rules"] == ctx.rules and old["op_sha256"] == op.sha256
        res["ka"], res["arms"] = old["ka"], old["arms"]
        Zo = np.load(ctx.work / f"arrays438-{held}-{fold_ops}.npz")
        for nm in Zo["names"]:
            arrays[str(nm)] = {k: Zo[f"{nm}__{k}"] for k in ("x", "a", "pds", "pdsl")}
        chk = res["arms"]["h1c1|K3+sh1.25"]["realised"]["16"]        # the resumed setup must reproduce the old rows
        log(f"RESUME: {len(res['arms'])} arms loaded")

    def rec(name, comp, key, q, scale=None, kmaxes=(16,), check=False, realised=True):
        bulk = FC.final(comp) if key is None else L.key_scale(FC.final(comp), W[key])[0]
        dr = L.realised_rows(bulk, P_DE, q=q, scale=scale, kmaxes=kmaxes if realised else ())
        sd, ad = L.row_score({c: dr[c]["design"] for c in SLOTS}, scored, truth, bet, nonpanel)
        if check and scale is None:
            s, _ = lib418.score(bulk, fold, truth, bet, q)
            dka = L.cmp(sd, s)
            assert dka < 1e-12, f"KA-REAL {name}: {dka}"
            res["ka"]["real_worst"] = max(res["ka"].get("real_worst", 0.0), dka)
        entry = dict(comp=comp, key=key, q=q, scale=scale, design=sd, realised={})
        arrays[f"{name}@D"] = ad
        for k in (kmaxes if realised else ()):
            sr, ar = L.row_score({c: dr[c]["realised"][k] for c in SLOTS}, scored, truth, bet, nonpanel)
            K = np.concatenate([dr[c]["stats"][k]["K"] for c in SLOTS])
            sr["frac_capped"], sr["K_mean"] = float(np.mean(K >= k)), float(np.mean(K))
            entry["realised"][str(k)] = sr
            arrays[f"{name}@R{k}"] = ar
        res["arms"][name] = entry
        r16 = entry["realised"].get("16")
        log(f"{name:30s} cos {sd['cos_agg']:.6f} pds@.5 {sd['pds_lam0.5']:.6f}"
            + (f" | R16 cos {r16['cos_agg']:.6f} pds@.5 {r16['pds_lam0.5']:.6f} capped {100 * r16['frac_capped']:.0f}%" if r16 else ""))
        return entry

    # ============================================================== KA-418 / KA-420 (A/B/C, held-out operators)
    if resume and res["arms"]:
        e = rec("h1c1|K3+sh1.25#resumecheck", "h1c1", "K3", 1.25)
        d = L.cmp(e["realised"]["16"], chk)
        del res["arms"]["h1c1|K3+sh1.25#resumecheck"], arrays["h1c1|K3+sh1.25#resumecheck@D"], arrays["h1c1|K3+sh1.25#resumecheck@R16"]
        assert d == 0.0, f"RESUME: the re-built setup does not reproduce the stored 417a analogue: {d}"
        log("RESUME: re-built setup reproduces the stored 417a-analogue row exactly")
    elif True:
        rec("h1c1|none", "h1c1", None, None, realised=False, check=abc)
    if abc and not resume:
        c418 = json.loads((W418 / f"comp418-{held}.json").read_text())["arms"]["PP|none"]
        d = L.cmp(res["arms"]["h1c1|none"]["design"], c418)
        assert d == 0.0, f"KA-418 PP|none: {d}"
        c420 = json.loads((W420 / f"comp420-{held}.json").read_text())["arms"]
        worst = 0.0
        for r420, mine in ALIAS420.items():
            comp, key, q = L.parse(mine)
            e = rec(mine, comp, key, q, kmaxes=L.KMAXES, check=True)
            d = L.cmp(e["design"], c420[r420]["design"])
            for k in L.KMAXES:
                d = max(d, L.cmp(e["realised"][str(k)], c420[r420]["realised"][str(k)]),
                        abs(e["realised"][str(k)]["frac_capped"] - c420[r420]["realised"][str(k)]["frac_capped"]))
            worst = max(worst, d)
            log(f"KA-420 {r420:22s} (design + realised K16/24/32) vs comp420: |d| {d:.1e}")
            assert d == 0.0, f"KA-420 FAILED {held} {r420}: {d}"
        res["ka"]["ka420_worst"] = worst
        res["ka"]["ka418_pp_none"] = 0.0
        log("ALL A/B/C KNOWN ANSWERS PASS")
    elif not resume:
        rec("h1c1|K3", "h1c1", "K3", None)

    # ============================================================== the grid
    for name, comp, key, q in L.grid():
        if name not in res["arms"]:
            rec(name, comp, key, q)
    rng = np.random.default_rng(420)
    W["V1shuf"] = W["V1"][rng.permutation(len(W["V1"]))]
    if "CTRL:h1c1|V1shuf+sh1.25" not in res["arms"]:
        rec("CTRL:h1c1|V1shuf+sh1.25", "h1c1", "V1shuf", 1.25)
    res["peak_rss_gib"] = E062.peak_rss_gib()
    out = ctx.work / f"comp438-{held}-{fold_ops}.json"
    out.write_text(json.dumps(res, indent=1, default=float))

    # ============================================================== re-levelled variants (factors from full438)
    if level_file:
        lv = json.loads(open(level_file).read())
        res["level_file"] = str(level_file)
        for name, comp, key, q in L.grid():
            if name not in lv["factors"] or f"{name}:L" in res["arms"]:
                continue
            fac = {c: float(lv["factors"][name][c]) for c in SLOTS}
            rec(f"{name}:L", comp, key, q, scale=fac)
        res["peak_rss_gib"] = E062.peak_rss_gib()
    out.write_text(json.dumps(res, indent=1, default=float))
    np.savez_compressed(ctx.work / f"arrays438-{held}-{fold_ops}.npz", names=np.asarray(list(arrays)),
                        d=np.asarray((truth * truth).sum(1)), **{f"{n}__{k}": v for n, a in arrays.items() for k, v in a.items()})
    log(f"done; {len(res['arms'])} arms; peak RSS {res['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("panel")
    ap.add_argument("fold")
    ap.add_argument("--fold-ops", default="held", choices=("held", "prod"))
    ap.add_argument("--level", default=None)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--resume", action="store_true")
    a = ap.parse_args()
    run(a.panel, a.fold, a.fold_ops, a.level, a.tag, a.resume)
