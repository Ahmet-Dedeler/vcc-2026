"""414 — re-price the BULK-only magnitude arms on PLOLO, now that ballast unwelds the DE channel (410b).

Under ballast the bulk (pseudobulk -> pds, mse) and the carriers (DE table -> fid/reach/jac/nmae) are separate
objects, so a bulk-only arm is priced on pds + mse alone. This script prices, on 162's PLOLO folds, with 228's
PRODUCTION composition and the a23 render read-back (plolo.deliver, minnorm), in the challenge frames A/B/C:

    pre --(x ALPHA)--> Operator(gamma 0.65) --> whitened [--(x m_gain, median-row-norm preserve)--]   (228's way, "gi")
    whitened + merged(CD4 @ 0.34, H1 union @ 0.28) --> final, panel coords zeroed
    final --(x key w_p * c, Frobenius preserve)--> bulk                                            (410's way, "K3"/"Kinf")
    bulk [--(x m_gain per gene, Frobenius preserve)--]                                               (410d's way, "gf")
    bulk [--(169 sharpening sign|h|^q, row renorm, per frame gate)--]                                ("sh"/"shL")
    deliver(a23 minnorm) in A/B/C  ->  pds, macro cos, sse -> cos_agg = sqrt(1 - sse), per frame, then mean

KNOWN-ANSWER ROWS (abort on failure), BEFORE any new number:
  KA-1  410 check410.py at effect level: key clip(0.3,3) and clip(0.3,inf) gains == check410.json (1e-9)
  KA-2  228 BASE (t = 0) and gain_reg@{0.1,0.25,0.5,1} rows == arms228-<fold>.json cos_agg / pds (1e-12)
  KA-3  169 a23 stack, sharpening q = 1.00 / 1.15 (h space, row) == plolo169.json pds / sse per frame (1e-9)

    uv run python arms414.py --folds k562_gwps orion_hct116 orion_hek293t kolf2
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
for p in ("experiments/162-retrodicting-validator", "experiments/218-independent-donors",
          "experiments/228-theta-gain-plolo", "experiments/410-pertscale-revisit",
          "experiments/169-sharpening-six-member"):
    sys.path.insert(0, str(ROOT / p))

import plolo  # noqa: E402
import arms228  # noqa: E402  (sets plolo.WORK to 228's dir on import; reset below)
from arms228 import CD4_F, UNION, UNION_G, apply_gain, multiplier  # noqa: E402
from compose218 import merged, solve_fraction_carrying, union_rows  # noqa: E402
from resid218 import ALPHA, GAMMA, LAM, source_residual  # noqa: E402
from key410 import shared_energy_key, weights  # noqa: E402
from sharp import sharpen, sharpen_lfc  # noqa: E402
from run_arms import CONFIGS, effects_for  # noqa: E402

plolo.WORK = WORK / "plolo"          # the copied, read-only fold caches live here
from plolo import E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for, deliver, endpoints, solve_fraction  # noqa: E402

TS = (0.1, 0.25, 0.5, 1.0)
QS = (1.05, 1.10, 1.15, 1.25)
FRAMES = {}


def frames():
    if not FRAMES:
        for ctx in "ABC":
            b = plolo.challenge_bctrl(ctx)
            expr = np.expm1(b)
            cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
            FRAMES[ctx] = (b, (cpm > plolo.MIN_CPM) & (b > 0))
    return FRAMES


def frob_preserve(E, M):
    """Scale M so its Frobenius norm equals E's (panel coords are zero in both, so this is the non-panel norm)."""
    return M * (np.linalg.norm(E) / max(np.linalg.norm(M), 1e-300))


def key_scale(final, w):
    """410's builder `apply_pert_scale`: rows x w, one constant restoring the non-panel Frobenius norm."""
    e2 = np.einsum("ij,ij->i", final, final)
    c = float(np.sqrt(e2.sum() / float((w * w * e2).sum())))
    return final * (w * c)[:, None], w * c


def gene_gain_final(bulk, m):
    """410d's builder `apply_bulk_gene_gain`: columns x m, then restore the non-panel Frobenius norm."""
    return frob_preserve(bulk, bulk * m[None, :])


def score(bulk, fold, truth, sharp=None):
    """a23 read-back in the three challenge frames. Returns summary + per-target arrays for bootstrap."""
    scored = fold["scored"]
    nonpanel = np.ones(bulk.shape[1], dtype=bool)
    nonpanel[fold["panel_positions"]] = False
    out, X, Aa, pdsp, cosp = {}, [], [], [], []
    AB = {}
    for ctx, (b, gate) in frames().items():
        F = bulk
        if sharp is not None:
            kind, q = sharp
            F = sharpen(bulk, q, gate, mode="row") if kind == "h" else sharpen_lfc(bulk, q, b, gate, mode="row")
        dl = deliver(F, b, gate, "minnorm")
        P = dl[:, scored]
        e = endpoints(P, truth)
        e["cos_agg"] = float(np.sqrt(max(0.0, 1.0 - e["sse"])))
        out[ctx] = e
        X.append((P * truth).sum(1))
        Aa.append((P * P).sum(1))
        Pu = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
        Tu = truth / np.maximum(np.linalg.norm(truth, axis=1, keepdims=True), 1e-300)
        S = Pu @ Tu.T
        own = np.diag(S).copy()
        n = S.shape[0]
        pdsp.append(((S < own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1)) / (n - 1))
        cosp.append(own)
        H = dl[:, nonpanel]
        a_ = float((H * H).sum(1).mean())
        Hc = H - H.mean(0, keepdims=True)
        b_ = float((Hc * Hc).sum(1).mean())
        AB[ctx] = dict(A=a_, B=b_, BoverA=b_ / a_)
    summ = {k: float(np.mean([out[c][k] for c in "ABC"])) for k in ("pds", "cos", "sse", "cos_agg", "sgn_w")}
    summ["per_frame"] = {c: {k: out[c][k] for k in ("pds", "cos", "sse", "cos_agg")} for c in "ABC"}
    summ["AB"] = AB
    arrays = dict(x=np.asarray(X), a=np.asarray(Aa), pds=np.asarray(pdsp), cos=np.asarray(cosp))
    return summ, arrays


def ka_check410(held, fold, op, truth_full):
    """KA-1: check410.py's effect-level numbers, recomputed here."""
    ref = json.loads((ROOT / "experiments/410-pertscale-revisit/work/check410.json").read_text())[held]
    pre = np.asarray(fold["pre"], float)
    scored, wcols = np.asarray(fold["scored"]), np.asarray(fold["wcols"])
    P = op.apply(pre, wcols, GAMMA, 0.0)[:, scored]
    T = truth_full[:, scored]
    live = (np.abs(T).sum(1) > 0) & (np.abs(P).sum(1) > 0)
    P, T = P[live], T[live]
    tg = [t for t, k in zip(fold["targets"], live) if k]
    key, carried, _ = shared_energy_key(tg, fold["sources"], op, wcols, cols=scored)
    pt, pp = (P * T).sum(1), (P * P).sum(1)
    agg = lambda w: (w * pt).sum() / np.sqrt((w ** 2 * pp).sum() * (T ** 2).sum())
    base = agg(np.ones(len(tg)))
    got = {}
    for lo, hi, tag in ((0.3, 3.0, "lo0.3_hi3.0"), (0.3, 1e9, "lo0.3_hi1000000000.0")):
        w, _ = weights(key, carried, lo, hi)
        g = float(agg(w) / base - 1)
        exp = ref["rules"][tag]["gain"]
        got[tag] = g
        print(f"[{held}] KA-1 check410 {tag:22s} gain {g:+.6f} (ref {exp:+.6f}, d {abs(g - exp):.1e})", flush=True)
        assert abs(g - exp) < 1e-9, f"KA-1 FAILED {held} {tag}"
    return got


def ka_169(held, fold, op, truth):
    """KA-3: 169's a23-stack sharpening rows (plolo169.json), challenge frames, q = 1.00 and 1.15."""
    ref = {(r["fold"], r["q"], r["space"]): r for r in
           json.loads((ROOT / "experiments/169-sharpening-six-member/work/plolo169.json").read_text())
           if r["arm"] == "a23"}
    cd4_old = cd4_for(fold["targets"])
    final, _ = effects_for(fold, CONFIGS["a23"], op, cd4_old)
    worst = 0.0
    for q in (1.00, 1.15):
        for ctx, (b, gate) in frames().items():
            dl = deliver(sharpen(final, q, gate, mode="row"), b, gate, "minnorm")
            e = endpoints(dl[:, fold["scored"]], truth)
            r = ref[(held, q, f"chal{ctx}/delivered")]
            d = max(abs(e["pds"] - r["pds"]), abs(e["sse"] - r["sse"]))
            worst = max(worst, d)
    print(f"[{held}] KA-3 169 sharpening a23 q 1.00/1.15 worst |d pds|,|d sse| {worst:.1e}", flush=True)
    assert worst < 1e-9, f"KA-3 FAILED {held}"
    return worst


def run_fold(held: str, args) -> dict:
    t0 = time.time()
    status = plolo.fold_cache_status(held)
    assert status["status"] in ("hit", "legacy"), f"fold cache not usable read-only: {status}"
    fold = build_fold(held)
    fold["scored"] = np.asarray(fold["scored"])
    wcols = np.asarray(fold["wcols"])
    truth_full = np.asarray(fold["truth"], dtype=np.float64)
    truth = truth_full[:, fold["scored"]]
    pre = np.asarray(fold["pre"], dtype=np.float64)
    targets = list(fold["targets"])
    op_path = (PROD_OPS / "whitener_union_pruned.npz" if FOLDS[held]["op_kind"] == "production"
               else Path(FOLDS[held]["op"]))
    op = Operator(op_path)
    res = dict(fold=held, n_targets=len(targets), ka={}, arms={})
    arrays = {}

    # ---------------- KA-1 (410) and KA-3 (169)
    res["ka"]["check410"] = ka_check410(held, fold, op, truth_full)
    res["ka"]["plolo169_worst"] = ka_169(held, fold, op, truth)

    # ---------------- production composition (228) and its gain
    cd4 = cd4_for(targets)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    z = np.load(ROOT / f"experiments/228-theta-gain-plolo/work/gain-{held}.npz", allow_pickle=True)
    live = np.asarray(z["live"], dtype=bool)
    greg = np.asarray(z["gain_reg"], dtype=np.float64)
    ref228 = {(r["arm"], r["a"]): r for r in
              json.loads((ROOT / f"experiments/228-theta-gain-plolo/work/arms228-{held}.json").read_text())}

    def compose(m_pre=None):
        whitened = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
        if m_pre is not None:
            whitened, _ = apply_gain(whitened, m_pre, wcols)
        f_cd4, _ = solve_fraction(whitened, cd4, CD4_F)
        f_u, _ = solve_fraction_carrying(whitened, h1, UNION_G)
        final, _ = E062.add_post_whiten_residual(whitened, merged(cd4, f_cd4, h1, f_u), 1.0, wcols)
        final[:, fold["panel_positions"]] = 0.0
        return final

    def record(name, bulk, sharp=None, note=""):
        s, a = score(bulk, fold, truth, sharp)
        s["note"] = note
        s["row_norm_q"] = np.quantile(np.linalg.norm(bulk, axis=1), [0.1, 0.5, 0.9]).tolist()
        res["arms"][name] = s
        arrays[name] = a
        print(f"[{held}] {name:28s} cos_agg {s['cos_agg']:.6f} pds {s['pds']:.6f} "
              f"B/A {np.mean([s['AB'][c]['BoverA'] for c in 'ABC']):.5f} ({time.time() - t0:.0f}s)", flush=True)
        return s

    # KA-2: 228's rows, t = 0 and gain_reg@t, before the residuals, NO key
    base = compose()
    s = record("base", base, note="228 base222 (t=0), KA-2")
    r = ref228[("BASE", 1.0)]
    assert abs(s["cos_agg"] - r["cos_agg"]) < 1e-12 and abs(s["pds"] - r["pds"]) < 1e-12, "KA-2 FAILED base"
    finals_gi = {}
    for t in TS:
        finals_gi[t] = compose(multiplier(greg, live, t))
        s = record(f"gi@{t:g}", finals_gi[t], note="228 gain_reg before residuals, KA-2")
        r = ref228[(f"gain_reg@{t:g}", 1.0)]
        d = max(abs(s["cos_agg"] - r["cos_agg"]), abs(s["pds"] - r["pds"]))
        print(f"[{held}] KA-2 gain_reg@{t:g} |d| {d:.1e}", flush=True)
        assert d < 1e-12, f"KA-2 FAILED gain_reg@{t}"
    res["ka"]["arms228"] = "PASS"

    # ---------------- the production key on the final row (cols = None, the builder's convention)
    key, carried, nd = shared_energy_key(targets, fold["sources"], op, wcols, cols=None)
    key_s, carried_s, _ = shared_energy_key(targets, fold["sources"], op, wcols, cols=fold["scored"])
    W = {"K3": weights(key, carried, 0.3, 3.0)[0], "Kinf": weights(key, carried, 0.3, 1e9)[0],
         "K3s": weights(key_s, carried_s, 0.3, 3.0)[0], "Kinfs": weights(key_s, carried_s, 0.3, 1e9)[0]}
    rng = np.random.default_rng(20260922)
    perm = rng.permutation(len(targets))
    res["key"] = dict(carried=int(carried.sum()), n_donor_hist={int(k): int((nd == k).sum()) for k in np.unique(nd)},
                      w_q={k: np.quantile(v, [0, .1, .25, .5, .75, .9, 1]).tolist() for k, v in W.items()})
    bulks = {}
    for k in ("K3", "Kinf"):
        bulks[k], fs = key_scale(base, W[k])
        res["key"][f"{k}_final_scale_q"] = np.quantile(fs, [0, .1, .25, .5, .75, .9, 1]).tolist()
        record(k, bulks[k], note=f"410 key {k} on the final row, energy-preserving")
    for k in ("K3s", "Kinfs"):
        record(k, key_scale(base, W[k])[0], note="key with inner products on scored cols (check410 convention)")
    record("Kinf_shuf", key_scale(base, W["Kinf"][perm])[0], note="CONTROL: same weights, permuted over targets")
    record("K3_shuf", key_scale(base, W["K3"][perm])[0], note="CONTROL: same weights, permuted over targets")

    # ---------------- (a) 228's gain on top of the key, (i) before residuals, (ii) on the final row
    for t in TS:
        record(f"gi@{t:g}+Kinf", key_scale(finals_gi[t], W["Kinf"])[0], note="(i) gain before residuals, key on final")
        m = multiplier(greg, live, t)
        record(f"Kinf+gf@{t:g}", gene_gain_final(bulks["Kinf"], m), note="(ii) key, then gain on the final row (410d)")
    record("gi@0.25+K3", key_scale(finals_gi[0.25], W["K3"])[0], note="(i) with the capped key")
    record("K3+gf@0.25", gene_gain_final(bulks["K3"], multiplier(greg, live, 0.25)), note="(ii) with the capped key")
    record("gf@0.25", gene_gain_final(base, multiplier(greg, live, 0.25)), note="final-row gain, no key")
    gperm = rng.permutation(len(greg))
    record("Kinf+gf_shuf@0.25", gene_gain_final(bulks["Kinf"], multiplier(greg[gperm], live[gperm], 0.25)),
           note="CONTROL: gene-permuted gain on the final row")

    # ---------------- (b) 169's sharpening on the bulk (h space, row renorm, per-frame gate), and lfc space
    record("sh@1.15", base, sharp=("h", 1.15), note="169 sharpening, no key")
    for q in QS:
        record(f"Kinf+sh@{q:g}", bulks["Kinf"], sharp=("h", q), note="169 sharpening on the key bulk")
    record("Kinf+shL@1.15", bulks["Kinf"], sharp=("L", 1.15), note="169 lfc-space sharpening on the key bulk")
    record("K3+sh@1.15", bulks["K3"], sharp=("h", 1.15), note="sharpening on the capped key bulk")
    record("Kinf+gf@0.25+sh@1.15", gene_gain_final(bulks["Kinf"], multiplier(greg, live, 0.25)), sharp=("h", 1.15),
           note="gain + sharpening on the key bulk")

    # ---------------- (c) global amplitude on the key bulk: A, B, cos_agg response
    for s_ in (0.97, 1.03):
        record(f"Kinf@x{s_:g}", bulks["Kinf"] * s_, note="global amplitude probe on the key bulk")

    np.savez_compressed(WORK / f"arrays414-{held}.npz", names=np.asarray(list(arrays)),
                        d=np.asarray((truth * truth).sum(1)),
                        **{f"{n}__{k}": v for n, a in arrays.items() for k, v in a.items()})
    (WORK / f"arms414-{held}.json").write_text(json.dumps(res, indent=1))
    print(f"[{held}] done in {time.time() - t0:.0f}s", flush=True)
    del op
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", nargs="*", default=["k562_gwps", "orion_hct116", "orion_hek293t", "kolf2"])
    args = ap.parse_args()
    for held in args.folds:
        if (WORK / f"arms414-{held}.json").exists():
            print(f"[{held}] exists, skipping", flush=True)
            continue
        run_fold(held, args)


if __name__ == "__main__":
    main()
