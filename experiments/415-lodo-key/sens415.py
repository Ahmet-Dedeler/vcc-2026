"""415 sensitivity — does the LODO key care that the SHIPPED bulk row is not 228's composition?

The shipped row (410b meta, builder 40772dd1) differs from 228's `build_final` in four ways:
  (1) normrel cap -1.5: per-target donor weights max(1, |row_j|/median_j)^-1.5, row norm-matched to the uniform row
      (norms over the donor's whole official coverage, panel columns included), PM50 correction then added at the
      uniform K562 share;
  (2) the base keeps its panel columns (the knockdown coordinates), so the all-rows sizing medians include them;
  (3) the CD4 residual is SVD-5 deflated before it is sized (VCC_CD4_DEFLATE_K=5, 246/249);
  (4) CD4's partial-condition rows are added at g 0.5 (carrying sizing), merged with CD4 + H1 union, one preserve.
"PP" below is that composition. KA-PP (A/B/C, aborts on failure): the PP constants reproduce the 410b build meta
(before/after whitening medians, operator preserve, CD4/H1/partial fractions, post-add preserve, deflation energy and
U_k sha256).

Then, on each PLOLO fold, three keys are scored on the PP prediction: K410, LODO-V1 built on 228's composition (the
written artifact's construction) and LODO-V1 built on PP; and on A/B/C the artifact's key is compared with the PP key.

    uv run python experiments/415-lodo-key/sens415.py
"""
from __future__ import annotations

import hashlib
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
for p in ("experiments/410-pertscale-revisit", "experiments/162-retrodicting-validator",
          "experiments/218-independent-donors", str(HERE.relative_to(ROOT))):
    sys.path.insert(0, str(ROOT / p))
import plolo  # noqa: E402
from lodo415 import (ALPHA, GAMMA, LAM, CD4_F, UNION_G, UNION, NBOOT, PREP_NAME, RULES, agg, agg_parts,  # noqa: E402
                     compose, spearman, whiten_split)
from resid218 import source_residual  # noqa: E402
from compose218 import solve_fraction_carrying, union_rows  # noqa: E402

plolo.WORK = WORK
from plolo import (CD4_NAMES, CD4_ROOT, DONOR_PATHS, E062, E064, PM_PATH, PM_BLEND, FOLDS, OFFICIAL, PANEL, PANEL_SET, cd4_for,  # noqa: E402
                   solve_fraction)
from key410 import raw_whiten, shared_energy_key, weights  # noqa: E402
from size410 import donor_rows  # noqa: E402

PARTIAL_G = 0.5
DEFLATE_K = 5
NORMREL_EXP = -1.5
MAIN = ("k562_gwps", "orion_hct116", "orion_hek293t")
META = json.loads((ROOT / "experiments/410-pertscale-revisit/work/build_b/"
                          "a-spreadsheet-achieved-enlightenment.h5ad.meta.json").read_text())["effects"]


def array_sha256(a):
    return hashlib.sha256(np.ascontiguousarray(a, dtype="<f8").tobytes()).hexdigest()


# ------------------------------------------------------------------ donor rows WITH panel columns (normrel's norms)
def full_rows(name: str):
    cache = WORK / f"full-{name}.npz"
    if cache.exists():
        z = np.load(cache)
        return z["effects"], z["available"]
    bank = E064.load_perturbmatch(PM_PATH) if name == "perturbmatch" else E062.Bank.load(name, DONOR_PATHS[name])
    mp = E062.map_bank_to_axis(bank, OFFICIAL)
    mean_all = E062.stream_column_mean(bank, mp.all_source_columns, chunk_rows=256, max_rss_gib=8.0)
    raw_median, _, _ = E062.centered_statistics(bank, mp.all_source_columns, mean_all, need_gene_scale=False,
                                                chunk_rows=256, max_rss_gib=8.0)
    at = {int(s): i for i, s in enumerate(mp.all_source_columns.tolist())}
    umean = np.asarray([mean_all[at[int(s)]] for s in mp.unique_source_columns])
    delta = np.load(bank.delta_path, mmap_mode="r")
    eff = np.zeros((len(PANEL), len(OFFICIAL)))
    avail = np.zeros(len(PANEL), bool)
    for i, t in enumerate(PANEL):
        r = bank.pert_pos.get(t)
        if r is None:
            continue
        eff[i, mp.unique_target_positions] = (np.asarray(delta[r, mp.unique_source_columns], float) - umean) \
            * (E062.EFFECT_SCALE / raw_median)
        avail[i] = True
    np.savez(cache, effects=eff, available=avail)
    return eff, avail


def cd4_frozen_and_partial(targets):
    """cd4_for (production's `complete mean3 only` freeze) and extract_cd4_partial, from ONE prepare_cd4_mean call."""
    banks = [E062.Bank.load(n, CD4_ROOT / n, require_cd4_scale=True) for n in CD4_NAMES]
    prep = E062.prepare_cd4_mean(banks, targets, OFFICIAL, PANEL_SET, "gene_norm", chunk_rows=256, max_rss_gib=8.0)
    count = np.asarray(prep.target_condition_count, dtype=np.int64)
    complete = count == len(CD4_NAMES)
    eff = np.asarray(prep.mean_effect, float).copy()
    votes = np.asarray(prep.gene_votes, np.uint8).copy()
    eff[~complete] = 0.0
    votes[~complete] = 0
    frozen = E062.CD4Prepared(mode=prep.mode, mean_effect=eff, gene_votes=votes,
                              target_condition_count=count.astype(np.uint8).copy(),
                              condition_details=prep.condition_details)
    partial = (count > 0) & ~complete
    pe = np.asarray(prep.mean_effect, float).copy()
    pv = np.asarray(prep.gene_votes, np.uint8).copy()
    pe[~partial] = 0.0
    pv[~partial] = 0
    part = E062.CD4Prepared(mode="cd4_partial", mean_effect=pe, gene_votes=pv,
                            target_condition_count=np.where(partial, count, 0).astype(np.uint8), condition_details=[])
    return frozen, part


def deflate(cd4, k):
    """make_submission_pm_cd4.deflate_cd4_residual (basis=None), verbatim algebra."""
    eff = np.asarray(cd4.mean_effect, float).copy()
    sup = np.asarray(cd4.gene_votes) > 0
    eff[~sup] = 0.0
    rows = np.flatnonzero(np.asarray(cd4.target_condition_count) > 0)
    R = eff[rows]
    e0 = float((R * R).sum())
    Rc = R - R.mean(axis=0, keepdims=True)
    Rc[~sup[rows]] = 0.0
    w, U = np.linalg.eigh(Rc @ Rc.T)
    Uk = U[:, np.argsort(w)[::-1][:k]]
    top = float(np.sort(w)[::-1][:k].sum() / max(w.sum(), 1e-300))
    R = Rc - Uk @ (Uk.T @ Rc)
    R[~sup[rows]] = 0.0
    eff[rows] = R
    audit = dict(k=k, n_carried=int(rows.size), energy_removed_frac=float(1 - (R * R).sum() / e0),
                 top_k_energy_frac_of_centred=top, basis_u_k_sha256=array_sha256(Uk))
    return E062.CD4Prepared(mode=cd4.mode, mean_effect=eff, gene_votes=np.asarray(cd4.gene_votes).copy(),
                            target_condition_count=np.asarray(cd4.target_condition_count).copy(),
                            condition_details=cd4.condition_details), audit


class Donors:
    """Per-target rows for a target list: full rows (panel columns kept) and the PM50 correction for the K562 slot."""

    def __init__(self, targets, sources):
        idx = {t: i for i, t in enumerate(PANEL)}
        ii = np.asarray([idx[t] for t in targets])
        self.sources = list(sources)
        self.F, self.A = {}, {}
        for s in sources:
            f, a = full_rows(s)
            self.F[s], self.A[s] = f[ii], a[ii]
        zk = np.load(WORK / "prep-k562_gwps.npz")
        zp = np.load(WORK / "prep-k562_gwps_pm50.npz")
        self.pmcorr = {"e064": (zp["effects"] - zk["effects"])[ii]}   # 162/228: shared non-panel K562&PM columns
        # production (apply_perturbmatch_k562_replacement, 061's broad policy): every official gene PM observes,
        # K562's slot 0 where K562 lacks the gene, targets carried by BOTH K562 and PM; uniform K562 share applied in pre()
        fpm, apm = full_rows("perturbmatch")
        fk, ak = full_rows("k562_gwps")
        both = (apm & ak)[ii]
        obs = np.zeros(fpm.shape[1], bool)
        bank = E064.load_perturbmatch(PM_PATH)
        obs[E062.map_bank_to_axis(bank, OFFICIAL).unique_target_positions] = True
        corr = PM_BLEND * (fpm[ii] - fk[ii])
        corr[:, ~obs] = 0.0
        corr[~both] = 0.0
        self.pmcorr["prod"] = corr
        # KA: panel-zeroed full rows == prep415 rows (prepare_main_source)
        pos = {g: i for i, g in enumerate(OFFICIAL)}
        pp = np.asarray([pos[g] for g in PANEL])
        for s in sources:
            z = np.load(WORK / f"prep-{s}.npz")
            f = self.F[s].copy()
            f[:, pp] = 0.0
            assert np.array_equal(self.A[s], z["available"][ii])
            assert np.abs(f - z["effects"][ii]).max() < 1e-12, f"full rows != prepared rows for {s}"
        # normrel medians: over the targets each donor carries in this list (production: the panel)
        self.med = {s: float(np.median(np.linalg.norm(self.F[s][self.A[s]], axis=1))) for s in sources}

    def pre(self, subset, normrel: bool, fallback_rows=None, pm: str | None = None):
        pm = pm or ("prod" if normrel else "e064")
        n, G = next(iter(self.F.values())).shape
        out = np.zeros((n, G))
        for t in range(n):
            have = [s for s in subset if self.A[s][t]]
            if not have:
                if fallback_rows is not None:
                    out[t] = fallback_rows[t]
                continue
            rows = np.stack([self.F[s][t] for s in have])
            u = rows.mean(0)
            if normrel and len(have) > 1:
                w = np.asarray([max(1.0, float(np.linalg.norm(self.F[s][t])) / self.med[s]) ** NORMREL_EXP
                                for s in have])
                cand = (w[:, None] * rows).sum(0) / w.sum()
                nu, nc = float(np.linalg.norm(u)), float(np.linalg.norm(cand))
                if nu > 0 and nc > 0:
                    cand = cand * (nu / nc)
            else:
                cand = u
            if "k562_gwps" in have:
                cand = cand + self.pmcorr[pm][t] / len(have)
            out[t] = cand
        return out


def compose_pp(op, pre_full, cd4_def, partial, h1, wcols, panel_positions):
    """apply_v9_whitening + apply_frozen_cd4_residual with CD4 (all-rows, deflated) + H1 union + partial, one preserve."""
    whitened = op.apply(pre_full * ALPHA, wcols, GAMMA, LAM)
    f_cd4, _ = solve_fraction(whitened, cd4_def, CD4_F)
    merged = np.where(np.asarray(cd4_def.gene_votes) > 0, np.asarray(cd4_def.mean_effect), 0.0) * f_cd4
    support = np.asarray(cd4_def.gene_votes) > 0
    fr = {}
    for tag, res, g in (("h1_union", h1, UNION_G), ("cd4_partial", partial, PARTIAL_G)):
        f, _ = solve_fraction_carrying(whitened, res, g)
        v = np.asarray(res.gene_votes) > 0
        merged += np.where(v, np.asarray(res.mean_effect), 0.0) * f
        support |= v
        fr[tag] = float(f)
    added = E062.CD4Prepared(mode="merged", mean_effect=merged, gene_votes=support.astype(np.uint8),
                             target_condition_count=np.asarray(cd4_def.target_condition_count), condition_details=[])
    final, info = E062.add_post_whiten_residual(whitened, added, 1.0, wcols)
    radd = np.where(support, merged, 0.0)
    before_op = float(np.median(np.linalg.norm((pre_full * ALPHA)[:, wcols], axis=1)))
    c, _ = whiten_split(op, pre_full * ALPHA, wcols)
    const = dict(before_operator_median_norm=before_op, operator_preserve_scale=c,
                 before_all_median_norm=float(np.median(np.linalg.norm(pre_full * ALPHA, axis=1))),
                 after_all_median_norm=float(np.median(np.linalg.norm(whitened, axis=1))),
                 frac_cd4=float(f_cd4), frac_h1=fr["h1_union"], frac_partial=fr["cd4_partial"],
                 s=float(info["post_add_global_preserve_scale"]))
    unzeroed = final.copy()
    final[:, panel_positions] = 0.0
    return final, unzeroed, const, radd, c


def loo_rows(op, dn, sources, normrel, c, s_add, radd, wcols, panel_positions, cols, fallback, pm=None):
    out = {}
    for j in sources:
        pj = dn.pre([i for i in sources if i != j], normrel, fallback_rows=None, pm=pm)
        bj = pj * ALPHA
        bj[:, wcols] = c * raw_whiten(op, bj, wcols)[:, wcols]
        Pj = s_add * (bj + radd)
        Pj[:, panel_positions] = 0.0
        out[j] = Pj if cols is None else Pj[:, cols]
    return out


def lodo_key(num_rows, d_rows, pp, A, carried, sources):
    key = np.full(len(pp), np.nan)
    for p in np.flatnonzero(carried):
        key[p] = float(np.mean([float(num_rows[j][p] @ d_rows[j][p]) for j in sources if A[j][p]])) / pp[p]
    return key


def abc():
    """KA-PP against the 410b meta, and the artifact key vs the PP key on the real panel."""
    targets = list(PANEL)
    pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([pos[g] for g in targets], dtype=np.int64)
    wcols = np.asarray(np.load(ROOT / "experiments/162-retrodicting-validator/work/fold-kolf2.npz")["wcols"])
    op = plolo.Operator(ROOT / "data/donors/whitener_union_pruned.npz", name="union_pruned")
    dn = Donors(targets, MAIN)
    frozen, partial = cd4_frozen_and_partial(targets)
    ref = cd4_for(targets)
    assert np.array_equal(ref.mean_effect, frozen.mean_effect) and np.array_equal(ref.gene_votes, frozen.gene_votes)
    cd4_def, daud = deflate(frozen, DEFLATE_K)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    out = {}
    for label, normrel, pm in (("pp", True, "prod"), ("pp_e064pm", True, "e064"), ("pp_no_normrel", False, "prod")):
        pre = dn.pre(MAIN, normrel, pm=pm)
        final, unz, const, radd, c = compose_pp(op, pre, cd4_def, partial, h1, wcols, panel_positions)
        prod = dict(before_operator_median_norm=META["whitening"]["before_operator_median_norm"],
                    operator_preserve_scale=META["whitening"]["operator_preserve_scale"],
                    before_all_median_norm=META["whitening"]["before_all_median_norm"],
                    after_all_median_norm=META["whitening"]["after_all_median_norm"],
                    frac_cd4=META["cd4"]["fraction"],
                    frac_h1=META["h1_union"]["residual_frame"]["applied_fraction"],
                    frac_partial=META["cd4_partial"]["residual_frame"]["applied_fraction"],
                    s=META["cd4"]["post_add_preserve"]["post_add_global_preserve_scale"])
        rel = {k: const[k] / prod[k] - 1 for k in prod}
        out[label] = dict(const=const, rel_to_410b_meta=rel, final_sha256=array_sha256(unz),
                          meta_final_sha256=META["final_sha256"])
        print(f"[A/B/C {label}] " + " ".join(f"{k} {const[k]:.6f} ({rel[k]:+.1e})" for k in prod), flush=True)
        if label == "pp":
            out["deflation"] = daud
            print(f"[A/B/C] deflation {json.dumps(daud)}  meta {json.dumps(META['cd4']['deflation'])}", flush=True)
            assert abs(daud["energy_removed_frac"] - META["cd4"]["deflation"]["energy_removed_frac"]) < 1e-9
            assert daud["n_carried"] == META["cd4"]["deflation"]["n_carried"]
            worst = max(abs(v) for v in rel.values())
            out["KA_PP_max_rel"] = worst
            print(f"[A/B/C] KA-PP worst relative constant error vs 410b meta {worst:.2e}; final sha "
                  f"{'MATCHES' if out['pp']['final_sha256'] == META['final_sha256'] else 'differs from'} meta", flush=True)
            assert worst < 1e-5, "KA-PP failed: PP does not reproduce the shipped composition constants"
            pp = (final * final).sum(1)
            carried = sum(dn.A[s].astype(int) for s in MAIN) >= 2
            num = loo_rows(op, dn, MAIN, True, c, const["s"], radd, wcols, panel_positions, None, None)
            d = {s: donor_rows(s, targets)[0] for s in MAIN}
            kpp = lodo_key(num, d, pp, dn.A, carried, MAIN)
            art = np.load(WORK / "pert_scale_415_lo0.3_hi3.npz")
            a410 = np.load(ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.3_hi3.npz")
            kart, k410 = art["key"], a410["key"]
            res = dict(spearman_pp_vs_artifact=spearman(kpp[carried], kart[carried]),
                       spearman_pp_vs_410=spearman(kpp[carried], k410[carried]),
                       spearman_artifact_vs_410=spearman(kart[carried], k410[carried]))
            for lo, hi in RULES:
                wpp, _ = weights(kpp, carried, lo, hi)
                wa, _ = weights(np.where(carried, kart, np.nan), carried, lo, hi)
                w4, _ = weights(np.where(carried, k410, np.nan), carried, lo, hi)
                res[f"lo{lo}_hi{hi:g}"] = dict(
                    frac_w_changed_gt_25pct_pp_vs_artifact=float((np.abs(np.log(wpp / wa)) > np.log(1.25)).mean()),
                    frac_w_changed_gt_25pct_pp_vs_410=float((np.abs(np.log(wpp / w4)) > np.log(1.25)).mean()),
                    # the aggregate 'cosine' of the artifact weights with the PP-key weights on the PP row energies
                    weight_cos_pp_artifact=float((wpp * wa * pp).sum() / np.sqrt((wpp ** 2 * pp).sum() * (wa ** 2 * pp).sum())),
                    weight_cos_pp_410=float((wpp * w4 * pp).sum() / np.sqrt((wpp ** 2 * pp).sum() * (w4 ** 2 * pp).sum())))
            out["keys"] = res
            np.savez(WORK / "abc_pp_key.npz", targets=np.asarray(targets), key=np.nan_to_num(kpp, nan=0.0),
                     carried=carried)
            print(f"[A/B/C] keys: {json.dumps(res)}", flush=True)
    del op
    return out


def plolo_fold(held, rng):
    fold = FOLDS[held]
    z = np.load(ROOT / f"experiments/162-retrodicting-validator/work/fold-{held}.npz", allow_pickle=True)
    truth = np.asarray(z["truth"], float)
    scored = np.asarray(z["scored"])
    wcols = np.asarray(z["wcols"])
    targets = [str(t) for t in z["targets"]]
    sources = list(fold["sources"])
    pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([pos[g] for g in PANEL if g in pos], dtype=np.int64)
    op = plolo.Operator(fold["op"])
    dn = Donors(targets, sources)
    # 228 composition (the artifact's construction), exactly as lodo415
    pre228 = dn.pre(sources, False)
    pre228[:, panel_positions] = 0.0
    cached = np.asarray(z["pre"], float)
    nocov = sum(dn.A[s].astype(int) for s in sources) == 0
    pre228[nocov] = cached[nocov]                                   # build_baseline's fallback rows
    assert np.abs(pre228 - cached).max() / np.abs(cached).max() < 1e-5, "228 pre rebuild"
    frozen, partial = cd4_frozen_and_partial(targets)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    f228, c228, radd228, wh228 = compose(op, pre228, frozen, h1, wcols, panel_positions)
    cc228, _ = whiten_split(op, pre228 * ALPHA, wcols)
    # PP composition
    cd4_def, daud = deflate(frozen, DEFLATE_K)
    fb = pre228.copy()
    prepp = dn.pre(sources, True, fallback_rows=fb)
    fpp, _, cpp, raddpp, ccpp = compose_pp(op, prepp, cd4_def, partial, h1, wcols, panel_positions)

    P = fpp[:, scored]
    T = truth[:, scored]
    live = (np.abs(T).sum(1) > 0) & (np.abs(P).sum(1) > 0)
    tg = [t for t, k in zip(targets, live) if k]
    n = len(tg)
    P, T = P[live], T[live]
    pt, pp, tt = agg_parts(P, T)
    sstar = pt / pp
    base = agg(np.ones(n), pt, pp, tt)
    A = {s: dn.A[s][live] for s in sources}
    carried = sum(A[s].astype(int) for s in sources) >= 2
    d = {s: donor_rows(s, tg)[0][:, scored] for s in sources}
    keys = {}
    k410, car, _ = shared_energy_key(tg, sources, op, wcols, cols=scored)
    assert np.array_equal(car, carried)
    keys["K410"] = k410
    # LODO-V1 on 228's composition (its own P for the denominator), then scored on PP
    P228 = f228[live][:, scored]
    num228 = loo_rows(op, dn, sources, False, cc228, c228["s"], radd228, wcols, panel_positions, scored, None)
    num228 = {j: v[live] for j, v in num228.items()}
    keys["V1_228comp"] = lodo_key(num228, d, (P228 * P228).sum(1), A, carried, sources)
    numpp = loo_rows(op, dn, sources, True, ccpp, cpp["s"], raddpp, wcols, panel_positions, scored, None)
    numpp = {j: v[live] for j, v in numpp.items()}
    keys["V1_PP"] = lodo_key(numpp, d, pp, A, carried, sources)
    out = dict(n=n, carried=int(carried.sum()), base_agg=base, constants_pp=cpp, deflation=daud,
               oracle_gain=float(np.sqrt((np.maximum(pt, 0) ** 2 / pp).sum() / tt) / base - 1),
               cos_P228_vs_PP_median=float(np.median((P228 * P).sum(1) / np.sqrt((P228 ** 2).sum(1) * pp))),
               keys={})
    names = list(keys)
    for lo, hi in RULES:
        tag = f"lo{lo}_hi{hi:g}"
        W = {k: weights(np.nan_to_num(keys[k], nan=0.0), carried, lo, hi)[0] for k in names}
        g = {k: agg(W[k], pt, pp, tt) / base - 1 for k in names}
        bs = {k: [] for k in names}
        for _ in range(NBOOT):
            idx = rng.integers(0, n, n)
            b = agg(np.ones(n), pt, pp, tt, idx)
            for k in names:
                bs[k].append(agg(W[k], pt, pp, tt, idx) / b - 1)
        for k in names:
            e = out["keys"].setdefault(k, dict(spearman_sstar=spearman(keys[k][carried], sstar[carried])))
            e[tag] = dict(gain=g[k], ci=[float(np.quantile(bs[k], .025)), float(np.quantile(bs[k], .975))])
            if k != "K410":
                dd = np.asarray(bs[k]) - np.asarray(bs["K410"])
                e[tag].update(diff_vs_410=g[k] - g["K410"],
                              diff_ci=[float(np.quantile(dd, .025)), float(np.quantile(dd, .975))])
        dd = np.asarray(bs["V1_228comp"]) - np.asarray(bs["V1_PP"])
        out.setdefault("v1_228_minus_v1_pp", {})[tag] = dict(
            diff=g["V1_228comp"] - g["V1_PP"], ci=[float(np.quantile(dd, .025)), float(np.quantile(dd, .975))])
    out["spearman_V1_228_vs_V1_PP"] = spearman(keys["V1_228comp"][carried], keys["V1_PP"][carried])
    del op
    return out


def main():
    t0 = time.time()
    res = {"abc": abc()}
    rng = np.random.default_rng(4150)
    for held in FOLDS:
        r = plolo_fold(held, rng)
        res[held] = r
        print(f"[{held}] PP n {r['n']} base {r['base_agg']:.5f} oracle {r['oracle_gain']:+.3f} median cos(P228,PP) "
              f"{r['cos_P228_vs_PP_median']:.4f} rho(V1_228,V1_PP) {r['spearman_V1_228_vs_V1_PP']:+.3f}", flush=True)
        for k, e in r["keys"].items():
            line = f"   {k:12s} rho {e['spearman_sstar']:+.3f}"
            for lo, hi in RULES:
                tag = f"lo{lo}_hi{hi:g}"
                line += f" | {tag} {e[tag]['gain']:+.3f} [{e[tag]['ci'][0]:+.3f},{e[tag]['ci'][1]:+.3f}]"
                if "diff_vs_410" in e[tag]:
                    line += f" d410 {e[tag]['diff_vs_410']:+.4f} [{e[tag]['diff_ci'][0]:+.4f},{e[tag]['diff_ci'][1]:+.4f}]"
            print(line, flush=True)
        print(f"   V1_228 - V1_PP: {json.dumps(r['v1_228_minus_v1_pp'])}", flush=True)
        (WORK / "sens415.json").write_text(json.dumps(res, indent=1))
    print(f"done {time.time() - t0:.0f}s, peak RSS {E062.peak_rss_gib():.2f} GiB")


if __name__ == "__main__":
    main()
