"""500: reproduce bulk-channel (pds, mse) suspicions on synthetic inputs.

Runs production functions. Does not need data/challenge. Prints a JSON record to stdout
and writes work/repro.json.

    uv run python experiments/500-cloud-bulk-review/repro_bulk.py
"""
from __future__ import annotations

import importlib.util
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

# The builder reads VCC_* at import. Keep the import side-effect free, then pin the
# production key only for the apply_pert_scale call below.
for k in list(os.environ):
    if k.startswith("VCC_"):
        os.environ.pop(k)

TS = 50_000.0
G = 4096
N_CELLS = 400
N_PERT = 8
N_PANEL = 64
Q_SHARP = 1.75
SEED = 500


def load_builder():
    os.environ["VCC_PERT_SCALE"] = "lodo415_lo0.5_hi3"
    os.environ["VCC_PERT_SCALE_BALLAST"] = "1"
    path = ROOT / "scripts" / "make_submission_pm_cd4.py"
    name = "make_submission_pm_cd4"
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    if na == 0.0 or nb == 0.0:
        return float("nan")
    return float(a @ b / (na * nb))


def effect_of(profile: np.ndarray, b: np.ndarray, mask: np.ndarray) -> np.ndarray:
    return (np.log1p(TS * profile) - b)[mask]


def make_control(rng: np.random.Generator, g: int) -> tuple[np.ndarray, np.ndarray]:
    raw = rng.lognormal(mean=-1.2, sigma=1.8, size=g)
    p = raw / raw.sum()
    b = np.log1p(TS * p)
    return p, b


def make_effect(rng: np.random.Generator, g: int, panel: np.ndarray, target: int) -> np.ndarray:
    h = rng.normal(0.0, 0.12, size=g)
    h += rng.standard_t(df=3, size=g) * 0.04
    h[panel] *= 0.25
    h[target] = -1.6
    return h


def readback(profile: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.log1p(TS * profile) - b


def profile_from_gated(b: np.ndarray, gated: np.ndarray) -> np.ndarray:
    profile = np.clip(np.expm1(b + gated), 0.0, None) / TS
    return profile / profile.sum()


def main() -> None:
    from cell_eval2.metrics.delta import _across_pert_budget, _row_weights
    from cell_eval2.moments import jackknife_correction

    M = load_builder()
    rng = np.random.default_rng(SEED)
    out: dict = {"seed": SEED, "G": G, "n_pert": N_PERT, "sharpen_q": Q_SHARP}

    # ---------------------------------------------------------------- jackknife
    p_ctrl, b_ctrl = make_control(rng, G)
    libs = np.clip(rng.lognormal(np.log(20109), 0.35, size=N_CELLS), 3275, 52420).astype(np.int64)
    equal = np.vstack([rng.multinomial(int(libs[i]), p_ctrl) for i in range(N_CELLS)])
    k_demo = 16
    # ballast-shaped: first 384 at control depth, last 16 near the ballast floor, same profile
    ball_libs = libs.copy()
    ball_libs[-k_demo:] = M.BALLAST_LB
    ballast_block = np.vstack([rng.multinomial(int(ball_libs[i]), p_ctrl) for i in range(N_CELLS)])

    def both(block: np.ndarray) -> dict:
        own = M.block_jackknife(block)
        ce = float(jackknife_correction(sparse.csr_matrix(block), np.zeros(block.shape[0], dtype=np.intp), 1, TS)[0])
        return {"builder": own, "cell_eval2": ce, "rel": (own - ce) / ce if ce else float("nan")}

    ref400 = both(equal)
    bal = both(ballast_block)
    # column of exact zeros: both should ignore it
    zblock = equal.copy()
    zblock[:, 0] = 0
    zero_gene = both(zblock)
    out["jackknife"] = {
        "equal_depth_control": ref400,
        "ballast_shaped_control_profile": bal,
        "zero_gene_column": zero_gene,
        "Chat_over_C400": bal["cell_eval2"] / ref400["cell_eval2"],
        "note": "C400 is 400 multinomial cells from the same profile at control depth; Chat is the ballast-shaped block",
    }
    print("jackknife", json.dumps(out["jackknife"]), flush=True)

    # ------------------------------------------------------- sharpen x ballast
    panel = np.zeros(G, dtype=bool)
    panel[:N_PANEL] = True
    genes = np.array([f"g{i}" for i in range(G)])
    # pretend the first N_PANEL genes are the panel
    pert_names = [f"g{i}" for i in range(N_PERT)]

    rows = []
    intended = {}  # (w, i) -> effect in non-panel log space
    achieved = {}
    for i in range(N_PERT):
        h = make_effect(rng, G, panel, i)
        de = M.profile_from_effect(h, b_ctrl)
        # carriers: multinomial from the UNSCALED design, using the first n-K libraries is decided per w
        base_draws = [rng.multinomial(int(libs[c]), de) for c in range(N_CELLS)]
        base = np.vstack(base_draws)
        for w in (0.5, 1.0, 3.0):
            bulk = M.profile_from_effect(h * w, b_ctrl, sharpen_q=Q_SHARP, protect=panel)
            # K from the production sizer, then the carriers are the first n-K cells (builder contract)
            k, ratio_min = M.ballast_count(de, bulk, libs, N_CELLS)
            carriers = base[: N_CELLS - k]
            full, diag = M.ballast_rows(carriers, de, bulk, k, seed=1000 + i)
            tot = full.sum(0).astype(np.float64)
            got = tot / tot.sum()
            mask = ~panel
            eb = effect_of(bulk, b_ctrl, mask)
            ed = effect_of(de, b_ctrl, mask)
            ea = effect_of(got, b_ctrl, mask)
            # mse space: drop only the target gene
            mse_mask = np.ones(G, dtype=bool)
            mse_mask[i] = False
            move = eb - ed
            got_move = ea - ed
            denom = float(move @ move)
            retain = float(got_move @ move / denom) if denom > 0 else float("nan")
            # column sums must survive the jackknife boost
            col_ok = bool(np.array_equal(full.sum(0), tot.astype(np.int64)))
            rec = {
                "pert": i, "w": w, "K": k, "design_min_ratio": ratio_min,
                "clip_genes": diag["ballast_clip_genes"], "clip_mass": diag["ballast_clip_mass"],
                "trim": diag["ballast_trim"], "jk_rho": diag["ballast_jk_rho"],
                "jk_ratio": diag["ballast_jk_ratio"],
                "cos_pds_vs_bulk": cosine(ea, eb),
                "cos_pds_vs_de": cosine(ea, ed),
                "cos_pds_de_vs_bulk": cosine(ed, eb),
                "sharpen_retain": retain,
                "mse_sse_vs_bulk": float(np.sum((effect_of(got, b_ctrl, mse_mask) - effect_of(bulk, b_ctrl, mse_mask)) ** 2)),
                "mse_sse_vs_de": float(np.sum((effect_of(got, b_ctrl, mse_mask) - effect_of(de, b_ctrl, mse_mask)) ** 2)),
                "colsum_preserved": col_ok,
                "extra_counts_vs_alloc": int(tot.sum() - (carriers.sum() + k * M.BALLAST_LB)),
            }
            rows.append(rec)
            intended[(w, i)] = eb
            achieved[(w, i)] = ea
            print(f"ballast w={w} pert={i} K={k} clip_mass={diag['ballast_clip_mass']} "
                  f"cos_bulk={rec['cos_pds_vs_bulk']:.6f} retain={retain:.4f}", flush=True)

    # self-match in the pds space: does the achieved vector still rank its own intended bulk first?
    rank = {}
    for w in (0.5, 1.0, 3.0):
        flips = 0
        margins = []
        for i in range(N_PERT):
            a = achieved[(w, i)]
            sims = [cosine(a, intended[(w, j)]) for j in range(N_PERT)]
            winner = int(np.argmax(sims))
            own = sims[i]
            other = max(sims[j] for j in range(N_PERT) if j != i)
            margins.append(own - other)
            flips += int(winner != i)
        sub = [r for r in rows if r["w"] == w]
        rank[str(w)] = {
            "self_match_flips": flips,
            "margin_min": float(min(margins)),
            "margin_median": float(np.median(margins)),
            "cos_vs_bulk_min": float(min(r["cos_pds_vs_bulk"] for r in sub)),
            "cos_vs_bulk_median": float(np.median([r["cos_pds_vs_bulk"] for r in sub])),
            "retain_median": float(np.median([r["sharpen_retain"] for r in sub])),
            "clip_mass_max": int(max(r["clip_mass"] for r in sub)),
            "clip_genes_median": float(np.median([r["clip_genes"] for r in sub])),
            "K_median": float(np.median([r["K"] for r in sub])),
            "jk_ratio_min": float(min(r["jk_ratio"] for r in sub)),
        }
    out["ballast"] = {"per_row": rows, "summary": rank}
    print("ballast summary", json.dumps(rank), flush=True)

    # K sized from the first 384 libraries, not from the n-K carriers that will actually ship.
    # Pick a design ratio that sits just under the K=1 -> K=2 boundary on 384 cells and over it on 399.
    s384 = float(libs[: N_CELLS - M.BALLAST_KMAX].sum())
    # need = s * (1/r - 1) = BALLAST_LB  => r = 1 / (1 + LB/s)
    r_star = 1.0 / (1.0 + M.BALLAST_LB / s384)
    de_u = np.full(G, 1.0 / G)
    bulk_u = de_u.copy()
    bulk_u[0] = de_u[0] * r_star * 0.999  # just below the K=2 edge on the 384-cell total
    bulk_u /= bulk_u.sum()
    k_edge, ratio_edge = M.ballast_count(de_u, bulk_u, libs, N_CELLS)
    s_actual = float(libs[: N_CELLS - k_edge].sum())
    need_actual = s_actual * (1.0 / max(ratio_edge, 1e-6) - 1.0)
    k_fixed = int(min(M.BALLAST_KMAX, max(1, math.ceil(need_actual / M.BALLAST_LB))))
    out["k_sizer"] = {
        "s384": s384, "design_min_ratio": ratio_edge, "K_returned": k_edge,
        "carrier_total_at_that_K": s_actual, "need_at_that_total": need_actual,
        "K_fixed_point": k_fixed,
        "short_by_one": bool(k_fixed > k_edge),
    }
    print("k_sizer", json.dumps(out["k_sizer"]), flush=True)

    # --------------------------------------------------------------- read-back
    rb = []
    for i in range(N_PERT):
        h = make_effect(np.random.default_rng(SEED + 17 + i), G, panel, i)
        gate = (1e6 * np.expm1(b_ctrl) / np.expm1(b_ctrl).sum() > M.MIN_CPM) & (b_ctrl > 0)
        for w in (0.5, 1.0, 3.0):
            gated = np.where(gate, h * w, 0.0)
            sharp = M.bulk_sharpen(gated, Q_SHARP, gate, panel)
            prof = profile_from_gated(b_ctrl, sharp)
            got = readback(prof, b_ctrl)
            on = gate & ~panel
            fixed = M.mass_consistent_effect(sharp, b_ctrl, gate, mode="minnorm")
            prof_f = profile_from_gated(b_ctrl, fixed)
            got_f = readback(prof_f, b_ctrl)
            clipped = int(np.sum(on & (b_ctrl + sharp <= 0)))
            rb.append({
                "pert": i, "w": w,
                "l2_ratio_sharp_over_pre": float(np.linalg.norm(sharp[on]) / max(np.linalg.norm(gated[on]), 1e-300)),
                "cos_readback_vs_designed": cosine(got[on], sharp[on]),
                "cos_readback_after_massfix": cosine(got_f[on], fixed[on]),
                "n_clipped_nonpanel": clipped,
                "rel_l2_readback_over_designed": float(np.linalg.norm(got[on]) / max(np.linalg.norm(sharp[on]), 1e-300)),
            })
    out["readback"] = rb
    print("readback", json.dumps({
        "cos_min": min(r["cos_readback_vs_designed"] for r in rb),
        "cos_median": float(np.median([r["cos_readback_vs_designed"] for r in rb])),
        "cos_massfix_min": min(r["cos_readback_after_massfix"] for r in rb),
        "l2_restore_max_abs_err": max(abs(r["l2_ratio_sharp_over_pre"] - 1.0) for r in rb),
    }), flush=True)

    # ------------------------------------------------- renderer exponent vs rho*
    # Synthetic B(f) through profile_from_effect + the scorer's own budget.
    effects = np.vstack([make_effect(rng, G, panel, i % N_PANEL) for i in range(N_PERT)])
    scales = [0.85, 1.0, 1.15, 1.30]
    budgets = []
    mus = []
    for f in scales:
        profs = np.vstack([
            M.profile_from_effect(effects[i] * f, b_ctrl, sharpen_q=Q_SHARP, protect=panel)
            for i in range(N_PERT)
        ])
        mu = np.log1p(TS * profs)
        mus.append(mu)
        cols = {i: i % N_PANEL for i in range(N_PERT)}  # each row drops a panel gene
        wts = _row_weights(N_PERT, G, cols)
        budgets.append(float(_across_pert_budget(mu, cols, wts)))
    exponents = []
    for a, b in zip(scales, scales[1:]):
        ia, ib = scales.index(a), scales.index(b)
        e = math.log(budgets[ib] / budgets[ia]) / math.log(b / a)
        exponents.append({"f0": a, "f1": b, "e": e, "B0": budgets[ia], "B1": budgets[ib]})
    q = 0.31 / (0.31 + 1.78)
    z = float(norm.ppf(q))
    sd = 0.01
    rho_star_code = math.exp(2.0 * sd * z)

    def expected_raw(mu: float) -> float:
        xs = np.linspace(mu - 10 * sd, mu + 10 * sd, 8001)
        pdf = norm.pdf(xs, loc=mu, scale=sd)
        loss = np.where(xs < 0.0, 0.31 * (-xs), 1.78 * xs)
        return float(np.trapezoid(loss * pdf, xs))

    mu_opt = sd * z  # ln(s / s_cliff), slopes are raw mse per unit (ln) scale
    loss_opt = expected_raw(mu_opt)
    # 448's measured local exponents (RESULTS.md): D/E/F between f=1 and lvl1, and F near the top
    measured = {"D_1_to_lvl1": 1.89, "E_1_to_lvl1": 1.90, "F_1_to_lvl1": 2.30, "F_local_top": 2.71}
    hedge = {}
    for name, e in measured.items():
        mu_code = (2.0 / e) * sd * z  # where relevel_iter lands: rho*=exp(2*sd*z), rho∝f^e
        rho_ev = math.exp(e * sd * z)
        hedge[name] = {
            "e": e,
            "rho_star_in_relevel_iter": rho_star_code,
            "rho_ev": rho_ev,
            "ln_scale_code": mu_code,
            "ln_scale_ev": mu_opt,
            "scale_ratio_code_over_ev": math.exp(mu_code - mu_opt),
            "excess_expected_raw_mse": expected_raw(mu_code) - loss_opt,
        }
    out["relevel"] = {
        "synthetic_exponents": exponents,
        "synthetic_B": budgets,
        "q": q, "z": z, "rho_star_hardcoded": rho_star_code,
        "hedge_at_448_exponents": hedge,
        "slope_unit": "raw mse per unit ln-scale, the unit ntclevel417.py states for the 0.31/1.78 slopes",
    }
    print("relevel", json.dumps({"exponents": exponents, "rho_star": rho_star_code, "hedge": hedge}), flush=True)

    # ------------------------------------------------------------- pert-scale energy
    zpath = ROOT / "experiments/420-bulk-stack/work/pert_scale_415pp_lo0.5_hi3.npz"
    zz = np.load(zpath, allow_pickle=False)
    targets = [str(t) for t in zz["targets"]]
    wkey = np.asarray(zz["w"], dtype=np.float64)
    genes_ps = targets + [f"X{i}" for i in range(1500)]
    E = rng.normal(0.0, 0.1, size=(len(targets), len(genes_ps)))
    for i, t in enumerate(targets):
        E[i, i] = -1.5
    scaled, audit = M.apply_pert_scale(E, targets, genes_ps)
    panel_set = set(targets)
    nonpanel = np.array([g not in panel_set for g in genes_ps])
    frob_before = float(np.linalg.norm(E[:, nonpanel]))
    frob_after = float(np.linalg.norm(scaled[:, nonpanel]))
    # mse keeps the other panel genes: drop only column i for row i, compare row energy ratios
    mse_ratios = []
    for i in range(0, len(targets), 30):
        m = np.ones(len(genes_ps), dtype=bool)
        m[i] = False
        mse_ratios.append(float(np.linalg.norm(scaled[i, m]) / np.linalg.norm(E[i, m])))
    out["pert_scale"] = {
        "key": "lodo415_lo0.5_hi3",
        "sha256": audit["sha256"],
        "energy_constant": audit["energy_constant"],
        "w_min": float(wkey.min()), "w_max": float(wkey.max()),
        "frac_at_floor_0.5": float(np.mean(np.isclose(wkey, 0.5))),
        "frac_at_cap_3": float(np.mean(np.isclose(wkey, 3.0))),
        "nonpanel_frobenius_ratio": frob_after / frob_before,
        "mse_subspace_row_energy_ratio_sample": mse_ratios,
    }
    print("pert_scale", json.dumps({k: out["pert_scale"][k] for k in
                                    ("energy_constant", "nonpanel_frobenius_ratio", "frac_at_floor_0.5", "frac_at_cap_3",
                                     "w_min", "w_max")}), flush=True)

    # ---------------------------------------------------------------- zero repair
    from vcc2026.zero_repair import zero_repair
    block = base.copy()
    gated = (1e6 * p_ctrl > 5.0)
    # force one gated gene to zero totals
    gidx = np.flatnonzero(gated)
    victim = int(gidx[len(gidx) // 2])
    block[:, victim] = 0
    cpm = 1e6 * p_ctrl[gidx]
    repaired, diag = zero_repair(block, n_car=block.shape[0], gated=gated, ctrl_cpm_gated=cpm,
                                 target_col=-1, rng=np.random.default_rng(1), floor_log2=-5.0)
    delta = repaired.sum(0).astype(np.float64) - block.sum(0)
    moved_genes = np.flatnonzero(delta != 0)
    out["zero_repair"] = {
        "diag": diag,
        "genes_whose_total_changed": int(moved_genes.size),
        "l1_count_change": float(np.abs(delta).sum()),
        "max_abs_count_change": float(np.abs(delta).max()),
        "block_total": int(block.sum()),
    }
    print("zero_repair", json.dumps(out["zero_repair"]), flush=True)

    dest = Path(__file__).resolve().parent / "work" / "repro.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=1))
    print(f"wrote {dest}", flush=True)


if __name__ == "__main__":
    main()
