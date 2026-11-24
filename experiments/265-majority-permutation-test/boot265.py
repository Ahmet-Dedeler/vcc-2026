"""265 -- the proper statistical check of 263's `kolf2` majority result (and 264's cancer-fold construction, and 263's
two-donor cells): per-target additivity (KA), a paired target-level bootstrap, and a permutation test over drop
assignments.  ARITHMETIC ON SAVED ARRAYS ONLY -- nothing is staged, no fold is rebuilt, no operator is loaded.

Inputs (read-only):
  experiments/263-donor-selection-arbiters/work/{sweep263,pertarget263,picks263}-<fold>.*
  experiments/264-majority-rule-cancer-folds/work/{sweep264[m],pertarget264[m],picks264[m]}-<fold>.*
Outputs: work/boot265.json, work/boot265.log, work/tables265.md.

The statistic.  A target's `pds` (resid218._per_target, tag a23) is, per challenge frame A/B/C, the fraction of the OTHER
truths whose cosine to the rendered prediction row is below the row's cosine to its own truth (ties at half), averaged
over the three frames.  It is a function of the RENDERED prediction row t and the (fixed) truth matrix only.  Whether
the PIPELINE is row-local is a separate question, measured here (H-2, H-3): `plolo.Operator.apply` rescales every
whitened row by ONE scalar (the un-whitened median row norm over the un-preserved whitened median), `E062.
add_post_whiten_residual` by another, and `plolo.deliver` is nonlinear in the row amplitude -- so a different set of
chosen rows moves the un-preserved median, hence every row's amplitude, hence (through the render) every row's cosine
matrix row.  MEASURED OUTCOME: strict additivity FAILS (untouched rows move by up to 0.012 `pds` on `kolf2`), but the
failure cancels on the TOTALS to <= 2.6e-4 over the 13 norm-matched `kolf2` arms (0.04 of a null SD), so the permutation
test is run on the ADDITIVE PROXY (an assignment's total = the sum of its per-target `pds` read off the always-drop-j
arms), labelled as such, with the proxy error reported for every real arm and the p-value also given for the real total.

    PYTHONDONTWRITEBYTECODE=1 uv run python boot265.py
"""
from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
E263 = ROOT / "experiments/263-donor-selection-arbiters/work"
E264 = ROOT / "experiments/264-majority-rule-cancer-folds/work"
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
B_BOOT = 20_000
B_PERM = 100_000
SEED = 20260905
CHUNK = 10_000
TOL_KA = 1e-12
# the outside re-analysis (GPT, 2026-09-05) quoted in the mandate, reproduced as a cross-check
GPT = {"majority - incumbent": (+0.01209, -0.0034, +0.0280), "majority - mean(random drop)": (+0.0053, -0.0115, +0.0228)}

t0 = time.time()
LINES: list[str] = []
OUT: dict = dict(seed=SEED, b_boot=B_BOOT, b_perm=B_PERM, ka=[], failures=[])


def P(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def check(label: str, got: float, want: float, tol: float, *, fatal: bool = True) -> bool:
    """fatal=True: a regression check (a mismatch is a bug, the run fails).  fatal=False: a MEASURED hypothesis
    (the mandate's additivity test) -- its verdict is reported, not enforced."""
    d = abs(float(got) - float(want))
    ok = d <= tol
    OUT["ka"].append(dict(label=label, got=float(got), want=float(want), abs_diff=d, tol=tol, ok=bool(ok), fatal=fatal))
    if not ok and fatal:
        OUT["failures"].append(f"{label}: got {got!r} want {want!r} |d| {d:.3e}")
    P(f"  {'KA' if fatal else 'H '} {'PASS' if ok else 'FAIL'} {label}: |d| {d:.1e}")
    return ok


def load(work: Path, tag: str, fold: str):
    sw = json.loads((work / f"sweep{tag}-{fold}.json").read_text())
    rows = {r["arm"]: r for r in sw["rows"]}
    pt = {k: np.asarray(v, dtype=np.float64) for k, v in np.load(work / f"pertarget{tag}-{fold}.npz").items()}
    pk = {k: np.asarray(v, dtype=np.int64) for k, v in np.load(work / f"picks{tag}-{fold}.npz").items()}
    return sw["meta"], rows, pt, pk


# ------------------------------------------------------------------------------------------------ statistics
def boot(diff: np.ndarray, rng: np.random.Generator, B: int = B_BOOT) -> dict:
    """Paired target bootstrap of the fold-mean difference: resample the n targets with replacement (differences are
    exactly 0 on the targets an arm does not touch, which are part of the fold mean)."""
    n = diff.size
    idx = rng.integers(0, n, size=(B, n))
    means = diff[idx].mean(1)
    lo, hi = np.percentile(means, [2.5, 97.5])
    return dict(point=float(diff.mean()), lo=float(lo), hi=float(hi), sd=float(means.std(ddof=1)),
                frac_le0=float((means <= 0).mean()), n=int(n), n_nonzero=int((diff != 0).sum()))


def boot_folds(diffs: list[np.ndarray], rng: np.random.Generator, B: int = B_BOOT) -> dict:
    """The K-fold mean of fold means, each fold resampled within itself (stratified)."""
    acc = np.zeros(B)
    for d in diffs:
        n = d.size
        idx = rng.integers(0, n, size=(B, n))
        acc += d[idx].mean(1)
    means = acc / len(diffs)
    lo, hi = np.percentile(means, [2.5, 97.5])
    return dict(point=float(np.mean([d.mean() for d in diffs])), lo=float(lo), hi=float(hi), sd=float(means.std(ddof=1)),
                frac_le0=float((means <= 0).mean()), per_fold=[float(d.mean()) for d in diffs])


def perm_nulls(Pm: np.ndarray, d_obs: np.ndarray, rest: float, n: int, rng: np.random.Generator, B: int = B_PERM) -> dict:
    """Pm: (m, J) per-target pds by dropped donor on the m gated targets; d_obs: the arm's drop vector on them; rest:
    the sum of the base per-target pds over the n - m untouched targets.  Returns the two nulls (fold-mean pds)."""
    m, J = Pm.shape
    ar = np.arange(m)
    uni = np.empty(B)
    frq = np.empty(B)
    for s in range(0, B, CHUNK):
        e = min(B, s + CHUNK)
        b = e - s
        idx = rng.integers(0, J, size=(b, m))
        uni[s:e] = (Pm[ar[None, :], idx].sum(1) + rest) / n
        order = np.argsort(rng.random((b, m)), axis=1)              # b independent uniform permutations of the targets
        dp = d_obs[order]
        frq[s:e] = (Pm[ar[None, :], dp].sum(1) + rest) / n
    return dict(uniform=uni, freq=frq)


def position(null: np.ndarray, obs: float) -> dict:
    B = null.size
    ge = int((null >= obs).sum())
    le = int((null <= obs).sum())
    lo, hi = np.percentile(null, [2.5, 97.5])
    sd = float(null.std(ddof=1))
    return dict(obs=float(obs), p_upper=(ge + 1) / (B + 1), p_lower=(le + 1) / (B + 1), null_mean=float(null.mean()),
                null_sd=sd, null_lo=float(lo), null_hi=float(hi), null_max=float(null.max()), null_min=float(null.min()),
                z=float((obs - null.mean()) / sd) if sd > 0 else float("nan"), B=int(B))


def fmt_b(b: dict) -> str:
    return f"{b['point']:+.5f} [{b['lo']:+.5f}, {b['hi']:+.5f}] sd {b['sd']:.5f} P(<=0) {b['frac_le0']:.4f}"


def fmt_p(p: dict) -> str:
    return (f"obs {p['obs']:.6f}  null mean {p['null_mean']:.6f} sd {p['null_sd']:.6f} 95% [{p['null_lo']:.6f}, {p['null_hi']:.6f}] "
            f"max {p['null_max']:.6f}  p_upper {p['p_upper']:.5f}  z {p['z']:+.2f}")


# ------------------------------------------------------------------------------------------------ one fold
def analyse_fold(label: str, meta: dict, rows: dict, pt: dict, pk: dict, *, base_arm: str, drop_arms: list[str],
                 donors: list[str], compare_to: list[str], rng: np.random.Generator, extra_boot: list[tuple[str, str, str]]) -> dict:
    """base_arm: the arm whose rows the drop arms leave untouched (kolf2: cd4svd5; cancer folds: mean3 = (a))."""
    n = meta["n_targets"]
    res: dict = dict(label=label, n=n, donors=donors)
    P(f"\n=== {label}: n {n}, donors {donors}, gates {meta['gates']}")

    # KA-1: every saved per-target array averages to its JSON total
    worst = 0.0
    for arm, r in rows.items():
        k = f"{arm}|pds"
        if k in pt:
            worst = max(worst, abs(pt[k].mean() - r["pds"]))
    check(f"[{label}] KA-1 mean(per-target pds) == JSON pds, all {len(rows)} rows", worst, 0.0, TOL_KA)

    base = pt[f"{base_arm}|pds"]
    Pfull = np.stack([pt[f"{a}|pds"] for a in drop_arms], 1)                 # (n, J)
    gate_all = pk[drop_arms[0]] >= 0
    for j, a in enumerate(drop_arms):
        assert np.array_equal(pk[a] >= 0, gate_all) and np.all(pk[a][gate_all] == j), a
    m = int(gate_all.sum())
    # H-2 (measured): do the drop arms leave the untouched targets at the base value?  (row-locality of the pipeline)
    unt = {}
    for j, a in enumerate(drop_arms):
        dd = Pfull[~gate_all, j] - base[~gate_all]
        cc = pt[f"{a}|cos"][~gate_all] - pt[f"{base_arm}|cos"][~gate_all]
        unt[a] = dict(n_untouched=int((~gate_all).sum()), n_off=int((dd != 0).sum()), max_abs_pds=float(np.abs(dd).max()),
                      mean_pds=float(dd.mean()), sum_pds=float(dd.sum()), max_abs_own_cos=float(np.abs(cc).max()))
    res["untouched"] = unt
    worst = max(v["max_abs_pds"] for v in unt.values())
    check(f"[{label}] H-2 drop-j arms == {base_arm} on the {n - m} untouched targets (row-locality)", worst, 0.0, 0.0, fatal=False)
    for a, v in unt.items():
        P(f"     {a:28s} untouched {v['n_untouched']}: off {v['n_off']:2d}, max|d pds| {v['max_abs_pds']:.4f}, mean {v['mean_pds']:+.5f}, sum {v['sum_pds']:+.4f}; max|d own cos| {v['max_abs_own_cos']:.1e}")

    # KA-3: ADDITIVITY -- every drop-assignment arm equals its per-target choices read off the drop-j arms
    add = {}
    for arm in [a for a in pk if a not in drop_arms and f"{arm}|pds" in pt and not a.startswith("only:")]:
        d = pk[arm]
        if d.max() < 0:
            continue
        sel = d >= 0
        if not np.all(gate_all[sel]):
            continue                                    # (picks outside the three-donor gate: not a drop arm)
        pred = base.copy()
        pred[sel] = Pfull[np.flatnonzero(sel), d[sel]]
        diff = pt[f"{arm}|pds"] - pred
        add[arm] = dict(n_selected=int(sel.sum()), max_abs=float(np.abs(diff).max()), n_targets_off=int((diff != 0).sum()),
                        total_diff=float(diff.mean()), total_json=float(rows[arm]["pds"]), total_pred=float(pred.mean()))
    res["additivity"] = add
    worst_arm = max(add, key=lambda a: add[a]["max_abs"])
    for arm, v in add.items():
        P(f"  additivity {arm:26s} sel {v['n_selected']:3d}  max|d| per target {v['max_abs']:.1e}  targets off {v['n_targets_off']:3d}  "
          f"total(JSON) {v['total_json']:.6f} vs sum of choices {v['total_pred']:.6f}  |d| {abs(v['total_diff']):.1e}")
    check(f"[{label}] H-3 additivity, worst arm {worst_arm}, per-target max|d|", add[worst_arm]["max_abs"], 0.0, 0.0, fatal=False)
    additive = add[worst_arm]["max_abs"] == 0.0
    res["additive"] = bool(additive)
    # the size of the failure ON THE TOTALS, over the norm-matched arms (the ones that share the base's global fractions)
    nm_tot = {a: v["total_diff"] for a, v in add.items() if a.endswith(":nm")}
    res["proxy_total_error_nm"] = dict(max_abs=float(max(abs(x) for x in nm_tot.values())), per_arm=nm_tot)
    P(f"  proxy error on the TOTAL over the {len(nm_tot)} norm-matched arms: max |real - sum of choices| {res['proxy_total_error_nm']['max_abs']:.1e}")

    # ---------------------------------------------------------------- bootstrap
    bt = {}
    for arm, ref, name in extra_boot:
        if arm.startswith("mean("):
            seeds = [s.strip() for s in arm[5:-1].split(",")]
            x = np.mean([pt[f"{s}|pds"] for s in seeds], 0)
        else:
            x = pt[f"{arm}|pds"]
        if ref.startswith("mean("):
            seeds = [s.strip() for s in ref[5:-1].split(",")]
            y = np.mean([pt[f"{s}|pds"] for s in seeds], 0)
        else:
            y = pt[f"{ref}|pds"]
        bt[name] = boot(x - y, rng)
        P(f"  boot {name:52s} {fmt_b(bt[name])}")
    res["boot"] = bt

    # ---------------------------------------------------------------- permutation over drop assignments
    # Run on the ADDITIVE PROXY: the total of an assignment = sum over gated targets of the per-target pds read off the
    # always-drop-j arm + the base's sum over the untouched targets.  Where additivity holds this is the arm's total
    # exactly; where it does not, the proxy's error on the totals is measured above (res['proxy_total_error_nm']) and
    # every real arm on the same gate is positioned in the nulls by BOTH its real and its proxy total.
    mode = "exact" if additive else "additive proxy"
    pm = {}
    for arm in compare_to:
        d = pk[arm]
        sel = d >= 0
        sub = np.flatnonzero(sel)
        Pm = Pfull[sub]
        rest = float(base[~sel].sum())
        nulls = perm_nulls(Pm, d[sub], rest, n, rng)
        obs_real = float(pt[f"{arm}|pds"].mean())
        obs_proxy = float((Pm[np.arange(sub.size), d[sub]].sum() + rest) / n)
        counts = {donors[j]: int((d[sub] == j).sum()) for j in range(len(donors))}
        e_uniform = float((Pm.mean(1).sum() + rest) / n)
        e_freq = float((Pm @ (np.bincount(d[sub], minlength=len(donors)) / sub.size)).sum() / n + rest / n)
        entry = dict(mode=mode, n_selected=int(sub.size), drop_counts=counts, expected_uniform=e_uniform, expected_freq=e_freq,
                     obs_real=obs_real, obs_proxy=obs_proxy, proxy_error=obs_real - obs_proxy,
                     uniform=position(nulls["uniform"], obs_proxy), freq=position(nulls["freq"], obs_proxy),
                     uniform_real=position(nulls["uniform"], obs_real), freq_real=position(nulls["freq"], obs_real))
        # the other arms on the same gate, positioned in the same nulls (real and proxy totals)
        others = {}
        for o in pk:
            if o == arm or f"{o}|pds" not in pt or o.startswith("only:") or o in drop_arms:
                continue
            do = pk[o]
            if not np.array_equal(do >= 0, sel):
                continue
            oo = float(pt[f"{o}|pds"].mean())
            op_ = float((Pm[np.arange(sub.size), do[sub]].sum() + rest) / n)
            others[o] = dict(obs_real=oo, obs_proxy=op_, drop_counts={donors[j]: int((do[sub] == j).sum()) for j in range(len(donors))},
                             p_upper_uniform=position(nulls["uniform"], op_)["p_upper"], p_upper_freq=position(nulls["freq"], op_)["p_upper"],
                             p_upper_uniform_real=position(nulls["uniform"], oo)["p_upper"], p_upper_freq_real=position(nulls["freq"], oo)["p_upper"],
                             z_uniform=float((op_ - nulls["uniform"].mean()) / nulls["uniform"].std(ddof=1)),
                             z_freq=float((op_ - nulls["freq"].mean()) / nulls["freq"].std(ddof=1)))
        entry["others"] = others
        pm[arm] = entry
        P(f"  perm [{mode}] {arm} (sel {sub.size}; drops {counts}; E[uniform] {e_uniform:.6f}, E[freq] {e_freq:.6f}; "
          f"obs real {obs_real:.6f} proxy {obs_proxy:.6f} (error {obs_real - obs_proxy:+.1e}))")
        P(f"     uniform: {fmt_p(entry['uniform'])}   [p_upper with the REAL total {entry['uniform_real']['p_upper']:.5f}]")
        P(f"     freq   : {fmt_p(entry['freq'])}   [p_upper with the REAL total {entry['freq_real']['p_upper']:.5f}]")
        for o, v in others.items():
            P(f"     {o:28s} real {v['obs_real']:.6f} proxy {v['obs_proxy']:.6f}  p_upper uniform {v['p_upper_uniform']:.5f} (real {v['p_upper_uniform_real']:.5f}; z {v['z_uniform']:+.2f})  "
              f"freq {v['p_upper_freq']:.5f} (real {v['p_upper_freq_real']:.5f}; z {v['z_freq']:+.2f})")
    res["perm"] = pm
    return res


# ================================================================================================= tables
def write_tables() -> None:
    L: list[str] = ["# 265 tables (generated by boot265.py from work/boot265.json)", ""]

    def brow(name, b):
        return f"| {name} | {b['point']:+.5f} | [{b['lo']:+.5f}, {b['hi']:+.5f}] | {b['sd']:.5f} | {b['frac_le0']:.4f} |"

    k = OUT["kolf2"]
    L += ["## T1 kolf2 -- paired target bootstrap (20,000 resamples of the 273 targets; exact, real arm evaluations)", "",
          "| comparison | point | 95 % interval | boot sd | P(<= 0) |", "|---|---:|---|---:|---:|"]
    L += [brow(n, b) for n, b in k["boot"].items()]
    L += ["", "## T2 kolf2 -- permutation over drop assignments (100,000 draws per null; ADDITIVE PROXY, see H-3)", "",
          "| arm | gate n | drop counts (k562 / hct116 / hek293t) | null | null mean (Δ vs incumbent) | null sd | null 95 % band (Δ) | obs proxy (Δ) | obs real (Δ) | p_upper proxy | p_upper real | z |",
          "|---|---:|---|---|---:|---:|---|---:|---:|---:|---:|---:|"]
    inc = float(OUT["kolf2_incumbent"])
    for arm, e in k["perm"].items():
        c = e["drop_counts"]
        cs = " / ".join(str(v) for v in c.values())
        for nm_, key, keyr in (("uniform", "uniform", "uniform_real"), ("freq-preserving", "freq", "freq_real")):
            p = e[key]
            L.append(f"| {arm} | {e['n_selected']} | {cs} | {nm_} | {p['null_mean'] - inc:+.5f} | {p['null_sd']:.5f} | [{p['null_lo'] - inc:+.5f}, {p['null_hi'] - inc:+.5f}] | "
                     f"{e['obs_proxy'] - inc:+.5f} | {e['obs_real'] - inc:+.5f} | {p['p_upper']:.4f} | {e[keyr]['p_upper']:.4f} | {p['z']:+.2f} |")
    L += ["", "### T2b the other arms on the same gates, positioned in the majority's nulls (proxy total; p with the real total in brackets)", "",
          "| gate arm | other arm | drop counts | obs real (Δ) | p_upper uniform | p_upper freq | z uniform | z freq |", "|---|---|---|---:|---:|---:|---:|---:|"]
    for arm, e in k["perm"].items():
        for o, v in e["others"].items():
            cs = " / ".join(str(x) for x in v["drop_counts"].values())
            L.append(f"| {arm} | {o} | {cs} | {v['obs_real'] - inc:+.5f} | {v['p_upper_uniform']:.4f} [{v['p_upper_uniform_real']:.4f}] | {v['p_upper_freq']:.4f} [{v['p_upper_freq_real']:.4f}] | {v['z_uniform']:+.2f} | {v['z_freq']:+.2f} |")
    L += ["", "### T3 kolf2 -- the real random-drop seeds against the uniform null", "",
          "| gate | seed | drop counts | agreement with the majority's drop | with the oracle's drop | real Δ | proxy Δ | p_upper (uniform) |", "|---|---|---|---:|---:|---:|---:|---:|"]
    for gate, sd in k["seeds"].items():
        for s, v in sd["seeds"].items():
            L.append(f"| {gate} | {s} | {' / '.join(map(str, v['counts']))} | {v['agree_majority']:.3f} | {v['agree_oracle_wht']:.3f} | {v['real'] - inc:+.5f} | {v['proxy'] - inc:+.5f} | {v['p_upper_uniform']:.4f} |")
        L.append(f"| {gate} | P(mean of 3 null draws >= the seeds' mean) | | (chance agreement 0.333, sd {sd['agreement_sd_expected']:.3f}) | majority vs oracle {sd['majority_agree_oracle_wht']:.3f} | | | **{sd['p_mean_of_three_ge_obs']:.4f}** |")
    for key, title in (("cancer_264m", "T4 cancer folds, 264 MATCHED construction"), ("cancer_264", "T5 cancer folds, 264 LETTER construction (equal three-donor mean)")):
        C = OUT[key]
        L += ["", f"## {title} -- bootstrap per fold and 3-fold mean (stratified)", "",
              "| comparison | k562_gwps | orion_hct116 | orion_hek293t | 3-fold mean [95 %] | P(<= 0) |", "|---|---|---|---|---|---:|"]
        names = list(next(iter(C["folds"].values()))["boot"].keys())
        for nme in names:
            cells = []
            for fold in CANCER:
                b = C["folds"][fold]["boot"][nme]
                cells.append(f"{b['point']:+.4f} [{b['lo']:+.4f}, {b['hi']:+.4f}]")
            if nme in C["pooled"]:
                pb = C["pooled"][nme]
                pooled = f"**{pb['point']:+.5f}** [{pb['lo']:+.5f}, {pb['hi']:+.5f}]"
                pl = f"{pb['frac_le0']:.4f}"
            else:
                pooled, pl = "", ""
            L.append(f"| {nme} | " + " | ".join(cells) + f" | {pooled} | {pl} |")
        L += ["", f"### {title[:2]}b permutation per fold (proxy; Δ vs (a) = mean3)", "",
              "| fold | arm | drop counts (donor order as listed) | proxy total error max over :nm arms | null | null mean Δa | null sd | obs proxy Δa | obs real Δa | p_upper proxy | p_upper real | z | always-drop-kolf2 p_upper freq | oracle (wht) z freq |",
              "|---|---|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
        for fold in CANCER:
            r = C["folds"][fold]
            a0 = float(OUT[f"{key}_mean3"][fold])
            for arm, e in r["perm"].items():
                cs = " / ".join(f"{kk} {v}" for kk, v in e["drop_counts"].items())
                for nm_, kk, kr in (("uniform", "uniform", "uniform_real"), ("freq-preserving", "freq", "freq_real")):
                    p = e[kk]
                    dk = e["others"].get("drop:kolf2:all", {})
                    orc = e["others"].get("oracle_drop:wht:all:nm", {})
                    L.append(f"| {fold} | {arm} | {cs} | {r['proxy_total_error_nm']['max_abs']:.1e} | {nm_} | {p['null_mean'] - a0:+.5f} | {p['null_sd']:.5f} | {e['obs_proxy'] - a0:+.5f} | {e['obs_real'] - a0:+.5f} | "
                             f"{p['p_upper']:.4f} | {e[kr]['p_upper']:.4f} | {p['z']:+.2f} | {dk.get('p_upper_freq', float('nan')):.4f} | {orc.get('z_freq', float('nan')):+.2f} |")
    L += ["", "## T6 263's two-donor cells on the cancer folds -- paired bootstrap vs cd4svd5", "",
          "| cell | k562_gwps | orion_hct116 | orion_hek293t | 3-fold mean [95 %] | P(<= 0) |", "|---|---|---|---|---|---:|"]
    for cell, v in OUT["two_donor_263"].items():
        cells = [f"{v['per_fold'][f]['point']:+.4f} [{v['per_fold'][f]['lo']:+.4f}, {v['per_fold'][f]['hi']:+.4f}]" for f in CANCER]
        pb = v["pooled"]
        L.append(f"| {cell} | " + " | ".join(cells) + f" | **{pb['point']:+.5f}** [{pb['lo']:+.5f}, {pb['hi']:+.5f}] | {pb['frac_le0']:.4f} |")
    L += ["", "## T7 additivity (H-2 / H-3) per fold", "",
          "| fold | construction | untouched targets | max |Δ pds| on an untouched target | max |Δ own cosine| untouched | max per-target |Δ| among :nm arms | max TOTAL error among :nm arms | null sd (uniform, gate all) | error / null sd |",
          "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for key, cons in (("kolf2", "production donors"), ("cancer_264m", "matched"), ("cancer_264", "letter")):
        folds = {"kolf2": OUT["kolf2"]} if key == "kolf2" else OUT[key]["folds"]
        for fold, r in folds.items():
            u = r["untouched"]
            nu = next(iter(u.values()))["n_untouched"]
            mx = max(v["max_abs_pds"] for v in u.values())
            mc = max(v["max_abs_own_cos"] for v in u.values())
            mt = max(v["max_abs"] for a, v in r["additivity"].items() if a.endswith(":nm"))
            te = r["proxy_total_error_nm"]["max_abs"]
            sd = r["perm"]["majority:raw:all:nm"]["uniform"]["null_sd"]
            L.append(f"| {fold} | {cons} | {nu} | {mx:.4f} | {mc:.1e} | {mt:.3f} | {te:.1e} | {sd:.5f} | {te / sd:.2f} |")
    L += ["", f"Known-answer (regression) checks {OUT['ka_pass']}/{OUT['ka_total']}; measured hypotheses (H-2 row-locality, H-3 additivity) {OUT['h_pass']}/{OUT['h_total']} hold. "
          f"Seed {SEED}; {OUT['seconds']} s; peak RSS {OUT['peak_rss_gib']:.2f} GiB."]
    (WORK / "tables265.md").write_text("\n".join(L) + "\n")


# ================================================================================================= kolf2
def main() -> None:
    rng = np.random.default_rng(SEED)
    P(f"265 boot/permutation -- seed {SEED}, B_boot {B_BOOT}, B_perm {B_PERM}; numpy {np.__version__}")

    meta4, rows4, pt4, pk4 = load(E264, "264", "kolf2")
    meta3, rows3, pt3, pk3 = load(E263, "263", "kolf2")
    donors = meta4["donors"]
    P(f"\n### kolf2 (264's arrays; 263's for the cross-check): donors {donors}")
    # KA-0: 264's kolf2 rows are 263's (R2) -- per-target arrays bit-identical
    for a3, a4 in (("cd4svd5", "cd4svd5"), ("majority:raw:all:nm", "majority:raw:all:nm"), ("majority:raw:tau0:nm", "majority:raw:tau0:nm"),
                   ("majority:raw:all:rev:nm", "reverse:raw:all:nm")):
        check(f"[kolf2] KA-0 264 {a4} per-target pds == 263 {a3}", float(np.abs(pt4[f"{a4}|pds"] - pt3[f"{a3}|pds"]).max()), 0.0, 0.0)
    check("[kolf2] KA-0 263 majority:raw:all:nm drop vector == 264's", float(np.abs(pk3["majority:raw:all:nm"] - pk4["majority:raw:all:nm"]).max()), 0.0, 0.0)
    check("[kolf2] KA-0 263 majority Δpds == +0.012089 (R263 §7)", rows4["majority:raw:all:nm"]["pds"] - rows4["cd4svd5"]["pds"], 0.012089, 5e-7)

    drop_arms = [f"drop:{d}:all:nm" for d in donors]
    rnd = ["random:all:s0:nm", "random:all:s1:nm", "random:all:s2:nm"]
    rnd0 = ["random:tau0:s0:nm", "random:tau0:s1:nm", "random:tau0:s2:nm"]
    extra = [("majority:raw:all:nm", "cd4svd5", "majority - incumbent"),
             *[(f"majority:raw:all:nm", s, f"majority - {s}") for s in rnd],
             ("majority:raw:all:nm", f"mean({','.join(rnd)})", "majority - mean(random drop, 3 seeds)"),
             *[(f"majority:raw:all:nm", a, f"majority - {a}") for a in drop_arms],
             ("majority:raw:all:nm", "reverse:raw:all:nm", "majority - reverse"),
             *[(s, "cd4svd5", f"{s} - incumbent") for s in rnd],
             (f"mean({','.join(rnd)})", "cd4svd5", "mean(random drop) - incumbent"),
             *[(a, "cd4svd5", f"{a} - incumbent") for a in drop_arms],
             ("reverse:raw:all:nm", "cd4svd5", "reverse - incumbent"),
             ("majority:raw:tau0:nm", "cd4svd5", "majority tau0 - incumbent"),
             *[("majority:raw:tau0:nm", s, f"majority tau0 - {s}") for s in rnd0],
             ("majority:raw:tau0:nm", f"mean({','.join(rnd0)})", "majority tau0 - mean(random drop tau0)"),
             ("majority:wht:all:nm", "cd4svd5", "majority (whitened space) - incumbent"),
             ("oracle_drop:raw:all:nm", "cd4svd5", "ORACLE drop raw - incumbent"),
             ("oracle_drop:wht:all:nm", "cd4svd5", "ORACLE drop whitened - incumbent")]
    rk = analyse_fold("kolf2", meta4, rows4, pt4, pk4, base_arm="cd4svd5", drop_arms=drop_arms, donors=donors,
                      compare_to=["majority:raw:all:nm", "majority:raw:tau0:nm"], rng=rng, extra_boot=extra)
    OUT["kolf2"] = rk
    OUT["kolf2_incumbent"] = float(rows4["cd4svd5"]["pds"])
    # the bootstrap's own known answers: the point estimates are the sweep's deltas
    check("[kolf2] boot point majority - incumbent == JSON delta", rk["boot"]["majority - incumbent"]["point"],
          rows4["majority:raw:all:nm"]["pds"] - rows4["cd4svd5"]["pds"], 1e-12)
    check("[kolf2] boot point mean(random drop) - incumbent == +0.0068 (R264 §6, three seeds)", rk["boot"]["mean(random drop) - incumbent"]["point"], 0.006831, 5e-5)
    for name, (pt_g, lo_g, hi_g) in GPT.items():
        key = "majority - incumbent" if name.startswith("majority - incumbent") else "majority - mean(random drop, 3 seeds)"
        b = rk["boot"][key]
        P(f"  cross-check vs the outside re-analysis: {name}: mine {b['point']:+.5f} [{b['lo']:+.4f}, {b['hi']:+.4f}]  GPT {pt_g:+.5f} [{lo_g:+.4f}, {hi_g:+.4f}]")
        OUT.setdefault("gpt_crosscheck", {})[name] = dict(mine=(b["point"], b["lo"], b["hi"]), gpt=(pt_g, lo_g, hi_g))
    # the identity E[uniform null] == mean over j of the always-drop-j totals, corrected for the untouched-row coupling
    # (the proxy keeps the base's values on the 33 untouched targets; the real drop-j arms move them by sum_pds_j)
    e_id = float(np.mean([rows4[a]["pds"] - rk["untouched"][a]["sum_pds"] / meta4["n_targets"] for a in drop_arms]))
    check("[kolf2] E[uniform drop null] == mean of the three always-drop-j totals (untouched-row coupling removed)",
          rk["perm"]["majority:raw:all:nm"]["expected_uniform"], e_id, 1e-12)
    # ---- the three real random-drop seeds per gate: are they typical draws of the null?  (counts, agreement, joint p)
    seeds_diag = {}
    for gate, arm, seeds in (("all", "majority:raw:all:nm", rnd), ("tau0", "majority:raw:tau0:nm", rnd0)):
        d = pk4[arm]
        sel = d >= 0
        sub = np.flatnonzero(sel)
        Pm = np.stack([pt4[f"{a}|pds"] for a in drop_arms], 1)[sub]
        rest = float(pt4["cd4svd5|pds"][~sel].sum())
        n = meta4["n_targets"]
        nulls = perm_nulls(Pm, d[sub], rest, n, rng, B=B_PERM)
        uni = nulls["uniform"]
        orc = pk4["oracle_drop:wht:all:nm"]
        rows_ = {}
        prox = []
        for s in seeds:
            ds = pk4[s]
            p_ = float((Pm[np.arange(sub.size), ds[sub]].sum() + rest) / n)
            prox.append(p_)
            rows_[s] = dict(counts=np.bincount(ds[sub], minlength=3).tolist(), agree_majority=float((ds[sub] == d[sub]).mean()),
                            agree_oracle_wht=float((ds[sub] == orc[sub]).mean()), proxy=p_, real=float(pt4[f"{s}|pds"].mean()),
                            p_upper_uniform=float(((uni >= p_).sum() + 1) / (uni.size + 1)))
        trip = uni[: uni.size // 3 * 3].reshape(-1, 3).mean(1)
        obs = float(np.mean(prox))
        seeds_diag[gate] = dict(seeds=rows_, seeds_mean_proxy=obs, p_mean_of_three_ge_obs=float(((trip >= obs).sum() + 1) / (trip.size + 1)),
                                agreement_sd_expected=float(np.sqrt(sub.size / 9 * 2) / sub.size), majority_agree_oracle_wht=float((d[sub] == orc[sub]).mean()))
        P(f"  seeds [{gate}]: " + "; ".join(f"{s} counts {v['counts']} agree maj {v['agree_majority']:.3f} orc {v['agree_oracle_wht']:.3f} p_upper {v['p_upper_uniform']:.4f}" for s, v in rows_.items())
          + f"; P(mean of 3 null draws >= seeds' mean) {seeds_diag[gate]['p_mean_of_three_ge_obs']:.4f}; majority agrees with the oracle on {seeds_diag[gate]['majority_agree_oracle_wht']:.3f}")
    OUT["kolf2"]["seeds"] = seeds_diag
    check("[kolf2] uniform-null mean (MC, 100k) == that identity within 4 MC-se",
          rk["perm"]["majority:raw:all:nm"]["uniform"]["null_mean"], e_id, 4 * rk["perm"]["majority:raw:all:nm"]["uniform"]["null_sd"] / np.sqrt(B_PERM))
    check("[kolf2] freq-null mean (MC, 100k) == its closed form within 4 MC-se",
          rk["perm"]["majority:raw:all:nm"]["freq"]["null_mean"], rk["perm"]["majority:raw:all:nm"]["expected_freq"],
          4 * rk["perm"]["majority:raw:all:nm"]["freq"]["null_sd"] / np.sqrt(B_PERM))
    # the split of the majority's per-target gain by dropped donor (263 kolfctrl: 110 / 65 / 65)
    d = pk4["majority:raw:all:nm"]
    dd = pt4["majority:raw:all:nm|pds"] - pt4["cd4svd5|pds"]
    OUT["kolf2"]["by_dropped"] = {donors[j]: dict(n=int((d == j).sum()), sum=float(dd[d == j].sum()), mean=float(dd[d == j].mean())) for j in range(3)}
    P("  majority per-target Δpds by dropped donor: " + ", ".join(f"{k} n {v['n']} sum {v['sum']:+.3f}" for k, v in OUT["kolf2"]["by_dropped"].items()))
    check("[kolf2] majority drop counts == 263 kolfctrl 110 / 65 / 65", sum(abs(OUT["kolf2"]["by_dropped"][k]["n"] - w) for k, w in zip(donors, (110, 65, 65))), 0, 0)

    # ============================================================================================= cancer folds (264)
    for tag, cname in (("264m", "matched"), ("264", "letter (equal three-donor mean)")):
        P(f"\n### cancer folds, 264 {cname} construction (tag {tag})")
        folds = {}
        diffs_ba, diffs_bc, diffs_bk, diffs_brm, diffs_ac = [], [], [], [], []
        for fold in CANCER:
            meta, rows, pt, pk = load(E264, tag, fold)
            donors = meta["donors"]
            drop_arms = [f"drop:{d}:all:nm" for d in donors]
            rnd = ["random:all:s0:nm", "random:all:s1:nm", "random:all:s2:nm"]
            extra = [("majority:raw:all:nm", "mean3", "(b) - (a)"), ("majority:raw:all:nm", "cd4svd5", "(b) - (c)"),
                     ("mean3", "cd4svd5", "(a) - (c)"),
                     ("majority:raw:all:nm", f"mean({','.join(rnd)})", "(b) - mean(random drop)"),
                     *[("majority:raw:all:nm", s, f"(b) - {s}") for s in rnd],
                     ("majority:raw:all:nm", "drop:kolf2:all:nm", "(b) - always-drop-kolf2"),
                     ("drop:kolf2:all:nm", "mean3", "always-drop-kolf2 - (a)"),
                     ("majority:raw:all:nm", "reverse:raw:all:nm", "(b) - reverse"),
                     ("majority:raw:tau0:nm", "mean3", "(b) tau0 - (a)"), ("majority:raw:tau0:nm", "cd4svd5", "(b) tau0 - (c)"),
                     ("majority:wht:all:nm", "mean3", "majority whitened - (a)"),
                     ("oracle_drop:wht:all:nm", "mean3", "ORACLE drop whitened - (a)")]
            r = analyse_fold(f"{fold} [{tag}]", meta, rows, pt, pk, base_arm="mean3", drop_arms=drop_arms, donors=donors,
                             compare_to=["majority:raw:all:nm", "majority:raw:tau0:nm"], rng=rng, extra_boot=extra)
            folds[fold] = r
            diffs_ba.append(pt["majority:raw:all:nm|pds"] - pt["mean3|pds"])
            diffs_bc.append(pt["majority:raw:all:nm|pds"] - pt["cd4svd5|pds"])
            diffs_bk.append(pt["majority:raw:all:nm|pds"] - pt["drop:kolf2:all:nm|pds"])
            diffs_brm.append(pt["majority:raw:all:nm|pds"] - np.mean([pt[f"{s}|pds"] for s in rnd], 0))
            diffs_ac.append(pt["mean3|pds"] - pt["cd4svd5|pds"])
        pooled = {"(b) - (a)": boot_folds(diffs_ba, rng), "(b) - (c)": boot_folds(diffs_bc, rng), "(a) - (c)": boot_folds(diffs_ac, rng),
                  "(b) - always-drop-kolf2": boot_folds(diffs_bk, rng), "(b) - mean(random drop)": boot_folds(diffs_brm, rng)}
        for k, v in pooled.items():
            P(f"  3-fold mean {k:28s} {fmt_b(v)}  per fold {', '.join(f'{x:+.4f}' for x in v['per_fold'])}")
        OUT[f"cancer_{tag}"] = dict(folds=folds, pooled=pooled)
        OUT[f"cancer_{tag}_mean3"] = {fold: float(load(E264, tag, fold)[1]["mean3"]["pds"]) for fold in CANCER}
        want = {"264m": (-0.00509, -0.01073), "264": (0.03628, -0.01945)}[tag]
        check(f"[{tag}] 3-fold (b) - (a) == R264 {want[0]:+.5f}", pooled["(b) - (a)"]["point"], want[0], 1e-5)
        check(f"[{tag}] 3-fold (b) - (c) == R264 {want[1]:+.5f}", pooled["(b) - (c)"]["point"], want[1], 1e-5)

    # ============================================================================================= 263's two-donor cells
    P("\n### 263's two-donor cells on the cancer folds (paired bootstrap vs cd4svd5; 3-fold mean stratified)")
    cells = ["cd4d:raw:tau0:rev:nm", "union:raw:tau0:rev", "union:wht:tau0:nm", "union:raw:tau0:rev:nm", "basal:raw:tau0:nm",
             "random:tau0:s0:nm", "oracle_wht:tau0:nm", "oracle_wht:all:nm"]
    two = {}
    for cell in cells:
        diffs, per = [], {}
        for fold in CANCER:
            meta, rows, pt, pk = load(E263, "263", fold)
            d = pt[f"{cell}|pds"] - pt["cd4svd5|pds"]
            diffs.append(d)
            per[fold] = boot(d, rng)
        pooled = boot_folds(diffs, rng)
        two[cell] = dict(per_fold=per, pooled=pooled)
        P(f"  {cell:24s} 3-fold {fmt_b(pooled)}  | per fold " + " | ".join(f"{f[:5]} {fmt_b(per[f])}" for f in CANCER))
    OUT["two_donor_263"] = two
    check("[263] union:wht:tau0:nm 3-fold point == R263 -0.0004", two["union:wht:tau0:nm"]["pooled"]["point"], -0.000405, 5e-5)
    check("[263] oracle_wht:all:nm 3-fold point == R263 +0.0459", two["oracle_wht:all:nm"]["pooled"]["point"], 0.045860, 5e-5)

    # the null's calibration where it can be checked: every real random-drop seed's position in its fold's uniform null
    cal = []
    for key in ("cancer_264m", "cancer_264"):
        for fold, r in OUT[key]["folds"].items():
            for arm, e in r["perm"].items():
                for o, v in e["others"].items():
                    if o.startswith("random:") and o.endswith(":nm"):
                        cal.append(dict(construction=key, fold=fold, gate=arm, seed=o, p_upper_uniform_proxy=v["p_upper_uniform"], p_upper_uniform_real=v["p_upper_uniform_real"]))
    OUT["seed_calibration_cancer"] = cal
    pp = np.array([c["p_upper_uniform_proxy"] for c in cal if c["construction"] == "cancer_264m"])
    P(f"\n  null calibration on the cancer folds (matched): {pp.size} real random-drop seeds at p_upper mean {pp.mean():.3f}, min {pp.min():.3f}, max {pp.max():.3f}, "
      f"n <= 0.05: {(pp <= 0.05).sum()}, n >= 0.95: {(pp >= 0.95).sum()}")

    OUT["seconds"] = round(time.time() - t0, 1)
    OUT["peak_rss_gib"] = peak_gib()
    OUT["ka_pass"] = sum(1 for k in OUT["ka"] if k["ok"] and k["fatal"])
    OUT["ka_total"] = sum(1 for k in OUT["ka"] if k["fatal"])
    OUT["h_pass"] = sum(1 for k in OUT["ka"] if k["ok"] and not k["fatal"])
    OUT["h_total"] = sum(1 for k in OUT["ka"] if not k["fatal"])
    (WORK / "boot265.json").write_text(json.dumps(OUT, indent=1, default=float))
    write_tables()
    P(f"\nknown-answer {OUT['ka_pass']}/{OUT['ka_total']} (regression checks), measured hypotheses {OUT['h_pass']}/{OUT['h_total']} hold, "
      f"{len(OUT['failures'])} failures, {OUT['seconds']}s, peak RSS {OUT['peak_rss_gib']:.2f} GiB")
    for f in OUT["failures"]:
        P(f"  *** {f}")
    (WORK / "boot265.log").write_text("\n".join(LINES) + "\n")
    if OUT["failures"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
