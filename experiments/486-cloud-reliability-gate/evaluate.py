"""Leave-one-line-out test of a learned per-(target, gene) reliability gate.

Lines are distinct cell lines. K562 is the merged gwps + essential profile, so the
held-out line is never a source. For each held-out line the gate is fit only on
leave-one-out folds of the remaining lines.

Rows
    equal        equal-weight mean of the other lines
    fixed        w = 0.15 + 0.85 * z^2/(z^2+8) * n_src/(n_src+4)
    learned      Nelder-Mead gate, nested
    perm_refit   same fit, but features permuted across genes (one permutation)
    shuffle      learned weights, source targets shuffled at test time
    fixed_ncells sensitivity: fixed formula with n := mean source cell count

Cosine excludes the queried target's own gene. The rank metric also drops every
evaluated target's gene (the pds panel exclusion).
"""

from __future__ import annotations

import ctypes
import gc
import json
import time
from pathlib import Path

import numpy as np

from gate import (
    FEATURE_NAMES,
    N_FEATURES,
    bootstrap_summary,
    fit_theta,
    fixed_weight,
    fixed_weight_ncells,
    gate_weight,
    normalised_rank,
    paired_bootstrap,
    prepare_task,
    score_weight,
)

PRIMARY = ("equal", "fixed", "learned", "perm_refit", "shuffle")
CONTRASTS = [
    ("learned", "equal"),
    ("learned", "fixed"),
    ("fixed", "equal"),
    ("perm_refit", "equal"),
    ("shuffle", "equal"),
    ("fixed_ncells", "equal"),
]


def load_line(path: Path, name: str) -> dict:
    z = np.load(path, allow_pickle=False)
    genes = np.array([str(g).upper() for g in z["genes"]])
    perts = np.array([str(p).upper() for p in z["perts"]])
    if len(set(genes.tolist())) != len(genes):
        raise RuntimeError(f"{name}: gene symbols collide after uppercasing")
    if len(set(perts.tolist())) != len(perts):
        raise RuntimeError(f"{name}: perturbation symbols collide after uppercasing")
    return {
        "name": name,
        "genes": genes,
        "perts": perts,
        "lfc": np.asarray(z["lfc"], dtype=np.float64),
        "var": np.asarray(z["var"], dtype=np.float64),
        "n_cells": np.asarray(z["n_cells"], dtype=np.float64),
        "ctrl_log": np.asarray(z["ctrl_log"], dtype=np.float64),
    }


def align(lines: list[dict]) -> dict:
    common = sorted(set.intersection(*[set(l["genes"].tolist()) for l in lines]))
    if len(common) < 1000:
        raise RuntimeError(f"only {len(common)} genes in the intersection")
    from collections import Counter

    cnt = Counter()
    for line in lines:
        cnt.update(line["perts"].tolist())
    perts = sorted(p for p, k in cnt.items() if k >= 3)
    if len(perts) < 50:
        raise RuntimeError(f"only {len(perts)} perturbations are in >= 3 lines")
    g_common = {g: i for i, g in enumerate(common)}
    p_common = {p: i for i, p in enumerate(perts)}
    L, P, G = len(lines), len(perts), len(common)
    lfc = np.full((L, P, G), np.nan, dtype=np.float64)
    var = np.full((L, P, G), np.nan, dtype=np.float64)
    n_cells = np.full((L, P), np.nan, dtype=np.float64)
    ctrl = np.full((L, G), np.nan, dtype=np.float64)
    for li, line in enumerate(lines):
        g_src = np.array([g_common.get(g, -1) for g in line["genes"].tolist()])
        p_src = {p: i for i, p in enumerate(line["perts"].tolist())}
        keep_g = np.flatnonzero(g_src >= 0)
        dest_g = g_src[keep_g]
        ctrl[li, dest_g] = line["ctrl_log"][keep_g]
        for p, pi in p_src.items():
            dest_p = p_common.get(p)
            if dest_p is None:
                continue
            lfc[li, dest_p, dest_g] = line["lfc"][pi, keep_g]
            var[li, dest_p, dest_g] = line["var"][pi, keep_g]
            n_cells[li, dest_p] = line["n_cells"][pi]
    own_col = np.array([g_common.get(p, -1) for p in perts], dtype=np.int64)
    return {
        "lines": [l["name"] for l in lines],
        "genes": np.array(common),
        "perts": np.array(perts),
        "lfc": lfc,
        "var": var,
        "n_cells": n_cells,
        "ctrl": ctrl,
        "own_col": own_col,
    }


def _sources_have(lfc: np.ndarray, ids: list[int]) -> np.ndarray:
    """(n_pert,) count of source lines with a finite response row."""
    if not ids:
        return np.zeros(lfc.shape[1], dtype=np.int64)
    finite = np.isfinite(lfc[ids]).any(axis=2)
    return finite.sum(axis=0)


def make_task(bundle: dict, source_ids: list[int], held_id: int, target_ok: np.ndarray) -> dict:
    return prepare_task(
        bundle["lfc"][source_ids],
        bundle["var"][source_ids],
        bundle["n_cells"][source_ids],
        bundle["ctrl"][source_ids],
        bundle["lfc"][held_id],
        bundle["own_col"],
        target_ok,
    )


def permute_features(task: dict, perm: np.ndarray) -> dict:
    out = dict(task)
    out["Z"] = np.ascontiguousarray(task["Z"][:, perm, :])
    return out


def permute_targets(task: dict, perm: np.ndarray) -> dict:
    """Move source profiles (and the features that describe them) onto the wrong targets.

    The held-out truth and its cosine mask stay put.
    """
    out = dict(task)
    out["mu"] = task["mu"][perm]
    out["Z"] = task["Z"][perm]
    out["z2"] = task["z2"][perm]
    out["n_src"] = task["n_src"][perm]
    out["n_cells_mean"] = task["n_cells_mean"][perm]
    # Eligibility stays with the held-out target, not the source row that moved.
    return out


def _subset(scored: dict, rank: np.ndarray, top1: np.ndarray) -> dict[str, np.ndarray]:
    ok = scored["ok"]
    return {
        "cos": scored["cos"][ok],
        "energy": scored["energy"][ok],
        "rank": rank,
        "top1": top1.astype(np.float64),
    }


def panel_keep(genes: np.ndarray, perts: np.ndarray, use: np.ndarray) -> np.ndarray:
    """True = gene stays. Drop every evaluated target's own gene."""
    keep = np.ones(len(genes), dtype=bool)
    gindex = {g: i for i, g in enumerate(genes.tolist())}
    for p in perts[use].tolist():
        j = gindex.get(p)
        if j is not None:
            keep[j] = False
    if keep.sum() < 100:
        raise RuntimeError("panel exclusion removed almost every gene")
    return keep


def point_from_subset(sub: dict) -> dict[str, float]:
    e = sub["energy"]
    c = sub["cos"]
    return {
        "uw": float(c.mean()),
        "ew": float(np.sum(e * c) / np.sum(e)),
        "rank": float(np.nanmean(sub["rank"])),
        "top1": float(np.mean(sub["top1"])),
    }


def run(seed: int = 486, B: int = 4000, n_restarts: int = 2, maxiter: int = 400, maxfev: int = 500) -> dict:
    root = Path(__file__).resolve().parent
    resp = root / "work" / "responses"
    names = ["k562", "rpe1", "hepg2", "jurkat"]
    lines = [load_line(resp / f"{n}.npz", n) for n in names]
    bundle = align(lines)
    del lines
    gc.collect()
    print(
        f"aligned lines={bundle['lines']} perts={len(bundle['perts'])} genes={len(bundle['genes'])}",
        flush=True,
    )
    out_dir = root / "work" / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    report = {
        "lines": bundle["lines"],
        "n_genes_common": int(len(bundle["genes"])),
        "n_perts_ge3": int(len(bundle["perts"])),
        "seed": seed,
        "bootstrap_draws": B,
        "n_restarts": n_restarts,
        "maxiter": maxiter,
        "maxfev": maxfev,
        "folds": [],
    }
    per_target = {}

    for held_i, held_name in enumerate(bundle["lines"]):
        t0 = time.time()
        train_ids = [i for i in range(len(bundle["lines"])) if i != held_i]
        # Inner tasks: each training line held out, predicted from the other training lines.
        inner = []
        for h in train_ids:
            sources = [i for i in train_ids if i != h]
            n_src = _sources_have(bundle["lfc"], sources)
            held_has = np.isfinite(bundle["lfc"][h]).any(axis=1)
            # "at least two others" inside the training set.
            ok = held_has & (n_src >= 2)
            inner.append(make_task(bundle, sources, h, ok))
        print(
            f"fold {held_name}: inner targets {[int(t['target_ok'].sum()) for t in inner]}",
            flush=True,
        )
        fit = fit_theta(inner, seed=seed, n_restarts=n_restarts, maxiter=maxiter, maxfev=maxfev)
        gc.collect()
        try:
            ctypes.CDLL("libc.so.6").malloc_trim(0)
        except OSError:
            pass
        float64_train = float(
            np.mean([score_weight(t, gate_weight(fit["theta"], t["Z"]))["ew"] for t in inner])
        )
        if abs(float64_train - fit["train_ew"]) > 1e-3:
            raise RuntimeError(
                f"{held_name}: float32 objective {fit['train_ew']} != float64 rescore {float64_train}"
            )
        print(
            f"  fit {held_name}: train_ew={fit['train_ew']:.4f} float64={float64_train:.4f} "
            f"nfev={fit['nfev']} coef={fit['coef']}",
            flush=True,
        )
        # Known answer: theta = 0 reproduces the equal-weight cosine on an inner task.
        base = score_weight(inner[0], np.ones_like(inner[0]["mu"]))
        zero = score_weight(inner[0], gate_weight(np.zeros(N_FEATURES + 1), inner[0]["Z"]))
        if not np.isclose(base["ew"], zero["ew"], atol=1e-6):
            raise RuntimeError(f"theta=0 cosine {zero['ew']} != equal {base['ew']}")

        # Feature-permutation null: one gene permutation, shared by every task in the fold.
        perm_g = rng.permutation(bundle["lfc"].shape[2])
        perm_tasks = [permute_features(t, perm_g) for t in inner]
        fit_perm = fit_theta(
            perm_tasks,
            seed=seed + 1,
            n_restarts=n_restarts,
            maxiter=maxiter,
            maxfev=maxfev,
        )
        del perm_tasks, inner
        gc.collect()
        try:
            ctypes.CDLL("libc.so.6").malloc_trim(0)
        except OSError:
            pass

        n_src = _sources_have(bundle["lfc"], train_ids)
        held_has = np.isfinite(bundle["lfc"][held_i]).any(axis=1)
        ok = held_has & (n_src >= 2)
        task = make_task(bundle, train_ids, held_i, ok)
        use = task["target_ok"]
        if int(use.sum()) < 30:
            raise RuntimeError(f"{held_name}: only {int(use.sum())} evaluation targets")

        w_fixed = fixed_weight(task["z2"], task["n_src"])
        w_ncells = fixed_weight_ncells(task["z2"], task["n_cells_mean"][:, None])
        w_learned = gate_weight(fit["theta"], task["Z"])
        task_perm = permute_features(task, perm_g)
        w_perm = gate_weight(fit_perm["theta"], task_perm["Z"])
        # Shuffle only among eligible targets so each one receives another eligible profile.
        perm_t = np.arange(task["mu"].shape[0])
        eligible = np.flatnonzero(use)
        perm_t[eligible] = rng.permutation(eligible)
        task_shuf = permute_targets(task, perm_t)
        w_shuf = gate_weight(fit["theta"], task_shuf["Z"])

        methods_w = {
            "equal": (task, np.ones_like(task["mu"])),
            "fixed": (task, w_fixed),
            "fixed_ncells": (task, w_ncells),
            "learned": (task, w_learned),
            "perm_refit": (task_perm, w_perm),
            "shuffle": (task_shuf, w_shuf),
        }
        scored_all = {name: score_weight(tk, w) for name, (tk, w) in methods_w.items()}
        # Paired targets: every method must have a finite cosine and positive truth energy.
        common = use.copy()
        for scored in scored_all.values():
            common &= scored["ok"] & np.isfinite(scored["cos"]) & (scored["energy"] > 0)
        if int(common.sum()) < 30:
            raise RuntimeError(f"{held_name}: only {int(common.sum())} paired evaluation targets")
        keep = panel_keep(bundle["genes"], bundle["perts"], common)
        subsets = {}
        points = {}
        ablation = []
        for name, (tk, w) in methods_w.items():
            scored = scored_all[name]
            rank, top1 = normalised_rank(scored["pred"][common], tk["truth"][common], keep)
            sub = {
                "cos": scored["cos"][common],
                "energy": scored["energy"][common],
                "rank": rank,
                "top1": top1.astype(np.float64),
            }
            if not np.isfinite(rank).all():
                raise RuntimeError(f"{held_name} {name}: rank has non-finite rows")
            subsets[name] = sub
            points[name] = point_from_subset(sub)
            print(
                f"  {held_name} {name}: n={len(sub['cos'])} uw={points[name]['uw']:.4f} "
                f"ew={points[name]['ew']:.4f} rank={points[name]['rank']:.4f} top1={points[name]['top1']:.4f}",
                flush=True,
            )

        # Ablation: zero one coefficient of the fitted gate, score the outer fold.
        for j, fname in enumerate(FEATURE_NAMES):
            theta = fit["theta"].copy()
            theta[j + 1] = 0.0
            scored = score_weight(task, gate_weight(theta, task["Z"]))
            e = scored["energy"][common]
            c = scored["cos"][common]
            ew = float(np.sum(e * c) / np.sum(e))
            ablation.append(
                {
                    "feature": fname,
                    "coef": float(fit["theta"][j + 1]),
                    "ew": ew,
                    "drop_vs_learned": float(points["learned"]["ew"] - ew),
                }
            )

        draws = paired_bootstrap(subsets, B=B, seed=seed)
        summary = bootstrap_summary(points, draws, CONTRASTS)
        fold = {
            "held_out": held_name,
            "train_lines": [bundle["lines"][i] for i in train_ids],
            "n_targets": int(len(subsets["equal"]["cos"])),
            "n_panel_genes_dropped": int((~keep).sum()),
            "n_genes_in_rank": int(keep.sum()),
            "fit": {
                "train_ew": fit["train_ew"],
                "train_ew_float64": float64_train,
                "nfev": fit["nfev"],
                "intercept": fit["intercept"],
                "coef": fit["coef"],
                "restarts": fit["restarts"],
            },
            "fit_perm": {
                "train_ew": fit_perm["train_ew"],
                "intercept": fit_perm["intercept"],
                "coef": fit_perm["coef"],
            },
            "metrics": summary,
            "ablation": ablation,
            "seconds": time.time() - t0,
        }
        report["folds"].append(fold)
        per_target[held_name] = {name: {k: v.astype(np.float32) for k, v in sub.items()} for name, sub in subsets.items()}
        (out_dir / f"fold_{held_name}.json").write_text(json.dumps(fold) + "\n")
        gc.collect()
        try:
            ctypes.CDLL("libc.so.6").malloc_trim(0)
        except OSError:
            pass
        print(f"fold {held_name} done in {fold['seconds']:.0f}s", flush=True)

    # How many held-out lines the learned gate beats equal and fixed on, by energy-weighted cosine.
    def beats(gain_key: str, metric: str = "ew") -> int:
        n = 0
        for fold in report["folds"]:
            g = fold["metrics"]["gains"][gain_key][metric]
            if g["ci95"][0] > 0:
                n += 1
        return n

    report["n_lines_learned_beats_equal_ew"] = beats("learned_minus_equal")
    report["n_lines_learned_beats_fixed_ew"] = beats("learned_minus_fixed")
    report["n_lines_fixed_beats_equal_ew"] = beats("fixed_minus_equal")
    report["n_lines"] = len(report["folds"])

    (out_dir / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    # Flat tables for the lead to re-read without the nested JSON.
    rows = ["held_out,method,n_targets,uw,ew,rank,top1,uw_lo,uw_hi,ew_lo,ew_hi,rank_lo,rank_hi,top1_lo,top1_hi"]
    for fold in report["folds"]:
        for method, block in fold["metrics"]["methods"].items():
            rows.append(
                ",".join(
                    [
                        fold["held_out"],
                        method,
                        str(fold["n_targets"]),
                        f"{block['uw']['point']:.6f}",
                        f"{block['ew']['point']:.6f}",
                        f"{block['rank']['point']:.6f}",
                        f"{block['top1']['point']:.6f}",
                        f"{block['uw']['ci95'][0]:.6f}",
                        f"{block['uw']['ci95'][1]:.6f}",
                        f"{block['ew']['ci95'][0]:.6f}",
                        f"{block['ew']['ci95'][1]:.6f}",
                        f"{block['rank']['ci95'][0]:.6f}",
                        f"{block['rank']['ci95'][1]:.6f}",
                        f"{block['top1']['ci95'][0]:.6f}",
                        f"{block['top1']['ci95'][1]:.6f}",
                    ]
                )
            )
    (out_dir / "scores_by_line.csv").write_text("\n".join(rows) + "\n")
    grows = ["held_out,contrast,metric,point,lo,hi,boot_sd,frac_positive"]
    for fold in report["folds"]:
        for contrast, metrics in fold["metrics"]["gains"].items():
            for metric, g in metrics.items():
                grows.append(
                    f"{fold['held_out']},{contrast},{metric},{g['point']:.6f},{g['ci95'][0]:.6f},"
                    f"{g['ci95'][1]:.6f},{g['boot_sd']:.6f},{g['frac_positive']:.4f}"
                )
    (out_dir / "gains.csv").write_text("\n".join(grows) + "\n")
    crows = ["held_out,term,value"]
    arows = ["held_out,feature,coef,ew,drop_vs_learned"]
    for fold in report["folds"]:
        crows.append(f"{fold['held_out']},intercept,{fold['fit']['intercept']:.6f}")
        for feat, val in fold["fit"]["coef"].items():
            crows.append(f"{fold['held_out']},{feat},{val:.6f}")
        for ab in fold["ablation"]:
            arows.append(
                f"{fold['held_out']},{ab['feature']},{ab['coef']:.6f},{ab['ew']:.6f},{ab['drop_vs_learned']:.6f}"
            )
    (out_dir / "coefficients.csv").write_text("\n".join(crows) + "\n")
    (out_dir / "ablation.csv").write_text("\n".join(arows) + "\n")
    np.savez_compressed(out_dir / "per_target.npz", **_flatten_per_target(per_target))
    print(json.dumps({k: report[k] for k in report if k != "folds"}, indent=2), flush=True)
    return report


def _flatten_per_target(per_target: dict) -> dict:
    flat = {}
    for line, methods in per_target.items():
        for method, arrays in methods.items():
            for key, val in arrays.items():
                flat[f"{line}__{method}__{key}"] = val
    return flat


if __name__ == "__main__":
    run()
