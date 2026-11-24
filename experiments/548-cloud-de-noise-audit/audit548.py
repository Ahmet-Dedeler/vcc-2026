#!/usr/bin/env python3
"""Adversarial audit of the DE-channel noise figure (experiment 548).

Runs from the repo root on tracked files only:

    uv run python experiments/548-cloud-de-noise-audit/audit548.py

experiments/547-de-noise-check and experiments/537-arm-pipeline are not in
this checkout (origin/master at the commit this branch was cut from). This
script therefore does not replay 547's pair list. It rebuilds the quantities
the claim is about — rms and quadrature of DE-member deltas on pairs that
share a DE design, the 230 sd, and same-design seed dispersion — and it
records what each questionable inclusion does to the rms.

Outputs land in experiments/548-cloud-de-noise-audit/work/.
"""
from __future__ import annotations

import csv
import glob
import json
import math
import os
import re
from collections import defaultdict

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(ROOT, "experiments", "548-cloud-de-noise-audit", "work")

# Board scaling slopes for a DELTA. scaled = slope * (raw - b), so
# Δscaled = slope * Δraw. From docs/BOARD-LEDGER.md:162-163 (fitted on 612
# entries) and experiments/429-zs-production/RESULTS.md:46.
SLOPE = {
    "pds": 2.2138,
    "mse": -1.051,
    "nmae": -1.6532,
    "fid": 3.3645,
    "reach": 1.1292,
    "jac": 2.7263,
}
# 467's H1-like spans (RESULTS.md:43-44): overall-equiv = Σ Δraw / span / 6.
SPAN_467 = {
    "pds": 0.45,
    "mse": -0.955,
    "nmae": -0.60,
    "fid": 0.30,
    "reach": 0.90,
    "jac": 0.37,
}
MEMBERS = ("pds", "mse", "fid", "reach", "jac", "nmae")
DE_MEMBERS = ("fid", "reach", "jac", "nmae")
SCALED_KEY = {
    "pds": "score_pds",
    "mse": "score_mse",
    "fid": "score_fid",
    "reach": "score_reach",
    "jac": "score_jac",
    "nmae": "score_nmae",
}

# Config keys that change the DE channel's construction or its random draws.
# Bulk profile knobs (sharpen, CD4 dose, union g, normrel, whitener, alpha,
# pert_scale weights, bulk_load) are deliberately absent: a pair that matches
# on this fingerprint and differs only there is a bulk-only arm.
DE_PLAIN = (
    "emitter",
    "emission",
    "emission_seed",
    "overdispersion",
    "dispersion_argument",
    "spread",
    "spread_hi",
    "spread_lo",
    "spread_map",
    "tail_spread",
    "resample_x",
    "resample_invert",
    "resample_keep_top",
    "seed",
    "library_pick_seed",
    "rc_jitter_seed_base",
    "rc_callset_cpm",
    "rc_callset_head",
    "rc_callset_z_cap",
    "rc_key",
    "rc_sign_policy",
    "rc_scorer_min_cpm",
    "de_sign_prior",
    "de_sign_prior_sha256",
    "de_sign_protect_q",
    "de_sign_shrink_floor",
    "de_sign_classes",
    "de_sign_emitter_sha256",
    "de_sign_prior_context_A",
    "min_cpm",
    "gate_unexpressed",
    "keep_top",
    "keep_top_frac",
    "keep_spread",
    "keep_key",
    "mass_fix",
    "sign_repair",
    "carrier_gamma",
)
SEED_KEYS = ("emission_seed", "library_pick_seed", "seed", "rc_jitter_seed_base")
# Dropping any one of these is a "questionable inclusion" sensitivity.
SENSITIVITY_DROPS = (
    "sameline",
    "ctx_a_prior",
    "ballast",
    "de_sign_emitter_sha256",
    "zero_repair",
    "de_mag",
    "de_mde",
    "callset",
    "seeds",
)


def canon(obj):
    return json.dumps(obj, sort_keys=True, default=str)


def pick(d, keys):
    if not isinstance(d, dict):
        return d
    return {k: d.get(k) for k in keys}


def fingerprint(cfg, drop=()):
    """DE-design identity. Missing keys are None, which does not equal an
    explicit value. `drop` names components to omit (sensitivity)."""
    drop = set(drop)
    fp = {}
    if "seeds" not in drop:
        for k in DE_PLAIN:
            if k in SEED_KEYS:
                fp[k] = cfg.get(k, None)
    for k in DE_PLAIN:
        if k in SEED_KEYS:
            continue
        if k == "de_sign_emitter_sha256" and "de_sign_emitter_sha256" in drop:
            continue
        fp[k] = cfg.get(k, None)
    if "de_mag" not in drop:
        fp["de_mag"] = pick(cfg.get("de_mag"), ("variant", "module_sha256")) if cfg.get("de_mag") else None
    if "de_mde" not in drop:
        fp["de_mde"] = (
            pick(cfg.get("de_mde"), ("c", "cap", "coarse", "min_count", "scope"))
            if cfg.get("de_mde")
            else None
        )
    if "zero_repair" not in drop:
        zr = cfg.get("de_zero_repair")
        fp["de_zero_repair"] = pick(zr, ("floor_log2",)) if isinstance(zr, dict) else zr
    if "callset" in drop:
        for k in (
            "rc_callset_cpm",
            "rc_callset_head",
            "rc_callset_z_cap",
            "rc_key",
            "rc_sign_policy",
            "rc_scorer_min_cpm",
            "keep_top",
            "keep_top_frac",
            "keep_spread",
            "keep_key",
        ):
            fp.pop(k, None)
    if "ballast" not in drop:
        b = cfg.get("pert_scale_ballast")
        fp["pert_scale_ballast"] = pick(b, ("kmax", "lb")) if isinstance(b, dict) else b
    if "sameline" not in drop:
        s = cfg.get("sameline")
        if not s:
            fp["sameline"] = None
        else:
            ctx = {}
            for c, d in sorted((s.get("contexts") or {}).items()):
                ctx[c] = {
                    "programme_default": d.get("programme_default"),
                    "w": d.get("w"),
                    "pt": d.get("pt"),
                }
            fp["sameline"] = {
                "sha256": s.get("sha256"),
                "prog": s.get("prog"),
                "pt": s.get("pt"),
                "w": s.get("w"),
                "contexts": ctx,
            }
    g = cfg.get("de_prior_gate")
    if g and "sameline" not in drop:
        # q/tau/weak/strength are the knob. Per-context counts are outcomes.
        fp["de_prior_gate"] = pick(g, ("q", "tau", "weak", "strength"))
    elif "sameline" in drop:
        fp.pop("de_prior_gate", None)
    if "ctx_a_prior" not in drop:
        keys = (
            "ctx_a_prior_mode",
            "ctx_a_prior_beta",
            "ctx_a_prior_context",
            "ctx_a_prior_sha256",
            "ctx_a_prior_targets",
        )
        if any(k in cfg for k in keys):
            fp["ctx_a_prior"] = {k: cfg.get(k) for k in keys}
        else:
            fp["ctx_a_prior"] = None
    return fp


def flatten(obj, prefix=""):
    out = {}
    if isinstance(obj, dict):
        for k, v in obj.items():
            p = f"{prefix}.{k}" if prefix else k
            if isinstance(v, dict):
                out.update(flatten(v, p))
            else:
                out[p] = canon(v)
    return out


def rms(xs):
    if not xs:
        return None
    return math.sqrt(sum(x * x for x in xs) / len(xs))


def sample_sd(xs):
    if len(xs) < 2:
        return None
    m = sum(xs) / len(xs)
    var = sum((x - m) ** 2 for x in xs) / (len(xs) - 1)
    if var < 0.0:
        var = 0.0
    return math.sqrt(var)


def mean(xs):
    return sum(xs) / len(xs) if xs else None


def load_board():
    rows = []
    for path in sorted(glob.glob(os.path.join(ROOT, "data", "board-entries", "*.json"))):
        blob = json.load(open(path))
        e = blob["entry"]
        pkg = blob.get("package") or ""
        meta = None
        if pkg.endswith(".vcc"):
            cand = pkg.replace(".vcc", ".h5ad.meta.json")
            if os.path.exists(os.path.join(ROOT, cand)):
                meta = cand
        scores = {m: e[SCALED_KEY[m]] for m in MEMBERS}
        rows.append(
            {
                "entry_id": e["entry_id"],
                "date": e.get("submission_date"),
                "note": blob.get("note"),
                "model_name": e.get("model_name"),
                "overall": e.get("score_avg"),
                "package": pkg,
                "meta": meta,
                "sha256": blob.get("sha256"),
                "scores": scores,
                "capture": os.path.relpath(path, ROOT),
            }
        )
    return rows


def load_meta(rel):
    blob = json.load(open(os.path.join(ROOT, rel)))
    cfg = blob.get("config") or {}
    impl = ((blob.get("hashes") or {}).get("implementation") or {})
    return cfg, impl


def arm_name(package):
    base = os.path.basename(package)
    base = re.sub(r"\.vcc$", "", base)
    # short label: the token before the first hyphenated recipe, else the file
    m = re.match(r"^([A-Za-z0-9]+)", base)
    return m.group(1) if m else base


def pair_record(a, b):
    rec = {
        "a": a["label"],
        "b": b["label"],
        "a_entry": a["entry_id"],
        "b_entry": b["entry_id"],
        "a_date": a["date"],
        "b_date": b["date"],
        "a_meta": a["meta"],
        "b_meta": b["meta"],
    }
    deltas = {}
    for m in MEMBERS:
        deltas[m] = b["scores"][m] - a["scores"][m]
    rec["delta_scaled"] = deltas
    rec["delta_overall"] = b["overall"] - a["overall"]
    rec["delta_de_overall"] = sum(deltas[m] for m in DE_MEMBERS) / 6.0
    rec["profile_pinned"] = abs(deltas["pds"]) < 1e-9 and abs(deltas["mse"]) < 1e-9
    rec["de_abs_max"] = max(abs(deltas[m]) for m in DE_MEMBERS)
    # config diffs outside the fingerprint, at top-level keys
    da, db = a["flat"], b["flat"]
    keys = set(da) | set(db)
    differ = sorted(k for k in keys if da.get(k) != db.get(k))
    rec["n_config_paths_differ"] = len(differ)
    rec["config_paths_differ"] = differ
    top = sorted({k.split(".")[0] for k in differ})
    rec["config_top_differ"] = top
    rec["donor_transfer_sha_a"] = a["impl"].get("src/vcc2026/donor_transfer.py")
    rec["donor_transfer_sha_b"] = b["impl"].get("src/vcc2026/donor_transfer.py")
    rec["emitter_code_differs"] = rec["donor_transfer_sha_a"] != rec["donor_transfer_sha_b"]
    rec["builder_sha_a"] = a["impl"].get("scripts/make_submission_pm_cd4.py")
    rec["builder_sha_b"] = b["impl"].get("scripts/make_submission_pm_cd4.py")
    rec["builder_code_differs"] = rec["builder_sha_a"] != rec["builder_sha_b"]
    return rec


DEISH = re.compile(
    r"de_|rc_|emit|seed|sameline|ballast|callset|prior|mde|jitter|"
    r"dispersion|spread|overdisp|resample|sign_repair|carrier_gamma|"
    r"min_cpm|gate_unexpressed|keep_top|keep_spread|keep_key|mass_fix"
)


def summarise(pairs, tag):
    if not pairs:
        return {"tag": tag, "n": 0}
    de = [p["delta_de_overall"] for p in pairs]
    per = {m: [p["delta_scaled"][m] for p in pairs] for m in DE_MEMBERS}
    member_rms = {m: rms(per[m]) for m in DE_MEMBERS}
    quad = math.sqrt(sum(v * v for v in member_rms.values())) / 6.0
    ordering = [
        (p["delta_scaled"]["fid"] + p["delta_scaled"]["reach"] + p["delta_scaled"]["jac"]) / 6.0
        for p in pairs
    ]
    nmae_ov = [p["delta_scaled"]["nmae"] / 6.0 for p in pairs]
    pinned = [p for p in pairs if p["profile_pinned"]]
    moved = [p for p in pairs if not p["profile_pinned"]]
    return {
        "tag": tag,
        "n_pairs": len(pairs),
        "n_profile_pinned": len(pinned),
        "n_profile_moved": len(moved),
        "de_overall_rms": rms(de),
        "de_overall_mean": mean(de),
        "de_overall_sd": sample_sd(de),
        "de_overall_max_abs": max(abs(x) for x in de),
        "member_rms_scaled": member_rms,
        "quadrature_over_6": quad,
        "ordering_fid_reach_jac_rms": rms(ordering),
        "nmae_over_6_rms": rms(nmae_ov),
        "steps": [
            {
                "a": p["a"],
                "b": p["b"],
                "delta_de_overall": p["delta_de_overall"],
                "delta_ordering_overall": (
                    p["delta_scaled"]["fid"] + p["delta_scaled"]["reach"] + p["delta_scaled"]["jac"]
                )
                / 6.0,
                "delta_nmae_overall": p["delta_scaled"]["nmae"] / 6.0,
                "delta_scaled": p["delta_scaled"],
                "config_top_differ": p["config_top_differ"],
            }
            for p in pairs
        ],
        "pinned_de_overall_rms": rms([p["delta_de_overall"] for p in pinned]) if pinned else None,
        "moved_de_overall_rms": rms([p["delta_de_overall"] for p in moved]) if moved else None,
        "pinned_quadrature_over_6": (
            math.sqrt(
                sum(rms([p["delta_scaled"][m] for p in pinned]) ** 2 for m in DE_MEMBERS)
            )
            / 6.0
            if pinned
            else None
        ),
    }


def leave_one_out(pairs):
    rows = []
    if len(pairs) < 2:
        return rows
    for i, p in enumerate(pairs):
        rest = pairs[:i] + pairs[i + 1 :]
        rows.append(
            {
                "dropped": f"{p['a']}->{p['b']}",
                "dropped_de_overall": p["delta_de_overall"],
                "dropped_profile_pinned": p["profile_pinned"],
                "rms_without": rms([q["delta_de_overall"] for q in rest]),
                "quadrature_without": summarise(rest, "loo")["quadrature_over_6"],
            }
        )
    return rows


def trace_230():
    fit_path = os.path.join(ROOT, "experiments", "230-toll-magnitude", "work", "fit230.json")
    fit = json.load(open(fit_path))
    identical = [t for t in fit["table"] if t["digest_same"] == 900 and t["relL1"] == 0.0]
    ys = [t["d_reach"] for t in identical]
    sd = sample_sd(ys)
    # overall deltas transcribed from the board record, same six arms.
    # Sources are cited in RESULTS.md; stored here so the sd is recomputed.
    overall_delta = {
        "a8": 0.1543607 - 0.1503084,          # BOARD-LEDGER.md:170-171
        "a10": 0.1541841 - 0.1543607,         # BOARD-LEDGER.md:172-171
        "a21": -0.021458,                     # docs/LOG.md:2867
        "a22": -0.002518,                     # docs/LOG.md:2867
        "rc_hybrid": 0.0093325,               # BOARD-LEDGER.md:241
        "rc_callset": 0.0084448,              # BOARD-LEDGER.md:242
    }
    # reach cross-check from the scaled-member table
    reach_check = {
        "a8": 0.124349 - 0.102951,            # BOARD-LEDGER.md:171-170
        "a10": 0.123595 - 0.124349,           # BOARD-LEDGER.md:172-171
        "rc_hybrid": 0.170308 - 0.115368,     # BOARD-LEDGER.md:253, a15 reach :177
        "rc_callset": 0.192898 - 0.170308,    # BOARD-LEDGER.md:253
    }
    ovs = [overall_delta[t["arm"]] for t in identical]
    def without(names):
        sub = [t["d_reach"] for t in identical if t["arm"] not in names]
        return {
            "arms": [t["arm"] for t in identical if t["arm"] not in names],
            "sd_reach": sample_sd(sub) if len(sub) > 1 else None,
            "sd_reach_over_6": (sample_sd(sub) / 6.0) if len(sub) > 1 else None,
        }
    return {
        "source": "experiments/230-toll-magnitude/work/fit230.json falsifier2.zero_sd",
        "recorded_zero_sd": fit["falsifier2"]["zero_sd"],
        "recomputed_sd_reach_ddof1": sd,
        "recomputed_sd_reach_over_6": sd / 6.0,
        "n": len(ys),
        "arms": [
            {
                "arm": t["arm"],
                "incumbent": t["incumbent"],
                "d_reach_scaled": t["d_reach"],
                "what": t["what"],
                "lib_cv": t["lib_cv"],
                "digest_same": t["digest_same"],
                "reach_check_from_ledger": reach_check.get(t["arm"]),
                "reach_check_abs_err": (
                    abs(reach_check[t["arm"]] - t["d_reach"]) if t["arm"] in reach_check else None
                ),
                "overall_delta_transcribed": overall_delta[t["arm"]],
            }
            for t in identical
        ],
        "sd_of_transcribed_overall_deltas": sample_sd(ovs),
        "rms_of_transcribed_overall_deltas": rms(ovs),
        "designs_shared": False,
        "note": (
            "These six numbers are treatment effects of different emission laws "
            "against different incumbents. They are not re-rolls of one design."
        ),
        "sensitivity": {
            "drop_rc_hybrid": without({"rc_hybrid"}),
            "drop_rc_hybrid_and_a8": without({"rc_hybrid", "a8"}),
            "drop_two_largest_abs": without({"rc_hybrid", "a8"}),
            "only_a10_a22": without({"a8", "a21", "rc_hybrid", "rc_callset"}),
        },
    }


def h1_seed_groups():
    roots = [
        "experiments/435-zs-on-coarse",
        "experiments/517-flex-programme",
        "experiments/528-callprop-block",
        "experiments/480-nmae-instrument",
        "experiments/490-de-emitter-review",
        "experiments/467-strength-gated-programme",
        "experiments/429-zs-production",
    ]
    groups = defaultdict(list)
    n_files = 0
    for rel_root in roots:
        base = os.path.join(ROOT, rel_root)
        if not os.path.isdir(base):
            continue
        for dirpath, _, files in os.walk(base):
            for fn in files:
                if not fn.startswith("metrics") or not fn.endswith(".json"):
                    continue
                path = os.path.join(dirpath, fn)
                try:
                    blob = json.load(open(path))
                except Exception:
                    continue
                raw = blob.get("raw")
                if not isinstance(raw, dict) or "seed" not in blob:
                    continue
                if not all(m in raw for m in MEMBERS):
                    continue
                n_files += 1
                tag = str(blob.get("tag") or fn)
                tag_key = re.sub(r"_s\d+_", "_s*_", tag)
                # directory distinguishes designs that share a tag
                key = os.path.relpath(dirpath, ROOT) + "|" + tag_key
                groups[key].append(
                    {
                        "path": os.path.relpath(path, ROOT),
                        "seed": blob["seed"],
                        "tag": tag,
                        "raw": {m: float(raw[m]) for m in MEMBERS},
                    }
                )
    summaries = []
    for key, items in sorted(groups.items()):
        # one row per seed; if duplicates, keep the first path
        by_seed = {}
        for it in items:
            by_seed.setdefault(it["seed"], it)
        if len(by_seed) < 2:
            continue
        seeds = [by_seed[s] for s in sorted(by_seed)]
        raw_sd = {}
        for m in MEMBERS:
            raw_sd[m] = sample_sd([s["raw"][m] for s in seeds])
        def overall_of(raw, slopes):
            return sum(slopes[m] * raw[m] for m in MEMBERS) / 6.0
        ov429 = [overall_of(s["raw"], SLOPE) for s in seeds]
        # 467 spans: slope = 1/span
        slope467 = {m: 1.0 / SPAN_467[m] for m in MEMBERS}
        ov467 = [overall_of(s["raw"], slope467) for s in seeds]
        de429 = [
            sum(SLOPE[m] * s["raw"][m] for m in DE_MEMBERS) / 6.0 for s in seeds
        ]
        pds_const = len({s["raw"]["pds"] for s in seeds}) == 1
        mse_const = len({s["raw"]["mse"] for s in seeds}) == 1
        summaries.append(
            {
                "group": key,
                "n_seeds": len(seeds),
                "seeds": [s["seed"] for s in seeds],
                "paths": [s["path"] for s in seeds],
                "raw_sd": raw_sd,
                "raw_by_seed": {str(s["seed"]): s["raw"] for s in seeds},
                "overall_sd_slope429": sample_sd(ov429),
                "overall_values_slope429": ov429,
                "de_only_sd_slope429": sample_sd(de429),
                "overall_sd_span467": sample_sd(ov467),
                "pds_identical_across_seeds": pds_const,
                "mse_identical_across_seeds": mse_const,
                "scaled_to_900_sqrt": (
                    sample_sd(ov429) * math.sqrt(85 / 900) if len(seeds) > 1 else None
                ),
                "scaled_to_300_sqrt": (
                    sample_sd(ov429) * math.sqrt(85 / 300) if len(seeds) > 1 else None
                ),
            }
        )
    return {"n_metrics_files_with_raw_and_seed": n_files, "groups_with_ge2_seeds": summaries}


def per_target_concentration():
    """Does one H1 re-roll's panel-mean shift come from a few targets?"""
    directory = os.path.join(ROOT, "experiments", "435-zs-on-coarse", "work", "h1")
    stems = {
        0: "per-target-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.csv",
        1: "per-target-base_s1_bps_prodem_mde0.1cap0.1min0.3coarse.csv",
        2: "per-target-base_s2_bps_prodem_mde0.1cap0.1min0.3coarse.csv",
    }
    tables = {}
    for seed, fn in stems.items():
        path = os.path.join(directory, fn)
        if not os.path.exists(path):
            return {"available": False, "missing": fn}
        acc = defaultdict(dict)
        with open(path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                acc[row["metric"]][row["perturbation"]] = float(row["value"])
        tables[seed] = acc
    # metric name -> board member
    name_map = {
        "pds_cosine": "pds",
        "expr_mse_unbiased_capped_norm": "mse",
        "de_wilcoxon_direction_fidelity_yield_raw": "fid",
        "de_wilcoxon_direction_reach_raw": "reach",
        "de_wilcoxon_sig_jaccard": "jac",
        "de_wilcoxon_lfc_nmae": "nmae",
    }
    out = {"available": True, "metrics": {}}
    for metric, member in name_map.items():
        if metric not in tables[0]:
            continue
        def finite(seed, t):
            v = tables[seed][metric].get(t)
            return v is not None and math.isfinite(v)

        targets = sorted(
            t
            for t in set(tables[0][metric]) & set(tables[1][metric]) & set(tables[2][metric])
            if finite(0, t) and finite(1, t) and finite(2, t)
        )
        # pairwise mean shifts and concentration of s0-s1
        pair_stats = []
        for i, j in ((0, 1), (0, 2), (1, 2)):
            diffs = [tables[j][metric][t] - tables[i][metric][t] for t in targets]
            abs_sorted = sorted((abs(d), d) for d in diffs)
            abs_sum = sum(a for a, _ in abs_sorted) or 1.0
            top5 = abs_sorted[-5:]
            top10 = abs_sorted[-max(1, len(diffs) // 10) :]
            pair_stats.append(
                {
                    "seeds": [i, j],
                    "n_targets": len(diffs),
                    "mean_diff": mean(diffs),
                    "rms_per_target": rms(diffs),
                    "sd_mean_if_iid": rms(diffs) / math.sqrt(len(diffs)),
                    "top5_share_of_abs": sum(a for a, _ in top5) / abs_sum,
                    "top_decile_share_of_abs": sum(a for a, _ in top10) / abs_sum,
                }
            )
        # observed sd of the three panel means
        means = []
        for s in (0, 1, 2):
            vals = []
            for t in targets:
                v = tables[s][metric][t]
                if v == v:  # finite
                    vals.append(v)
            means.append(mean(vals) if vals else float("nan"))
        out["metrics"][member] = {
            "metric": metric,
            "n_targets": len(targets),
            "panel_means": means,
            "panel_mean_sd": sample_sd(means),
            "pairs": pair_stats,
        }
    return out


def seed_key_census(metas):
    """Every meta, not only scored ones: do any two differ only in a seed?"""
    census = {k: defaultdict(int) for k in SEED_KEYS}
    present = {k: 0 for k in SEED_KEYS}
    for _rel, cfg, _impl, _fp in metas:
        for k in SEED_KEYS:
            if k in cfg:
                present[k] += 1
                census[k][canon(cfg.get(k))] += 1
    # full-config identity ignoring seed keys
    buckets = defaultdict(list)
    for rel, cfg, _impl, _fp in metas:
        flat = flatten(cfg)
        for k in SEED_KEYS:
            flat.pop(k, None)
            # nested copies
            flat = {p: v for p, v in flat.items() if not p.startswith(k + ".") and p != k}
        buckets[canon(flat)].append(rel)
    seed_only_groups = []
    for key, rels in buckets.items():
        if len(rels) < 2:
            continue
        # confirm the dropped keys actually differ inside the group
        seed_vals = defaultdict(set)
        for rel, cfg, _impl, _fp in metas:
            if rel not in rels:
                continue
            for k in SEED_KEYS:
                if k in cfg:
                    seed_vals[k].add(canon(cfg.get(k)))
        varying = {k: sorted(v) for k, v in seed_vals.items() if len(v) > 1}
        if varying:
            seed_only_groups.append({"n": len(rels), "metas": rels, "seed_keys_that_vary": varying})
    return {
        "n_metas": len(metas),
        "seed_key_present": present,
        "seed_key_values": {k: dict(v) for k, v in census.items()},
        "groups_identical_except_seed_keys": seed_only_groups,
    }


def main():
    os.makedirs(OUT, exist_ok=True)
    trace = trace_230()
    json.dump(trace, open(os.path.join(OUT, "trace230.json"), "w"), indent=2)

    board = load_board()
    joined = []
    unmatched = []
    for row in board:
        if not row["meta"]:
            unmatched.append({"entry_id": row["entry_id"], "package": row["package"], "note": row["note"]})
            continue
        cfg, impl = load_meta(row["meta"])
        row["cfg"] = cfg
        row["impl"] = impl
        row["fp"] = fingerprint(cfg)
        row["fp_key"] = canon(row["fp"])
        row["flat"] = flatten(cfg)
        row["label"] = arm_name(row["package"])
        joined.append(row)

    # all metas, for the seed census (includes unscored packages)
    all_metas = []
    for path in sorted(glob.glob(os.path.join(ROOT, "submissions", "*.meta.json"))):
        rel = os.path.relpath(path, ROOT)
        cfg, impl = load_meta(rel)
        all_metas.append((rel, cfg, impl, fingerprint(cfg)))
    census = seed_key_census(all_metas)
    json.dump(census, open(os.path.join(OUT, "seed_census.json"), "w"), indent=2)

    # cliques of scored arms with an identical strict fingerprint
    cliques = defaultdict(list)
    for row in joined:
        cliques[row["fp_key"]].append(row)
    for members in cliques.values():
        members.sort(key=lambda r: r["date"] or "")

    all_pairs = []
    consecutive = []
    clique_report = []
    for key, members in cliques.items():
        if len(members) < 2:
            continue
        labels = [m["label"] for m in members]
        de_vals = [sum(m["scores"][k] for k in DE_MEMBERS) / 6.0 for m in members]
        clique_report.append(
            {
                "labels": labels,
                "entries": [m["entry_id"] for m in members],
                "n": len(members),
                "de_overall_values": de_vals,
                "de_overall_sd": sample_sd(de_vals),
                "dates": [m["date"] for m in members],
            }
        )
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                all_pairs.append(pair_record(members[i], members[j]))
            if i + 1 < len(members):
                consecutive.append(pair_record(members[i], members[i + 1]))

    # adversarial: a differing path that looks DE-ish means the fingerprint missed a key
    missed = []
    for p in all_pairs:
        bad = [k for k in p["config_paths_differ"] if DEISH.search(k)]
        if bad:
            missed.append({"pair": f"{p['a']}->{p['b']}", "paths": bad})

    # single top-level-key diffs among strict pairs (closest to one-knob bulk arms)
    single_key = [p for p in all_pairs if len(p["config_top_differ"]) == 1]
    few_key = [p for p in all_pairs if len(p["config_top_differ"]) <= 3]

    # sensitivity: rebuild pairs after dropping one fingerprint component.
    # Only report pairs that the strict fingerprint separated and the drop joins,
    # plus the rms of the enlarged consecutive set.
    def pairs_for(drop):
        buckets = defaultdict(list)
        for row in joined:
            fp = fingerprint(row["cfg"], drop=drop)
            buckets[canon(fp)].append(row)
        cons, extra = [], []
        strict_keys = {}
        for row in joined:
            strict_keys[row["entry_id"]] = row["fp_key"]
        for members in buckets.values():
            members.sort(key=lambda r: r["date"] or "")
            if len(members) < 2:
                continue
            for i in range(len(members) - 1):
                rec = pair_record(members[i], members[i + 1])
                same_strict = strict_keys[members[i]["entry_id"]] == strict_keys[members[i + 1]["entry_id"]]
                rec["same_strict_fingerprint"] = same_strict
                cons.append(rec)
                if not same_strict:
                    extra.append(rec)
        return cons, extra

    sensitivity = {}
    for name in SENSITIVITY_DROPS:
        cons, extra = pairs_for((name,))
        sensitivity[name] = {
            "consecutive_summary": summarise(cons, f"drop:{name}"),
            "n_extra_pairs_the_drop_joins": len(extra),
            "extra_pairs": [
                {
                    "a": p["a"],
                    "b": p["b"],
                    "delta_de_overall": p["delta_de_overall"],
                    "delta_scaled": p["delta_scaled"],
                    "profile_pinned": p["profile_pinned"],
                    "config_top_differ": p["config_top_differ"],
                }
                for p in extra
            ],
        }

    # prior-only: strict fingerprints differ, and the ONLY fingerprint
    # components that differ are sameline / de_prior_gate / ctx_a_prior /
    # de_sign_prior*. Built by dropping those and requiring the rest to match,
    # then keeping pairs whose strict fingerprints differ.
    prior_cons, prior_extra = pairs_for(("sameline", "ctx_a_prior"))
    prior_only = [p for p in prior_extra]
    # filter to pairs whose config diff, outside prior keys, is empty or bulk-only.
    # We already know the non-prior fingerprint matches. Report them as effects.
    prior_effects = []
    for p in prior_only:
        prior_effects.append(
            {
                "a": p["a"],
                "b": p["b"],
                "delta_de_overall": p["delta_de_overall"],
                "delta_overall": p["delta_overall"],
                "delta_scaled": p["delta_scaled"],
                "profile_pinned": p["profile_pinned"],
                "config_top_differ": p["config_top_differ"],
            }
        )

    summaries = {
        "strict_all_unordered_pairs": summarise(all_pairs, "strict_all"),
        "strict_time_consecutive": summarise(consecutive, "strict_consecutive"),
        "strict_single_top_key": summarise(single_key, "strict_single_key"),
        "strict_at_most_3_top_keys": summarise(few_key, "strict_few_keys"),
        "strict_consecutive_profile_pinned": summarise(
            [p for p in consecutive if p["profile_pinned"]], "strict_consecutive_pinned"
        ),
    }

    h1 = h1_seed_groups()
    conc = per_target_concentration()

    payload = {
        "checkout_note": (
            "experiments/547-de-noise-check and experiments/537-arm-pipeline are absent "
            "from this checkout. Pair set is reconstructed, not replayed."
        ),
        "n_board_entries": len(board),
        "n_joined_to_meta": len(joined),
        "unmatched_board_entries": unmatched,
        "n_strict_cliques_of_size_ge2": len(clique_report),
        "cliques": clique_report,
        "summaries": summaries,
        "leave_one_out_consecutive": leave_one_out(consecutive),
        "missed_deish_paths_inside_strict_pairs": missed,
        "sensitivity_drop_one_component": {
            k: {
                "consecutive_summary": v["consecutive_summary"],
                "n_extra_pairs_the_drop_joins": v["n_extra_pairs_the_drop_joins"],
                "extra_pairs": v["extra_pairs"],
            }
            for k, v in sensitivity.items()
        },
        "prior_only_effects_not_noise": prior_effects,
        "consecutive_pairs": consecutive,
        "all_unordered_pairs": all_pairs,
    }
    json.dump(payload, open(os.path.join(OUT, "pairs548.json"), "w"), indent=2)

    # flat csv of consecutive strict pairs, the primary set
    csv_path = os.path.join(OUT, "pairs548.csv")
    fields = [
        "a",
        "b",
        "a_date",
        "b_date",
        "profile_pinned",
        "delta_overall",
        "delta_de_overall",
        "d_pds",
        "d_mse",
        "d_fid",
        "d_reach",
        "d_jac",
        "d_nmae",
        "n_config_paths_differ",
        "config_top_differ",
        "emitter_code_differs",
        "builder_code_differs",
    ]
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for p in consecutive:
            w.writerow(
                {
                    "a": p["a"],
                    "b": p["b"],
                    "a_date": p["a_date"],
                    "b_date": p["b_date"],
                    "profile_pinned": p["profile_pinned"],
                    "delta_overall": f"{p['delta_overall']:.8f}",
                    "delta_de_overall": f"{p['delta_de_overall']:.8f}",
                    "d_pds": f"{p['delta_scaled']['pds']:.8f}",
                    "d_mse": f"{p['delta_scaled']['mse']:.8f}",
                    "d_fid": f"{p['delta_scaled']['fid']:.8f}",
                    "d_reach": f"{p['delta_scaled']['reach']:.8f}",
                    "d_jac": f"{p['delta_scaled']['jac']:.8f}",
                    "d_nmae": f"{p['delta_scaled']['nmae']:.8f}",
                    "n_config_paths_differ": p["n_config_paths_differ"],
                    "config_top_differ": "|".join(p["config_top_differ"]),
                    "emitter_code_differs": p["emitter_code_differs"],
                    "builder_code_differs": p["builder_code_differs"],
                }
            )

    json.dump({"groups": h1, "per_target_435": conc}, open(os.path.join(OUT, "h1_seeds.json"), "w"), indent=2)

    def in_family(p, labels):
        return p["a"] in labels and p["b"] in labels

    modern = {"433h", "433k", "446b", "446bj2", "446bd2"}
    sharpen = {"410b", "417a", "417e"}
    early_pairs = [p for p in consecutive if p["a"] not in modern | sharpen and p["b"] not in modern | sharpen]
    modern_pairs = [p for p in consecutive if in_family(p, modern)]
    sharpen_pairs = [p for p in consecutive if in_family(p, sharpen)]
    # Programme arms opened by dropping sameline. extra_pairs are summary dicts.
    cont_de = [p["delta_de_overall"] for p in modern_pairs] + [
        p["delta_de_overall"]
        for p in sensitivity["sameline"]["extra_pairs"]
        if (p["a"], p["b"]) in {("446bd2", "458d"), ("458d", "458g")}
    ]

    def paired_effect(arm_glob, base_glob):
        arm = {}
        base = {}
        for path in glob.glob(arm_glob):
            blob = json.load(open(path))
            if "raw" in blob and "seed" in blob:
                arm[blob["seed"]] = blob["raw"]
        for path in glob.glob(base_glob):
            blob = json.load(open(path))
            if "raw" in blob and "seed" in blob and "coarse" in path and "dnk" not in path:
                base[blob["seed"]] = blob["raw"]
        seeds = sorted(set(arm) & set(base))
        diffs = []
        for s in seeds:
            d_raw = {m: arm[s][m] - base[s][m] for m in MEMBERS}
            d_ov = sum(SLOPE[m] * d_raw[m] for m in MEMBERS) / 6.0
            diffs.append({"seed": s, "d_raw": d_raw, "d_overall": d_ov})
        return {
            "seeds": seeds,
            "per_seed": diffs,
            "sd_overall": sample_sd([d["d_overall"] for d in diffs]),
            "mean_overall": mean([d["d_overall"] for d in diffs]),
        }

    paired = {
        "gated33_minus_435": paired_effect(
            os.path.join(ROOT, "experiments/467-strength-gated-programme/work/h1_gated33/metrics-base_s*.json"),
            os.path.join(ROOT, "experiments/435-zs-on-coarse/work/h1/metrics-base_s*_bps_prodem_mde0.1cap0.1min0.3coarse.json"),
        ),
        "gated67_minus_435": paired_effect(
            os.path.join(ROOT, "experiments/467-strength-gated-programme/work/h1_gated67/metrics-base_s*.json"),
            os.path.join(ROOT, "experiments/435-zs-on-coarse/work/h1/metrics-base_s*_bps_prodem_mde0.1cap0.1min0.3coarse.json"),
        ),
    }
    # 446b vs 446bja from the joined board rows (sign-prior arm, not in the strict set)
    by_label = {r["label"]: r for r in joined}
    sign_prior_pairs = []
    if "446b" in by_label and "446bja" in by_label:
        sign_prior_pairs.append(pair_record(by_label["446b"], by_label["446bja"]))

    h1_sds = [
        g["de_only_sd_slope429"]
        for g in h1["groups_with_ge2_seeds"]
        if g["n_seeds"] >= 3 and g["de_only_sd_slope429"] is not None
    ]
    h1_sds_sorted = sorted(h1_sds)
    families = {
        "modern_433h_to_446bd2": summarise(modern_pairs, "modern"),
        "sharpen_410b_to_417e": summarise(sharpen_pairs, "sharpen"),
        "early_222_to_319p2": summarise(early_pairs, "early"),
        "modern_plus_programme_arms": {
            "deltas": cont_de,
            "rms": rms(cont_de),
            "note": "modern steps plus 446bd2/458d/458g links opened by dropping sameline",
        },
        "sign_prior_446b_to_446bja": [
            {
                "delta_de_overall": p["delta_de_overall"],
                "delta_scaled": p["delta_scaled"],
                "profile_pinned": p["profile_pinned"],
                "config_top_differ": p["config_top_differ"],
            }
            for p in sign_prior_pairs
        ],
        "paired_h1": paired,
        "modern_leave_one_step_out": leave_one_out(modern_pairs),
        "h1_de_only_sd_n3": {
            "n_groups": len(h1_sds_sorted),
            "min": h1_sds_sorted[0] if h1_sds_sorted else None,
            "median": h1_sds_sorted[len(h1_sds_sorted) // 2] if h1_sds_sorted else None,
            "max": h1_sds_sorted[-1] if h1_sds_sorted else None,
        },
    }
    json.dump(families, open(os.path.join(OUT, "families548.json"), "w"), indent=2)

    # stdout: the numbers RESULTS.md will quote
    print("joined", len(joined), "of", len(board), "board entries;", len(unmatched), "unmatched")
    print("strict cliques", len(clique_report), "consecutive pairs", len(consecutive), "all pairs", len(all_pairs))
    print("230 sd", trace["recomputed_sd_reach_ddof1"], "over 6", trace["recomputed_sd_reach_over_6"])
    print("recorded", trace["recorded_zero_sd"])
    for name, s in summaries.items():
        print(
            f"{name}: n={s.get('n_pairs')} rms={s.get('de_overall_rms')} "
            f"quad={s.get('quadrature_over_6')} pinned_rms={s.get('pinned_de_overall_rms')}"
        )
    print("missed deish paths", len(missed))
    print("seed-only groups", len(census["groups_identical_except_seed_keys"]))
    print("h1 groups", len(h1["groups_with_ge2_seeds"]), "from", h1["n_metrics_files_with_raw_and_seed"], "files")
    print("prior-only extra consecutive", len(prior_effects))


if __name__ == "__main__":
    main()
