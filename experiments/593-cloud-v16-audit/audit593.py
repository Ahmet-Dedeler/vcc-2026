"""593 — recompute every changed number in manuscript Draft v16 from primary artifacts.

Reads tracked files. Draft v16, 537's noise json, and 590 live on origin/master
(this branch was cut before those commits); board entries, build metas, fit230 and
the H1 diag jsons are in the working tree. Nothing outside this experiment is written.

Usage: uv run python experiments/593-cloud-v16-audit/audit593.py
"""

from __future__ import annotations

import csv
import json
import math
import subprocess
from pathlib import Path

import numpy as np
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "work" / "audit593.out"
REF = "origin/master"

E = {
    "446b": "UDgjOKL9gI1udAexq6ue",
    "446bj2": "YvVru3DAF89Jiyr34olR",
    "446bj2z": "be90smZUDcIV7j5MgQ03",
    "446bd2": "R5uE5Rve3R9nWXGbeSLM",
    "453dr": "q3VuTImh4FRGKXH2waCp",
    "457dr": "1FfirZsGsAZvaXYTd39J",
    "460dr": "LrkhsIXWMjEsvRXSGToa",
    "446bj2zcp": "K1Q0y2OFYjTNRqh9jKfk",
    "458d": "b7UI3NRy1rdVNRRyZwmU",
    "458g": "JL8llNLHaTyDQzm798TV",
}
MEM = ["pds", "mse", "nmae", "fid", "reach", "jac"]
RAW = {
    "pds": "pds_cosine",
    "mse": "expr_mse_unbiased_capped_norm",
    "nmae": "de_wilcoxon_lfc_nmae",
    "fid": "de_wilcoxon_direction_fidelity_yield_raw",
    "reach": "de_wilcoxon_direction_reach_raw",
    "jac": "de_wilcoxon_sig_jaccard",
}
DE = ["nmae", "fid", "reach", "jac"]
METAS = {
    "446b": "446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5",
    "458": "458-cap15-partial50-arch1prog-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5",
    "458d": "458d-cap15-partial50-arch1dprog-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5",
    "458g": "458g-cap15-partial50-arch1dprog-gate33-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5",
}
B1 = [
    ("417a", "410b", "T88GsCjp53Bmk8xQcJMS", "4WgVi6FIE4iVw9PmRUP3"),
    ("417e", "417a", "wPP8U70kj4mOIvtA83CO", "T88GsCjp53Bmk8xQcJMS"),
    ("433k", "433h", "T76JWD3tOajk5IYvzwUW", "wMoDAi1mca08FA2Bfqkd"),
    ("446b", "433k", "UDgjOKL9gI1udAexq6ue", "T76JWD3tOajk5IYvzwUW"),
    ("446bj2", "446b", "YvVru3DAF89Jiyr34olR", "UDgjOKL9gI1udAexq6ue"),
    ("446bd2", "446b", "R5uE5Rve3R9nWXGbeSLM", "UDgjOKL9gI1udAexq6ue"),
    ("457dr", "453dr", "1FfirZsGsAZvaXYTd39J", "q3VuTImh4FRGKXH2waCp"),
]
# Reach's (r-b) is not 1. Overall contribution of a raw delta is raw * INV / 6,
# with INV the reciprocal of the member span used by 537 (noise537.py). Recomputed
# here from the captured score and the raw member, not copied from that script.
lines: list[str] = []


def emit(s: str = "") -> None:
    lines.append(s)
    print(s)


def text_of(rel: str) -> str:
    p = ROOT / rel
    if p.is_file():
        return p.read_text()
    return subprocess.check_output(["git", "show", f"{REF}:{rel}"], cwd=ROOT, text=True)


def json_of(rel: str):
    return json.loads(text_of(rel))


def entry(key: str) -> dict:
    j = json_of(f"data/board-entries/{E[key]}.json")
    return j.get("entry", j)


def entry_id(eid: str) -> dict:
    j = json_of(f"data/board-entries/{eid}.json")
    return j.get("entry", j)


def delta(a: dict, b: dict) -> dict:
    raw = {m: a[RAW[m]] - b[RAW[m]] for m in MEM}
    scaled = {m: a["score_" + m] - b["score_" + m] for m in MEM}
    share = {m: scaled[m] / 6 for m in MEM}
    return {
        "overall": a["score_avg"] - b["score_avg"],
        "bulk": share["pds"] + share["mse"],
        "de": sum(share[m] for m in DE),
        "raw": raw,
        "scaled": scaled,
        "share": share,
    }


def fmt(x: float, n: int = 7) -> str:
    return f"{x:+.{n}f}"


def main() -> None:
    emit("593 audit — Draft v16 numbers recomputed from primary artifacts")
    emit(f"branch files from the work tree; missing paths via git show {REF}")
    emit()

    cache = {k: entry(k) for k in E}

    emit("=== 1. Board deltas (score_* / 6 = contribution to overall; raw = cell-eval2 member) ===")
    pairs = [
        ("460dr", "453dr"),
        ("446bj2", "446b"),
        ("460dr", "446bj2"),
        ("446bd2", "446b"),
        ("453dr", "446b"),
        ("457dr", "446b"),
        ("457dr", "453dr"),
        ("446bj2zcp", "446bj2z"),
        ("458d", "446b"),
        ("458g", "458d"),
        ("458g", "446b"),
        ("460dr", "446b"),
    ]
    D = {}
    for a, b in pairs:
        d = delta(cache[a], cache[b])
        D[(a, b)] = d
        emit(f"{a} - {b}: overall {fmt(d['overall'])}  bulk {fmt(d['bulk'])}  DE {fmt(d['de'])}")
        emit(
            "    raw  "
            + "  ".join(f"{m} {d['raw'][m]:+.6e}" for m in MEM)
        )
        emit(
            "    share "
            + "  ".join(f"{m} {d['share'][m]:+.7f}" for m in MEM)
        )
    d = D[("460dr", "453dr")]
    k = D[("446bj2", "446b")]
    emit(f"mse ratio 460dr/key = {d['raw']['mse'] / k['raw']['mse']:.6f}")
    emit(f"pds ratio 460dr/key = {d['raw']['pds'] / k['raw']['pds']:.6f}")
    fj = D[("446bj2zcp", "446bj2z")]
    emit(
        f"zcp d(fid+jac) overall = {fj['share']['fid'] + fj['share']['jac']:+.7f}"
    )
    h1_fid, h1_jac = -0.0013, 0.0091  # 528 RESULTS table, f 1, lcap-kept
    emit(
        f"fid ratio board/H1 = {fj['raw']['fid'] / h1_fid:.4f}; "
        f"jac ratio = {fj['raw']['jac'] / h1_jac:.4f}"
    )
    emit("3σ bands on 460dr-453dr raw: "
         f"|fid| {abs(d['raw']['fid']):.6e} ≤ 0.0001; "
         f"|jac| {abs(d['raw']['jac']):.6e} ≤ 0.00005; "
         f"|reach| {abs(d['raw']['reach']):.6e} ≤ 0.0040")
    emit()
    emit("DE shares vs 446b (the C-82 cluster):")
    for name in ("446bj2", "446bd2", "453dr", "457dr", "460dr"):
        emit(f"  {name:8} {D[(name, '446b')]['de']:+.7f}")

    emit()
    emit("=== 2. Delivered dose (build metas, output.rc_callset) ===")
    lams = {}
    for name, stem in METAS.items():
        m = json_of(f"submissions/{stem}.h5ad.meta.json")
        rc = m["output"]["rc_callset"]
        lams[name] = rc["prior_lambda"]
        emit(f"{name}: blocks {rc.get('blocks')}  prior_lambda {lams[name]}")
    med = lams["458d"]["median"]
    emit(f"458d median exact {med!r}  == 887/2048? {math.isclose(med, 887/2048)}")
    emit(f"10-step bisection grid quantum 1/1024 = {1/1024}")
    emit(f"largest k/1024 ≤ median: {math.floor(med * 1024)}/1024 = {math.floor(med * 1024)/1024}")
    strong, weak = 214 + 217 + 217, 900 - (214 + 217 + 217)
    emit(f"468 gate: strong {strong} weak {weak}; ≥ {450 - weak} strong blocks ≤ median")
    emit("consensus-prior metas whose pooled median is not 1:")
    n_meta = 0
    for p in sorted((ROOT / "submissions").glob("*.h5ad.meta.json")):
        m = json.loads(p.read_text())
        rc = (m.get("output") or {}).get("rc_callset") or {}
        pl = rc.get("prior_lambda")
        if not isinstance(pl, dict) or "median" not in pl:
            continue
        n_meta += 1
        cfg = m.get("config") or {}
        prior = str(cfg.get("de_sign_prior", ""))
        if abs(pl["median"] - 1.0) > 1e-9:
            emit(f"  {p.name.split('-')[0]:8} median {pl['median']:.6f} min {pl['min']:.4f} "
                 f"prior {prior!r} blocks {rc.get('blocks')}")
    emit(f"metas with prior_lambda: {n_meta}")

    emit()
    emit("=== 3. H1 programme prior_lambda (diag jsons; not the RESULTS sentence) ===")
    diags = [
        "experiments/462-broad16-lab-signature/work/h1_arcbank/diag-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.json",
        "experiments/462-broad16-lab-signature/work/h1_arcfull/diag-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.json",
        "experiments/467-strength-gated-programme/work/h1_gated33/diag-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.json",
        "experiments/517-flex-programme/work/h1_g33-none/diag-base_s0_bps_hookem_mde0.1cap0.1min0.3coarse.json",
        "experiments/517-flex-programme/work/h1_arc-none-q0.7/diag-base_s0_bps_hookem_mde0.1cap0.1min0.3coarse.json",
    ]
    for rel in diags:
        rows = json_of(rel)
        low = [(r.get("target") or r.get("gene"), r["prior_lambda"]) for r in rows if r["prior_lambda"] < 1 - 1e-12]
        emit(f"{rel.split('work/')[1]}: {len(rows) - len(low)}/{len(rows)} at 1; below: {low}")

    emit()
    emit("=== 4. DE noise: seven B1 pairs from board JSON, χ² interval, 230's sd ===")
    reach, fid, jac, nmae, de_ov = [], [], [], [], []
    for _a, _b, ida, idb in B1:
        dlt = delta(entry_id(ida), entry_id(idb))
        reach.append(dlt["raw"]["reach"])
        fid.append(dlt["raw"]["fid"])
        jac.append(dlt["raw"]["jac"])
        nmae.append(dlt["raw"]["nmae"])
        de_ov.append(dlt["de"])
        emit(f"  {_a}-{_b}: dreach {dlt['raw']['reach']:+.6f} dfid {dlt['raw']['fid']:+.6e} "
             f"djac {dlt['raw']['jac']:+.6e} dnmae {dlt['raw']['nmae']:+.6e} DE {dlt['de']:+.5f}")

    def rms(xs: list[float]) -> float:
        a = np.asarray(xs, float)
        return float(np.sqrt(np.mean(a ** 2)))

    def chi_ci(r: float, n: int) -> tuple[float, float]:
        # mean-zero normal: n * rms^2 / σ^2 ~ χ²_n
        return r * math.sqrt(n / chi2.ppf(0.975, n)), r * math.sqrt(n / chi2.ppf(0.025, n))

    rr, fr, jr = rms(reach), rms(fid), rms(jac)
    lo, hi = chi_ci(rr, 7)
    emit(f"reach rms {rr:.5f}  χ²95 [{lo:.5f}, {hi:.5f}]")
    emit(f"fid rms {fr:.6e}  jac rms {jr:.6e}  nmae rms {rms(nmae):.6e}  DE-overall rms {rms(de_ov):.5f}")
    # span implied by score_reach vs raw reach on one entry, to convert rms to overall
    e0 = cache["446b"]
    # not used: share already computed per pair. Reach's overall rms:
    reach_share = []
    for _a, _b, ida, idb in B1:
        dlt = delta(entry_id(ida), entry_id(idb))
        reach_share.append(dlt["share"]["reach"])
    emit(f"reach share rms {rms(reach_share):.5f}")

    fit = json_of("experiments/230-toll-magnitude/work/fit230.json")
    zero = fit["falsifier2"]["zero_sd"]
    emit(f"fit230 falsifier2.zero_sd {zero:.15f}  /6 = {zero/6:.15f}  n {fit['falsifier2']['zero_n']}")
    noise = json_of("experiments/537-arm-pipeline/work/noise537b.json")
    dreach = [p["d_reach"] for p in noise["ka230"]["pairs"]]
    emit(f"ka230 six d_reach {dreach}")
    emit(f"  sample sd ddof1 {float(np.std(dreach, ddof=1)):.15f}")
    emit(f"  population sd   {float(np.std(dreach, ddof=0)):.15f}")
    emit(f"  stored sd_scaled {noise['ka230']['sd_scaled']}  sd_overall {noise['ka230']['sd_overall']}")
    h1 = noise["h1"]["base_reroll_independent"]["1.0"]["reach"]
    emit(f"H1 f1 reach per_target_sd {h1['per_target_sd']}")
    emit(f"  h1_panel_sd {h1['h1_panel_sd']}  n_targets {h1['n_targets']}")
    pred = math.sqrt(2) * h1["per_target_sd"] / math.sqrt(900)
    naive = h1["h1_panel_sd"] * math.sqrt(85 / 900)
    emit(f"  √2 * per_target_sd / √900 = {pred:.6f}  (stored board_diff_sd_raw {h1['board_diff_sd_raw']:.6f})")
    emit(f"  h1_panel_sd * √(85/900) = {naive:.6f}   (this is NOT 0.00123)")

    redraw = json_of("experiments/537-arm-pipeline/work/redraw537.json")
    bulk = redraw["bulk-only step (ps bulk scale ^1.2)"]
    emit("bulk step pass-1 redrawn (1 - pass1_equal_after_frac) * 85:")
    for s, row in bulk["per_seed"].items():
        n_eq = row["pass1_equal_after_frac"] * 85
        emit(f"  seed {s}: K moved {row['n_K_changed']}/85  pass1 still equal {n_eq:.2f}/85  "
             f"same-K targets {row['sameK_after']} of which pass1 equal {row['sameK_after_pass1_equal']}")
    rho = bulk.get("rho_same_seed", {})
    emit(f"  rho_same_seed keys {list(rho)[:8]}")
    emit(f"  rho {json.dumps(rho)[:800]}")

    emit()
    emit("=== 4b. What 'scaled by √(85/900)' does and does not equal ===")
    sig = h1["per_target_sd"]
    pred_code = math.sqrt(2) * sig / math.sqrt(900)
    panel = h1["h1_panel_sd"]
    emit(f"board_diff_sd_raw stored {h1['board_diff_sd_raw']:.6f}")
    emit(f"√2 * per_target_sd / √900 = {pred_code:.6f}")
    emit(f"h1_panel_sd * √(85/900) = {panel * math.sqrt(85/900):.6f}  (not 0.00123)")
    emit(f"n_targets finite {h1['n_targets']} (the prose says 85)")

    emit()
    emit("=== 5. 590 identity, analyse590's cohort, recomputed from the CSVs and fold masks ===")
    cache_dir = Path(__file__).resolve().parent / "work" / "src590"
    cache_dir.mkdir(parents=True, exist_ok=True)
    blobs = {
        "plain.csv": "experiments/590-identity-validation/work/match590.csv",
        "flex0.csv": "experiments/590-identity-validation/work/match590_flexbias_fold0.csv",
        "flex1.csv": "experiments/590-identity-validation/work/match590_flexbias_fold1.csv",
        "fold0.npz": "experiments/590-identity-validation/work/flexbias590_fold0.npz",
        "fold1.npz": "experiments/590-identity-validation/work/flexbias590_fold1.npz",
    }
    for name, rel in blobs.items():
        dest = cache_dir / name
        if not dest.exists():
            dest.write_bytes(subprocess.check_output(["git", "show", f"{REF}:{rel}"], cwd=ROOT))
    REG = {
        "HCT116_LARGE_INTESTINE", "K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "DLD1_LARGE_INTESTINE",
        "JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "HEPG2_LIVER", "A549_LUNG", "HT29_LARGE_INTESTINE",
        "MCF7_BREAST", "BXPC3_PANCREAS", "HAP1_ENGINEERED", "CALU3_LUNG", "THP1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
    }
    UNREADABLE = {"sciplex3", "cui", "liang", "gehring"}

    def load_rows(name: str) -> list[dict]:
        with (cache_dir / name).open() as f:
            rows = list(csv.DictReader(f))
        return [r for r in rows if r["dataset"] not in UNREADABLE and r["query"] != "donor:kolf2"]

    def as_bool(v: str) -> bool:
        return v == "True"

    def flex_combined() -> list[dict]:
        out = []
        for k in (0, 1):
            d = load_rows(f"flex{k}.csv")
            z = np.load(cache_dir / f"fold{k}.npz", allow_pickle=True)
            applies = {str(q).split(":", 1)[1] for q in z["apply_to"].astype(str)}
            for r in d:
                mix = r["dataset"].startswith("mcfarland")
                lid = r["query"].split(":", 1)[1] if ":" in r["query"] else r["query"]
                keep = (mix and lid in applies) or ((not mix) if k == 0 else False)
                if keep:
                    out.append(r)
        return out

    flex_rows = flex_combined()
    flex_queries = {r["query"] for r in flex_rows}
    plain_rows = [r for r in load_rows("plain.csv")
                  if r["query"] in flex_queries or r["dataset"] in ("challenge", "anchor")]

    def report(label: str, rows: list[dict]) -> None:
        inn = [r for r in rows if as_bool(r["true_in_ccle"])]
        wrong = [r for r in inn if not as_bool(r["correct"])]
        una = [r for r in inn if as_bool(r["unanimous"])]
        una_ok = [r for r in una if as_bool(r["correct"])]
        corr = [r for r in inn if as_bool(r["correct"])]
        cu = [r for r in corr if as_bool(r["unanimous"])]
        reg_lin = {r["true_lineage"] for r in rows if r["true"] in REG and r["true_lineage"]}
        sister = [r for r in inn if r["true"] not in REG and r["true_lineage"] in reg_lin]
        sw = [r for r in rows if as_bool(r["switch_on"])]
        sw_fp = [r for r in sw if r["call"] != r["true"]]
        regrows = [r for r in rows if r["true"] in REG]
        reg_ok = [r for r in regrows if as_bool(r["switch_on"]) and as_bool(r["correct"])]
        reg_ok_gate = [r for r in reg_ok if float(r["min_s1"]) >= 0.65]
        emit(f"-- {label}: n {len(rows)}  in CCLE {len(inn)}  top1 {len(corr)/len(inn):.6f}  "
             f"acc|una {len(una_ok)/len(una):.6f}")
        emit(f"  wrong {len(wrong)}  sister {len(sister)} (lineages {sorted(reg_lin)})")
        emit(f"  switch_on {len(sw)}  false {len(sw_fp)}  registry correct switches {len(reg_ok)}  "
             f"still on at 0.65 {len(reg_ok_gate)}")
        for r in wrong:
            emit(f"    WRONG {r['query']}: {r['true']} -> {r['call']} una {r['unanimous']} "
                 f"min_s1 {float(r['min_s1']):.4f}")
        for r in sw_fp:
            emit(f"    FALSE SWITCH {r['query']} -> {r['call']} min_s1 {float(r['min_s1']):.4f}")
        q_u = float(np.quantile([float(r["min_s1"]) for r in cu], 0.01))
        q_c = float(np.quantile([float(r["min_s1"]) for r in corr], 0.01))
        emit(f"  q01 min_s1 correct-unanimous {q_u:.4f} (n {len(cu)}); all correct {q_c:.4f} (n {len(corr)})")
        emit(f"  fraction ≥0.68: correct {np.mean([float(r['min_s1'])>=0.68 for r in corr]):.4f}; "
             f"correct-unanimous {np.mean([float(r['min_s1'])>=0.68 for r in cu]):.4f}")
        rejected = [r for r in reg_ok if float(r["min_s1"]) < 0.65]
        emit("  gate rejects correct registry: " +
             ", ".join(f"{r['query']} {r['chemistry']} {float(r['min_s1']):.3f}" for r in rejected))

    report("plain", plain_rows)
    report("flex", flex_rows)
    abc = {r["query"]: float(r["min_s1"]) for r in plain_rows if r["query"] in ("A", "B", "C")}
    emit(f"A/B/C min_s1 { {k: round(v, 4) for k, v in abc.items()} }")
    emit(f"n_valid_all A {next(r['n_valid_all'] for r in plain_rows if r['query']=='A')}")
    from collections import Counter
    emit("plain dataset counts " + str(dict(Counter(r["dataset"] for r in plain_rows))))

    emit()
    emit("=== 6. 548 earlier-family figure is a different family (quoted, not re-derived) ===")
    emit("548 RESULTS: earlier family 222–319p2 DE-overall rms 0.00321; no emission_seed pair (99/99 are 1000003).")
    emit("537 §2.6 restates the seed count and treats bulk-only pairs as re-rolls via the stream shift.")
    OUT.write_text("\n".join(lines) + "\n")
    emit(f"wrote {OUT}")


if __name__ == "__main__":
    main()
