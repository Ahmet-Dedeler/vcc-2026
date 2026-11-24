#!/usr/bin/env python3
"""553 — the ACCESS: score the emitted arms of a surface against its truth, with the surface's own anchors and noise band.

THIS OPENS A TRUTH FILE. On r7b it runs only with `--owner-go <file>` (the owner's preregistered claim; its sha256 is
recorded in the output and re-verified by every phase) and never before the owner says so (RESERVE.md sections 4-6; the
ledger's R7 section). On h1 (not a reserve) it is the end-to-end test; on h1r7 it is the R7b-shaped dress rehearsal.

The orchestrator (no --phase) checks the inputs, writes plan.json, then runs each phase as its OWN subprocess (macOS keeps
a process's freed large allocations in its RSS, so phases release memory only by exiting; h1-test's phases peak ~3-6 GB):
  ctrl   the control row of the reference: exact integer sums, count, jackknife (control rows only)
  ref    the truth reference (evaluate323.reference generalised: any control count, any truth cells per target), cross-
         checked against the prediction side (same gate, same control sums, same targets); KA-R on h1
  score  every arm x seed's stats (written by run553.py before the access) scored with evaluate323.score's calls
         (lib553.score_stats); KA-S on h1 (== the recorded metrics json)
  split  the reference's split-half re-draws: per split k (seeds SeedSequence(553).generate_state(n)) each target's truth
         cells are permuted and cut into two disjoint halves (floor(n/2) each, file order within a half), each half's
         statistics against the FULL control block, every arm x seed scored against half A and half B as the truth
  abase  anchors323's full stats + dispersed seed-0 baseline (cell-eval2's own profile, emission and fast_score)
  rh0A … rh4B  anchors323's 5-split replicate, one half per process (only that half's rows are read)
  arep   the replicate's per-split scores from the halves; KA-A on h1 (== 323's h1-anchors.json)
  report per arm and member: the paired effect vs the incumbent `cpm` (same seed), its seed sd, its split-half sd (/sqrt 2
         as the full-truth estimate), a suggested band sd, the target-bootstrap CI of the overall-equivalent, and the
         overall-equivalent (mean of the six members scaled by cell-eval2's score_one at the surface's own anchors)
Outputs in work/access553_<surface>/ (access553.json, access553.md, the phases' files).

usage (h1, the test):
  uv run python experiments/553-r7b-harness/access553.py h1 --arms lcap-kept,... --seeds 0,1,2 --stats-map h1test
usage (r7b, OWNER ONLY, after the claim is preregistered and its sha256 is in the ledger):
  uv run python experiments/553-r7b-harness/access553.py r7b --owner-go <PREREGISTRATION file> --arms <arm>,... --seeds 0,1,2
usage (a floored R7b, r7bd<floor>, OWNER ONLY; needs 552's floored truth in
  experiments/552-r7-reserve/work/r7b_floor<floor>/{truth_r7b.h5ad, build_r7b.json}, refused without it):
  uv run python experiments/553-r7b-harness/access553.py r7bd8000 --owner-go <file> --exclude-targets BST2 --arms <arm>,... --seeds 0,1,2
readiness check that never opens the truth:
  uv run python experiments/553-r7b-harness/access553.py r7bd8000 --dry-run --exclude-targets BST2 --arms <arm>,... --seeds 0,1,2
(an option value that begins with a dash needs the `=` form: --tag=-flex)
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import lib553 as S  # noqa: E402

TAG = "base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse"
INCUMBENT = "cpm"
N_ANCHOR_SPLITS = 5
HALF_PHASES = tuple(f"rh{i}{h}" for i in range(N_ANCHOR_SPLITS) for h in "AB")
PHASES = ("ctrl", "ref", "score", "split", "csplit", "abase", *HALF_PHASES, "arep", "report")


# ---------------------------------------------------------------------------------------------------------------- helpers
def stats_path(surface: str, mode: str, seed: int, zs: bool, stats_map: str | None, tag: str = ""):
    """Where an arm's stats live. Default: this directory's work/<surface>_<mode><tag>/ (tag "-flex" for run553's flex class
    rule). `--stats-map h1test` also reads the recorded H1 real-emitter runs (435 base, 537, 528), which KA-H1 shows this
    harness reproduces bit for bit."""
    t = TAG.format(seed=seed) + ("_dnk250zs" if zs else "")
    cands = [HERE / "work" / f"{surface}_{mode}{tag}"]
    if stats_map == "h1test":
        if mode == "cpm":
            cands.append(ROOT / "experiments/435-zs-on-coarse/work/h1")
        cands += [ROOT / "experiments/537-arm-pipeline/work" / f"h1_{mode}",
                  ROOT / "experiments/528-callprop-block/work" / f"h1_{mode}"]
    for d in cands:
        p, m = d / f"stats-{t}.npz", d / f"metrics-{t}.json"
        if p.exists():
            return p, (m if m.exists() else None)
    return None, None


def load_stats(p: Path):
    z = np.load(p, allow_pickle=True)
    return dict(targets=z["targets"].astype(str), sums=np.asarray(z["sums"]), jk=np.asarray(z["jk"], float),
                p=np.asarray(z["p"], float), q=np.asarray(z["q"], float), lfc=np.asarray(z["lfc"], float))


def all_stats(plan):
    return {(arm, s): load_stats(ROOT / plan["stats"][f"{arm}|s{s}"]["path"]) for arm in plan["arms"] for s in plan["seeds"]}


def per_target(rows, targets):
    out = {m: np.full(len(targets), np.nan) for m in ("pds", "nmae", "fid", "reach", "jac")}
    ix = {str(t): i for i, t in enumerate(targets)}
    inv = {v: k for k, v in S.FULL.items()}
    for r in rows:
        m = inv.get(r["metric"])
        if m in out and str(r["perturbation"]) in ix:
            out[m][ix[str(r["perturbation"])]] = r["value"]
    return out


def load_ref(outdir: Path):
    import polars as pl
    with np.load(outdir / "reference_stats.npz", allow_pickle=False) as z:
        ref = {k: z[k] for k in z.files}
    return ref, pl.read_parquet(outdir / "reference_de.parquet")


def half_reference(T, targets, c, full_ref, seed_k: int):
    """One split: two disjoint half-truths (A, B) per target, each against the FULL control block (the prediction's DE was
    computed against the same controls). Returns [(refA, deA), (refB, deB)], each with the full reference's control row."""
    import evaluate323 as E
    rng = np.random.default_rng(seed_k)
    halves = ([], [])
    cells = ([], [])
    lab, xt, pos = T["lab"], T["xt"], T["pos"]
    for t in targets:
        rows = np.flatnonzero(lab == t)
        perm = rng.permutation(rows)
        h = perm.size // 2
        for j, part in enumerate((perm[:h], perm[h:2 * h])):
            halves[j].append(E.block_stats(xt[pos[np.sort(part)]], c))
            cells[j].append(int(part.size))
    out = []
    for j in (0, 1):
        s, jk, p, q, l = (list(v) for v in zip(*halves[j]))
        ref = dict(genes=full_ref["genes"], targets=full_ref["targets"], perts=full_ref["perts"],
                   sums=np.vstack([*s, full_ref["sums"][-1]]), counts=np.array([*cells[j], full_ref["counts"][-1]]),
                   jk=np.array([*jk, full_ref["jk"][-1]]))
        out.append((ref, E.de_frame(full_ref["targets"], full_ref["genes"][c.keep], p, q, l)))
    return out


def phase_record(outdir: Path, phase: str, rec: dict, t0: float) -> None:
    rec.update(phase=phase, seconds=round(time.monotonic() - t0, 1), peak_rss_gb=S.rss_gb())
    S.write_json(outdir / f"phase_{phase}.json", rec)
    S.log(f"phase {phase} done in {rec['seconds']} s, peak RSS {rec['peak_rss_gb']:.2f} GB")


# ---------------------------------------------------------------------------------------------------------------- phases
def run_phase(phase: str, a, spec: dict, outdir: Path, plan: dict) -> None:
    t0 = time.monotonic()
    if spec["sealed"]:
        if not a.owner_go or S.sha256(Path(a.owner_go)) != plan["owner_go"]["sha256"]:
            raise SystemExit("sealed surface: the phase's --owner-go file does not match the plan's")
        S.arm_access(spec["name"])
    cr = S.load_controlref(spec, need_cs=False)
    exclude = tuple(plan.get("exclude", []))
    targets = cr["targets"].astype(str)
    targets = targets[~np.isin(targets, list(exclude))]
    rec = {}
    if phase == "ctrl":
        T = S.truth_rows(spec, xt=False, exclude=exclude)
        cp, n, cj = S.control_moments(T, expect_ctrl_sums=cr["ctrl_sums"])
        np.savez(outdir / "ctrl_moments.npz", sums=cp, n=np.int64(n), jk=np.float64(cj))
        rec.update(n_controls=n, jk=cj)
    elif phase == "ref":
        T = S.truth_rows(spec, exclude=exclude)
        with np.load(outdir / "ctrl_moments.npz") as z:
            cm = (z["sums"], int(z["n"]), float(z["jk"]))
        ref, c, de = S.truth_reference(spec, T, cm, expect_keep=cr["keep"], expect_targets=targets)
        de.write_parquet(outdir / "reference_de.parquet")
        np.savez(outdir / "reference_stats.npz", **ref)
        n_real = ((np.asarray(ref["q"]) < .05) & ~np.equal.outer(ref["targets"], ref["genes"][c.keep])).sum(1)
        rec.update(n_targets=len(targets), n_controls=int(ref["counts"][-1]),
                   truth_cells_per_target=dict(zip(map(str, targets), map(int, ref["counts"][:-1]))),
                   n_real_quantiles=np.percentile(n_real, [0, 10, 25, 50, 75, 90, 100]).tolist(),
                   n_real_by_target=dict(zip(map(str, targets), map(int, n_real))))
        if spec["name"] == "h1":
            import polars as pl
            with np.load(S.E323 / "work/reference.npz", allow_pickle=False) as z:
                kar = {k: bool(np.array_equal(z[k], ref[k])) for k in ("sums", "counts", "jk", "p", "q", "lfc", "keep")}
            kar["reference_de"] = bool(pl.read_parquet(S.E323 / "work/reference_de.parquet").equals(de))
            rec["ka_reference"] = kar
            S.log(f"KA-R (truth reference == 323's reference.npz / reference_de.parquet): {kar}")
            if not all(kar.values()):
                phase_record(outdir, phase, rec, t0)
                raise SystemExit("KA-R FAILED")
    elif phase == "score":
        ref, de = load_ref(outdir)
        stats = all_stats(plan)
        raw, pt, ka_s = {}, {}, {}
        for (arm, s), st in stats.items():
            r_, rows, _ = S.score_stats(ref, ref["keep"], de, st["sums"], st["jk"], st["p"], st["q"], st["lfc"])
            key = f"{arm}|s{s}"
            raw[key], pt[key] = r_, per_target(rows, targets)
            mj = plan["stats"][key].get("metrics")
            if spec["name"] == "h1" and mj:
                rr = json.loads((ROOT / mj).read_text())["raw"]
                ka_s[key] = {m: abs(float(rr[m]) - r_[m]) for m in S.MEMBERS}
        S.write_json(outdir / "scores_full.json", raw)
        np.savez(outdir / "pertarget_full.npz", **{f"{k}|{m}": v[m] for k, v in pt.items() for m in v})
        if spec["name"] == "h1":
            mx = max(max(v.values()) for v in ka_s.values()) if ka_s else None
            rec.update(ka_scores_absdiff=ka_s, ka_scores_max_absdiff=mx)
            S.log(f"KA-S (scores == recorded metrics json): {len(ka_s)} files, max |d| {mx}")
            if mx is None or mx > 1e-12:
                phase_record(outdir, phase, rec, t0)
                raise SystemExit("KA-S FAILED")
    elif phase == "split":
        import exact_de323
        T = S.truth_rows(spec, exclude=exclude)
        c = exact_de323.control_block(T["xc"])
        T["xc"] = None                  # the halves read target rows only (h1-test: -0.94 GB in this phase)
        ref, _ = load_ref(outdir)
        if not np.array_equal(c.keep, ref["keep"]):
            raise AssertionError("control block differs from the reference's")
        stats = all_stats(plan)
        split_seeds = [int(v) for v in np.random.SeedSequence(553).generate_state(a.splits)]
        hraw = {}
        hi = 0
        for k_i, sk in enumerate(split_seeds):
            for j, (href, hde) in enumerate(half_reference(T, targets, c, ref, sk)):
                for (arm, s), st in stats.items():
                    r_, _, _ = S.score_stats(href, c.keep, hde, st["sums"], st["jk"], st["p"], st["q"], st["lfc"])
                    hraw[f"{hi}|{arm}|s{s}"] = r_
                S.log(f"split {k_i + 1}/{a.splits} half {'AB'[j]} scored ({time.monotonic() - t0:.1f}s, peak RSS {S.rss_gb():.2f} GB)")
                hi += 1
        S.write_json(outdir / "scores_halves.json", hraw)
        rec.update(split_seeds=split_seeds, n_half=hi)
    elif phase == "csplit":
        # 553 (10-01): the split-half re-draw above halves TARGET cells against the full, fixed control block, so it
        # cannot see the control sample's own noise -- which a second R7b-shaped rehearsal draw showed dominates (RESULTS).
        # --ctrl-splits N halves the CONTROLS too: per split, a random half of the controls (control block on the full
        # reference's gene universe, its own control row) and a random half of each target's cells; every arm x seed is
        # scored against both halves. The prediction side stays as emitted (its control-sample dependence is not re-drawn).
        if a.ctrl_splits <= 0:
            rec.update(skipped=True)
        else:
            import evaluate323 as E
            from cell_eval2.moments import jackknife_correction
            T = S.truth_rows(spec, exclude=exclude)
            ref, _ = load_ref(outdir)
            keep = np.asarray(ref["keep"], dtype=bool)
            stats = all_stats(plan)
            cseeds = [int(v) for v in np.random.SeedSequence(5531).generate_state(a.ctrl_splits)]
            craw = {}
            hi = 0
            lab, xt, pos, xc = T["lab"], T["xt"], T["pos"], T["xc"]
            for k_i, sk in enumerate(cseeds):
                rng = np.random.default_rng(sk)
                cperm = rng.permutation(xc.shape[0])
                hc = cperm.size // 2
                tparts = {}
                for t in targets:
                    rows = rng.permutation(np.flatnonzero(lab == t))
                    h = rows.size // 2
                    tparts[t] = (np.sort(rows[:h]), np.sort(rows[h:2 * h]))
                for j, cpart in enumerate((np.sort(cperm[:hc]), np.sort(cperm[hc:2 * hc]))):
                    xc_j = xc[cpart]
                    c_j = S.control_block_keep(xc_j, keep)
                    csum = S.int_colsums(xc_j).astype(np.float64)
                    cj = float(jackknife_correction(xc_j, np.zeros(xc_j.shape[0], dtype=np.int64), 1, S.TS, chunk=32)[0])
                    hs = [E.block_stats(xt[pos[tparts[t][j]]], c_j) for t in targets]
                    s_, jk_, p_, q_, l_ = (list(v) for v in zip(*hs))
                    href = dict(genes=ref["genes"], targets=ref["targets"], perts=ref["perts"],
                                sums=np.vstack([*s_, csum]), counts=np.array([*[tparts[t][j].size for t in targets], hc]),
                                jk=np.array([*jk_, cj]))
                    hde = E.de_frame(ref["targets"], ref["genes"][keep], p_, q_, l_)
                    for (arm, s), st in stats.items():
                        r_, _, _ = S.score_stats(href, keep, hde, st["sums"], st["jk"], st["p"], st["q"], st["lfc"])
                        craw[f"{hi}|{arm}|s{s}"] = r_
                    S.log(f"csplit {k_i + 1}/{a.ctrl_splits} half {'AB'[j]} scored ({hc} controls; "
                          f"{time.monotonic() - t0:.1f}s, peak RSS {S.rss_gb():.2f} GB)")
                    hi += 1
            S.write_json(outdir / "scores_csplit.json", craw)
            rec.update(ctrl_split_seeds=cseeds, n_csplit=hi)
    elif phase == "abase":
        import exact_de323
        T = S.truth_rows(spec, exclude=exclude)
        c = exact_de323.control_block(T["xc"])
        with np.load(outdir / "ctrl_moments.npz") as z:
            cm = (z["sums"], int(z["n"]), float(z["jk"]))
        base, full, bdiag = S.anchors_base(spec, T, ctrl_block=c, ctrl_moments=cm)
        full["de"].write_parquet(outdir / "full_de.parquet")
        rec.update(baseline=base, baseline_emission=bdiag["baseline_emission"])
    elif phase in HALF_PHASES:
        # one replicate half in its own process: only that half's rows are read (each half gates its own control block,
        # so its size differs from every other half's and macOS would not reuse the freed regions within one process)
        i, h = int(phase[2]), phase[3]
        ga, gb = S.replicate_halves(spec, i, n_splits=N_ANCHOR_SPLITS, exclude=exclude)
        st = S.replicate_half_stats(spec, ga if h == "A" else gb)
        st["de"].write_parquet(outdir / f"rhalf_{i}{h}_de.parquet")
        np.savez(outdir / f"rhalf_{i}{h}.npz", **{k: st[k] for k in ("perts", "genes", "sums", "counts", "jk")})
    elif phase == "arep":
        import polars as pl
        from cell_eval2.anchor import _derive_seeds
        full_de = pl.read_parquet(outdir / "full_de.parquet")
        splits = []
        for i in range(N_ANCHOR_SPLITS):
            half = {}
            for h in "AB":
                with np.load(outdir / f"rhalf_{i}{h}.npz", allow_pickle=False) as z:
                    half[h] = {k: z[k] for k in z.files}
                half[h]["de"] = pl.read_parquet(outdir / f"rhalf_{i}{h}_de.parquet")
            raw, cohort = S.replicate_score(half["A"], half["B"], full_de)
            splits.append({"split": i, "seed": int(_derive_seeds(0, N_ANCHOR_SPLITS)[i]), "raw": raw,
                           "nmae_full_gate_cohort": cohort})
            S.log(f"anchors: split {i + 1}/{N_ANCHOR_SPLITS} {json.dumps(raw)}")
        rep = {m: float(np.mean([sp_["raw"][m] for sp_ in splits])) for m in S.MEMBERS}
        rec.update(replicate=rep, replicate_splits=splits)
        if spec["name"] == "h1":
            base = json.loads((outdir / "phase_abase.json").read_text())["baseline"]
            want = json.loads((S.E323 / "work/h1-anchors.json").read_text())
            kaa = {f"baseline.{m}": abs(base[m] - want["baseline"][m]) for m in S.MEMBERS}
            kaa.update({f"replicate.{m}": abs(rep[m] - want["replicate"][m]) for m in S.MEMBERS})
            for sp_, sw in zip(splits, want["replicate_splits"]):
                kaa.update({f"split{sp_['split']}.{m}": abs(sp_["raw"][m] - sw["raw"][m]) for m in S.MEMBERS})
            rec.update(ka_anchors_absdiff=kaa, ka_anchors_max_absdiff=max(kaa.values()))
            S.log(f"KA-A (anchors == 323's h1-anchors.json): max |d| {rec['ka_anchors_max_absdiff']:.3e}")
            if rec["ka_anchors_max_absdiff"] > 1e-12:
                phase_record(outdir, phase, rec, t0)
                raise SystemExit("KA-A FAILED")
    elif phase == "report":
        report(a, spec, outdir, plan, targets)
    phase_record(outdir, phase, rec, t0)


def report(a, spec, outdir, plan, targets):
    raw_s = json.loads((outdir / "scores_full.json").read_text())
    hraw_s = json.loads((outdir / "scores_halves.json").read_text())
    base = json.loads((outdir / "phase_abase.json").read_text())["baseline"]
    rep = json.loads((outdir / "phase_arep.json").read_text())["replicate"]
    seeds, arms = plan["seeds"], plan["arms"]
    n_half = json.loads((outdir / "phase_split.json").read_text())["n_half"]
    with np.load(outdir / "pertarget_full.npz") as z:
        pt = {k: z[k] for k in z.files}
    raw = {(k.split("|")[0], int(k.split("|s")[1])): v for k, v in raw_s.items()}
    hraw = {(int(k.split("|")[0]), k.split("|")[1], int(k.split("|s")[1])): v for k, v in hraw_s.items()}
    scl = {k: S.scaled(r_, base, rep) for k, r_ in raw.items()}
    hscl = {k: S.scaled(r_, base, rep) for k, r_ in hraw.items()}
    cpath = outdir / "scores_csplit.json"
    cscl, n_c = {}, 0
    if cpath.exists():
        craw_s = json.loads(cpath.read_text())
        if craw_s:
            craw_ = {(int(k.split("|")[0]), k.split("|")[1], int(k.split("|s")[1])): v for k, v in craw_s.items()}
            cscl = {k: S.scaled(r_, base, rep) for k, r_ in craw_.items()}
            n_c = 1 + max(k[0] for k in craw_)
    rng = np.random.default_rng(55301)
    T = len(targets)
    boots = rng.integers(0, T, (a.boot, T))
    den = {m: rep[m] - base[m] for m in S.MEMBERS}

    def paired(arm, ref_arm):
        """The paired effect of `arm` vs `ref_arm` (same seed): per member raw and scaled, seed sd, split-half sd, the
        suggested band sd, and the target-bootstrap CI of the overall-equivalent."""
        out = {"effect": {}, "effect_scaled": {}}
        for m in (*S.MEMBERS, "overall"):
            if m != "overall":
                d = np.array([raw[(arm, s)][m] - raw[(ref_arm, s)][m] for s in seeds])
                hd = np.array([np.mean([hraw[(h, arm, s)][m] - hraw[(h, ref_arm, s)][m] for s in seeds])
                               for h in range(n_half)])
                wd = hd[0::2] - hd[1::2]
                out["effect"][m] = dict(mean=float(d.mean()), seed_values=d.tolist(),
                                        seed_sd=float(d.std(ddof=1)) if len(d) > 1 else None,
                                        half_values=hd.tolist(), split_half_sd=float(hd.std(ddof=1)),
                                        split_half_sd_within=float(np.sqrt(np.mean(wd ** 2) / 2)))
            ds = np.array([scl[(arm, s)][m] - scl[(ref_arm, s)][m] for s in seeds])
            hds = np.array([np.mean([hscl[(h, arm, s)][m] - hscl[(h, ref_arm, s)][m] for s in seeds])
                            for h in range(n_half)])
            sd_seed = float(ds.std(ddof=1)) if len(ds) > 1 else None
            sd_full = float(hds.std(ddof=1) / np.sqrt(2))
            out["effect_scaled"][m] = dict(mean=float(ds.mean()), seed_values=ds.tolist(), seed_sd=sd_seed,
                                           half_values=hds.tolist(), split_half_sd=float(hds.std(ddof=1)),
                                           split_half_sd_full_truth_est=sd_full,
                                           # a SUGGESTED band for the claim (not a reading rule): the sd of the
                                           # seed-averaged paired effect from both re-roll sources
                                           band_sd=(float(np.sqrt(sd_seed ** 2 / len(ds) + sd_full ** 2))
                                                    if sd_seed is not None else None))
            if n_c:
                # the cells + controls re-draw (csplit): the same paired effect on half-truths with half the controls
                cds = np.array([np.mean([cscl[(h, arm, s)][m] - cscl[(h, ref_arm, s)][m] for s in seeds])
                                for h in range(n_c)])
                csd_full = float(cds.std(ddof=1) / np.sqrt(2))
                out["effect_scaled"][m].update(
                    csplit_values=cds.tolist(), csplit_sd=float(cds.std(ddof=1)), csplit_sd_full_truth_est=csd_full,
                    band_sd_ctrl=(float(np.sqrt(sd_seed ** 2 / len(ds) + csd_full ** 2)) if sd_seed is not None else None))
        # target bootstrap of the overall-equivalent's per-target members (seed-averaged, paired); mse at its value
        dpt = {m: np.nanmean([pt[f"{arm}|s{s}|{m}"] - pt[f"{ref_arm}|s{s}|{m}"] for s in seeds], axis=0)
               for m in ("pds", "nmae", "fid", "reach", "jac")}
        bs = [np.mean([np.nanmean(dpt[m][bi]) / den[m] for m in ("pds", "nmae", "fid", "reach", "jac")]
                      + [out["effect_scaled"]["mse"]["mean"]]) for bi in boots]
        out["overall_boot_ci95"] = [float(np.quantile(bs, .025)), float(np.quantile(bs, .975))]
        out["targets_improved"] = {m: int(np.sum(np.nan_to_num(dpt[m]) * (-1 if m == "nmae" else 1) > 0)) for m in dpt}
        return out

    arms_out = {}
    for arm in arms:
        row = {"seeds": seeds, "raw_mean": {m: float(np.mean([raw[(arm, s)][m] for s in seeds])) for m in S.MEMBERS},
               "scaled_mean": {m: float(np.mean([scl[(arm, s)][m] for s in seeds])) for m in (*S.MEMBERS, "overall")},
               "effect": {}, "effect_scaled": {}}
        if arm != INCUMBENT:
            row.update(paired(arm, INCUMBENT))
        arms_out[arm] = row
    # --pairs A:B,...: switch-vs-switch effects (e.g. the hook on top of CP: lcap-kept+c4h:lcap-kept), same machinery;
    # the same bootstrap draws are reused, so a pair's CI is paired with the arms'
    pairs_out = {}
    for pr in [x for x in (getattr(a, "pairs", "") or "").split(",") if x]:
        A, B = pr.split(":")
        if A not in arms or B not in arms:
            raise SystemExit(f"--pairs {pr}: both arms must be in --arms")
        pairs_out[f"{A} vs {B}"] = paired(A, B)
    rec = dict(plan=plan, anchors=dict(baseline=base, replicate=rep),
               reference=json.loads((outdir / "phase_ref.json").read_text()),
               phases={p: json.loads((outdir / f"phase_{p}.json").read_text()) for p in PHASES[:-1]
                       if (outdir / f"phase_{p}.json").exists()},
               scaled_crosscheck_max_absdiff=max(abs(scl[k][m] - S.scaled_simple(raw[k], base, rep)[m])
                                                 for k in scl for m in (*S.MEMBERS, "overall")),
               arms=arms_out, n_half=n_half)
    if pairs_out:
        rec["pairs"] = pairs_out
    retro = {}
    if getattr(a, "retro", ""):
        # board-retrodiction rows (552's requirement): each raw member's effect here over the same effect in a reference
        # access (H1's price: work/access553_h1), for every arm vs `cpm` and every pair present in both
        rp = Path(a.retro)
        rp = rp if rp.is_absolute() else HERE / rp
        rref = json.loads((rp / "access553.json").read_text())
        both = {f"{arm} vs cpm": (arms_out[arm], rref["arms"][arm]) for arm in arms_out
                if arm != INCUMBENT and arm in rref["arms"] and arms_out[arm]["effect"]}
        both.update({k: (v, rref.get("pairs", {})[k]) for k, v in pairs_out.items() if k in rref.get("pairs", {})})
        for k, (mine, theirs) in both.items():
            retro[k] = {m: dict(here=mine["effect"][m]["mean"], reference=theirs["effect"][m]["mean"],
                                ratio=(mine["effect"][m]["mean"] / theirs["effect"][m]["mean"]
                                       if abs(theirs["effect"][m]["mean"]) > 1e-12 else None))
                        for m in ("nmae", "fid", "reach", "jac")}
        rec["retrodiction"] = dict(reference=str(rp), rows=retro)
    S.write_json(outdir / "access553.json", rec)
    lines = [f"# access553 — {spec['name']} ({plan['started']}), form {plan['form']}, seeds {seeds}, {n_half} half-truths",
             "", f"anchors: baseline {json.dumps({m: round(base[m], 5) for m in S.MEMBERS})}",
             f"         replicate {json.dumps({m: round(rep[m], 5) for m in S.MEMBERS})}", "",
             "| arm | member | raw effect vs cpm (seed mean) | seed sd | split-half sd | scaled effect | scaled seed sd | "
             "scaled split-half sd (half / full-est) | band sd (suggested) |" + (" band sd, cells + controls |" if n_c else ""),
             "|---|---|---:|---:|---:|---:|---:|---:|---:|" + ("---:|" if n_c else "")]
    for arm, row in arms_out.items():
        if arm == INCUMBENT:
            continue
        for m in (*S.MEMBERS, "overall"):
            es, e = row["effect_scaled"][m], row["effect"].get(m)
            lines.append(f"| {arm} | {m} | {'' if e is None else f'{e['mean']:+.5f}'} | "
                         f"{'' if e is None or e['seed_sd'] is None else f'{e['seed_sd']:.5f}'} | "
                         f"{'' if e is None else f'{e['split_half_sd']:.5f}'} | {es['mean']:+.5f} | "
                         f"{'' if es['seed_sd'] is None else f'{es['seed_sd']:.5f}'} | "
                         f"{es['split_half_sd']:.5f} / {es['split_half_sd_full_truth_est']:.5f} | "
                         f"{'' if es['band_sd'] is None else f'{es['band_sd']:.5f}'} |"
                         + (f" {es['band_sd_ctrl']:.5f} |" if n_c and es.get('band_sd_ctrl') is not None else (" |" if n_c else "")))
        lines.append(f"| {arm} | overall target-bootstrap 95% | | | | [{row['overall_boot_ci95'][0]:+.5f}, "
                     f"{row['overall_boot_ci95'][1]:+.5f}] | | | |" + (" |" if n_c else ""))
    if pairs_out:
        lines += ["", "| pair (switch on top of a recipe) | overall-equiv effect | seed sd | split-half sd (full-est) | band sd | "
                  "effect / band sd | target-bootstrap 95 % | raw nmae / fid / reach / jac |" + (" band sd, cells + controls |" if n_c else ""),
                  "|---|---:|---:|---:|---:|---:|---|---|" + ("---:|" if n_c else "")]
        for name, pr in pairs_out.items():
            o = pr["effect_scaled"]["overall"]
            # one seed (a control-half surface): no seed sd and no band
            seed_s = "" if o["seed_sd"] is None else f"{o['seed_sd']:.5f}"
            band_s = "" if o["band_sd"] is None else f"{o['band_sd']:.5f}"
            ratio_s = "" if o["band_sd"] is None else f"{o['mean'] / o['band_sd']:+.1f}"
            lines.append(f"| {name} | {o['mean']:+.5f} | {seed_s} | {o['split_half_sd_full_truth_est']:.5f} | "
                         f"{band_s} | {ratio_s} | [{pr['overall_boot_ci95'][0]:+.5f}, "
                         f"{pr['overall_boot_ci95'][1]:+.5f}] | "
                         + " / ".join(f"{pr['effect'][m]['mean']:+.4f}" for m in ("nmae", "fid", "reach", "jac")) + " |"
                         + (f" {o['band_sd_ctrl']:.5f} |" if n_c and o.get('band_sd_ctrl') is not None else (" |" if n_c else "")))
    if retro:
        lines += ["", f"retrodiction: raw member effect here / the same effect in {Path(a.retro).name} (H1's price)", "",
                  "| effect | nmae here / ref (ratio) | fid here / ref (ratio) | reach here / ref (ratio) | jac here / ref (ratio) |",
                  "|---|---|---|---|---|"]
        def fmt(d):
            ratio = "n/a" if d["ratio"] is None else f"{d['ratio']:+.2f}"
            return f"{d['here']:+.4f} / {d['reference']:+.4f} ({ratio})"

        for k, row in retro.items():
            lines.append(f"| {k} | " + " | ".join(fmt(row[m]) for m in ("nmae", "fid", "reach", "jac")) + " |")
    lines += ["", "| arm | " + " | ".join(f"{m} raw (scaled)" for m in S.MEMBERS) + " | overall |",
              "|---|" + "---:|" * (len(S.MEMBERS) + 1)]
    for arm, row in arms_out.items():
        lines.append(f"| {arm} | " + " | ".join(f"{row['raw_mean'][m]:.5f} ({row['scaled_mean'][m]:.4f})" for m in S.MEMBERS)
                     + f" | {row['scaled_mean']['overall']:.5f} |")
    lines += ["", "| phase | seconds | peak RSS GB |", "|---|---:|---:|"]
    lines += [f"| {p} | {v['seconds']} | {v['peak_rss_gb']:.2f} |" for p, v in rec["phases"].items()]
    (outdir / "access553.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)


# ---------------------------------------------------------------------------------------------------------------- main
LEDGER = ROOT / "docs/CONFIRMATION-LEDGER.md"
R7_HEADING = "## Reserve R7"


def r7_ledger_section(text: str) -> str:
    """The ledger's `## Reserve R7` section: from its heading to the next level-2 heading (a `### ` line does not end it).
    555 (2026-10-02): a sealed access arms only if the claim's sha256 is recorded HERE, not merely somewhere in the file."""
    text = "\n" + text                                       # a heading on the file's first line counts too
    i = text.find("\n" + R7_HEADING)
    if i < 0:
        return ""
    j = text.find("\n## ", i + 1)
    return text[i:j] if j >= 0 else text[i:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("surface", help="h1 | h1r7 | h1t7 | r7b | r7bd<floor> (R7b with a human-library floor) | fx")
    ap.add_argument("--owner-go", default=None, help="r7b only: the owner's preregistered claim file (its sha is recorded)")
    ap.add_argument("--arms", required=True, help="comma-separated MODEs; the incumbent `cpm` is always included")
    ap.add_argument("--seeds", default="0,1,2")
    ap.add_argument("--splits", type=int, default=5)
    ap.add_argument("--stats-map", default=None, choices=[None, "h1test"])
    ap.add_argument("--tag", default="", help="stats directory suffix, e.g. -flex (run553's flex class rule)")
    ap.add_argument("--exclude-targets", default="", help="targets not emitted (run553 --exclude-targets); their truth cells "
                                                          "are treated as absent; the tag gains -x<targets> as run553's")
    ap.add_argument("--form", default="zs", choices=["zs", "pre"], help="score the zs (446bj2z's DE form) or pre-zs stats")
    ap.add_argument("--boot", type=int, default=4000)
    ap.add_argument("--ctrl-splits", type=int, default=0,
                    help="csplit phase: N splits halving the CONTROLS as well as each target's cells (0 = off, the "
                         "mandate's band only); the report then adds a cells + controls band (band_sd_ctrl)")
    ap.add_argument("--retro", default="", help="report only: a reference access directory (e.g. work/access553_h1, H1's "
                                                "price); adds each raw member's effect ratio here / there (board retrodiction)")
    ap.add_argument("--halves", default="", help="after the access, the control-half accesses of splits k,... and band553 "
                                                "(e.g. 0,1,2; the halves' arms must be emitted at seed 0 beforehand)")
    ap.add_argument("--half-arms", default="lcap-kept,cpm+c4h,lcap-kept+c4h,cpm-n1600,lcap-kept-n1600+c4h")
    ap.add_argument("--half-pairs", default="lcap-kept+c4h:lcap-kept,lcap-kept-n1600+c4h:lcap-kept+c4h")
    ap.add_argument("--pairs", default="", help="report only: switch-vs-switch paired effects A:B,... (both in --arms), "
                                                "e.g. lcap-kept+c4h:lcap-kept (the hook on top of CP)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--dry-run", action="store_true",
                    help="readiness only: every stats file present and hashed, the bundle vs build_r7b.json (the truth by sha "
                         "only); the truth is NOT opened")
    ap.add_argument("--phase", default=None, choices=PHASES, help="internal: run one phase (the orchestrator calls this)")
    a = ap.parse_args()
    import re
    if not (a.surface in ("h1", "r7b", "fx") or a.surface.startswith(("h1r7", "h1t7")) or re.fullmatch(S.CTRL_HALF_RE, a.surface)
            or (a.surface.startswith("r7bd") and a.surface[4:].isdigit())):
        ap.error(f"unknown surface {a.surface!r}")
    spec = S.surface_spec(a.surface)
    xtag = ("-x" + "+".join(sorted(x for x in a.exclude_targets.split(",") if x))) if a.exclude_targets else ""
    outdir = Path(a.out) if a.out else HERE / "work" / f"access553_{a.surface}{a.tag}{xtag}"
    outdir.mkdir(parents=True, exist_ok=True)
    if a.phase:
        return run_phase(a.phase, a, spec, outdir, json.loads((outdir / "plan.json").read_text()))
    t_start = time.monotonic()
    seeds = [int(s) for s in a.seeds.split(",")]
    arms = [INCUMBENT] + [m for m in a.arms.split(",") if m and m != INCUMBENT]
    exclude = tuple(sorted(x for x in a.exclude_targets.split(",") if x))
    for pr in [x for x in a.pairs.split(",") if x]:            # 555: checked before any phase opens the truth
        A_, sep, B_ = pr.partition(":")
        if not sep or A_ not in arms or B_ not in arms or A_ == B_:
            raise SystemExit(f"--pairs {pr!r}: each pair is A:B, two different arms, both in --arms (cpm included)")
    plan = dict(surface=a.surface, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), argv=sys.argv[1:],
                arms=arms, seeds=seeds, splits=a.splits, form=a.form, tag=a.tag, exclude=list(exclude), stats={})
    if spec["sealed"] and not a.dry_run:
        if not a.owner_go or not Path(a.owner_go).is_file():
            raise SystemExit(f"{a.surface} is SEALED: the access needs --owner-go <the owner's preregistered claim file>")
        plan["owner_go"] = dict(path=str(Path(a.owner_go).resolve()), sha256=S.sha256(Path(a.owner_go)))
        if plan["owner_go"]["sha256"] not in r7_ledger_section(LEDGER.read_text()):
            raise SystemExit(f"{a.surface} is SEALED: the claim's sha256 {plan['owner_go']['sha256']} is not in the ledger's "
                             f"'{R7_HEADING}' section ({LEDGER}): record it there first (555)")
        plan["owner_go"]["in_ledger_r7_section"] = True
    if spec["sealed"]:
        plan["bundle_vs_build_r7b"] = S.check_r7b_outputs(check_truth_sha=True)
    if spec.get("depth_floor"):
        # a floored prediction side is scored only against 552's floored truth (this access never floors a truth)
        plan["floored_truth"] = dict(path=str(spec["truth_h5ad"]), built_by_552=bool(spec.get("floored_build_json")))
        if not spec.get("floored_build_json"):
            msg = (f"{a.surface}: 552 has not built the floored truth (experiments/552-r7-reserve/work/r7b_floor"
                   f"{spec['depth_floor']}/truth_r7b.h5ad + build_r7b.json); the access cannot run until it does")
            if not a.dry_run:
                raise SystemExit(msg)
            S.log("NOTE: " + msg)
    cr = S.load_controlref(spec, need_cs=False)
    targets = cr["targets"].astype(str)
    targets = targets[~np.isin(targets, list(exclude))]
    tag = a.tag + (("-x" + "+".join(exclude)) if exclude else "")
    if spec.get("floored_build_json"):
        # 552's floored surface: every output against its build record (the truth by sha256 only), its target list == the
        # prediction side's kept targets, its control bundle == the prediction side's controls (column sums and count;
        # control cells only)
        fdir = Path(spec["floored_build_json"]).parent
        fb = json.loads(Path(spec["floored_build_json"]).read_text())
        fchk = {rel: S.sha256(fdir / rel) == want for rel, want in fb["outputs"].items()}
        if not all(fchk.values()):
            raise SystemExit(f"552's floored build: sha256 mismatch {fchk}")
        if list(map(str, fb["targets"])) != list(targets):
            raise SystemExit("552's floored target list differs from the prediction side's kept targets")
        xb = S.read_rows(fdir / "bundle/context_R.h5ad", None)
        same_ctrl = (xb.shape[0] == int(cr["n_controls"])
                     and np.array_equal(S.int_colsums(xb).astype(np.float64), np.asarray(cr["ctrl_sums"], dtype=np.float64)))
        del xb
        if not same_ctrl:
            raise SystemExit("552's floored control bundle differs from the prediction side's controls")
        plan["floored_build"] = dict(json=str(Path(spec["floored_build_json"]).relative_to(ROOT)), outputs_sha_ok=fchk,
                                     targets_equal=True, controls_equal_prediction_side=True,
                                     cells_per_target_after=fb.get("floor", {}).get("cells_per_target_after"))
        S.log(f"552's floored build verified: {len(fchk)} outputs by sha256 (the truth hashed only), {len(targets)} targets, "
              f"{int(cr['n_controls'])} controls equal to the prediction side's")
    for arm in arms:
        for s in seeds:
            p, m = stats_path(a.surface, arm, s, a.form == "zs", a.stats_map, tag)
            if p is None:
                raise SystemExit(f"no stats for {arm} seed {s} ({a.form}); emit it with run553.py first")
            if list(np.load(p, allow_pickle=True)["targets"].astype(str)) != list(targets):
                raise SystemExit(f"{p}: targets differ from the surface's")
            plan["stats"][f"{arm}|s{s}"] = dict(path=str(p.relative_to(ROOT)), sha256=S.sha256(p),
                                                metrics=None if m is None else str(m.relative_to(ROOT)))
    if a.halves:
        # the control halves' readiness: their prediction side and their seed-0 stats must exist before the access
        missing = []
        for k in [int(x) for x in a.halves.split(",")]:
            for h in "AB":
                hs = f"{a.surface}c{k}{h}"
                hspec = S.surface_spec(hs)
                if not Path(hspec["controlref"]).exists():
                    missing.append(f"{hs}: controlref")
                for arm in [INCUMBENT] + [m for m in a.half_arms.split(",") if m and m != INCUMBENT]:
                    if stats_path(hs, arm, 0, a.form == "zs", None, tag)[0] is None:
                        missing.append(f"{hs}: {arm} s0")
        half_arms = [INCUMBENT] + [m for m in a.half_arms.split(",") if m and m != INCUMBENT]
        half_pairs = [x for x in a.half_pairs.split(",") if x]
        bad = [x for x in half_pairs if len(x.split(":")) != 2 or any(m not in half_arms for m in x.split(":"))]
        if bad:
            raise SystemExit(f"--half-pairs {bad}: both arms of a half pair must be in --half-arms")
        # 555: an effect of this access that the halves do not carry gets no band_total (report-only); listed here, in the
        # plan, before the truth opens, and again by band553 (which used to drop such effects silently)
        unbanded = ([f"{m} vs cpm" for m in arms if m != INCUMBENT and m not in half_arms]
                    + [f"{x.split(':')[0]} vs {x.split(':')[1]}" for x in a.pairs.split(",") if x and x not in half_pairs])
        plan["halves"] = dict(splits=a.halves, arms=a.half_arms, pairs=a.half_pairs, missing=missing, unbanded=unbanded)
        if missing:
            raise SystemExit(f"control halves not ready ({len(missing)} missing): {missing[:6]} ...")
        if unbanded:
            S.log(f"NOTE: no band_total for {unbanded} (not in --half-arms / --half-pairs): report-only effects")
        S.log(f"control halves ready: splits {a.halves}, {len(a.half_arms.split(',')) + 1} arms x seed 0 each")
    S.write_json(outdir / "plan.json", plan)
    S.log(f"plan: {len(plan['stats'])} stats files ({len(arms)} arms x {len(seeds)} seeds, form {a.form})")
    if a.dry_run:
        S.log(f"DRY RUN OK: {len(plan['stats'])} stats files hashed; {len(plan['stats']) * (1 + 2 * a.splits)} scorings "
              f"planned; the truth was not opened -> {outdir / 'plan.json'}")
        return
    if spec["sealed"]:
        S.log(f"ACCESS ARMED for {a.surface} by {a.owner_go} ({plan['owner_go']['sha256']})")
    base_cmd = [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]]
    for ph in PHASES:
        S.log(f"--- phase {ph}")
        r = subprocess.run([*base_cmd, "--phase", ph])
        if r.returncode != 0:
            raise SystemExit(f"phase {ph} failed (rc {r.returncode}); see the log above")
    S.log(f"done in {time.monotonic() - t_start:.1f} s -> {outdir}")
    if a.halves:
        # the control re-draw band (RESULTS, "the band"): the surface's control halves <surface>c<k><A|B> (prediction side
        # rebuilt and emitted beforehand at seed 0, `prep553.py chalf`), each accessed against the same target cells with
        # its half of the controls, then band553 combines them with this access's band
        for k in [int(x) for x in a.halves.split(",")]:
            for h in "AB":
                cmd = [sys.executable, str(Path(__file__).resolve()), f"{a.surface}c{k}{h}", "--arms", a.half_arms,
                       "--seeds", "0", "--pairs", a.half_pairs]
                if a.exclude_targets:
                    cmd += ["--exclude-targets", a.exclude_targets]
                if a.owner_go:
                    cmd += ["--owner-go", a.owner_go]
                S.log(f"--- control half {a.surface}c{k}{h}")
                if subprocess.run(cmd).returncode != 0:
                    raise SystemExit(f"control half {a.surface}c{k}{h} failed")
        r = subprocess.run([sys.executable, str(HERE / "band553.py"), a.surface, "--splits", a.halves, "--tag", tag])
        if r.returncode != 0:
            raise SystemExit("band553 failed")
        S.log(f"control re-draw band -> {HERE / 'work' / f'band553_{a.surface}{tag}.json'}")


if __name__ == "__main__":
    main()
