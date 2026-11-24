#!/usr/bin/env python3
"""553 — the prediction side of a surface, from its CONTROL cells and public donors only (RESERVE.md section 6).

Steps (each writes only under experiments/553-r7b-harness/work/):
  controlref h1 --ka   lib553.controlref on h1-test's 18,400 controls must equal 323's reference.npz control arrays bit for
                       bit (keep, cs, tie, ref_lin, emit_cs, emit_tie, emit_ref_lin, libraries, b_ctrl, the control row's
                       sums / jackknife / count). Nothing large is written (work/h1ka/controlref_ka.json).
  controlref r7b       the same function on context_R.h5ad -> work/r7b/controlref553.npz (+ .json: sha256 of the npz and of
                       every array, the control-block summary). Control-only statistics (allowed).
  callblock r7b        528's callblock528.py (the October form of the call-propensity key: log10 control CPM + 22 x excess
                       dropout, median-rescaled to H1's dropout scale) on the R7b bundle -> work/r7b/callblock_r7b.npz
  callblock h1 --ka    callblock528.py --ka (read-only: the H1 key from h1-test's controls gives keys528's lcap blocks)
  cd4 r7b              047's extract_cd4.py for the R7b panel (relaxed thresholds, as def_pipeline on a non-A/B/C panel)
                       -> work/r7b/cd4_compact_r7b
  cd4 h1 --ka          the same extraction for 273's H1 panel into work/h1ka/cd4_ka; every array must equal
                       273's cd4_compact_h1 bit for bit (then the copy is deleted)
  audit r7b            the prediction-side audit: bundle shas vs build_r7b.json (truth: sha only), donor <-> target
                       coverage, forbidden-source scan of every input path, operator / key / effects shas, control stats
usage: uv run python experiments/553-r7b-harness/prep553.py <step> <surface> [--ka]
"""
from __future__ import annotations

import json
import shutil
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

E323 = ROOT / "experiments/323-current-family-validation"
CD4_SOURCE = ROOT / "experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad"
CD4_SOURCE_BYTES = 16_786_240_107
CR_KEYS = ("keep", "cs", "tie", "ref_lin", "emit_cs", "emit_tie", "emit_ref_lin", "libraries", "b_ctrl")


def controlref(surface: str, ka: bool):
    spec = S.surface_spec(surface)
    t0 = time.time()
    if surface == "h1":
        assert ka, "h1's control reference is 323's; this step only checks it (--ka)"
        res = {}
        ref_npz = E323 / "work/reference.npz"

        def check(part, keys):
            with np.load(ref_npz, allow_pickle=False) as zz:
                for k in keys:
                    res[k] = bool(np.array_equal(zz[k], part[k]))
            print(f"KA-C partial {dict((k, res[k]) for k in keys)}; peak RSS {S.rss_gb():.2f} GB", flush=True)
        cr = S.controlref(Path(spec["control_h5ad"]), check=check)
        with np.load(ref_npz, allow_pickle=False) as z:
            for k in CR_KEYS:
                if k not in res:
                    res[k] = bool(np.array_equal(z[k], cr[k]))
            res["genes"] = bool(np.array_equal(z["genes"], cr["genes"]))
            res["ctrl_sums"] = bool(np.array_equal(z["sums"][-1], cr["ctrl_sums"]))
            res["ctrl_jk"] = bool(z["jk"][-1] == cr["ctrl_jk"])
            res["n_controls"] = bool(int(z["counts"][-1]) == int(cr["n_controls"]))
        classes = S.library_classes(cr["libraries"])
        out = dict(step="controlref", surface="h1", equal=res, pass_=all(res.values()), classes=classes,
                   seconds=round(time.time() - t0, 1), peak_rss_gb=S.rss_gb())
        S.write_json(HERE / "work/h1ka/controlref_ka.json", out)
        print(f"KA-C {'PASS' if out['pass_'] else 'FAIL'}: {res}", flush=True)
        sys.exit(0 if out["pass_"] else 1)
    if surface == "r7b" or surface.startswith("r7bd"):
        S.check_r7b_outputs(check_truth_sha=False)
        targets = S.surface_targets(spec)
    else:                       # h1r7: the target list is the file's own labels
        lab, _ = S.read_labels_genes(Path(spec["control_h5ad"]))
        targets = np.array(sorted(set(lab) - {S.CTRL}))
    cr = S.controlref(Path(spec["control_h5ad"]))
    out_npz = Path(spec["controlref"])
    out_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez(out_npz, targets=targets, **cr)
    lib = cr["libraries"]
    rec = dict(surface=surface, source=str(Path(spec["control_h5ad"]).relative_to(ROOT)),
               source_sha256=S.sha256(Path(spec["control_h5ad"])), statistics_sha256=S.sha256(out_npz),
               arrays_sha256={k: S.arr_sha(v) for k, v in cr.items() if isinstance(v, np.ndarray)},
               n_controls=int(cr["n_controls"]), n_genes=int(len(cr["genes"])), n_gated=int(cr["keep"].sum()),
               n_targets=int(len(targets)), library_quantiles={q: float(np.quantile(lib, q)) for q in (0, .1, .25, .5, .75, .9, 1)},
               classes=S.library_classes(lib), ctrl_jk=float(cr["ctrl_jk"]),
               zero_genes=int((cr["ctrl_sums"] == 0).sum()), b_ctrl_positive=int((cr["b_ctrl"] > 0).sum()),
               note="control cells only (RESERVE.md section 6): no perturbed cell was read",
               seconds=round(time.time() - t0, 1), peak_rss_gb=S.rss_gb())
    S.write_json(Path(spec["controlref_json"]), rec)
    print(json.dumps({k: v for k, v in rec.items() if k != "arrays_sha256"}), flush=True)


def callblock(surface: str, ka: bool):
    script = ROOT / "experiments/528-callprop-block/callblock528.py"
    if surface == "h1":
        assert ka
        r = subprocess.run([sys.executable, str(script), "--ka"], capture_output=True, text=True)
        print(r.stdout, r.stderr, flush=True)
        S.write_json(HERE / "work/h1ka/callblock_ka.json", dict(rc=r.returncode, stdout=r.stdout[-2000:]))
        sys.exit(r.returncode)
    if surface == "h1cb" or surface.startswith(("h1r7", "r7bd")):
        # the October key formula on h1-test's controls, on h1-test's own axis (callblock528.py itself needs a bundle on the
        # official axis): the same two functions it calls, feat528.features and key_from_features, artifact in its format
        import hashlib
        sys.path.insert(0, str(ROOT / "experiments/528-callprop-block"))
        from callblock528 import MEDIAN_X_H1, key_from_features, median_x_h1  # noqa: PLC0415
        from feat528 import features  # noqa: PLC0415
        med_h1, _ = median_x_h1()
        F = features(ROOT / "data/harness/h1-test.h5ad" if surface == "h1cb" else S.surface_spec(surface)["control_h5ad"])
        if surface.startswith("r7bd"):
            import pandas as pd
            off = list(pd.read_csv(S.R7B / "bundle/gene_names.csv")["gene_name"].astype(str))
            if list(map(str, F["genes"])) != off:
                raise SystemExit("the floored control file's gene order differs from gene_names.csv")
        k, med = key_from_features(F, med_h1)
        out = HERE / "work" / surface / f"callblock_{surface}.npz"
        out.parent.mkdir(parents=True, exist_ok=True)
        np.savez(out, genes=np.asarray(F["genes"]), key_H1=k, meta=np.array(json.dumps(dict(beta=22.0, median_x_h1=med_h1, median_x=med))))
        sha = hashlib.sha256(out.read_bytes()).hexdigest()
        Path(str(out) + ".json").write_text(json.dumps(dict(sha256=sha, path=str(out), median_x=med, scale=med_h1 / med,
                                                            n_universe=int(F["keep"].sum()), n_ctrl=int(F["n_ctrl"])), indent=1))
        print(f"{surface} key: median x {med} (H1 pinned {MEDIAN_X_H1}); sha256 {sha}")
        return
    spec = S.surface_spec(surface)
    S.check_r7b_outputs(check_truth_sha=False)
    out = HERE / "work" / surface / f"callblock_{surface}.npz"
    out.parent.mkdir(parents=True, exist_ok=True)
    r = subprocess.run([sys.executable, str(script), "--bundle", str(spec["bundle"]), "--contexts", "R", "--out", str(out)],
                       capture_output=True, text=True)
    print(r.stdout, r.stderr, flush=True)
    if r.returncode != 0:
        sys.exit(r.returncode)
    # the swap this key makes against the CPM block, per target head is not known here: the no-head approximation is in
    # the artifact's json (callblock528.py); the per-target swap is read from the emission diagnostics
    rec = json.loads(Path(str(out) + ".json").read_text())
    print(json.dumps(rec), flush=True)


def cd4(surface: str, ka: bool):
    if not CD4_SOURCE.is_file() or CD4_SOURCE.stat().st_size != CD4_SOURCE_BYTES:
        raise SystemExit(f"CD4 source missing or wrong size: {CD4_SOURCE}")
    script = ROOT / "experiments/047-cd4-direct/extract_cd4.py"
    if surface == "h1":
        assert ka
        ref = ROOT / "experiments/273-flex-reference-instrument/work/cd4/cd4_compact_h1"
        man = json.loads((ref / "manifest.json").read_text())
        out = HERE / "work/h1ka/cd4_ka"
        cmd = [sys.executable, str(script), "--input", str(CD4_SOURCE), "--output", str(out),
               "--panel-csv", man["panel_csv"], "--genes-csv", str(ROOT / "data/challenge/gene_names.csv"),
               "--overwrite", "--min-targets", "0", "--min-union-targets", "0"]
        subprocess.run(cmd, check=True)
        res = {}
        for cond in ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr"):
            for f in sorted((ref / cond).iterdir()):
                if f.suffix == ".npy":
                    res[f"{cond}/{f.name}"] = bool(np.array_equal(np.load(f), np.load(out / cond / f.name), equal_nan=True))
                elif f.suffix in (".txt", ".csv"):
                    res[f"{cond}/{f.name}"] = f.read_bytes() == (out / cond / f.name).read_bytes()
        ok = all(res.values())
        S.write_json(HERE / "work/h1ka/cd4_ka.json", dict(pass_=ok, equal=res))
        print(f"KA-CD4 {'PASS' if ok else 'FAIL'}: {sum(res.values())}/{len(res)} files equal", flush=True)
        if ok:
            shutil.rmtree(out)
        sys.exit(0 if ok else 1)
    spec = S.surface_spec(surface)
    S.check_r7b_outputs(check_truth_sha=False)
    out = Path(spec["cd4"])
    cmd = [sys.executable, str(script), "--input", str(CD4_SOURCE), "--output", str(out),
           "--panel-csv", str(Path(spec["bundle"]) / "pert_counts.csv"),
           "--genes-csv", str(Path(spec["bundle"]) / "gene_names.csv"),
           "--overwrite", "--min-targets", "0", "--min-union-targets", "0"]
    subprocess.run(cmd, check=True)
    man = json.loads((out / "manifest.json").read_text())
    print(json.dumps({k: man[k] for k in ("panel_union_coverage", "panel_union_missing", "mapped_official_genes")}), flush=True)


def audit(surface: str, ka: bool):
    """The R7b prediction side: what it read, what it covers, and that none of it is R7b's own line."""
    assert surface == "r7b"
    spec = S.surface_spec(surface)
    rec = {"bundle_vs_build_r7b": S.check_r7b_outputs(check_truth_sha=True)}
    build = S.r7b_build_record()
    targets = list(S.surface_targets(spec))
    cov = {}
    for n in S.DONORS:
        perts = {ln.strip() for ln in (ROOT / f"data/donors/{n}/perts.txt").read_text().splitlines() if ln.strip()}
        cov[n] = sorted(set(targets) & perts)
    ess = {ln.strip() for ln in (ROOT / "data/donors/k562_essential/perts.txt").read_text().splitlines() if ln.strip()}
    per_t = {t: [n for n in S.DONORS if t in cov[n]] for t in targets}
    # every file the prediction side reads (paths), scanned for R7b's own study / line
    inputs = [*(ROOT / f"data/donors/{n}" for n in S.DONORS), ROOT / "data/external/perturbmatch", CD4_SOURCE,
              ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz",
              ROOT / "experiments/528-callprop-block/callblock528.py", ROOT / "experiments/528-callprop-block/work/feat528_H1.npz",
              Path(spec["cd4"]), Path(spec["operator"]), Path(spec["effects"]), Path(spec["controlref"])]
    S.assert_clean_sources(inputs)
    cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    banks = list(map(str, cp["bank_names"]))
    rec.update(
        n_targets=len(targets), cells_per_target=build["cells_per_target"],
        cells_per_target_quantiles={q: float(np.quantile(list(build["cells_per_target"].values()), q)) for q in (0, .25, .5, .75, 1)},
        n_controls=build["n_controls"], genes=build["genes"],
        donor_coverage={n: len(v) for n, v in cov.items()},
        donors_per_target_hist={k: sum(1 for v in per_t.values() if len(v) == k) for k in (0, 1, 2, 3)},
        essential_targets=sorted(set(targets) & ess),
        clean4_banks_used=["k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2"], consensus_prior_banks=banks,
        forbidden_source_scan=dict(patterns=list(S.FORBIDDEN_SOURCES), inputs=[str(p.relative_to(ROOT)) for p in inputs],
                                   hits=0),
        shas={p.name: S.sha256(p) for p in (Path(spec["operator"]), Path(spec["effects"]), Path(spec["controlref"]),
                                             HERE / "work/r7b/callblock_r7b.npz") if p.exists()})
    # what the emitter reads on R7b's own axis: the CLEAN4 block prior's coverage, the hook's signed genes (tau 0.05),
    # the call-propensity key, the operator, the bulk key, the carriers
    cr = S.load_controlref(spec, need_cs=False)
    keep = np.asarray(cr["keep"], bool)
    genes = cr["genes"].astype(str)
    P = cp["P_banks"].astype(float)[[banks.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
    pos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
    col = np.array([pos.get(g, -1) for g in genes])
    Pf = np.full(len(genes), np.nan)
    Pf[col >= 0] = P[col[col >= 0]]
    hook = np.isfinite(Pf[keep]) & (np.abs(Pf[keep] - 0.5) >= 0.05)
    opj = json.loads(Path(spec["operator"]).with_suffix(".json").read_text())
    comp = json.loads(Path(spec["compose"]).read_text())
    cb = json.loads((HERE / "work/r7b/callblock_r7b.npz.json").read_text())
    with np.load(HERE / "work/r7b/bulk553.npz", allow_pickle=False) as z:
        bulk = dict(n_carried=int(np.sum(z["carried"])), median_key=float(z["median"]), w_min=float(z["w"].min()),
                    w_max=float(z["w"].max()), energy_c=float(z["c"]), operator_sha256=str(z["operator_sha256"]))
    own_code = sorted(str(f.name) for f in HERE.glob("*.py"))
    # the patterns are assembled so this scanner's own source does not match itself
    pats = ("sunshine" + "hein2023" + ".h5ad", "sunshine" + "_" + "calu3", "scperturb/" + "sunshine")
    hits = [f for f in own_code if any(w in (HERE / f).read_text().lower() for w in pats)]
    rec.update(
        gated_axis=int(keep.sum()), clean4_prior_coverage_gated=float(np.isfinite(Pf[keep]).mean()),
        clean4_prior_down_frac_gated=float((Pf[keep][np.isfinite(Pf[keep])] > 0.5).mean()),
        hook_tau005_signed_genes=int(hook.sum()), hook_down=int((hook & (Pf[keep] > 0.5)).sum()),
        callblock=cb["contexts"]["R"], operator=dict(path=str(Path(spec["operator"]).relative_to(ROOT)), sha256=opj["sha256"],
                                                    axis_sha256=opj["axis_sha256"], n_cols=opj["n_cols"], basis=opj["basis"],
                                                    dropped_in_panel=opj.get("dropped_in_panel")),
        carriers=dict(final_sha256=comp["arms"]["cap15"]["final_sha256"], effects_sha256=comp.get("effects_sha256"),
                      cd4=comp["arms"]["cap15"]["cd4"].get("deflation"), whitening=comp["arms"]["cap15"]["whitening"]),
        bulk_key=bulk, own_code_mentions_of_the_source_file=hits)
    S.write_json(HERE / "work/r7b/audit553.json", rec)
    print(json.dumps({k: v for k, v in rec.items() if k not in ("cells_per_target",)}, indent=1), flush=True)


def h1r7(surface: str, ka: bool):
    """`prep553.py h1r7 build`: the R7b-shaped dress rehearsal surface from h1-test (NOT a reserve): R7b's control count
    and its cells-per-target distribution (build_r7b.json METADATA: counts only) mapped onto h1-test's 85 targets
    (quantiles of R7b's counts in a seeded random order, capped at 400); cells drawn without replacement, seed 55307."""
    import anndata as ad
    import pandas as pd
    rec = S.r7b_build_record()
    # --controls N / --cells-scale f: the R7b shape to emulate (default: R7b as built, 3,536 controls, its counts); the
    # floored candidates keep 2,095 (r7bd8000) / 2,378 (r7bd5000) controls and, as an ESTIMATE (target cells' totals were
    # not read; the controls' kept fraction is assumed), ~0.59 / ~0.67 of each target's cells. The queued rehearsal
    # emulates r7bd8000 (--controls 2095 --cells-scale 0.59).
    name = "h1r7" if surface == "build" else surface                    # `h1r7 build` = the first draw (seed 55307)
    assert name.startswith("h1r7"), name
    seed = int(sys.argv[sys.argv.index("--seed") + 1]) if "--seed" in sys.argv else 55307
    n_ctrl = int(sys.argv[sys.argv.index("--controls") + 1]) if "--controls" in sys.argv else int(rec["n_controls"])
    cscale = float(sys.argv[sys.argv.index("--cells-scale") + 1]) if "--cells-scale" in sys.argv else 1.0
    counts = np.sort(np.asarray(list(rec["cells_per_target"].values()), dtype=float)) * cscale
    src = ROOT / "data/harness/h1-test.h5ad"
    lab, genes = S.read_labels_genes(src)
    targets = np.array(sorted(set(lab) - {S.CTRL}))
    rng = np.random.default_rng(seed)
    crow = np.flatnonzero(lab == S.CTRL)
    pick = [np.sort(rng.choice(crow, n_ctrl, replace=False))]
    qs = np.quantile(counts, (np.arange(len(targets)) + 0.5) / len(targets), method="nearest").astype(int)
    assigned = np.minimum(qs[rng.permutation(len(targets))], 400)
    per = {}
    for t, n in zip(targets, assigned):
        trow = np.flatnonzero(lab == t)
        pick.append(np.sort(rng.choice(trow, int(n), replace=False)))
        per[str(t)] = int(n)
    rows = np.sort(np.concatenate(pick))
    X = S.read_rows(src, rows)
    out = HERE / f"work/{name}/{name}.h5ad"
    out.parent.mkdir(parents=True, exist_ok=True)
    a = ad.AnnData(X=X, obs=pd.DataFrame({"target_gene": pd.Categorical(lab[rows])}, index=[f"h1r7_{i}" for i in rows]),
                   var=pd.DataFrame(index=pd.Index(genes)))
    a.write_h5ad(out)
    info = dict(source=str(src.relative_to(ROOT)), source_sha256=S.sha256(src), seed=seed, n_controls=n_ctrl,
                cells_scale=cscale,
                r7b_cells_per_target_quantiles={q: float(np.quantile(counts, q)) for q in (0, .25, .5, .75, 1)},
                cells_per_target=per, n_rows=int(len(rows)), sha256=S.sha256(out),
                note="h1-test (not a reserve) cut to R7b's shape; R7b metadata only (counts), no R7b expression")
    S.write_json(out.with_suffix(".json"), info)
    print(json.dumps({k: v for k, v in info.items() if k != "cells_per_target"}), flush=True)


def h1t7(surface: str, ka: bool):
    """`prep553.py h1t7 build`: h1r7's target cells (the same h1-test rows, read from h1r7.h5ad's obs names) plus ALL of
    h1-test's 18,400 controls, in h1-test's row order -- the rehearsal's truth shape with h1's control set, so the access
    can run h1's own emitted stats against it (an ablation; not a reserve, no R7b data)."""
    import anndata as ad
    import pandas as pd
    src = ROOT / "data/harness/h1-test.h5ad"
    r7name = "h1r7" if surface == "build" else surface                  # the rehearsal draw whose target cells to take
    assert r7name.startswith("h1r7"), r7name
    name = "h1t7" + r7name[len("h1r7"):]                                # h1r7 -> h1t7, h1r7s2 -> h1t7s2
    r7 = HERE / f"work/{r7name}/{r7name}.h5ad"
    lab, genes = S.read_labels_genes(src)
    obs7 = ad.read_h5ad(r7, backed="r").obs
    rows7 = np.array([int(i.split("_", 1)[1]) for i in obs7.index])
    assert np.array_equal(lab[rows7], obs7["target_gene"].astype(str).to_numpy()), "h1r7 rows do not map onto h1-test"
    trow = rows7[lab[rows7] != S.CTRL]
    rows = np.sort(np.concatenate([np.flatnonzero(lab == S.CTRL), trow]))
    X = S.read_rows(src, rows)
    out = HERE / f"work/{name}/{name}.h5ad"
    out.parent.mkdir(parents=True, exist_ok=True)
    a = ad.AnnData(X=X, obs=pd.DataFrame({"target_gene": pd.Categorical(lab[rows])}, index=[f"h1t7_{i}" for i in rows]),
                   var=pd.DataFrame(index=pd.Index(genes)))
    a.write_h5ad(out)
    info = dict(source=str(src.relative_to(ROOT)), source_sha256=S.sha256(src), rehearsal=r7name, h1r7_sha256=S.sha256(r7),
                n_controls=int((lab[rows] == S.CTRL).sum()), n_target_cells=int(trow.size), n_rows=int(rows.size),
                sha256=S.sha256(out), note="h1r7's target cells + all h1-test controls (ablation); no R7b data")
    S.write_json(out.with_suffix(".json"), info)
    print(json.dumps(info), flush=True)


def chalf(base: str, ka: bool):
    """`prep553.py chalf <base> <k>`: split k of the base surface's CONTROLS into two disjoint halves A / B (seeded:
    SeedSequence(5532).generate_state(16)[k]; floor(n/2) each, file order), one file per half: <base>c<k>A / B. On a
    rehearsal base (h1r7*) the half file carries ALL the base's target cells too (it is the half surface's truth); on an
    R7b base (r7bd<floor>) it holds the half's controls only (control cells: RESERVE section 6) and a json with their obs
    names, which the armed access uses to keep the same controls of the truth."""
    import anndata as ad
    k = int(sys.argv[3])
    spec = S.surface_spec(base)
    src = Path(spec["control_h5ad"])
    if base.startswith("r7bd"):
        S.check_r7b_outputs(check_truth_sha=False)
    lab, genes = S.read_labels_genes(src)
    crow = np.flatnonzero(lab == S.CTRL)
    trow = np.flatnonzero(lab != S.CTRL)
    seed = int(np.random.SeedSequence(5532).generate_state(16)[k])
    perm = np.random.default_rng(seed).permutation(crow)
    h = perm.size // 2
    a = ad.read_h5ad(src)
    for half, part in (("A", perm[:h]), ("B", perm[h:2 * h])):
        name = f"{base}c{k}{half}"
        rows = np.sort(np.concatenate([np.sort(part), trow]))
        out = HERE / "work" / name / f"{name}.h5ad"
        out.parent.mkdir(parents=True, exist_ok=True)
        a[rows].copy().write_h5ad(out)
        rec = dict(base=base, split=k, half=half, seed=seed, n_controls=int(part.size), n_target_cells=int(trow.size),
                   source=str(src.relative_to(ROOT)), source_sha256=S.sha256(src), sha256=S.sha256(out),
                   control_obs_names=[str(x) for x in a.obs_names[np.sort(part)]],
                   note="a control half of the base surface (controls only were selected; no truth was read)")
        S.write_json(out.with_suffix(".json"), rec)
        print(json.dumps({k_: v for k_, v in rec.items() if k_ != "control_obs_names"}), flush=True)


def subset(surface: str, ka: bool):
    """`prep553.py subset r7bd<floor>`: the R7b control cells whose official-axis library >= floor, as a control file for
    the prediction side (control cells only; RESERVE.md section 6)."""
    import anndata as ad
    import pandas as pd
    spec = S.surface_spec(surface)
    S.check_r7b_outputs(check_truth_sha=False)
    src = S.surface_spec("r7b")["control_h5ad"]
    lab, genes = S.read_labels_genes(src)
    x = S.read_rows(src, None)
    lib = np.asarray(x.sum(1)).ravel()
    keep = lib >= float(spec["depth_floor"])
    out = Path(spec["control_h5ad"])
    out.parent.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(src)
    a[keep].copy().write_h5ad(out)
    q = np.quantile(lib[keep], [0, .25, .5, .75, 1])
    rec = dict(surface=surface, floor=spec["depth_floor"], kept=int(keep.sum()), of=int(keep.size),
               library_quantiles=[float(v) for v in q], q75_over_q25=float(q[3] / q[1]), sha256=S.sha256(out),
               note="control cells only")
    S.write_json(out.with_suffix(".json"), rec)
    print(json.dumps(rec), flush=True)


def main():
    step, surface = sys.argv[1], sys.argv[2]
    ka = "--ka" in sys.argv
    if step == "h1r7":
        return h1r7(surface, ka)
    if step == "h1t7":
        return h1t7(surface, ka)
    if step == "chalf":
        return chalf(surface, ka)
    if step == "subset":
        return subset(surface, ka)
    {"controlref": controlref, "callblock": callblock, "cd4": cd4, "audit": audit}[step](surface, ka)


if __name__ == "__main__":
    main()
