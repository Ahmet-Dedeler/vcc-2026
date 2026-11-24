#!/usr/bin/env python3
"""376 — reusable same-line calibration artifact builder.

Generalizes 372's hand-adapted Jurkat script (which mirrors 346 fit346.py verbatim)
to any recipient bank: ridge map from donor-consensus directions to the recipient
line's own public CRISPRi responses, fitted on non-panel targets, evaluated on panel
targets with a full donor consensus.

Estimator primitives are imported verbatim from 341/334; nothing is copied/retuned.
The recipient's panel-target rows are never read during fitting; the poison replay
proves it by refitting with those rows replaced by noise and requiring byte identity.

Output matches the 373 `blend_rows` builder mode: `matrix` is the UNBLENDED ridge
prediction (372's `prediction_raw`); the builder applies the fixed-norm blend at beta
itself. Key set and save recipe match
experiments/375-ctx-a-arm/work/ctx_a_blend_rows_b035.npz exactly
(uncompressed np.savez, keys targets/genes/matrix/beta/mode/context/norm_rule).

Writes only the artifact npz + one JSON sidecar. Never touches scripts/src/builders.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, os, resource, sys, time
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for _p in (str(ROOT / "experiments/334-text-conditioned-transfer"),
           str(ROOT / "experiments/341-known-line-calibration")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import core334 as P
import core341 as C

SEED = 3460910  # 346's derangement seed, kept verbatim (null control + poison noise)
MODE = "blend_rows"
NORM_RULE = "fixed_norm_on_covered"


def names(p: Path) -> list[str]:
    return p.read_text().splitlines()


def sha_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def bank(path: Path, targets, genes):
    """346 fit346.py::bank shape: exact-name row/column selection + fixture."""
    perts, gene_names = names(path / "perts.txt"), names(path / "genes.txt")
    raw = np.load(path / "delta.npy", mmap_mode="r")
    if raw.shape != (len(perts), len(gene_names)):
        raise ValueError(f"{path}: delta shape {raw.shape} vs name lists")
    pi = {p: i for i, p in enumerate(perts)}
    gi = {g: i for i, g in enumerate(gene_names)}
    rows = [pi[t] for t in targets]
    cols = [gi[g] for g in genes]
    out = np.array(raw[np.ix_(rows, cols)], dtype=np.float64)
    basal = np.array(np.load(path / "b_ctrl.npy")[cols], dtype=np.float64)
    if not (np.isfinite(out).all() and np.isfinite(basal).all()):
        raise ValueError(f"{path}: non-finite values in the selected block")
    np.testing.assert_array_equal(
        out[:3, :3], [[raw[pi[t], gi[g]] for g in genes[:3]] for t in targets[:3]])
    return out, basal


def build_axes(recipient: Path, donors: list[Path], panel_csv: Path, axis_csv: Path):
    panel = [r["target_gene"] for r in csv.DictReader(panel_csv.open())]  # file order
    canonical = [r["gene_name"] for r in csv.DictReader(axis_csv.open())]
    panel_set = set(panel)
    donor_perts = [set(names(d / "perts.txt")) for d in donors]
    donor_genes = [set(names(d / "genes.txt")) for d in donors]
    r_perts = names(recipient / "perts.txt")
    r_genes = names(recipient / "genes.txt")
    r_pert_set = set(r_perts)
    # Training excludes ALL panel symbols (372 asserted the overlap was 0 for
    # Jurkat; a general recipient e.g. K562 carries most of the panel — exclusion,
    # not assertion, is the safe general rule).
    recip_panel_overlap = sorted(r_pert_set & panel_set)
    training = sorted((set.intersection(*donor_perts) & r_pert_set) - panel_set)
    common = set.intersection(*donor_genes) & set(r_genes) & set(canonical)
    # Optional measured mask (Pan H1 precedent, 346): intersect if present.
    mask_note = "absent"
    mp = recipient / "measured_mask.npy"
    if mp.exists():
        mask = np.load(mp, mmap_mode="r")
        if mask.shape[0] == len(canonical) and (names(recipient / "genes.txt") == canonical):
            row0 = np.asarray(mask[0], dtype=bool)
            measured = {g for g, k in zip(canonical, row0) if k}
            common = common & measured
            mask_note = f"applied ({len(measured)} measured)"
        else:
            mask_note = "present but not on canonical axis — IGNORED"
    excluded_genes = set(training) | panel_set
    covered = sorted(common - excluded_genes)
    test = [t for t in panel if all(t in s for s in donor_perts)]
    missing = [t for t in panel if t not in set(test)]
    panel_on_axis = sorted(common & panel_set)
    train_on_axis = sorted(common & set(training))
    return dict(panel=panel, canonical=canonical, training=training, covered=covered,
                test=test, missing=missing, panel_on_axis=panel_on_axis,
                train_on_axis=train_on_axis, common_n=len(common),
                r_perts=r_perts, r_genes=r_genes,
                recip_panel_overlap=recip_panel_overlap, mask_note=mask_note)


def assemble(axes, recipient: Path, donors: list[Path], poison: bool):
    order = axes["training"] + axes["test"]
    genes = axes["covered"]
    raws, basals = [], []
    for d in donors:
        a, b = bank(d, order, genes)
        raws.append(a)
        basals.append(b)
    r_set = set(axes["r_perts"])
    have_train = [t for t in axes["training"] if t in r_set]
    assert len(have_train) == len(axes["training"]), "training ⊆ recipient required"
    rec_train, rec_basal = bank(recipient, axes["training"], genes)
    rec = np.zeros((len(order), len(genes)), dtype=np.float64)
    rec[:len(axes["training"])] = rec_train
    del rec_train
    if poison:
        rng = np.random.default_rng(SEED + 77)
        rec[len(axes["training"]):] = 1e6 * rng.standard_normal(
            (len(axes["test"]), len(genes)))
    raws.append(rec)
    basals.append(rec_basal)
    raw = np.asarray(raws)
    del raws
    direction, rms = P.rms_direction(raw)
    del raw  # 372 held this through the poison pass (12.22 GiB); free it here.
    if poison:
        rms[len(donors), len(axes["training"]):] = 1e6
    return direction, np.asarray(basals), rms


def fit_all(axes, direction, basal, rms, n_donors: int):
    n_train = len(axes["training"])
    train = np.arange(n_train)
    test = np.arange(n_train, direction.shape[1])
    outer = n_donors
    donors = list(range(n_donors))
    baseline = direction[donors].mean(0)
    inner = P.source_blocks(direction, basal, rms, outer, train)
    cross_gain = P.diagonal(inner)
    x = np.concatenate([f[..., 0] for _, f, _ in inner])
    cross_y = np.concatenate([r for _, _, r in inner])
    same_y = np.tile(direction[outer, train], (len(donors), 1))
    permutation = P.derangement(n_train, SEED)
    shuffled_y = np.tile(direction[outer, train[permutation]], (len(donors), 1))
    models = {"crossline_full": C.fit(x, cross_y, gain=cross_gain),
              "same_line_full": C.fit(x, same_y),
              "labels_shuffled_full": C.fit(x, shuffled_y)}
    del cross_y, same_y, shuffled_y
    query = baseline[test]
    predictions = {"source_mean": query.copy(),
                   "crossline_diag": query * cross_gain,
                   "same_line_diag": query * models["same_line_full"]["gain"]}
    predictions.update({a: C.predict(m, query) for a, m in models.items()})
    subset = np.mean([C.predict(models["same_line_full"],
                                direction[[d for d in donors if d != omit]][:, test].mean(0))
                      for omit in donors], axis=0)
    subset_error = float(np.abs(subset - predictions["same_line_full"]).max())
    sc = models["same_line_full"]["scale"]
    beta = models["same_line_full"]["beta"]
    residual = float(np.abs(
        (x / sc).T @ ((x / sc) @ beta) + len(x) * beta
        - (x / sc).T @ (np.tile(direction[outer, train], (len(donors), 1))
                        - x * models["same_line_full"]["gain"])).max())
    return models, predictions, dict(cross_gain=cross_gain, permutation=permutation,
                                     subset_linear_max_difference=subset_error,
                                     normal_equation_max_residual=residual,
                                     n_rows=int(len(x)), n_features=int(x.shape[1]))


def heldout_two_split(axes, direction, basal, rms, n_donors: int, seed: int):
    """372's within-bank held-out sign protocol (diagnostic, not a score)."""
    n = len(axes["training"])
    idx = np.arange(n)
    halves = {"evenodd": (idx[0::2], idx[1::2])}
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    halves["random"] = (np.sort(perm[:n // 2]), np.sort(perm[n // 2:]))
    out = {}
    for name, (ia, ib) in halves.items():
        donors = list(range(n_donors))
        outer = n_donors
        inner = P.source_blocks(direction, basal, rms, outer, ia)
        cross_gain = P.diagonal(inner)
        x = np.concatenate([f[..., 0] for _, f, _ in inner])
        same_y = np.tile(direction[outer, ia], (len(donors), 1))
        permA = P.derangement(len(ia), SEED)
        shuf_y = np.tile(direction[outer, ia[permA]], (len(donors), 1))
        m_same = C.fit(x, same_y)
        m_shuf = C.fit(x, shuf_y)
        baseB = direction[donors][:, ib].mean(0)
        pred_same = C.predict(m_same, baseB)
        pred_diag = baseB * cross_gain
        pred_shuf = C.predict(m_shuf, baseB)
        truth = direction[outer, ib]
        s_same = C.sign_endpoints(pred_same, truth)
        s_diag = C.sign_endpoints(pred_diag, truth)
        s_shuf = C.sign_endpoints(pred_shuf, truth)
        d_same = s_same["weighted_sign"] - s_diag["weighted_sign"]
        d_shuf = s_shuf["weighted_sign"] - s_diag["weighted_sign"]
        rng2 = np.random.default_rng(seed + 1)
        bs = []
        for s in range(0, 10000, 100):
            nn = min(100, 10000 - s)
            ti = rng2.integers(len(ib), size=(nn, len(ib)))
            bs.extend(d_same[ti].mean(1))
        out[name] = dict(n_fit=int(len(ia)), n_eval=int(len(ib)),
                         mean_same=float(s_same["weighted_sign"].mean()),
                         mean_diag=float(s_diag["weighted_sign"].mean()),
                         mean_shuf=float(s_shuf["weighted_sign"].mean()),
                         gain_same_minus_diag=float(d_same.mean()),
                         ci95=np.quantile(bs, [.025, .975]).tolist(),
                         frac_improved=float((d_same > 0).mean()),
                         gain_shuf_minus_diag=float(d_shuf.mean()),
                         top64_gain=float((s_same["top64_weighted_sign"]
                                           - s_diag["top64_weighted_sign"]).mean()))
    return out


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Build a same-line calibration artifact.")
    ap.add_argument("--recipient", required=True, help="recipient bank dir")
    ap.add_argument("--donors", required=True, nargs="+", help="donor bank dirs")
    ap.add_argument("--panel", required=True, help="pert_counts.csv")
    ap.add_argument("--axis", required=True, help="gene_names.csv")
    ap.add_argument("--beta", required=True, type=float)
    ap.add_argument("--context", required=True)
    ap.add_argument("--out", required=True, help="output .npz path")
    ap.add_argument("--heldout-seed", type=int, default=None,
                    help="also run the two-split held-out sign check with this seed")
    return ap.parse_args(argv)


def main(argv=None) -> None:
    args = parse_args(argv)
    started, clock = datetime.now(timezone.utc).isoformat(), time.monotonic()
    if not (0.0 < args.beta <= 1.0):
        raise ValueError(f"beta must be in (0,1]; got {args.beta!r}")
    recip = Path(args.recipient).resolve()
    donors = [Path(d).resolve() for d in args.donors]
    dropped = None
    if recip in donors:
        donors = [d for d in donors if d != recip]
        dropped = str(recip)
        print(f"recipient is itself a donor — dropped; {len(donors)} donors remain",
              flush=True)
    if len(donors) < 2:
        raise ValueError(f"need >=2 donors after drop, have {len(donors)}")
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    def rel(p: Path) -> str:
        try:
            return str(p.relative_to(ROOT))
        except ValueError:
            return str(p)

    axes = build_axes(recip, donors, Path(args.panel), Path(args.axis))
    print(json.dumps({k: (len(v) if isinstance(v, list) else v)
                      for k, v in axes.items()
                      if k not in ("canonical", "r_perts", "r_genes")},
                     indent=2, sort_keys=True), flush=True)
    n_panel_overlap = len(axes["recip_panel_overlap"])
    print(f"recipient∩panel targets (excluded from training): {n_panel_overlap}",
          flush=True)
    direction, basal, rms = assemble(axes, recip, donors, poison=False)
    print(f"assemble done {time.monotonic() - clock:.1f}s dir{str(direction.shape)}",
          flush=True)
    models, predictions, receipt = fit_all(axes, direction, basal, rms, len(donors))
    print(f"fit done subset={receipt['subset_linear_max_difference']:.3e} "
          f"resid={receipt['normal_equation_max_residual']:.3e} "
          f"rows={receipt['n_rows']} feat={receipt['n_features']}", flush=True)
    # Poison replay: recipient query rows -> noise; everything fitted must be
    # byte-identical (panel rows never read).
    pd, pb, pr = assemble(axes, recip, donors, poison=True)
    np.testing.assert_array_equal(direction[:len(donors)], pd[:len(donors)])
    np.testing.assert_array_equal(rms[:len(donors)], pr[:len(donors)])
    np.testing.assert_array_equal(basal, pb)
    pm, pp, prec = fit_all(axes, pd, pb, pr, len(donors))
    for arm in models:
        for key in models[arm]:
            np.testing.assert_array_equal(models[arm][key], pm[arm][key])
    for arm in predictions:
        np.testing.assert_array_equal(predictions[arm], pp[arm])
    np.testing.assert_array_equal(receipt["cross_gain"], prec["cross_gain"])
    print("poison PASS byte-identical", flush=True)
    del pd, pb, pr, pm, pp

    ho = None
    if args.heldout_seed is not None:
        ho = heldout_two_split(axes, direction, basal, rms, len(donors),
                               args.heldout_seed)
        print(json.dumps(ho, indent=2, sort_keys=True), flush=True)

    mat = np.asarray(predictions["same_line_full"], dtype=np.float64)
    assert np.isfinite(mat).all()
    # Panel-gene exclusion gate: the covered axis must carry zero panel symbols
    # (the 373 builder refuses panel-coordinate support).
    assert len(set(axes["covered"]) & set(axes["panel"])) == 0
    np.savez(str(out),
             targets=np.array(axes["test"]),
             genes=np.array(axes["covered"]),
             matrix=mat,
             beta=np.array(float(args.beta)),
             mode=np.array(MODE),
             context=np.array(str(args.context)),
             norm_rule=np.array(NORM_RULE))
    bank_shas = {}
    for d in donors + [recip]:
        bank_shas[rel(d)] = dict(delta=sha_file(d / "delta.npy"),
                                 perts=sha_file(d / "perts.txt"),
                                 genes=sha_file(d / "genes.txt"),
                                 b_ctrl=sha_file(d / "b_ctrl.npy"))
    sidecar = dict(context=str(args.context), mode=MODE, beta=float(args.beta),
                   norm_rule=NORM_RULE, donor_banks=[rel(d) for d in donors],
                   recipient_bank=rel(recip),
                   dropped_recipient_from_donors=dropped,
                   n_training_targets=len(axes["training"]),
                   n_covered_genes=len(axes["covered"]),
                   n_query_targets=len(axes["test"]),
                   n_panel_total=len(axes["panel"]),
                   n_panel_missing_no_consensus=len(axes["missing"]),
                   missing_targets=axes["missing"],
                   recipient_panel_overlap=axes["recip_panel_overlap"],
                   n_recipient_panel_overlap=n_panel_overlap,
                   n_panel_genes_dropped_from_axis=len(axes["panel_on_axis"]),
                   panel_genes_dropped=axes["panel_on_axis"],
                   n_training_genes_dropped_from_axis=len(axes["train_on_axis"]),
                   common_genic_pre_exclusion=int(axes["common_n"]),
                   recipient_measured_mask=axes["mask_note"],
                   ridge_rows=receipt["n_rows"], ridge_features=receipt["n_features"],
                   ridge_lambda_rule="core341.fit: lambda = number of rows",
                   ridge_lambda=int(receipt["n_rows"]),
                   held_target_poison_pass=True,
                   subset_linear_max_difference=receipt["subset_linear_max_difference"],
                   normal_equation_max_residual=receipt["normal_equation_max_residual"],
                   covered_fraction_of_canonical=float(
                       len(axes["covered"]) / len(axes["canonical"])),
                   heldout_within_bank=ho,
                   bank_shas=bank_shas,
                   started_utc=started,
                   finished_utc=datetime.now(timezone.utc).isoformat(),
                   seconds=time.monotonic() - clock,
                   max_rss_bytes_macos=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    sidecar_path = out.with_suffix(".json")
    sidecar_path.write_text(json.dumps(sidecar, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: sidecar[k] for k in (
        "n_training_targets", "n_covered_genes", "n_query_targets", "ridge_lambda",
        "held_target_poison_pass", "dropped_recipient_from_donors")}, indent=2))
    print(f"wrote {out} + {sidecar_path}", flush=True)


if __name__ == "__main__":
    main()
