"""Build a VCC 2026 submission from the donor-transfer model.

Reads the official controls bundle, predicts all 300 perturbations in each of the three
contexts, and writes a format-valid .h5ad ready for `vcc prep`.

The model is `predicted effect = beta * m + alpha * s_g` (see src/vcc2026/donor_transfer.py).
Perturbations the donor does not carry get `beta * m` only — the mean response, which scores
about zero but avoids the catastrophic `fid` penalty for calling nothing.

    uv run python scripts/make_submission.py --alpha 0.1 --dispersion 0.1 \
        -o submissions/v1.h5ad
    vcc prep submissions/v1.h5ad -g data/challenge/gene_names.csv \
        --perts data/challenge/pert_counts.csv -o submissions/v1.vcc --dry-run
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_submission import (  # noqa: E402
    CELLS_PER_PERT, CONTEXT_COL, MAX_STORED_ENTRIES, N_GENES, PERT_COL,
    ChallengeSpec, _write_categorical, _write_str, _str_dtype,
)
from vcc2026.donor_transfer import Donor, TransferModel, emit_counts, TS  # noqa: E402
from vcc2026.multi_donor import MultiDonorModel  # noqa: E402
from vcc2026.whiten import Whitener  # noqa: E402

import h5py  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bundle", default="data/challenge")
    ap.add_argument("--donor", default="data/donors/k562_gwps")
    ap.add_argument("--donors", nargs="+", default=None,
                    help="name=path[:weight] for multi-donor mode, e.g. "
                         "k562_gwps=data/donors/k562_gwps:2 orion_hct116=data/donors/orion_hct116:1")
    ap.add_argument("--no-eb", action="store_true", help="disable empirical-Bayes shrinkage")
    ap.add_argument("-o", "--output", default="submissions/v1.h5ad")
    ap.add_argument("--alpha", type=float, default=0.1)
    ap.add_argument("--beta", type=float, default=0.0)
    ap.add_argument("--dispersion", type=float, default=0.1)
    ap.add_argument("--spread", type=float, default=0.0,
                    help="fraction of Poisson cell-to-cell noise for --emit quasi "
                         "(experiment 034). Each gene's panel total is preserved EXACTLY, so "
                         "the pseudobulk -- and therefore pds, mse and nmae -- is bit-"
                         "identical to the even emitter's; only the Wilcoxon's view of "
                         "cell-to-cell spread changes, which is what fid, reach and jac read.")
    ap.add_argument("--emit", default="even", choices=["nb", "even", "quasi"],
                    help="cell emission model. 'even' is minimum-variance allocation, which "
                         "maximises DE significance and therefore fid coverage (exp 005/007).")
    ap.add_argument("--min-cpm", type=float, default=5.0,
                    help="predict no change below this CPM in the target's controls. "
                         "5 matches the scorer's own DE filter; genes below it are never "
                         "tested, and gating higher makes us predict 0 on genes the "
                         "reference DOES call, which is a full nmae penalty.")
    ap.add_argument("--space", default="log1p", choices=["log1p", "foldchange"])
    ap.add_argument("--no-gate", action="store_true")
    ap.add_argument("--whiten-gamma", type=float, default=0.0,
                    help="spectral whitening exponent for the panel effect matrix "
                         "(experiment 017). 0 disables. The operator is fitted on donor "
                         "effects only and never sees the target context.")
    ap.add_argument("--whiten-lam", type=float, default=0.0)
    ap.add_argument("--whiten-basis-perts", type=int, default=6000,
                    help="perturbations sampled per donor to estimate the effect covariance")
    ap.add_argument("--whiten-cache", default="data/donors/whitener.npz")
    ap.add_argument("--whiten-axis", default="inter", choices=["inter", "union"],
                    help="which genes the operator covers. 'inter' is genes EVERY donor "
                         "measures (7,312) -- what v6/v8 shipped, and it leaves ~a third of "
                         "the tested universe untransformed under our largest lever. 'union' "
                         "is genes ANY donor measures (17,942); experiment 038 measures "
                         "+0.022 macro PDS for it, with the gamma optimum moving 0.85 -> 0.65.")
    ap.add_argument("--sharpen", type=float, default=1.0,
                    help="power transform on predicted effect magnitudes, "
                         "sign(e)*|e|^q, applied AFTER whitening (experiment 030: 'after' "
                         "beats 'before' at every exponent). 1.0 disables. Magnitude is "
                         "renormalised afterwards so alpha keeps its meaning.")
    ap.add_argument("--whiten-renorm", default="global", choices=["global", "preserve"],
                    help="how to restore magnitude after whitening. 'global' rescales the "
                         "whole effect matrix, which silently down-weights the ~23%% of "
                         "tested genes the operator axis does not cover by the same ~4.6x "
                         "the on-axis block grew. 'preserve' rescales only the on-axis "
                         "block back to its own prior norm, leaving off-axis genes alone "
                         "(experiment 025: +0.006 macro PDS).")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    spec = ChallengeSpec.load(args.bundle)
    weights = None
    if args.donors:
        dd, weights, cpp = {}, {}, {}
        for item in args.donors:
            name, rest = item.split("=", 1)
            path, _, wt = rest.partition(":")
            dd[name] = Donor.load(Path(path))
            weights[name] = float(wt) if wt else 1.0
            mp = Path(path) / "meta.json"
            cpp[name] = (json.loads(mp.read_text())["cells_per_pert"]["median"]
                         if mp.exists() else 150)
        model = MultiDonorModel(dd, spec.genes, gene_norm=False,
                                eb=not args.no_eb, cells_per_pert=cpp)
        cov = model.coverage(spec.perts)
        print(f"panel {len(spec.perts)} perturbations | donors {weights} | "
              f"coverage {({k: round(v, 3) for k, v in cov.items()})}")
    else:
        model = TransferModel(Donor.load(Path(args.donor)), spec.genes)
        covered = model.covered(spec.perts)
        print(f"panel {len(spec.perts)} perturbations, donor covers {covered.sum()} "
              f"({covered.mean():.0%}); gene-axis coverage {model.coverage:.1%}")
    print(f"alpha={args.alpha} beta={args.beta} dispersion={args.dispersion} "
          f"emit={args.emit} min_cpm={args.min_cpm} space={args.space} "
          f"whiten_gamma={args.whiten_gamma}\n")

    whitener = None
    wcols = None
    if args.whiten_gamma > 0:
        if weights is None:
            raise SystemExit("--whiten-gamma requires multi-donor mode (--donors)")
        # Operator axis: genes every donor measures, minus the panel target genes -- exactly
        # the axis experiment 017 validated on. Genes outside it keep their raw effect.
        if args.whiten_axis == "union":
            shared = set()
            for n in model.names:
                shared |= set(model.target_genes[model.tcols[n]])
        else:
            shared = set(model.target_genes[model.tcols[model.names[0]]])
            for n in model.names[1:]:
                shared &= set(model.target_genes[model.tcols[n]])
        shared -= set(spec.perts)
        tindex = {g: i for i, g in enumerate(model.target_genes)}
        wcols = np.array(sorted(tindex[g] for g in shared), dtype=np.int64)
        cache = Path(args.whiten_cache)
        if cache.exists() and json.loads(cache.with_suffix(".json").read_text()).get(
                "n_cols") == len(wcols):
            print(f"  whitener: loading cached operator ({len(wcols)} genes)", flush=True)
            whitener = Whitener.load(cache)
        else:
            rng0 = np.random.default_rng(0)
            rows = []
            for n in model.names:
                S = np.zeros((model.s[n].shape[0], len(model.target_genes)), np.float32)
                S[:, model.tcols[n]] = model.s[n][:, model.dcols[n]]
                if S.shape[0] > args.whiten_basis_perts:
                    S = S[np.sort(rng0.choice(S.shape[0], args.whiten_basis_perts, False))]
                rows.append(model._renorm[n] * S[:, wcols].astype(np.float64))
            B = np.vstack(rows)
            print(f"  whitener: fitting on donor basis {B.shape} ...", flush=True)
            whitener = Whitener.fit(B)
            whitener.save(cache)
            cache.with_suffix(".json").write_text(json.dumps(
                {"n_cols": int(len(wcols)), "basis": list(B.shape),
                 "donors": model.names}))
            del B
        print(f"  whitener ready: gamma={args.whiten_gamma} lam={args.whiten_lam}",
              flush=True)

    # Two independent streams. `pick_rng` chooses which control cells' library sizes to
    # reuse; `rng` drives emission. They must not share state: the even emitter consumes no
    # random numbers and the quasi emitter consumes many, so a single stream makes the two
    # builds diverge in which control cells they sample -- which shifts the emitted depth and
    # breaks the claim that the mean profile is untouched. Seeding `pick_rng` with the same
    # seed reproduces the even-emitter builds' picks exactly, so a quasi build differs from
    # its even counterpart in emission and in nothing else.
    pick_rng = np.random.default_rng(args.seed)
    rng = np.random.default_rng(args.seed + 1_000_003)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    n_rows = spec.n_rows

    obs_pert: list[str] = []
    obs_ctx: list[str] = []
    indptr = np.zeros(n_rows + 1, dtype=np.int64)
    nnz = 0
    row = 0

    with h5py.File(out, "w") as f:
        X = f.create_group("X")
        X.attrs["encoding-type"] = "csr_matrix"
        X.attrs["encoding-version"] = "0.1.0"
        X.attrs["shape"] = np.array([n_rows, N_GENES], dtype=np.int64)
        d_data = X.create_dataset("data", shape=(0,), maxshape=(None,),
                                  dtype=np.float32, chunks=(1 << 20,))
        d_idx = X.create_dataset("indices", shape=(0,), maxshape=(None,),
                                 dtype=np.int32, chunks=(1 << 20,))

        for ctx in spec.contexts:
            path = Path(args.bundle) / f"context_{ctx}.h5ad"
            a = ad.read_h5ad(path)
            Xc = a.X.tocsr() if not sp.isspmatrix_csr(a.X) else a.X
            if list(a.var_names) != list(spec.genes):
                raise SystemExit(f"context {ctx} gene order differs from gene_names.csv")
            P = np.asarray(Xc.sum(0)).ravel().astype(np.float64)
            b_ctrl = np.log1p(TS * P / P.sum())
            libs = np.asarray(Xc.sum(1)).ravel()
            del a, Xc
            print(f"  context {ctx}: b_ctrl from {len(libs):,} control cells, "
                  f"median lib {np.median(libs):,.0f}", flush=True)

            # With whitening on, the whole panel has to be transformed together -- the
            # operator is a statement about the *set* of 300 predictions, not about any one
            # of them -- so effects are materialised as a matrix before any cell is emitted.
            eff_mat = None
            if whitener is not None:
                eff_mat = np.vstack([model.effect(p, args.alpha, args.beta, weights=weights)
                                     for p in spec.perts]).astype(np.float64)
                before_all = float(np.median(np.linalg.norm(eff_mat, axis=1)))
                before_op = float(np.median(np.linalg.norm(eff_mat[:, wcols], axis=1)))
                eff_mat[:, wcols] = whitener.apply(eff_mat[:, wcols], args.whiten_gamma,
                                                   args.whiten_lam, renorm=False)
                after_op = float(np.median(np.linalg.norm(eff_mat[:, wcols], axis=1))) or 1.0
                if args.whiten_renorm == "preserve":
                    eff_mat[:, wcols] *= before_op / after_op
                else:
                    after_all = float(np.median(np.linalg.norm(eff_mat, axis=1))) or 1.0
                    eff_mat *= before_all / after_all
                med = float(np.median(np.linalg.norm(eff_mat, axis=1)))
                print(f"    whitened panel ({args.whiten_renorm}): on-axis norm "
                      f"{before_op:.3f} -> {after_op:.3f}; overall {before_all:.3f} -> "
                      f"{med:.3f}", flush=True)
            if args.sharpen != 1.0:
                if eff_mat is None:
                    eff_mat = np.vstack([model.effect(p, args.alpha, args.beta,
                                                      weights=weights)
                                         for p in spec.perts]).astype(np.float64)
                pre = float(np.median(np.linalg.norm(eff_mat, axis=1)))
                eff_mat = np.sign(eff_mat) * np.abs(eff_mat) ** args.sharpen
                post = float(np.median(np.linalg.norm(eff_mat, axis=1))) or 1.0
                eff_mat *= pre / post
                print(f"    sharpened q={args.sharpen}: median effect norm held at "
                      f"{pre:.3f}", flush=True)

            for pi, pert in enumerate(spec.perts):
                if eff_mat is not None:
                    prof = model.profile_from_effect(
                        eff_mat[pi], b_ctrl, min_cpm=args.min_cpm,
                        gate_unexpressed=not args.no_gate)
                elif weights is not None:
                    prof = model.profile(pert, b_ctrl, args.alpha, args.beta,
                                         weights=weights, min_cpm=args.min_cpm,
                                         gate_unexpressed=not args.no_gate)
                else:
                    prof = model.profile(pert, b_ctrl, args.alpha, args.beta,
                                         space=args.space,
                                         gate_unexpressed=not args.no_gate,
                                         min_cpm=args.min_cpm)
                pick = pick_rng.choice(len(libs), spec.cells_per_pert,
                                       replace=len(libs) < spec.cells_per_pert)
                if args.emit == "quasi":
                    data, idx, ip = emit_counts(
                        prof, libs[pick], args.dispersion, rng, mode="quasi",
                        spread=np.full(prof.size, args.spread))
                else:
                    data, idx, ip = emit_counts(prof, libs[pick], args.dispersion, rng,
                                                mode=args.emit)

                d_data.resize((nnz + data.size,)); d_data[nnz:] = data
                d_idx.resize((nnz + idx.size,));  d_idx[nnz:] = idx
                indptr[row + 1: row + 1 + spec.cells_per_pert] = nnz + ip[1:]
                nnz += data.size
                row += spec.cells_per_pert
                obs_pert.extend([pert] * spec.cells_per_pert)
                obs_ctx.extend([ctx] * spec.cells_per_pert)

                if nnz > MAX_STORED_ENTRIES:
                    raise SystemExit(f"density cap exceeded at row {row}: {nnz:,}")
            print(f"    -> {row:,} rows, {nnz:,} stored entries "
                  f"({nnz/row:.0f}/cell, {100*nnz/MAX_STORED_ENTRIES:.1f}% of cap)", flush=True)

        X.create_dataset("indptr", data=indptr, dtype=np.int64)

        obs = f.create_group("obs")
        obs.attrs["encoding-type"] = "dataframe"; obs.attrs["encoding-version"] = "0.2.0"
        obs.attrs["_index"] = "_index"
        obs.attrs.create("column-order", [PERT_COL, CONTEXT_COL], dtype=_str_dtype())
        _write_str(obs, "_index", np.array([f"cell_{i}" for i in range(n_rows)], dtype=object))
        _write_categorical(obs, PERT_COL, np.asarray(obs_pert, dtype=object))
        _write_categorical(obs, CONTEXT_COL, np.asarray(obs_ctx, dtype=object))

        var = f.create_group("var")
        var.attrs["encoding-type"] = "dataframe"; var.attrs["encoding-version"] = "0.2.0"
        var.attrs["_index"] = "_index"
        var.attrs.create("column-order", np.empty(0, dtype=object), dtype=_str_dtype())
        _write_str(var, "_index", spec.genes.astype(object))

        for name in ("layers", "obsm", "varm", "obsp", "varp", "uns"):
            g = f.create_group(name)
            g.attrs["encoding-type"] = "dict"; g.attrs["encoding-version"] = "0.1.0"
        f.attrs["encoding-type"] = "anndata"; f.attrs["encoding-version"] = "0.1.0"

    meta = dict(alpha=args.alpha, beta=args.beta, dispersion=args.dispersion,
                whiten_gamma=args.whiten_gamma, whiten_lam=args.whiten_lam,
                whiten_renorm=args.whiten_renorm, sharpen=args.sharpen,
                whiten_axis=args.whiten_axis,
                emit=args.emit, spread=args.spread, min_cpm=args.min_cpm,
                space=args.space, gate=not args.no_gate, seed=args.seed,
                donor=args.donor, donors=args.donors, eb=not args.no_eb,
                n_rows=n_rows, nnz=int(nnz),
                nnz_per_cell=float(nnz / n_rows),
                pct_of_density_cap=float(100 * nnz / MAX_STORED_ENTRIES),
                donor_pert_coverage=(float(model.coverage(spec.perts)["any"])
                                     if weights is not None else float(covered.mean())))
    Path(str(out) + ".meta.json").write_text(json.dumps(meta, indent=2))
    print(f"\nwrote {out}  ({out.stat().st_size/1e9:.2f} GB)")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
