"""Verify the two built donor tables against the contract, with a control that can fail.

Three checks, in increasing strength:

1. `Donor.load` round-trips and the shapes agree.

2. The experiment-161 mass invariant, `sum(expm1(b_ctrl)) == sum(expm1(b_ctrl + delta)) == TS`.
   **On its own this is close to vacuous**: the builder writes `delta = prof - b_ctrl`, so in
   exact arithmetic the identity holds by construction and the check only measures float32
   round-trip error. It is reported here *comparatively* -- against the same statistic on the
   three shipped production donors -- so the question it answers is "is this table's numerical
   error the same order as a table we already trust?", which is a real question, rather than
   "does subtraction invert addition?", which is not.

3. The check that carries the weight: an INDEPENDENT recomputation of sampled effect rows
   straight from the raw cells, through `anndata`'s own backed row-slicing rather than the
   builder's group-sum path. A wrong cell-to-group assignment, an off-by-one in the
   perturbation index, a mismatched gene collapse or a wrong control group all break this and
   none of them break check 2.

   And its NEGATIVE CONTROL: the same recomputation with the target labels permuted must
   *fail*. Without that, "cosine 1.000" only proves the two code paths agree on something --
   it does not prove they agree on the right thing. AGENTS.md 11.1 has two logged cases of a
   control that was invariant by construction; this is the guard against a third.

Writes `work/verify.json`.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO / "src"))

from vcc2026.donor_transfer import TS, Donor  # noqa: E402

SEED = 20260902

BUILT = {
    "nadig_hepg2": dict(
        dir=HERE / "donors/nadig_hepg2",
        raw=REPO / "data/external/nadig2025/hepg2_raw_singlecell.h5ad",
        label_col="gene", gene_col="gene_name", control="non-targeting",
        stratum_col=None, stratum=None,
    ),
    "frangieh_melanoma": dict(
        dir=HERE / "donors/frangieh_melanoma",
        raw=REPO / "data/external/scperturb/FrangiehIzar2021_RNA.h5ad",
        label_col="perturbation", gene_col=None, control="control",
        stratum_col="perturbation_2", stratum="Control",
    ),
}
SHIPPED = ["k562_gwps", "orion_hct116", "orion_hek293t"]


def mass(d: Donor) -> dict:
    """Total delivered mass implied by b_ctrl and by each b_ctrl + delta row."""
    m0 = float(np.expm1(d.b_ctrl.astype(np.float64)).sum())
    rows = np.expm1((d.b_ctrl[None, :] + d.delta).astype(np.float64)).sum(axis=1)
    return dict(
        ctrl_mass=m0, ctrl_rel_err=abs(m0 - TS) / TS,
        pert_mass_median=float(np.median(rows)),
        pert_rel_err_max=float(np.abs(rows - TS).max() / TS),
        pert_rel_err_median=float(np.median(np.abs(rows - TS)) / TS),
    )


def recompute(cfg: dict, targets: list[str], genes_out: np.ndarray,
              permute: bool = False, rng: np.random.Generator | None = None):
    """Recompute effect rows from raw cells via anndata backed slicing.

    Deliberately a different code path from the builders: anndata resolves the row index and
    the sparse/dense decoding itself, so a mistake shared with the builder would have to be a
    mistake in anndata.
    """
    a = ad.read_h5ad(cfg["raw"], backed="r")
    lab = a.obs[cfg["label_col"]].astype(str).to_numpy()
    if permute:
        lab = lab.copy()
        rng.shuffle(lab)
    keep = np.ones(a.n_obs, dtype=bool)
    if cfg["stratum_col"]:
        keep = a.obs[cfg["stratum_col"]].astype(str).to_numpy() == cfg["stratum"]

    gsrc = np.asarray(a.var[cfg["gene_col"]] if cfg["gene_col"] else a.var_names, dtype=str)
    order = {g: i for i, g in enumerate(genes_out)}
    col = np.fromiter((order.get(g, -1) for g in gsrc), dtype=np.int64, count=len(gsrc))

    def pseudobulk(mask):
        idx = np.flatnonzero(mask)
        X = a[idx].to_memory().X
        v = np.asarray(X.sum(axis=0)).ravel().astype(np.float64)
        out = np.zeros(len(genes_out))          # collapse duplicate symbols by summing
        np.add.at(out, col[col >= 0], v[col >= 0])
        return out

    Pc = pseudobulk(keep & (lab == cfg["control"]))
    b_ctrl = np.log1p(TS * Pc / Pc.sum())
    rows = {}
    for t in targets:
        P = pseudobulk(keep & (lab == t))
        if P.sum() <= 0:
            continue
        rows[t] = np.log1p(TS * P / P.sum()) - b_ctrl
    return b_ctrl, rows


def main() -> None:
    rng = np.random.default_rng(SEED)
    report: dict = {"seed": SEED, "TS": TS, "tables": {}, "shipped_reference": {}}

    for name in SHIPPED:
        d = Donor.load(REPO / "data/donors" / name)
        report["shipped_reference"][name] = mass(d)
        m = report["shipped_reference"][name]
        print(f"[shipped] {name:<16} ctrl_rel_err {m['ctrl_rel_err']:.3e}  "
              f"pert_rel_err median {m['pert_rel_err_median']:.3e} "
              f"max {m['pert_rel_err_max']:.3e}")
    print()

    panel = set(pd.read_csv(REPO / "data/challenge/pert_counts.csv")["target_gene"])

    for name, cfg in BUILT.items():
        print(f"=== {name}")
        d = Donor.load(cfg["dir"])
        assert d.delta.shape == (len(d.perts), len(d.genes)) == (d.delta.shape[0], d.b_ctrl.size)
        rec: dict = {
            "path": str(cfg["dir"].relative_to(REPO)),
            "shape": list(d.delta.shape),
            "n_perts": len(d.perts), "n_genes": len(d.genes),
            "panel_overlap": len(panel & set(d.perts.tolist())),
            "mass": mass(d),
        }
        m = rec["mass"]
        print(f"  load OK: delta {d.delta.shape}, panel {rec['panel_overlap']}/300")
        print(f"  mass: ctrl_rel_err {m['ctrl_rel_err']:.3e}  "
              f"pert_rel_err median {m['pert_rel_err_median']:.3e} "
              f"max {m['pert_rel_err_max']:.3e}")

        # --- check 3: independent recomputation of 5 sampled rows ---
        pick = rng.choice(len(d.perts), size=5, replace=False)
        targets = [str(d.perts[i]) for i in pick]
        b_ref, rows = recompute(cfg, targets, d.genes)
        cos, mad = [], []
        for t in targets:
            i = int(np.flatnonzero(d.perts == t)[0])
            a_, b_ = d.delta[i].astype(np.float64), rows[t]
            cos.append(float(a_ @ b_ / (np.linalg.norm(a_) * np.linalg.norm(b_))))
            mad.append(float(np.abs(a_ - b_).max()))
        rec["recompute"] = dict(
            targets=targets, cosine=cos, max_abs_diff=mad,
            b_ctrl_max_abs_diff=float(np.abs(d.b_ctrl.astype(np.float64) - b_ref).max()),
        )
        print(f"  recompute {targets}")
        print(f"    cosine  min {min(cos):.6f}  max {max(cos):.6f}")
        print(f"    max|diff| {max(mad):.3e}   b_ctrl max|diff| "
              f"{rec['recompute']['b_ctrl_max_abs_diff']:.3e}")

        # --- negative control: permuted labels must NOT reproduce the table ---
        _, prows = recompute(cfg, targets, d.genes, permute=True,
                             rng=np.random.default_rng(SEED + 1))
        pcos = []
        for t in targets:
            if t not in prows:
                continue
            i = int(np.flatnonzero(d.perts == t)[0])
            a_, b_ = d.delta[i].astype(np.float64), prows[t]
            n = np.linalg.norm(a_) * np.linalg.norm(b_)
            pcos.append(float(a_ @ b_ / n) if n > 0 else 0.0)
        rec["permuted_control_cosine"] = pcos
        print(f"    PERMUTED control cosine  min {min(pcos):.4f}  max {max(pcos):.4f}"
              f"   (must be far below 1 -- else the check is vacuous)")
        rec["verdict"] = ("PASS" if min(cos) > 0.999 and max(pcos) < 0.5
                          else "FAIL")
        print(f"  verdict: {rec['verdict']}\n")
        report["tables"][name] = rec

    (HERE / "work/verify.json").write_text(json.dumps(report, indent=2))
    print("wrote work/verify.json")


if __name__ == "__main__":
    main()
