"""216 Phase 1c — how much scored energy can the CD4 residual never reach?

Phase 1 established that the published release's `GWCD4i.DE_stats.h5ad` is byte-for-byte the
object our banks were built from, and that the ONE materially different processed object is
`GWCD4i.pseudobulk_merged.h5ad`: **17,771** genes on the official axis against the DE table's
**10,113**, a strict superset (0 DE genes missing from it).

That difference is only worth anything if the 7,658 extra genes carry scored energy.  This
measures it directly on the four PLOLO folds -- on the truth (what the metric reads), on the
delivered prediction, and under production's CPM > 10 gate, which is what actually decides
whether a gene can contribute at all.

    uv run python energy216.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
from plolo import FOLDS, OFFICIAL, build_fold, cd4_for  # noqa: E402

CD4_ROOT = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"


def main() -> None:
    import pandas as pd

    gm = pd.read_csv(CD4_ROOT / "gene_map.csv")
    cd4_genes = set(gm["gene_name"].astype(str))
    extra = set(json.loads((WORK / "extra_official_genes_pseudobulk.json").read_text()))
    pos = {g: i for i, g in enumerate(OFFICIAL)}
    cd4_mask = np.zeros(len(OFFICIAL), dtype=bool)
    for g in cd4_genes:
        if g in pos:
            cd4_mask[pos[g]] = True
    extra_mask = np.zeros(len(OFFICIAL), dtype=bool)
    for g in extra:
        if g in pos:
            extra_mask[pos[g]] = True

    rep: dict[str, object] = {
        "cd4_genes_on_official_axis": int(cd4_mask.sum()),
        "pseudobulk_extra_genes_on_official_axis": int(extra_mask.sum()),
        "official_axis": len(OFFICIAL),
        "overlap_cd4_and_extra": int((cd4_mask & extra_mask).sum()),
    }
    folds = {}
    for held in FOLDS:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)
        scored = np.asarray(fold["scored"])
        cd4 = cd4_for(fold["targets"])
        votes = np.asarray(cd4.gene_votes)

        T = truth[:, scored]
        e_tot = float((T * T).sum())
        m_cd4 = cd4_mask[scored]
        m_extra = extra_mask[scored]
        row = dict(
            n_scored=int(len(scored)),
            scored_in_cd4_axis=int(m_cd4.sum()),
            scored_in_pseudobulk_extra=int(m_extra.sum()),
            truth_energy_on_cd4_axis=float((T[:, m_cd4] ** 2).sum() / e_tot),
            truth_energy_on_pseudobulk_extra=float((T[:, m_extra] ** 2).sum() / e_tot),
            truth_energy_off_both=float(
                (T[:, ~(m_cd4 | m_extra)] ** 2).sum() / e_tot))

        # the residual's ACTUAL support: a gene only receives the residual where CD4 voted
        sup = votes[:, scored] > 0
        row["truth_energy_on_actual_cd4_support"] = float((T[sup] ** 2).sum() / e_tot)
        row["support_density"] = float(sup.mean())

        # under production's CPM gate, per challenge frame
        gated = {}
        for ctx in ("A", "B", "C"):
            b = plolo.challenge_bctrl(ctx)
            expr = np.expm1(b)
            cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
            g = (cpm > plolo.MIN_CPM) & (b > 0)
            gs = g[scored]
            eg = float((T[:, gs] ** 2).sum())
            gated[ctx] = dict(
                n_gated_scored=int(gs.sum()),
                gated_in_cd4_axis=int((gs & m_cd4).sum()),
                gated_in_pseudobulk_extra=int((gs & m_extra).sum()),
                truth_energy_gated_frac=float(eg / e_tot),
                cd4_share_of_gated_energy=float((T[:, gs & m_cd4] ** 2).sum() / eg),
                extra_share_of_gated_energy=float((T[:, gs & m_extra] ** 2).sum() / eg))
        row["gated"] = gated
        row["gated_mean_cd4_share"] = float(np.mean(
            [gated[c]["cd4_share_of_gated_energy"] for c in "ABC"]))
        row["gated_mean_extra_share"] = float(np.mean(
            [gated[c]["extra_share_of_gated_energy"] for c in "ABC"]))
        folds[held] = row
        print(f"[{held}] truth energy on CD4 axis {row['truth_energy_on_cd4_axis']:.4f}, "
              f"on pseudobulk-extra {row['truth_energy_on_pseudobulk_extra']:.4f}; "
              f"gated: cd4 {row['gated_mean_cd4_share']:.4f} "
              f"extra {row['gated_mean_extra_share']:.4f}", flush=True)
        del fold, truth, T, cd4
    rep["folds"] = folds
    rep["macro"] = {
        k: float(np.mean([folds[h][k] for h in FOLDS]))
        for k in ("truth_energy_on_cd4_axis", "truth_energy_on_pseudobulk_extra",
                  "truth_energy_off_both", "truth_energy_on_actual_cd4_support",
                  "gated_mean_cd4_share", "gated_mean_extra_share", "support_density")}
    # the headroom a pseudobulk re-extraction could address, as a multiplier on reachable energy
    m = rep["macro"]
    rep["headroom"] = dict(
        reachable_now=m["gated_mean_cd4_share"],
        reachable_with_pseudobulk=m["gated_mean_cd4_share"] + m["gated_mean_extra_share"],
        multiplier=(m["gated_mean_cd4_share"] + m["gated_mean_extra_share"])
        / max(m["gated_mean_cd4_share"], 1e-12),
        note="upper bound: assumes the extra genes carry effect as informative per unit of "
             "energy as the ones already covered, which is optimistic since they are the "
             "genes DESeq2's own expression filter dropped")
    (WORK / "energy216.json").write_text(json.dumps(rep, indent=2))
    print(json.dumps({k: rep[k] for k in ("cd4_genes_on_official_axis",
                                          "pseudobulk_extra_genes_on_official_axis",
                                          "macro", "headroom")}, indent=2))


if __name__ == "__main__":
    main()
