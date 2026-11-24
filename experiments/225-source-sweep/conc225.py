"""225 — INBOX-CONCENTRATION: leave-top-k-out retention on BOTH cosine and `pds`, separately.

222 found that on its `h1union` arm the two channels have opposite fragility: Δ`pds`'s top-3
targets per fold carry 99.8% of the mean and leave-top-3-out retains **14.6%**, while the cosine
gain's leave-top-3-out retains **72.7%**.  A gain carried by three targets out of ~290 is a
different object from a gain spread over the panel, and only the second survives a new panel.

The INBOX asks specifically that **223's DLD-1 (+0.00228) gets this test before it is recommended
for a slot** — its two-library falsifier established that the gain is real, not that it is
distributed, and those are different claims.

Definition is 222's `conc222.py` verbatim: per-target arrays straight out of 218's own
`score_effects`, leave-top-k-out computed by DELETING those targets from the average. Three cancer
folds. Nothing new is defined here.

    uv run python conc225.py
"""
from __future__ import annotations

import json

import numpy as np

from harness225 import (CANCER, WORK, Fold, merged, solve_fraction_carrying,  # noqa: E402
                        source_residual, union_rows)

ARMS = {
    "dld1_lowmoi":            (("dld1_lowmoi",), 0.35),
    "CUM4":                   (("dld1_lowmoi", "jiang_5line_noK562",
                                "igvf_cardiomyocyte", "myllia_a549"), 0.28),
    "ALLSMALL":               (("dld1_lowmoi", "myllia_a549", "igvf_cardiomyocyte",
                                "tian_ipsc_neuron", "sunshine_calu3", "xucao_hek293",
                                "jiang_5line_noK562", "frangieh_melanoma", "norman_k562",
                                "tian_ipsc_neuron_crispra", "wessels_thp1",
                                "gasperini_k562", "tian2019_ipsc"), 0.35),
    "igvf_cardiomyocyte":     (("igvf_cardiomyocyte",), 0.10),
    "jiang_5line_noK562":     (("jiang_5line_noK562",), 0.28),
}
KS = (1, 3, 5, 10)


def main() -> None:
    per_fold: dict[str, dict[str, dict]] = {a: {} for a in ARMS}
    for held in CANCER:
        F = Fold(held)
        base = F.scored["base222"]
        b_pds = np.asarray(base["_pds_per_target_a23"], dtype=np.float64)
        b_cos = np.asarray(base["_cos_per_target_a23"], dtype=np.float64)
        for tag, (names, g) in ARMS.items():
            live = [n for n in names if n != held]
            R = union_rows([source_residual(n, F.targets, "gene_norm") for n in live])
            fr, _ = solve_fraction_carrying(F.whitened, R, g)
            assert fr > 0.0, f"{held}/{tag} sized to zero"
            m = merged(F.cd4, F.frac_cd4, F.h1, F.frac_union)
            e = F.score(merged(m, 1.0, R, fr), 1.0)
            d_pds = np.asarray(e["_pds_per_target_a23"], dtype=np.float64) - b_pds
            d_cos = np.asarray(e["_cos_per_target_a23"], dtype=np.float64) - b_cos
            carry = np.asarray(R.target_condition_count) > 0
            per_fold[tag][held] = dict(
                targets=list(F.targets), carry=carry.tolist(),
                d_pds=d_pds.tolist(), d_cos=d_cos.tolist(),
                fold_d_pds=float(e["pds"] - base["pds"]),
                fold_d_cos=float(e["cos"] - base["cos"]))
            print(f"[{held}/{tag}] n={len(live)} g={g} fold Δpds "
                  f"{e['pds'] - base['pds']:+.6f} Δcos {e['cos'] - base['cos']:+.6f}",
                  flush=True)
        del F
        import resid218
        resid218._PREP_CACHE.clear()

    out = {}
    print("\n" + "=" * 100)
    print("LEAVE-TOP-k-OUT RETENTION — pooled over the three cancer folds, BOTH channels")
    print("=" * 100)
    print(f"{'arm':22s}{'k':>3s}{'Δpds full':>11s}{'Δpds LOO':>11s}{'pds ret':>9s}"
          f"{'top-k share':>12s}{'Δcos full':>11s}{'Δcos LOO':>11s}{'cos ret':>9s}")
    for tag in ARMS:
        full_p = float(np.mean([per_fold[tag][h]["fold_d_pds"] for h in CANCER]))
        full_c = float(np.mean([per_fold[tag][h]["fold_d_cos"] for h in CANCER]))
        rec = dict(full_d_pds=full_p, full_d_cos=full_c, k={})
        for k in KS:
            lp, lc, sh_p, sh_c, top = [], [], [], [], {}
            for h in CANCER:
                f = per_fold[tag][h]
                dp = np.asarray(f["d_pds"]); dc = np.asarray(f["d_cos"])
                # each channel is ranked by ITS OWN per-target delta (222's definition)
                for arr, L, S in ((dp, lp, sh_p), (dc, lc, sh_c)):
                    o = np.argsort(-arr)[:k]
                    keep = np.ones(arr.size, bool); keep[o] = False
                    L.append(float(arr[keep].mean()))
                    S.append(float(arr[o].sum() / arr.sum()) if arr.sum() else np.nan)
                if k == 3:
                    o = np.argsort(-dp)[:3]
                    top[h] = [f["targets"][i] for i in o]
            mp, mc = float(np.mean(lp)), float(np.mean(lc))
            rec["k"][str(k)] = dict(
                loo_d_pds=mp, loo_d_cos=mc,
                pds_retained=(mp / full_p if full_p else np.nan),
                cos_retained=(mc / full_c if full_c else np.nan),
                top_share_pds=float(np.mean(sh_p)), top_share_cos=float(np.mean(sh_c)),
                **({"top3_targets_by_pds": top} if k == 3 else {}))
            r = rec["k"][str(k)]
            print(f"{tag:22s}{k:>3d}{full_p:>+11.5f}{mp:>+11.5f}"
                  f"{100 * r['pds_retained']:>8.1f}%{100 * r['top_share_pds']:>11.1f}%"
                  f"{full_c:>+11.6f}{mc:>+11.6f}{100 * r['cos_retained']:>8.1f}%")
        out[tag] = rec
        t3 = rec["k"]["3"]["top3_targets_by_pds"]
        print(f"{'':22s}    top-3 Δpds targets: "
              + " | ".join(f"{h[:8]}: {','.join(v)}" for h, v in t3.items()))
    (WORK / "conc225.json").write_text(json.dumps(
        dict(summary=out, per_fold=per_fold), indent=2))
    print(f"\n-> {WORK / 'conc225.json'}")


if __name__ == "__main__":
    main()
