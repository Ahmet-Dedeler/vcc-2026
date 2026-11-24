"""805 stage 0 statistics (EC2): known answers first, then the gate on P, the secondary on S (PREREGISTRATION-stage0.md).

    uv run python experiments/805-cd4-flex-surface/s0stats805.py ka   <ka_ref_dir>                 -> exit 1 on a KA fail
    uv run python experiments/805-cd4-flex-surface/s0stats805.py gate <condition> <outdir>         -> <outdir>/gate805.json

`ka`  : KA-1 = ref614.py's DE parquet sha on h1-test equals 614's ref_h1ka (1bb8c76d…); KA-2 = blob615.old_surface on
        h1-test reproduces 615's split-half 0.0487 / plain 0.0376 (|d| <= 0.002).
`gate`: n_real per target from <outdir>/ref_{P,S}/reference.json; knockdown of the own gene; blob split-half share of the
        P truth; Spearman of n_real vs the published n_downstream; the verdict.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/615-orion-hek-harness"))
import blob615  # noqa: E402  (read-only: old_surface, the 615 statistic)

KA_DE_SHA = "1bb8c76d1177047569e156af5a7c40e48f8fedce0ce3cf85b7d82998c5485e99"
KA_BLOB = dict(share_sh=0.0487, share=0.0376)
SUP = ROOT / "data/external/cd4_cell2026/suppl_tables_2026-08-24"


def ka(refdir: Path) -> None:
    """KA-1 (Amendment 1): n_real identical to 614's ref_h1ka for all 85 targets; the DE frame equal to 323's within 1e-9
    on p / p_adj / lfc with the same rows and the same call set. KA-2: the blob statistic on h1-test."""
    import polars as pl
    js = json.loads((refdir / "reference.json").read_text())
    want = json.loads((ROOT / "experiments/614-k562-real-emitter/work/ref_h1ka/reference.json").read_text())
    nr = np.array(list(js["n_real_by_target"].values()), float)
    import hashlib
    ref_frame = ROOT / "experiments/323-current-family-validation/work/reference_de.parquet"
    assert hashlib.sha256(ref_frame.read_bytes()).hexdigest() == KA_DE_SHA, "323's DE frame is not 1bb8c76d"
    a = pl.read_parquet(ref_frame)
    b = pl.read_parquet(refdir / "reference_de.parquet")
    keys = ["target", "feature"]
    same_rows = a.height == b.height and a.select(keys).sort(keys).equals(b.select(keys).sort(keys))
    diff = {}
    calls_same = False
    p_ok = False
    if same_rows:
        j = a.join(b, on=keys, suffix="_b")
        diff["log2_fold_change_abs"] = float((j["log2_fold_change"] - j["log2_fold_change_b"]).abs().max())
        p_ok = True
        for c in ("p_value", "p_adj"):
            x, y = j[c].to_numpy(), j[c + "_b"].to_numpy()
            rel = np.abs(x - y) / np.maximum(np.abs(x), 1e-300)
            diff[c + "_rel_max"] = float(rel.max())
            diff[c + "_share_rel_le_1e-9"] = float((rel <= 1e-9).mean())
            p_ok = p_ok and rel.max() <= 1e-3 and (rel <= 1e-9).mean() >= 0.999
        calls_same = bool(((j["p_adj"] < 0.05) == (j["p_adj_b"] < 0.05)).all())
    ka1 = bool(js["n_real_by_target"] == want["n_real_by_target"] and same_rows and calls_same and p_ok
               and diff.get("log2_fold_change_abs", 1.0) <= 1e-9)
    genes = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).to_numpy()
    r, _ = blob615.old_surface("data/harness/h1-test.h5ad", None, genes)
    out = dict(ka1_de_sha256=js["de_sha256"], ka1_mac_sha256=KA_DE_SHA,
               ka1_n_real_identical=js["n_real_by_target"] == want["n_real_by_target"], ka1_rows_same=same_rows,
               ka1_max_abs_diff=diff, ka1_calls_same=calls_same, ka1_pass=ka1,
               ka1_n_real_median=float(np.median(nr)), ka1_targets_ge10=int((nr >= 10).sum()),
               ka2=r, ka2_pass=bool(abs(r["share_sh"] - KA_BLOB["share_sh"]) <= 0.002 and abs(r["share"] - KA_BLOB["share"]) <= 0.002))
    out["pass"] = bool(out["ka1_pass"] and out["ka2_pass"])
    (refdir.parent / "ka805.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out))
    sys.exit(0 if out["pass"] else 1)


def own_gene(refdir: Path) -> dict:
    """Own-gene lfc / q from the reference (if in the DE universe) and the CPM 'remaining' ratio from the summed counts."""
    z = np.load(refdir / "reference.npz", allow_pickle=False)
    genes, targets = z["genes"].astype(str), z["targets"].astype(str)
    sums, counts = z["sums"].astype(np.float64), z["counts"].astype(np.float64)
    keep = z["keep"].astype(bool)
    kg = genes[keep]
    kpos = {g: i for i, g in enumerate(kg)}
    gpos = {g: i for i, g in enumerate(genes)}
    ctrl_cpm = 1e6 * sums[-1] / sums[-1].sum()
    out = {}
    for i, t in enumerate(targets):
        cpm = 1e6 * sums[i] / sums[i].sum()
        g = gpos.get(t)
        rem = float(cpm[g] / ctrl_cpm[g]) if g is not None and ctrl_cpm[g] > 0 else None
        k = kpos.get(t)
        out[t] = dict(remaining_pseudobulk=rem,
                      lfc=float(z["lfc"][i][k]) if k is not None else None,
                      q=float(z["q"][i][k]) if k is not None else None)
    return out


def summarise(refdir: Path, cond: str) -> dict:
    js = json.loads((refdir / "reference.json").read_text())
    nr = pd.Series(js["n_real_by_target"], dtype=float)
    z = np.load(refdir / "reference.npz", allow_pickle=False)
    q, lfc = np.asarray(z["q"], float), np.asarray(z["lfc"], float)
    sig = np.isfinite(q) & (q < 0.05)
    de = pd.read_csv(SUP / "DE_stats.suppl_table.csv")
    de = de[de["culture_condition"] == cond].set_index("target_contrast_gene_name")
    pub = de.reindex(nr.index)["n_downstream"].astype(float)
    kd = own_gene(refdir)
    rem = np.array([v["remaining_pseudobulk"] for v in kd.values() if v["remaining_pseudobulk"] is not None])
    lf = np.array([v["lfc"] for v in kd.values() if v["lfc"] is not None])
    qs = np.array([v["q"] for v in kd.values() if v["q"] is not None])
    return dict(n_targets=int(len(nr)), n_real_median=float(nr.median()), frac_ge10=float((nr >= 10).mean()),
                frac_ge20=float((nr >= 20).mean()),
                n_real_quantiles=dict(zip(["0", "10", "25", "50", "75", "90", "100"],
                                          np.percentile(nr, [0, 10, 25, 50, 75, 90, 100]).tolist())),
                de_calls_total=int(sig.sum()), frac_down=float((lfc[sig] < 0).mean()) if sig.any() else None,
                spearman_vs_published_n_downstream=float(spearmanr(nr, pub, nan_policy="omit")[0]),
                published_n_downstream_median=float(pub.median()),
                kd_remaining_quantiles=np.percentile(rem, [10, 25, 50, 75, 90]).tolist() if len(rem) else None,
                kd_own_lfc_quantiles=np.percentile(lf, [10, 25, 50, 75, 90]).tolist() if len(lf) else None,
                kd_own_in_universe=int(len(lf)), kd_own_significant=int((qs < 0.05).sum()),
                n_keep_genes=js["n_keep_genes"], median_control_library=js["median_control_library"],
                cells_per_target=js["cells_per_target"], n_real_by_target=js["n_real_by_target"], own_gene=kd)


def gate(cond: str, outdir: Path) -> None:
    res = dict(condition=cond, P=summarise(outdir / "ref_P", cond), S=summarise(outdir / "ref_S", cond))
    genes = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).to_numpy()
    rel = str((outdir / "truth_P.h5ad").relative_to(ROOT))
    blob, _ = blob615.old_surface(rel, None, genes)
    res["blob_P"] = blob
    p = res["P"]
    power = p["n_real_median"] >= 20 and p["frac_ge10"] >= 0.60
    res["gate"] = dict(power_pass=bool(power), blob_pass=bool(blob["share_sh"] <= 0.10),
                       verdict="PASS" if (power and blob["share_sh"] <= 0.10) else "FAIL",
                       rule="median n_real >= 20 AND >= 60% targets n_real >= 10 (P); blob split-half share <= 0.10")
    res["S_would_pass_power"] = bool(res["S"]["n_real_median"] >= 20 and res["S"]["frac_ge10"] >= 0.60)
    (outdir / "gate805.json").write_text(json.dumps(res, indent=1))
    short = {k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items()
                                                    if kk not in ("n_real_by_target", "own_gene")}) for k, v in res.items()}
    print(json.dumps(short, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "ka":
        ka(Path(sys.argv[2]).resolve())
    else:
        gate(sys.argv[2], Path(sys.argv[3]).resolve())
