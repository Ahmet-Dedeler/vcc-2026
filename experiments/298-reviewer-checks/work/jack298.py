"""298 item 1 — leave-one-file-out jackknife of claim C5.

Reads ONLY `experiments/295-reserve-r3-claim/work/reserve295.csv` (295's SAVED per-gene outcome
table).  No reserve file's `X` is opened.

Statistical code (spearman, gene bootstrap, file-block bootstrap) COPIED verbatim from
`experiments/295-reserve-r3-claim/work/stats295.py` lines 36-64, per that file's own
"copied, never imported across experiments" convention.

Writes work/jackknife298.csv and work/jackknife298.json.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("~/Code/virtual-cell-challenge-2026")
SRC = ROOT / "experiments/295-reserve-r3-claim/work/reserve295.csv"
WORK = ROOT / "experiments/298-reviewer-checks/work"

SEED = 20260906          # 295's seed, so the KA row reproduces bit-for-bit
N_BOOT = 2000
ADAMSON = "AdamsonWeissman2016_GSM2406681_10X010"


# ---- stats295.py:36-42 (== lib291.py:73-78) ---------------------------------------------------
def spearman(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 3:
        return float("nan")
    return float(stats.spearmanr(a[m], b[m]).statistic)


# ---- stats295.py:45-54 ------------------------------------------------------------------------
def boot_rho(a, b, groups, rng, n: int = N_BOOT) -> dict:
    groups = np.asarray(groups)
    uniq = np.array(sorted(set(groups.tolist())))
    idx_of = {g: np.flatnonzero(groups == g) for g in uniq}
    vals = [spearman(a[sel], b[sel]) for sel in
            (np.concatenate([idx_of[uniq[p]] for p in rng.integers(0, len(uniq), size=len(uniq))])
             for _ in range(n))]
    return dict(point=spearman(a, b), lo=float(np.nanpercentile(vals, 2.5)),
                hi=float(np.nanpercentile(vals, 97.5)), sd=float(np.nanstd(vals)),
                n=int((np.isfinite(a) & np.isfinite(b)).sum()))


# ---- stats295.py:57-64 ------------------------------------------------------------------------
def boot_rho_rows(a, b, rng, n: int = N_BOOT) -> dict:
    """plain row (gene) bootstrap."""
    m = np.isfinite(a) & np.isfinite(b)
    a, b = np.asarray(a)[m], np.asarray(b)[m]
    vals = [spearman(a[i], b[i]) for i in
            (rng.integers(0, a.size, size=a.size) for _ in range(n))]
    return dict(point=spearman(a, b), lo=float(np.nanpercentile(vals, 2.5)),
                hi=float(np.nanpercentile(vals, 97.5)), sd=float(np.nanstd(vals)), n=int(a.size))


def fmt(d) -> str:
    return f"{d['point']:+.4f} [{d['lo']:+.4f}, {d['hi']:+.4f}]"


def cell(sub: pd.DataFrame, xcol: str, rng, blocks: bool) -> dict:
    a = sub[xcol].to_numpy(float)
    b = sub["A"].to_numpy(float)
    gb = boot_rho_rows(a, b, rng)
    out = dict(rho=gb["point"], n=gb["n"],
               gene_lo=gb["lo"], gene_hi=gb["hi"], gene_sd=gb["sd"])
    if blocks and sub["file"].nunique() >= 2:
        fb = boot_rho(a, b, sub["file"].to_numpy(), rng)
        out.update(file_lo=fb["lo"], file_hi=fb["hi"], file_sd=fb["sd"],
                   n_files=int(sub["file"].nunique()))
    else:
        out.update(file_lo=np.nan, file_hi=np.nan, file_sd=np.nan,
                   n_files=int(sub["file"].nunique()))
    return out


def main() -> None:
    df = pd.read_csv(SRC)
    P = df[df["primary"]].reset_index(drop=True)
    files = sorted(P["file"].unique())

    # =========================== KA row: 295's pooled C5a / C5b ================================
    rng = np.random.default_rng(SEED)
    ka_a = cell(P, "C_held", rng, blocks=True)
    ka_b = cell(P, "gen_pc20", rng, blocks=True)
    print("=== KA (295 pooled, 174 primary genes, 8 files) ===", flush=True)
    print(f"  C5a rho {ka_a['rho']:+.6f}  (295 reported +0.5016)  "
          f"delta {abs(ka_a['rho'] - 0.5016):.2e}", flush=True)
    print(f"  C5b rho {ka_b['rho']:+.6f}  (295 reported +0.3709)  "
          f"delta {abs(ka_b['rho'] - 0.3709):.2e}", flush=True)
    ka = dict(n_primary=int(len(P)), n_files=len(files),
              C5a_rho=ka_a["rho"], C5a_ref=0.5016, C5a_delta=abs(ka_a["rho"] - 0.5016),
              C5b_rho=ka_b["rho"], C5b_ref=0.3709, C5b_delta=abs(ka_b["rho"] - 0.3709),
              C5a_file_block=[ka_a["file_lo"], ka_a["file_hi"]],
              C5b_file_block=[ka_b["file_lo"], ka_b["file_hi"]])
    assert ka["C5a_delta"] < 5e-4 and ka["C5b_delta"] < 5e-4, "KA row does not reproduce"

    # =========================== the jackknife ==================================================
    rows = []

    def add(label: str, sub: pd.DataFrame, kind: str, blocks: bool = True) -> None:
        r = np.random.default_rng(SEED)
        a = cell(sub, "C_held", r, blocks)
        r = np.random.default_rng(SEED)
        b = cell(sub, "gen_pc20", r, blocks)
        rec = dict(set=label, kind=kind, n_genes=int(len(sub)), n_files=a["n_files"],
                   C5a_rho=a["rho"], C5a_gene_lo=a["gene_lo"], C5a_gene_hi=a["gene_hi"],
                   C5a_file_lo=a["file_lo"], C5a_file_hi=a["file_hi"],
                   C5b_rho=b["rho"], C5b_gene_lo=b["gene_lo"], C5b_gene_hi=b["gene_hi"],
                   C5b_file_lo=b["file_lo"], C5b_file_hi=b["file_hi"])
        # 295's frozen rules: C5a rho >= 0.15 and lower bound > 0 ; C5b lower bound > 0
        rec["C5a_pass_gene"] = bool(a["rho"] >= 0.15 and a["gene_lo"] > 0)
        rec["C5b_pass_gene"] = bool(b["gene_lo"] > 0)
        rec["C5a_pass_file"] = (bool(a["rho"] >= 0.15 and a["file_lo"] > 0)
                                if np.isfinite(a["file_lo"]) else None)
        rec["C5b_pass_file"] = (bool(b["file_lo"] > 0)
                                if np.isfinite(b["file_lo"]) else None)
        rows.append(rec)

    add("ALL (8 files)", P, "pooled")
    for f in files:
        add(f"-{f}", P[P["file"] != f], "leave-one-file-out")
    add(f"ONLY {ADAMSON}", P[P["file"] == ADAMSON], "single file", blocks=False)
    add("ALL BUT Adamson", P[P["file"] != ADAMSON], "leave-one-file-out (Adamson)")

    out = pd.DataFrame(rows)
    out.to_csv(WORK / "jackknife298.csv", index=False)

    print("\n=== leave-one-file-out ===", flush=True)
    hdr = (f"  {'set':<46} {'n':>4} {'C5a rho':>9} {'gene 95%':>20} {'file 95%':>20} "
           f"{'C5b rho':>9} {'gene 95%':>20} {'file 95%':>20}")
    print(hdr, flush=True)
    for _, r in out.iterrows():
        def iv(lo, hi):
            return "        —          " if not np.isfinite(lo) else f"[{lo:+.4f}, {hi:+.4f}]"
        print(f"  {r['set']:<46} {r['n_genes']:>4d} {r['C5a_rho']:>+9.4f} "
              f"{iv(r['C5a_gene_lo'], r['C5a_gene_hi']):>20} "
              f"{iv(r['C5a_file_lo'], r['C5a_file_hi']):>20} "
              f"{r['C5b_rho']:>+9.4f} {iv(r['C5b_gene_lo'], r['C5b_gene_hi']):>20} "
              f"{iv(r['C5b_file_lo'], r['C5b_file_hi']):>20}", flush=True)

    # =========================== per-file composition ==========================================
    comp = []
    for f in files:
        sub = P[P["file"] == f]
        comp.append(dict(file=f, n_primary=int(len(sub)),
                         share=float(len(sub) / len(P)),
                         line=sub["line"].iloc[0], modality=sub["modality"].iloc[0],
                         median_A=float(sub["A"].median()),
                         median_C_held=float(sub["C_held"].median()),
                         median_r_split=float(sub["r_split"].median())))
    print("\n=== composition ===", flush=True)
    for c in comp:
        print(f"  {c['file']:<44} n {c['n_primary']:3d} ({c['share']:5.1%})  "
              f"{c['line']:<10} {c['modality']:<9} medA {c['median_A']:+.4f} "
              f"medC {c['median_C_held']:+.4f}", flush=True)

    ADA = out[out["set"] == "ALL BUT Adamson"].iloc[0]
    verdict = dict(
        C5a_rho=float(ADA["C5a_rho"]),
        C5a_gene=[float(ADA["C5a_gene_lo"]), float(ADA["C5a_gene_hi"])],
        C5a_file=[float(ADA["C5a_file_lo"]), float(ADA["C5a_file_hi"])],
        C5a_holds_gene=bool(ADA["C5a_pass_gene"]), C5a_holds_file=bool(ADA["C5a_pass_file"]),
        C5b_rho=float(ADA["C5b_rho"]),
        C5b_gene=[float(ADA["C5b_gene_lo"]), float(ADA["C5b_gene_hi"])],
        C5b_file=[float(ADA["C5b_file_lo"]), float(ADA["C5b_file_hi"])],
        C5b_holds_gene=bool(ADA["C5b_pass_gene"]), C5b_holds_file=bool(ADA["C5b_pass_file"]),
    )
    print("\n=== WITHOUT ADAMSON ===", flush=True)
    print(f"  C5a {verdict['C5a_rho']:+.4f} gene {verdict['C5a_gene']} file {verdict['C5a_file']} "
          f"-> gene {'HOLDS' if verdict['C5a_holds_gene'] else 'FAILS'} / "
          f"file {'HOLDS' if verdict['C5a_holds_file'] else 'FAILS'}", flush=True)
    print(f"  C5b {verdict['C5b_rho']:+.4f} gene {verdict['C5b_gene']} file {verdict['C5b_file']} "
          f"-> gene {'HOLDS' if verdict['C5b_holds_gene'] else 'FAILS'} / "
          f"file {'HOLDS' if verdict['C5b_holds_file'] else 'FAILS'}", flush=True)

    (WORK / "jackknife298.json").write_text(json.dumps(
        dict(ka=ka, composition=comp, jackknife=rows, without_adamson=verdict,
             seed=SEED, n_boot=N_BOOT, source=str(SRC)), indent=1, default=float))


if __name__ == "__main__":
    main()
