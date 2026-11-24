"""522 — 518's robust response STRENGTH for any panel (the "src" form's input), plus the per-target descriptors the
tail-origin analysis needs. Light (< 1 GB): only each screen's panel rows are read.

strength518.robust_strength with the panel as an argument (518 hard-codes A/B/C's panel):
    genes   = the screen's genes > 10 CPM in its control, minus THIS panel's genes;
    scale_g = 1.4826 MAD over the screen's covered panel targets;  z = delta / scale;
    noise_p = median_g |z_pg| / 0.6745;  strength_p = log1p #{g : |z_pg| > 4 noise_p}.
Screens: the three production sources (K562-GWPS, Orion HCT116, Orion HEK293T; data/donors) and CD4 T cells (the panel's
own cd4_compact extraction, 3 conditions, mean over those carrying the target: 518's cd4_robust).
Descriptors: the raw count and noise per screen; DepMap Chronos (17 lines, 360's subset) mean/min; membership of the K562
essential list (238's `k562_essential/perts.txt`, the list 238 excluded from efp/efr).

KNOWN ANSWER (abort): on A/B/C (data/challenge) the four arrays equal experiments/518-bulk-mse-gap/work/strength518.npz's
`k562_gwps`, `orion_hct116`, `orion_hek293t`, `cd4_robust` exactly (NaN positions included).

    uv run python experiments/522-october-bulk-key/strength522.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
import re  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
R233 = ROOT / "experiments/233-october-retention/work"
R238 = ROOT / "experiments/238-essential-free-panel/work"
PANELS = {
    "abc": dict(bundle=ROOT / "data/challenge", cd4=ROOT / "data/external/cd4t/cloud-compact/cd4_compact"),
    "efp": dict(bundle=R238 / "ef-bundle", cd4=R238 / "cd4/efp/cd4_compact"),
    "efr": dict(bundle=R238 / "ef-bundle-r", cd4=R238 / "cd4/efr/cd4_compact"),
    "real": dict(bundle=ROOT / "experiments/199-def-runbook/work/def-bundle-r", cd4=R233 / "cd4/real/cd4_compact"),
}
BANKS = {"k562_gwps": ROOT / "data/donors/k562_gwps", "orion_hct116": ROOT / "data/donors/orion_hct116",
         "orion_hek293t": ROOT / "data/donors/orion_hek293t"}
CONDS = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
DEP = ROOT / "experiments/360-depmap-three-lines/work/CRISPRGeneEffect_subset.csv"
ESS = ROOT / "data/donors/k562_essential/perts.txt"


def panel_of(bundle: Path):
    lines = [ln for ln in (bundle / "pert_counts.csv").read_text().splitlines() if ln.strip()]
    hdr = [c.strip() for c in lines[0].split(",")]
    j = hdr.index("target_gene")
    return [ln.split(",")[j].strip() for ln in lines[1:]]


def robust_strength(path: Path, targets, pset):
    """strength518.robust_strength with the panel (targets, its gene set) as arguments. Returns strength, count, noise."""
    perts = [ln.strip() for ln in (path / "perts.txt").read_text().splitlines() if ln.strip()]
    genes = [ln.strip() for ln in (path / "genes.txt").read_text().splitlines() if ln.strip()]
    b = np.load(path / "b_ctrl.npy").astype(np.float64)
    ex = np.expm1(b)
    cpm = 1e6 * ex / ex.sum()
    keep = (cpm > 10.0) & np.asarray([g not in pset for g in genes])
    pi = {p: i for i, p in enumerate(perts)}
    rows = [pi[t] for t in targets if t in pi]
    cov = [t for t in targets if t in pi]
    delta = np.load(path / "delta.npy", mmap_mode="r")
    R = np.asarray(delta[np.asarray(sorted(rows))], np.float64)
    order = np.argsort(np.argsort(rows))
    R = R[order][:, keep]
    med = np.median(R, 0)
    scale = 1.4826 * np.median(np.abs(R - med), 0)
    ok = scale > 0
    Z = R[:, ok] / scale[ok]
    noise = np.median(np.abs(Z), 1) / 0.6745
    cnt = (np.abs(Z) > 4.0 * noise[:, None]).sum(1)
    so, co, no = dict(zip(cov, np.log1p(cnt).tolist())), dict(zip(cov, cnt.tolist())), dict(zip(cov, noise.tolist()))
    f = lambda d: np.asarray([d.get(t, np.nan) for t in targets], float)
    return f(so), f(co), f(no)


def depmap(targets):
    import csv
    with open(DEP) as fh:
        r = csv.reader(fh)
        hdr = next(r)
        rows = [row for row in r]
    syms = [re.sub(r"\s*\(\d+\)$", "", c) for c in hdr[1:]]
    col = {}
    for j, s in enumerate(syms):
        col.setdefault(s, j)
    M = np.asarray([[float(x) if x != "" else np.nan for x in row[1:]] for row in rows])
    mean, mn = np.full(len(targets), np.nan), np.full(len(targets), np.nan)
    for i, t in enumerate(targets):
        if t in col:
            v = M[:, col[t]]
            if np.isfinite(v).any():
                mean[i], mn[i] = np.nanmean(v), np.nanmin(v)
    return mean, mn, len(rows)


def main():
    ess = set(ln.strip() for ln in ESS.read_text().splitlines() if ln.strip())
    out, summ = {}, {}
    for name, spec in PANELS.items():
        targets = panel_of(spec["bundle"])
        assert len(targets) == 300 and len(set(targets)) == 300
        pset = set(targets)
        out[f"{name}__targets"] = np.asarray(targets)
        for bk, path in BANKS.items():
            s, c, n = robust_strength(path, targets, pset)
            out[f"{name}__{bk}"], out[f"{name}__{bk}__count"], out[f"{name}__{bk}__noise"] = s, c, n
        vals, cnts = [], []
        for cnd in CONDS:
            s, c, _ = robust_strength(spec["cd4"] / cnd, targets, pset)
            vals.append(s)
            cnts.append(c)
        with np.errstate(all="ignore"):
            out[f"{name}__cd4_robust"] = np.nanmean(np.vstack(vals), 0)
        dm, dmin, nlines = depmap(targets)
        out[f"{name}__dep_mean"], out[f"{name}__dep_min"] = dm, dmin
        out[f"{name}__k562_essential"] = np.asarray([t in ess for t in targets])
        src = np.vstack([out[f"{name}__{bk}"] for bk in BANKS])
        with np.errstate(all="ignore"):
            smean = np.nanmean(src, 0)
        summ[name] = dict(
            n=len(targets), coverage={bk: int(np.isfinite(out[f"{name}__{bk}"]).sum()) for bk in BANKS},
            cd4_coverage=int(np.isfinite(out[f"{name}__cd4_robust"]).sum()),
            src_mean_strength_quantiles=np.nanquantile(smean, [.1, .25, .5, .75, .9]).tolist(),
            frac_src_count_le_10=float(np.nanmean(np.nanmean(np.vstack([out[f"{name}__{bk}__count"] for bk in BANKS]), 0) <= 10)),
            k562_essential=int(out[f"{name}__k562_essential"].sum()),
            depmap_covered=int(np.isfinite(dm).sum()), depmap_lines=nlines,
            frac_dep_mean_lt_m0_5=float(np.nanmean(dm < -0.5)), frac_dep_mean_lt_m1=float(np.nanmean(dm < -1.0)))
        print(f"{name}: {json.dumps(summ[name])}", flush=True)
    # ---- KNOWN ANSWER: A/B/C == 518's strength518.npz
    ref = np.load(ROOT / "experiments/518-bulk-mse-gap/work/strength518.npz")
    assert [str(t) for t in ref["panel"]] == [str(t) for t in out["abc__targets"]], "KA panel order"
    for k in ("k562_gwps", "orion_hct116", "orion_hek293t", "cd4_robust"):
        a, b = out[f"abc__{k}"], np.asarray(ref[k], float)
        same = np.array_equal(np.isnan(a), np.isnan(b)) and np.array_equal(a[~np.isnan(a)], b[~np.isnan(b)])
        print(f"KA abc {k:14s} == strength518.npz: {'PASS' if same else 'FAIL'}", flush=True)
        assert same, f"KA FAILED {k}"
    summ["ka"] = "PASS: abc k562_gwps / orion_hct116 / orion_hek293t / cd4_robust == strength518.npz (bit-identical)"
    np.savez_compressed(WORK / "strength522.npz", **out)
    (WORK / "strength522.json").write_text(json.dumps(summ, indent=1, default=float))


if __name__ == "__main__":
    main()
