#!/usr/bin/env python3
"""
294 — figures for PAPER.md (revised in experiment 533, 2026-09-29/30: F-12 and F-13 added; F-9's caption and
F-10's source line corrected).

Every panel re-plots numbers that a finished experiment already wrote to disk.
NOTHING here recomputes a statistic from raw data: the script reads saved CSV /
JSON / Markdown tables only, and the only arithmetic it does is the arithmetic
already reported in the source (e.g. reading a Spearman-Brown value that
`analyse285.json` stores, not computing one).

Run:  uv run python make_figs.py
Out:  figs/F-1.{png,svg} ... figs/F-13.{png,svg}, figs/MANIFEST.md
"""

from __future__ import annotations

import csv
import json
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# 533's copy lives one level deeper (experiments/533-science-lead/paper/) than 294's (experiments/294-writeup/).
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
assert os.path.isdir(os.path.join(REPO, "data", "board-entries")), REPO
HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "figs")
os.makedirs(FIGS, exist_ok=True)

E285 = os.path.join(REPO, "experiments", "285-reliability-benchmark", "work")
E284 = os.path.join(REPO, "experiments", "284-conservation-family-split", "work")
E290 = os.path.join(REPO, "experiments", "290-remainder-response-features", "work")
E291 = os.path.join(REPO, "experiments", "291-f13-independent-surfaces", "work")
E276 = os.path.join(REPO, "experiments", "276-conservation-map", "work")
E272 = os.path.join(REPO, "experiments", "272-noise-corrected-operator", "work")
E286 = os.path.join(REPO, "experiments", "286-crossproduct-estimator", "work")
E289 = os.path.join(REPO, "experiments", "289-broad16-state-features", "work")
E273 = os.path.join(REPO, "experiments", "273-flex-reference-instrument", "work")

plt.rcParams.update(
    {
        "figure.dpi": 140,
        "savefig.dpi": 200,
        "font.size": 8.5,
        "axes.titlesize": 9.0,
        "axes.labelsize": 8.5,
        "legend.fontsize": 7.2,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.5,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)

CB = {  # colour-blind-safe
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "red": "#D55E00",
    "purple": "#CC79A7",
    "sky": "#56B4E9",
    "yellow": "#F0E442",
    "grey": "#7F7F7F",
    "black": "#222222",
}

MANIFEST: list[tuple[str, str, str]] = []


def save(fig, name: str, caption: str, sources: str) -> None:
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(FIGS, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    MANIFEST.append((name, caption, sources))
    print(f"  wrote figs/{name}.png + .svg")


def read_json(path: str):
    with open(path) as fh:
        return json.load(fh)


def read_csv(path: str) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh))


def f(x, default=np.nan):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


# ----------------------------------------------------------------------------
# F-1  the cell ladder against Spearman-Brown, and cells vs UMIs
# ----------------------------------------------------------------------------
def fig1() -> None:
    bench = read_csv(os.path.join(E285, "bench285.csv"))
    an = read_json(os.path.join(E285, "analyse285.json"))

    ladders = {
        "h1-test (H1 hESC, Flex)": ("h1-test", CB["blue"], "h1-test"),
        "K562-GWPS (panel targets)": ("k562-gwps", CB["orange"], "K562-GWPS panel"),
        "Nadig Jurkat": ("nadig-jurkat", CB["green"], "Nadig Jurkat"),
        "Nadig HepG2": ("nadig-hepg2", CB["red"], "Nadig HepG2"),
    }

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.7))

    # (a) measured cell ladders + the Spearman-Brown extrapolation from n = 25
    ax = axes[0]
    for src, (key, colour, label) in ladders.items():
        rows = [
            r
            for r in bench
            if r["source"] == src
            and r["unit"].startswith("a cells")
            and "one construct" not in r["unit"]
            and (r["note"] in ("", "fixed cohort"))
        ]
        rows = sorted({int(r["n_per_side"]): r for r in rows}.items())
        if not rows:
            continue
        ns = [n for n, _ in rows]
        med = [f(r["median_cos_resid"]) for _, r in rows]
        lo = [f(r["lo"]) for _, r in rows]
        hi = [f(r["hi"]) for _, r in rows]
        ax.errorbar(
            ns,
            med,
            yerr=[np.array(med) - np.array(lo), np.array(hi) - np.array(med)],
            marker="o",
            ms=4,
            lw=1.5,
            capsize=2,
            color=colour,
            label=label,
        )
        dl = an["depth_law"].get(key, {}).get("from_smallest", {})
        if dl:
            sb_n = sorted(int(k) for k in dl)
            sb_v = [dl[str(n)]["sb_median"] for n in sb_n]
            ax.plot(
                [ns[0]] + sb_n,
                [med[0]] + sb_v,
                ls="--",
                lw=1.1,
                color=colour,
                alpha=0.75,
            )
    ax.plot([], [], ls="--", lw=1.1, color=CB["grey"], label="Spearman–Brown from n = 25")
    ax.set_xscale("log")
    ax.set_xticks([25, 50, 100, 200])
    ax.set_xticklabels(["25", "50", "100", "200"])
    ax.minorticks_off()
    ax.set_xlabel("cells per side (same construct, random cell halves)")
    ax.set_ylabel("median split-half residual cosine")
    ax.set_title("(a) the cell ladder is a Spearman–Brown ladder")
    ax.legend(frameon=False, loc="upper left")

    # (b) cells and UMI depth are one axis: thinned ladders on the cell ladder
    ax = axes[1]
    for src, colour, label, ns_keep in [
        ("h1-test (H1 hESC, Flex)", CB["blue"], "h1-test", (25, 50, 100, 200)),
        ("K562-GWPS (panel targets)", CB["orange"], "K562-GWPS panel", (25, 50, 100)),
    ]:
        cells = sorted(
            {
                int(r["n_per_side"]): f(r["median_cos_resid"])
                for r in bench
                if r["source"] == src
                and r["unit"].startswith("a cells")
                and "one construct" not in r["unit"]
                and r["note"] in ("", "fixed cohort")
                and int(r["n_per_side"]) in ns_keep
            }.items()
        )
        ax.plot(
            [n for n, _ in cells],
            [v for _, v in cells],
            marker="o",
            ms=4.5,
            lw=1.6,
            color=colour,
            label=f"{label}: full depth, n cells",
            zorder=3,
        )
        thin = [r for r in bench if r["source"] == src and r["unit"].startswith("a' thinned")]
        xs, ys = [], []
        for r in thin:
            frac = float(re.search(r"f=([0-9.]+)", r["unit"]).group(1))
            xs.append(int(r["n_per_side"]) * frac)  # equal-total-UMI abscissa
            ys.append(f(r["median_cos_resid"]))
        ax.scatter(xs, ys, s=18, facecolor="none", edgecolor=colour, lw=1.0, zorder=4)
    ax.scatter([], [], s=18, facecolor="none", edgecolor=CB["grey"], label="counts thinned, plotted at n × f")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("cells per side × depth fraction  (equal total UMIs per side)")
    ax.set_ylabel("median split-half residual cosine")
    ax.set_title("(b) cells and UMI depth are one axis")
    ax.legend(frameon=False, loc="upper left")

    save(
        fig,
        "F-1",
        "Reference reliability is one removable axis. (a) Median split-half residual cosine of a "
        "perturbation response against cells per side, same guide construct, random cell halves; "
        "unit = targets (n = 85 h1-test, 95 K562-GWPS panel, 336 Nadig Jurkat, 176 Nadig HepG2); "
        "bars are 95 % paired target bootstraps from `bench285.csv`; the dashed line is the "
        "Spearman–Brown prediction from the n = 25 rung stored in `analyse285.json` "
        "(`depth_law.*.from_smallest.sb_median`). (b) The same medians with counts thinned at a "
        "fixed cell number, plotted at equal effective total UMIs: the thinned points fall on the "
        "cell ladders, i.e. cells and depth trade one for one.",
        "experiments/285-reliability-benchmark/work/bench285.csv; .../analyse285.json",
    )


# ----------------------------------------------------------------------------
# F-2  the unit decomposition: what more measurement cannot remove
# ----------------------------------------------------------------------------
def fig2() -> None:
    """Rows are drawn from the two files 285 wrote: `analyse285.json` (`rho`,
    which carries a bootstrap where 285 computed one) and `equiv285.json`
    (`frac`, the same ratio on the best-powered two-construct cohorts).  Each
    row below names its file and key so the plot is traceable line by line."""
    rho = read_json(os.path.join(E285, "analyse285.json"))["rho"]
    eq = read_json(os.path.join(E285, "equiv285.json"))

    def R(key):  # analyse285.json rho block
        v = rho[key]
        lo, hi = v.get("lo"), v.get("hi")
        if lo is None or hi is None or not np.isfinite(lo) or abs(hi - lo) > 2.0:
            lo = hi = None
        return v["rho"], lo, hi, v.get("n")

    def Q(key, n):  # equiv285.json frac (no interval stored)
        return eq[key]["frac"], None, None, n

    groups = [
        (
            "guide construct",
            CB["red"],
            [
                ("K562-GWPS, 748 two-construct (CRISPRi)", Q("K562-GWPS, 748 two-construct perturbations @ 25", 748)),
                ("Nadig HepG2 (CRISPRi)", R("Nadig HepG2  b GUIDE construct / same-construct cells")),
                ("h1-test H1 hESC (CRISPRi)", R("h1-test  b GUIDE construct / same-construct cells")),
                ("Nadig Jurkat (CRISPRi)", R("Nadig Jurkat  b GUIDE construct / same-construct cells")),
                ("CD4 Stim48hr (CRISPRi)", R("CD4 Stim48hr  b GUIDE construct / cells (count split − 0.044) (T=1e+05)")),
                ("Broad-16, 12 lines (Cas9 KO)", R("Broad-16  b GUIDE construct / cells (T=2e+05)")),
            ],
        ),
        (
            "biological donor",
            CB["purple"],
            [
                ("CD4 Rest", R("CD4 Rest  c different DONOR / cells (count split − 0.044) (T=1e+05)")),
                ("CD4 Stim8hr", R("CD4 Stim8hr  c different DONOR / cells (count split − 0.044) (T=1e+05)")),
                ("CD4 Stim48hr", R("CD4 Stim48hr  c different DONOR / cells (count split − 0.044) (T=1e+05)")),
            ],
        ),
        (
            "cell line",
            CB["green"],
            [
                ("Broad-16, 66 line pairs (T = 5e5)", R("Broad-16  f different LINE / cells (T=5e+05)")),
                ("Broad-16, 66 line pairs (T = 2e5)", R("Broad-16  f different LINE / cells (T=2e+05)")),
            ],
        ),
        (
            "activation state (same donor, same construct)",
            CB["orange"],
            [
                ("CD4 Rest", R("CD4 Rest  s different STATE / cells (count split − 0.044) (T=1e+05)")),
                ("CD4 Stim48hr", R("CD4 Stim48hr  s different STATE / cells (count split − 0.044) (T=1e+05)")),
            ],
        ),
        (
            "different study, same line",
            CB["sky"],
            [("K562-GWPS vs the K562 essential screen", R("K562 GWPS-vs-essential  e different STUDY / within-study cells"))],
        ),
        (
            "batch / lane / gemgroup",
            CB["blue"],
            [
                ("267 K562-GWPS batch", R("267 k562_gwps  d BATCH / random cells")),
                ("Nadig Jurkat gemgroup", R("Nadig Jurkat  d GEMGROUP / random cells")),
                ("h1-test lane", R("h1-test  d LANE / random cells")),
                ("Nadig HepG2 gemgroup", R("Nadig HepG2  d GEMGROUP / random cells")),
                ("267 Orion HCT116 batch", R("267 orion_hct116  d BATCH / random cells")),
                ("K562-GWPS gemgroup", R("K562-GWPS  d GEMGROUP / random cells")),
                ("267 Orion HEK293T batch  (the one exception)", R("267 orion_hek293t  d BATCH / random cells")),
            ],
        ),
    ]

    labels, pts, colours, headers = [], [], [], []
    for gname, colour, items in groups:
        headers.append((gname, colour, len(labels)))
        for lab, val in items:
            labels.append(lab)
            pts.append(val)
            colours.append(colour)

    y = np.arange(len(labels))[::-1].astype(float)
    # push each group down a little so the header has room
    shift = np.zeros_like(y)
    for gi, (_, _, start) in enumerate(headers):
        shift[start:] -= 0.85
    y = y + shift

    fig, ax = plt.subplots(figsize=(8.6, 7.4))
    for yi, (v, lo, hi, n), c in zip(y, pts, colours):
        if lo is not None:
            ax.plot([lo, hi], [yi, yi], color=c, lw=1.3, alpha=0.55, zorder=2)
        ax.scatter([v], [yi], s=32, color=c, zorder=3)
        ax.text(1.24, yi, f"{v:+.2f}", va="center", ha="right", fontsize=7.0, color="#333333")
        ax.text(1.44, yi, f"n {n}", va="center", ha="right", fontsize=6.6, color=CB["grey"])
    ax.axvline(0.0, color=CB["black"], lw=0.9)
    ax.axvline(1.0, color=CB["grey"], lw=0.9, ls=":")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlim(-0.15, 1.45)
    ax.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel(r"$\rho$ = (agreement across the unit) / (agreement across cells of one construct, matched budget)")
    ax.set_title("What more measurement of the same construct cannot remove")
    for gname, colour, start in headers:
        ax.text(-0.14, y[start] + 0.62, gname.upper(), fontsize=7.4, color=colour, fontweight="bold", va="bottom")
    ax.text(0.0, y[-1] - 1.5, "0 = irreducible", fontsize=7.2, color=CB["black"], ha="center")
    ax.text(1.0, y[-1] - 1.5, "1.0 = free (more cells remove it)", fontsize=7.2, color=CB["grey"], ha="center")
    ax.set_ylim(y[-1] - 2.1, y[0] + 1.7)
    ax.grid(axis="y", alpha=0.0)

    save(
        fig,
        "F-2",
        "The unit decomposition. Each row is one dataset and one experimental unit; ρ is that unit's "
        "agreement divided by the agreement of two halves of ONE construct's cells at the same budget — "
        "1.0 means the unit costs nothing that more cells cannot buy back, 0 means the disagreement is "
        "irreducible. n is the number of targets (or target × line pairs) behind the ratio. Bars are the "
        "95 % paired target bootstrap where `analyse285.json` (`rho`) stores one; the K562 guide row is "
        "the 748-perturbation cohort from `equiv285.json` (`frac`), which carries no stored interval but "
        "is the best-powered construct contrast on disk (the tight statistic there is the CELLS − GUIDE "
        "difference, +0.0058 [+0.0045, +0.0069]). The guide construct, not the batch, is the dominant "
        "irreducible term on CRISPRi libraries; the Cas9-knockout row is a modality hypothesis, not a "
        "controlled contrast.",
        "experiments/285-reliability-benchmark/work/analyse285.json (`rho`); .../equiv285.json (`frac`)",
    )


# ----------------------------------------------------------------------------
# F-3  the conservation map: reproducibility and off-map transfer
# ----------------------------------------------------------------------------
def fig3() -> None:
    an276 = read_json(os.path.join(E276, "analyse276.json"))
    rel276 = read_json(os.path.join(E276, "rel276.json"))
    pred276 = read_json(os.path.join(E276, "predict276.json"))
    sa = read_json(os.path.join(E291, "surfA291.json"))
    sb = read_json(os.path.join(E291, "surfB291.json"))

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2))

    # (a) reproducibility
    ax = axes[0]
    rnd = an276["disjoint_line_split"]["random"]["values"]
    ax.hist(rnd, bins=14, color=CB["sky"], alpha=0.55, edgecolor="white", label="50 random line splits")
    dec = an276["disjoint_line_split"]["declared"]
    ax.axvline(dec["rho"], color=CB["blue"], lw=2.0, label=f"declared disjoint line sets  {dec['rho']:+.3f}")
    ax.axvspan(dec["ci"][0], dec["ci"][1], color=CB["blue"], alpha=0.12)
    st = rel276["strict"]
    ax.axvline(st["rho"], color=CB["red"], lw=2.0, label=f"independent cells, 4 sources  {st['rho']:+.3f}")
    ax.axvspan(st["ci"][0], st["ci"][1], color=CB["red"], alpha=0.12)
    ax.axvline(0.0, color=CB["black"], lw=0.8)
    ax.set_xlabel(r"Spearman $\rho$ between two independent estimates of a target's conservation")
    ax.set_ylabel("random splits")
    ax.set_title("(a) the map reproduces (n = 289 / 295 targets)")
    ax.legend(frameon=False, loc="upper left")

    # (b) prediction and off-map transfer
    ax = axes[1]
    rows = []
    tbl = {r["key"]: r for r in pred276["table"]}
    lolo = tbl["c_lolo"]["rho|pooled3"]
    rows.append(("PLOLO per-target `pds`, held line out\n(3 cancer folds, n 853)", lolo["rho"], lolo["ci"], CB["blue"]))
    h1 = pred276["h1_test"]["C_lolo (H1 sources removed)"]["rho_pds"]
    rows.append(("Flex-400 per-target `pds`, H1 sources out\n(85 targets)", h1["rho"], h1["ci"], CB["blue"]))
    ch = sa["primary"]["C_held"]
    rows.append(("Broad-16 12-line Cas9 conservation\n(82 targets, off-map)", ch["point"], [ch["lo"], ch["hi"]], CB["green"]))
    cb4 = sb["conditions"]["Stim48hr"]["boot"]["C_held"]
    rows.append(
        (
            "GWCD4i authors' CD4-vs-K562 statistic\n(1,773 genes, off-map)",
            cb4["point"],
            [cb4["lo"], cb4["hi"]],
            CB["green"],
        )
    )
    dd = tbl["dd"]["rho|pooled3"]
    rows.append(("baseline: the two-donor agreement the\npipeline already has (PLOLO, n 853)", dd["rho"], dd["ci"], CB["grey"]))

    y = np.arange(len(rows))[::-1]
    for yi, (lab, v, ci, c) in zip(y, rows):
        ax.plot([ci[0], ci[1]], [yi, yi], color=c, lw=2.0, alpha=0.6)
        ax.scatter([v], [yi], s=40, color=c, zorder=3)
        ax.text(v, yi + 0.22, f"{v:+.3f}", ha="center", fontsize=7.2, color="#333333")
    ax.axvline(0.0, color=CB["black"], lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows])
    ax.set_xlabel(r"Spearman $\rho$ with the outcome")
    ax.set_title("(b) the score predicts transfer, on and off the map")
    ax.set_ylim(y[-1] - 0.7, y[0] + 0.7)

    save(
        fig,
        "F-3",
        "Cross-line conservation is a reproducible target property that predicts transfer. "
        "(a) Unit = targets. Two independent estimates of the same target's conservation agree at "
        "ρ +0.435 [+0.340, +0.522] across the declared disjoint line sets (n = 289; `analyse276.json`) "
        "and +0.273 [+0.158, +0.388] across independent cells of four sources (n = 295; `rel276.json`); "
        "the histogram is 50 random line splits, all positive. (b) Unit = targets (a, b, c) or genes (d). "
        "The leave-the-held-line-out score predicts the incumbent's per-target `pds` on PLOLO "
        "(`predict276.json`) and on the Flex-400 reference, and the map's own per-gene score transfers "
        "to two surfaces outside the 28 banks (`surfA291.json`, `surfB291.json`). Intervals: 95 % paired "
        "target bootstrap (a, b PLOLO/Flex/Broad-16) and component-block bootstrap (CD4-vs-K562).",
        "experiments/276-conservation-map/work/{analyse276,rel276,predict276}.json; "
        "experiments/291-f13-independent-surfaces/work/{surfA291,surfB291}.json",
    )


# ----------------------------------------------------------------------------
# F-4  the family axis under strict folds; the remainder's cross-half reproducibility
# ----------------------------------------------------------------------------
def fig4() -> None:
    fam = read_json(os.path.join(E284, "fam284.json"))
    rem = read_json(os.path.join(E284, "rem284.json"))

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))

    ax = axes[0]
    tab = fam["family_axis_all"]["table"]
    obs = np.array([t["observed"] for t in tab])
    prd = np.array([t["predicted"] for t in tab])
    ns = np.array([t["n"] for t in tab])
    ax.scatter(prd, obs, s=np.clip(ns, 8, 160) * 0.55, color=CB["blue"], alpha=0.55, edgecolor="white", lw=0.5)
    for t in tab:
        if t["family"] in (
            "Zinc fingers C2H2-type",
            "Proteasome",
            "Mediator complex",
            "Homeoboxes",
            "Spliceosomal C complex",
        ):
            ax.annotate(
                t["family"],
                (t["predicted"], t["observed"]),
                textcoords="offset points",
                xytext=(6, 3),
                fontsize=6.4,
                color="#444444",
            )
    lim = [min(prd.min(), obs.min()) - 0.01, max(prd.max(), obs.max()) + 0.01]
    ax.plot(lim, lim, ls=":", color=CB["grey"], lw=1.0)
    ax.axhline(0, color=CB["black"], lw=0.7)
    ax.axvline(0, color=CB["black"], lw=0.7)
    r = fam["family_axis_all"]
    ax.set_title(
        f"(a) family axis, strict folds: ρ {r['rho']:+.3f} [{r['lo']:+.3f}, {r['hi']:+.3f}], {r['n_families']} families"
    )
    ax.set_xlabel("out-of-fold PREDICTED family median conservation")
    ax.set_ylabel("OBSERVED family median conservation")

    ax = axes[1]
    parts = rem["partitions"]
    ceil = np.array([p["ceiling"] for p in parts])
    remv = np.array([p["rho_remainder"] for p in parts])
    ax.scatter(ceil, remv, s=34, color=CB["orange"], alpha=0.8, edgecolor="white", lw=0.5, label="20 admissible partitions")
    pr = rem["primary"]
    ax.scatter([pr["ceiling"]], [pr["rho_remainder"]], s=95, marker="*", color=CB["red"], zorder=4, label="primary partition #1")
    ax.errorbar(
        [pr["ceiling"]],
        [pr["rho_remainder"]],
        xerr=[[pr["ceiling"] - pr["boot_ceiling"]["lo"]], [pr["boot_ceiling"]["hi"] - pr["ceiling"]]],
        yerr=[
            [pr["rho_remainder"] - pr["boot_remainder"]["lo"]],
            [pr["boot_remainder"]["hi"] - pr["rho_remainder"]],
        ],
        fmt="none",
        ecolor=CB["red"],
        lw=1.3,
        capsize=2,
    )
    nl = rem["primary"]["null"]
    ax.axhspan(nl["mean"] - 2 * nl["sd"], nl["mean"] + 2 * nl["sd"], color=CB["grey"], alpha=0.22)
    ax.text(0.02, nl["mean"] + 2 * nl["sd"] + 0.012, "within-component permutation null (0/200 reach the observed)", fontsize=6.6, color="#555555")
    lim = [0.0, 0.60]
    ax.plot(lim, lim, ls=":", color=CB["grey"], lw=1.0)
    ax.set_xlim(*lim)
    ax.set_ylim(-0.02, 0.60)
    ax.set_xlabel("ceiling: ρ between the two halves' RAW maps")
    ax.set_ylabel("ρ between the two halves' gene-level REMAINDERS")
    ax.set_title("(b) the per-gene remainder reproduces across disjoint source sets")
    ax.legend(frameon=False, loc="upper left")

    save(
        fig,
        "F-4",
        "The family axis survives strict folds and the per-gene remainder is real. (a) Unit = HGNC "
        "families; out-of-fold predicted vs observed family median conservation under folds in which a "
        "test gene shares NO family with any training row; n = 78 families, ρ +0.805 [+0.679, +0.883] "
        "(component block bootstrap; `fam284.json`), marker area ∝ genes in the family. (b) Unit = genes. "
        "For each of the 20 admissible partitions of the 28 sources into two halves sharing no source and "
        "no cell line, the cross-half agreement of the raw maps (the ceiling) against the cross-half "
        "agreement of the gene-level remainder after the 21 annotation features' out-of-fold prediction "
        "is removed; the star is the primary partition (n = 2,848 genes, remainder +0.328 [+0.241, +0.373] "
        "against a ceiling of +0.382), bars are its component bootstraps, and the grey band is the "
        "within-component permutation null (`rem284.json`).",
        "experiments/284-conservation-family-split/work/{fam284,rem284}.json",
    )


# ----------------------------------------------------------------------------
# F-5  response geometry: the sign inside and outside the map
# ----------------------------------------------------------------------------
def fig5() -> None:
    fit = read_json(os.path.join(E290, "fit290.json"))
    sa = read_json(os.path.join(E291, "surfA291.json"))
    sb = read_json(os.path.join(E291, "surfB291.json"))

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.0))

    # (a) inside the map: per-feature associations with the held half's remainder
    ax = axes[0]
    pf = fit["per_feature"]
    picks = [
        ("gen_pc20  (energy in the bank's top-20 PCs)", "gen_pc20"),
        ("gen_pc5   (top-5 PCs)", "gen_pc5"),
        ("pe_prolif_abs  (proliferation engagement)", "pe_prolif_abs"),
        ("breadth", "breadth"),
        ("gen_abscos  (|cos| to the bank mean row)", "gen_abscos"),
        ("logmag  (response magnitude)", "logmag"),
        ("gen_cos  (signed cos to the bank mean row)", "gen_cos"),
        ("ontarget  (on-target knockdown, less negative = weaker)", "ontarget"),
    ]
    labs = [p[0] for p in picks if p[1] in pf]
    vals = [pf[p[1]]["point"] for p in picks if p[1] in pf]
    los = [pf[p[1]]["lo"] for p in picks if p[1] in pf]
    his = [pf[p[1]]["hi"] for p in picks if p[1] in pf]
    y = np.arange(len(labs))[::-1]
    cols = [CB["green"] if v > 0 else CB["red"] for v in vals]
    for yi, v, lo, hi, c in zip(y, vals, los, his, cols):
        ax.plot([lo, hi], [yi, yi], color=c, lw=1.8, alpha=0.6)
        ax.scatter([v], [yi], s=36, color=c, zorder=3)
    ax.axvline(0, color=CB["black"], lw=0.9)
    # mark the two features whose preregistered sign was NEGATIVE
    prereg_neg = {"gen_pc5   (top-5 PCs)", "pe_prolif_abs  (proliferation engagement)"}
    for yi, lab, v in zip(y, labs, vals):
        if lab in prereg_neg:
            ax.text(v + 0.017, yi + 0.13, "*", fontsize=13, color=CB["red"], ha="left", va="center")
    ax.text(-0.31, y[-1] - 0.62, "* preregistered sign was NEGATIVE", fontsize=6.9, color=CB["red"])
    ax.set_yticks(y)
    ax.set_yticklabels(labs)
    ax.set_ylim(y[-1] - 0.95, y[0] + 0.6)
    ax.set_xlabel(r"Spearman $\rho$ with the OTHER half's conservation remainder (n = 2,848 genes)")
    ax.set_title("(a) inside the map: the preregistered mechanism inverted")
    ax.set_xlim(-0.33, 0.30)

    # (b) outside the map
    ax = axes[1]
    A = sa["primary"]
    B = sb["conditions"]["Stim48hr"]["boot"]
    items = [
        ("gen_pc20", A["gen_pc20"], B["gen_pc20"]),
        ("breadth", A["breadth"], B["breadth"]),
        ("logmag", A["logmag"], B["logmag"]),
        ("M_resp (geometry ridge)", A["M_resp"], B["M_resp"]),
        ("M_ann (annotation ridge)", A["M_ann"], B["M_ann"]),
        ("C_held (the map's own score)", A["C_held"], B["C_held"]),
    ]
    y = np.arange(len(items))[::-1]
    off = 0.17
    for yi, (lab, a, b) in zip(y, items):
        ax.plot([a["lo"], a["hi"]], [yi + off] * 2, color=CB["green"], lw=1.7, alpha=0.6)
        ax.scatter([a["point"]], [yi + off], s=32, color=CB["green"], zorder=3)
        ax.plot([b["lo"], b["hi"]], [yi - off] * 2, color=CB["purple"], lw=1.7, alpha=0.6)
        ax.scatter([b["point"]], [yi - off], s=32, color=CB["purple"], zorder=3)
    ax.scatter([], [], s=32, color=CB["green"], label="Broad-16, 82 Cas9 targets × 12 lines (n = 82)")
    ax.scatter([], [], s=32, color=CB["purple"], label="GWCD4i CD4-vs-K562, Stim48hr (n = 1,773 genes)")
    ax.axvline(0, color=CB["black"], lw=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([i[0] for i in items])
    ax.set_xlabel(r"Spearman $\rho$ with the surface's own conservation statistic")
    ax.set_title("(b) outside the map: the sign holds, the ranking does not")
    ax.legend(frameon=True, framealpha=0.92, edgecolor="none", loc="upper left")
    ax.set_ylim(y[-1] - 0.6, y[0] + 1.15)

    save(
        fig,
        "F-5",
        "Response geometry predicts conservation, with an inverted sign, and the sign holds off the map. "
        "(a) Unit = genes (n = 2,848 with ≥ 1 cross-line pair in both halves). Spearman of each "
        "response-derived feature measured on one half of the 28 sources against the OTHER half's "
        "conservation remainder; intervals are component block bootstraps (`fit290.json`, `per_feature`). "
        "`gen_pc20`, `gen_pc5` and proliferation engagement are POSITIVE — the preregistered prediction "
        "was negative for the latter two. (b) The same features and two ridge models applied to two "
        "surfaces outside the 28 banks: unit = 82 Broad-16 targets (paired target bootstrap) and 1,773 "
        "GWCD4i genes (component block bootstrap). `gen_pc20` clears zero on both; the geometry ridge "
        "does not beat the annotation ridge on either, and the map's own score `C_held` matches any "
        "model fitted on it.",
        "experiments/290-remainder-response-features/work/fit290.json; "
        "experiments/291-f13-independent-surfaces/work/{surfA291,surfB291}.json",
    )


# ----------------------------------------------------------------------------
# F-6  the whitener's noise directions and the retrieval loss
# ----------------------------------------------------------------------------
def fig6() -> None:
    a272 = read_json(os.path.join(E272, "analysis272.json"))
    a286 = read_json(os.path.join(E286, "analysisB286.json"))

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.0))

    ax = axes[0]
    ops = [
        ("held K562", "k562_gwps"),
        ("held HCT116", "orion_hct116"),
        ("held HEK293T", "orion_hek293t"),
        ("production", "production"),
    ]
    x = np.arange(len(ops))
    w = 0.26
    d_share = [a272["spectra"][k]["share_s_zero"] * 100 for _, k in ops]
    e_nd = [a272["spectra"][k]["whitened_energy_share_noise_dominated"] * 100 for _, k in ops]
    e_sz = [a272["spectra"][k]["whitened_energy_share_s_zero"] * 100 for _, k in ops]
    ax.bar(x - w, d_share, w, color=CB["grey"], label="directions with noise ≥ fitted eigenvalue")
    ax.bar(x, e_nd, w, color=CB["blue"], label="whitened prediction energy in the noise-dominated set")
    ax.bar(x + w, e_sz, w, color=CB["sky"], label="… of it, in the noise ≥ eigenvalue set")
    for xi, v in zip(x - w, d_share):
        ax.text(xi, v + 1.2, f"{v:.1f}", ha="center", fontsize=6.6)
    for xi, v in zip(x, e_nd):
        ax.text(xi, v + 1.2, f"{v:.1f}", ha="center", fontsize=6.6)
    for xi, v in zip(x + w, e_sz):
        ax.text(xi, v + 1.2, f"{v:.1f}", ha="center", fontsize=6.6)
    ax.set_xticks(x)
    ax.set_xticklabels([o[0] for o in ops], fontsize=7.6)
    ax.set_ylabel("per cent")
    ax.set_ylim(0, 128)
    ax.set_title("(a) where the whitened prediction's energy sits (p ≈ 17,900)")
    ax.legend(frameon=False, loc="upper center", fontsize=6.5, ncol=1, handlelength=1.2)

    ax = axes[1]
    arms = [
        ("hard cut (noise directions removed)", a272["arms"]["hard:0.5"]["pooled"], CB["red"]),
        ("Wiener damping κ = 4", a272["arms"]["wiener:k4"]["pooled"], CB["red"]),
        ("Wiener damping κ = 1", a272["arms"]["wiener:k1"]["pooled"], CB["red"]),
        ("Wiener κ = 1, shuffled control", a272["arms"]["wiener:k1:shuf"]["pooled"], CB["grey"]),
        ("noise floor c = 1", a272["arms"]["floor:c1"]["pooled"], CB["orange"]),
        ("noise floor c = 1, shuffled control", a272["arms"]["floor:c1:shuf"]["pooled"], CB["grey"]),
        ("286 signal-cov. Wiener φ = 0.25", a286["arms"]["e1:phi0.25:nm"]["pooled"], CB["purple"]),
        ("286 φ = 0.25, random-direction control", a286["arms"]["e1:phi0.25:rand:nm"]["pooled"], CB["grey"]),
        ("286 low-rank r = 100", a286["arms"]["e2:r100:nm"]["pooled"], CB["purple"]),
    ]
    y = np.arange(len(arms))[::-1]
    for yi, (lab, p, c) in zip(y, arms):
        ax.plot([p["lo"], p["hi"]], [yi, yi], color=c, lw=1.8, alpha=0.6)
        ax.scatter([p["point"]], [yi], s=34, color=c, zorder=3)
        ax.text(p["point"], yi + 0.24, f"{p['point']:+.4f}", ha="center", fontsize=6.6, color="#333333")
    ax.axvline(0, color=CB["black"], lw=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([a[0] for a in arms], fontsize=7.0)
    ax.set_xlabel(r"$\Delta$`pds` vs the incumbent, PLOLO, 3 cancer folds (n = 853 targets)")
    ax.set_title("(b) damping those directions loses retrieval")
    ax.set_ylim(y[-1] - 0.7, y[0] + 0.7)

    save(
        fig,
        "F-6",
        "Rank retrieval is paid in the whitener's noise-dominated directions. (a) Unit = eigen-directions "
        "of each fitted operator (p = 17,806–17,942). Grey: the share of directions on which an "
        "independent estimate of the donors' sampling noise exceeds the fitted eigenvalue; blue: the "
        "share of the incumbent's whitened prediction energy in the noise-dominated set (defined in "
        "R272 §1.1 as the directions whose estimated signal share is below one half); light blue: the "
        "share in the stricter noise ≥ eigenvalue set (`analysis272.json`, `spectra`). (b) Unit = targets; "
        "Δ`pds` against the incumbent pooled over three PLOLO cancer folds, 95 % paired target bootstrap, "
        "from `analysis272.json` (`arms`) and `analysisB286.json` (`arms`). Every damping of those "
        "directions loses retrieval; the shuffled and random-direction controls lose more, so the loss is "
        "not a generic amplitude effect.",
        "experiments/272-noise-corrected-operator/work/analysis272.json; "
        "experiments/286-crossproduct-estimator/work/analysisB286.json",
    )


# ----------------------------------------------------------------------------
# F-7  the per-target p53 mechanism on Broad-16
# ----------------------------------------------------------------------------
def fig7() -> None:
    rows = read_csv(os.path.join(E289, "p1_scores289.csv"))
    hpv = {"SKGII", "C4I"}
    rows = sorted(rows, key=lambda r: -f(r["p53_UBE3A"]))
    labs = [r["line"] for r in rows]
    p53 = [f(r["p53_UBE3A"]) for r in rows]
    e2f = [f(r["e2f_UBE3A"]) for r in rows]
    x = np.arange(len(labs))
    w = 0.38

    fig, ax = plt.subplots(figsize=(7.6, 3.6))
    c53 = [CB["red"] if l in hpv else CB["grey"] for l in labs]
    ax.bar(x - w / 2, p53, w, color=c53, label="P53_PATHWAY score of the line-specific residual")
    ax.bar(x + w / 2, e2f, w, color=CB["sky"], label="E2F_TARGETS score (specificity check)")
    ax.axhline(0, color=CB["black"], lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(labs, rotation=45, ha="right")
    ax.set_ylabel("programme score of the residual after UBE3A knockout")
    ax.set_title("HPV-E6 lines (red) are the two with a p53 response — found blind by the LOLO")
    ax.legend(frameon=False, loc="lower left")
    for xi, l in zip(x, labs):
        if l in hpv:
            ax.text(xi - w / 2, max(p53) * 1.04, "HPV-E6", ha="center", fontsize=6.6, color=CB["red"])

    save(
        fig,
        "F-7",
        "The source paper's per-target mechanism reproduces in our residual object. Unit = cell lines "
        "(n = 12 working Broad-16 lines). After UBE3A knockout, the P53_PATHWAY score of the "
        "line-specific residual is highest in the two HPV-E6 lines (SKGII +0.10, C4I +0.07; every other "
        "line negative); the E2F_TARGETS direction shows the opposite ordering. Exact permutation over "
        "the 55 line pairs gives P = 0.018 for the p53 direction and P = 1.0 for E2F "
        "(`p1p2_289.json`); values re-plotted from `p1_scores289.csv`.",
        "experiments/289-broad16-state-features/work/p1_scores289.csv",
    )


# ----------------------------------------------------------------------------
# F-8  the knee: pds as a threshold in the per-target transfer cosine
# ----------------------------------------------------------------------------
def fig8() -> None:
    """Parse section C of 273's saved table (the only file that carries both
    references' bin means side by side)."""
    txt = open(os.path.join(E273, "tables273.md")).read()
    block = txt.split("## C --")[1].split("\n## ")[0]  # to the next section header
    rows = []
    for line in block.splitlines():
        if not line.startswith("| ("):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        lab = cells[0]
        n_inc, pds_inc, n_base, pds_base, pds_plolo = cells[1:6]
        rows.append((lab, f(n_inc), f(pds_inc), f(n_base), f(pds_base), f(pds_plolo)))
    rows = [r for r in rows if np.isfinite(r[5])]
    assert len(rows) >= 6, f"F-8: parsed only {len(rows)} bins from tables273.md section C"

    labs = [
        r[0].replace("(-inf, ", "≤ ").replace("(0.05, inf]", "> 0.05")
        for r in rows
    ]
    x = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(7.8, 3.7))
    ax.plot(x, [r[5] for r in rows], marker="s", ms=5, lw=1.8, color=CB["blue"], label="PLOLO, incumbent (n = 853; 290 decoys)")
    flex = [r[2] for r in rows]
    ax.plot(x, flex, marker="o", ms=5, lw=1.8, color=CB["orange"], label="Flex-400 h1-test, incumbent (n = 85; 85 decoys)")
    for xi, r in zip(x, rows):
        if np.isfinite(r[1]) and r[1] > 0 and np.isfinite(r[2]):
            ax.text(xi, r[2] + 0.045, f"n {int(r[1])}", ha="center", fontsize=6.2, color=CB["orange"])
    knee = next(i for i, r in enumerate(rows) if r[0].startswith("(0.02"))
    ax.axvspan(knee - 0.5, knee + 0.5, color=CB["grey"], alpha=0.16)
    ax.text(knee, 0.20, "the knee:\nfirst bin ≥ 0.97\non both references", ha="center", fontsize=7, color="#444444")
    ax.axhline(0.5, ls=":", lw=1.0, color=CB["black"])
    ax.text(x[0], 0.52, "chance", fontsize=6.8, color="#444444", ha="left")
    ax.set_xticks(x)
    ax.set_xticklabels(labs, rotation=30, ha="right")
    ax.set_xlabel("per-target own cosine of the prediction to its own truth (bin)")
    ax.set_ylabel("mean per-target `pds`")
    ax.set_ylim(0, 1.10)
    ax.set_title("`pds` is a threshold in the transfer cosine at ≈ 0.02, on both references")
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(0.02, 0.72))

    save(
        fig,
        "F-8",
        "The scorer's rank statistic saturates. Unit = targets; mean per-target `pds` by bin of the "
        "prediction's own cosine to its own truth, on PLOLO (three cancer folds pooled, n = 853, 290 "
        "decoys) and on the 400-cell Flex reference (n = 85, 85 decoys). Both cross 0.97 in the "
        "(0.02, 0.03] bin; the slope inside (0, 0.05] is 9.5 vs 6.0 `pds` per unit cosine, and that "
        "difference is the decoy count, not biology. Counts above the orange curve are the Flex "
        "targets per bin. Values re-plotted from `tables273.md` section C (no interval is stored per bin).",
        "experiments/273-flex-reference-instrument/work/tables273.md (section C)",
    )


# ----------------------------------------------------------------------------
# F-9  two functionals of the same counts: ballast cells decouple the members
# ----------------------------------------------------------------------------
def fig9() -> None:
    """Board entries only (data/board-entries/*.json, the scored raw members)."""
    def entry(i):
        d = json.load(open(os.path.join(REPO, "data", "board-entries", f"{i}.json")))
        def find(x):
            if isinstance(x, dict):
                if "score_avg" in x:
                    return x
                for v in x.values():
                    y = find(v)
                    if y:
                        return y
            return None
        return find(d)
    raw = [("pds", "pds_cosine", 1), ("mse", "expr_mse_unbiased_capped_norm", -1), ("nmae", "de_wilcoxon_lfc_nmae", -1),
           ("fid", "de_wilcoxon_direction_fidelity_yield_raw", 1), ("reach", "de_wilcoxon_direction_reach_raw", 1),
           ("jac", "de_wilcoxon_sig_jaccard", 1)]
    pairs = [("410, bulk key rendered into every cell (vs 395)", "EETok0iqIbzQFg4CTYMQ", "RzcGA26xZEigSuXgUp5t", CB["red"]),
             ("410b, the same bulk on ballast cells (vs 395)", "4WgVi6FIE4iVw9PmRUP3", "RzcGA26xZEigSuXgUp5t", CB["blue"]),
             ("433a, a DE-design-only change (vs 417e)", "FCbBcvGxF1HOFPxPISUC", "wPP8U70kj4mOIvtA83CO", CB["green"])]
    x = np.arange(len(raw)); wdt = 0.26
    fig, ax = plt.subplots(figsize=(7.8, 3.8))
    for j, (lab, a, b, col) in enumerate(pairs):
        ea, eb = entry(a), entry(b)
        v = [sgn * (ea[k] - eb[k]) for _, k, sgn in raw]
        ax.bar(x + (j - 1) * wdt, v, wdt, color=col, label=lab)
        for xi, vi in zip(x + (j - 1) * wdt, v):
            if abs(vi) < 1e-12:
                ax.text(xi, 0.002, "0", ha="center", fontsize=6.5, color=col)
    ax.axhline(0, lw=0.8, color=CB["black"])
    ax.axvspan(-0.5, 1.5, color=CB["grey"], alpha=0.10)
    ax.text(0.5, -0.062, "pseudobulk of summed counts", ha="center", fontsize=7, color="#444444")
    ax.text(3.5, -0.062, "per-cell CPM means and Wilcoxon ranks", ha="center", fontsize=7, color="#444444")
    ax.set_xticks(x); ax.set_xticklabels([r[0] for r in raw])
    ax.set_ylabel("raw member change (improvement > 0)")
    ax.set_ylim(-0.07, 0.095)
    ax.set_title("The scorer reads two functionals of the same counts; a submission can set them apart")
    ax.legend(frameon=False, fontsize=6.8, loc="upper right")
    save(
        fig,
        "F-9",
        "Unit = scored board entries (A/B/C, 900 perturbations each); raw members oriented so that up is better (mse and "
        "nmae sign-flipped). Rendered into every cell, a per-target bulk magnitude key bought mse 0.080 and cost the "
        "fold-change member 0.054 (red). Carried by a few ~995k-count ballast cells per perturbation, the same bulk kept "
        "its mse gain and left all four DE members within 0.0007 (the largest move nmae −0.000643; blue). Conversely, a change to the DE design alone left "
        "pds and mse bit-identical (green, the zeros). Exact board values, no intervals.",
        "data/board-entries/{RzcGA26xZEigSuXgUp5t,EETok0iqIbzQFg4CTYMQ,4WgVi6FIE4iVw9PmRUP3,wPP8U70kj4mOIvtA83CO,FCbBcvGxF1HOFPxPISUC}.json",
    )


# ----------------------------------------------------------------------------
# F-10  the shared programme belongs to the lab (fifteen banks, four labs; matched perturbation sets)
# ----------------------------------------------------------------------------
def fig10() -> None:
    rows = json.load(open(os.path.join(REPO, "experiments", "450-same-line-programme", "work", "labline450_v3.json")))
    order = [("k562_gwps", "K562 · W"), ("k562_essential", "K562-ess · W"), ("nadig_jurkat", "Jurkat · W"),
             ("nadig_hepg2", "HepG2 · W"), ("rpe1", "RPE1 · W"), ("orion_hct116", "HCT116 · X"),
             ("orion_hek293t", "HEK293T · X"), ("jiang_K562", "K562 · J"), ("jiang_A549", "A549 · J"),
             ("jiang_MCF7", "MCF7 · J"), ("jiang_HT29", "HT29 · J"), ("jiang_HAP1", "HAP1 · J"), ("jiang_BXPC3", "BxPC3 · J"),
             ("viperturb_k562", "K562 · S"), ("xucao_hek293", "HEK293 · C")]
    idx = {k: i for i, (k, _) in enumerate(order)}
    M = np.full((len(order), len(order)), np.nan)
    np.fill_diagonal(M, 1.0)
    for r in rows:
        i, j = idx[r["a"]], idx[r["b"]]
        M[i, j] = M[j, i] = r["r"]
    fig, ax = plt.subplots(figsize=(8.4, 7.2))
    im = ax.imshow(M, cmap="RdBu_r", vmin=-0.8, vmax=0.8)
    for i in range(len(order)):
        for j in range(len(order)):
            if i != j and np.isfinite(M[i, j]):
                ax.text(j, i, f"{M[i, j]:+.2f}", ha="center", va="center", fontsize=5.2,
                        color="white" if abs(M[i, j]) > 0.45 else "#222222")
    ax.set_xticks(range(len(order))); ax.set_yticks(range(len(order)))
    ax.set_xticklabels([o[1] for o in order], rotation=40, ha="right", fontsize=7)
    ax.set_yticklabels([o[1] for o in order], fontsize=7)
    for b in (4.5, 6.5, 12.5, 13.5):
        ax.axhline(b, color="black", lw=1.0); ax.axvline(b, color="black", lw=1.0)
    ax.set_title("Programme agreement follows the lab more than the line (matched perturbation sets)\n(W Weissman, X Xaira, J Jiang 2025, S Satija, C a second HEK293 lab)", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="r of (P_down − 0.5)")
    save(
        fig,
        "F-10",
        "Unit = pairs of public CRISPRi Perturb-seq banks; per-gene P_down (the fraction of a bank's perturbations that lower "
        "a gene), with BOTH programmes of a pair computed over only the perturbations the two banks share (the challenge panel "
        "excluded; blank = fewer than 100 shared). On pairs sharing >= 1,000 perturbations, different lines within a lab agree at "
        "mean r 0.54 and the same line across labs at 0.21; the line adds ~0.2 inside the cross-lab K562 screen's row. Small "
        "pathway screens (J, ~200 shared) are noisy. On the benchmark's Flex H1 reference only a same-line screen carries the "
        "recipient's programme (§3(j)). Point estimates, no intervals.",
        "experiments/450-same-line-programme/work/labline450_v3.json",
    )


def fig11() -> None:
    w = os.path.join(REPO, "experiments", "462-broad16-lab-signature", "work")
    hv = json.load(open(os.path.join(w, "halves462.json")))
    h1 = json.load(open(os.path.join(w, "h1halves462.json")))
    d8 = json.load(open(os.path.join(REPO, "experiments", "458-arc-lab-programme", "work", "denoise458d.json")))
    fl = json.load(open(os.path.join(w, "fidlike462.json")))
    bl = json.load(open(os.path.join(w, "blend462.json")))
    fig, (a, b, c3) = plt.subplots(1, 3, figsize=(13.2, 4.2), gridspec_kw={"width_ratios": [1.3, 1, 0.62]})
    # (a) reliability with a shared vs disjoint control
    lines = hv["lines"]
    sh = [hv["within_line_shared_control"][l] for l in lines]; dj = [hv["within_line_disjoint_controls"][l] for l in lines]
    xs = np.arange(len(lines))
    a.scatter(xs, sh, color=CB["orange"], s=22, label="halves share one control")
    a.scatter(xs, dj, color=CB["blue"], s=22, label="halves on disjoint controls")
    for x, u, v in zip(xs, sh, dj):
        a.plot([x, x], [v, u], color=CB["grey"], lw=0.7)
    k = len(lines) + 1.2
    pairs = [("Arc H1 bank\n(400 cells, 18k ctrl rows)", h1["expressed_cpm5"]["shared_control_r"], h1["expressed_cpm5"]["disjoint_control_r"]),
             ("Arc H1 all cells\n(~1,071 cells, 38k ctrl)", d8["shared_control_r"], d8["disjoint_control_r"])]
    for i, (lab, u, v) in enumerate(pairs):
        x = k + 1.6 * i
        a.scatter([x], [u], color=CB["orange"], s=44, marker="D"); a.scatter([x], [v], color=CB["blue"], s=44, marker="D")
        a.plot([x, x], [v, u], color=CB["grey"], lw=0.9)
    a.set_xticks(list(xs) + [k, k + 1.6])
    a.set_xticklabels(lines + [p[0] for p in pairs], rotation=60, ha="right", fontsize=6.3)
    a.set_ylim(0, 1.12); a.set_ylabel("split-half r of per-gene P_down")
    a.set_title("(a) a single control sample is part of the 'programme'\n(Broad-16: 12 lines, Cas9 KO; Arc H1: CRISPRi Flex)", fontsize=8.6)
    a.legend(fontsize=7, loc="upper left", ncol=2, frameon=False)
    # (b) fid-like accuracy on H1
    names = [("CLEAN4 (4 Weissman banks)", fl["clean4"]["acc"], CB["grey"]), ("Arc H1 bank", fl["bank_h1_2025"]["acc"], CB["sky"]),
             ("Arc H1 all cells (458d)", fl["full_h1_2025_full"]["acc"], CB["blue"]), ("LOTO oracle", fl["oracle_loto"]["acc"], CB["black"])]
    for i, (lab, v, c) in enumerate(names):
        b.barh(i, v - 0.5, left=0.5, color=c, height=0.62)
        b.text(v + 0.0008, i, f"{v:.4f}", va="center", fontsize=7)
    b.set_yticks(range(len(names))); b.set_yticklabels([n[0] for n in names], fontsize=7.4)
    ws = [0.0, 0.25, 0.5, 0.75, 1.0]
    vals = [bl[f"w_arc={x}"]["fidlike"] for x in ws]
    c3.plot(ws, vals, "-o", color=CB["blue"], ms=4, lw=1.2)
    c3.set_xlabel("weight on Arc (vs CLEAN4)"); c3.set_ylabel("fid-like sign agreement")
    c3.set_title("(c) CLEAN4 + Arc blends on H1", fontsize=8.6)
    b.set_xlim(0.5, 0.568); b.set_xlabel("sign agreement, all block genes (fid-like)")
    b.set_title("(b) on H1's own line, the clean Arc programme\nreaches the recipient oracle", fontsize=8.6)
    save(
        fig,
        "F-11",
        "Unit = genes (P_down per gene over a screen's perturbations). (a) Split-half reliability of a per-gene programme "
        "when both halves are measured against the SAME control sample (orange) vs against disjoint halves of the controls "
        "(blue). Broad-16 lines (about 500 control cells, 43-129 per target, 82 matched targets): 0.80 vs 0.12. Arc's H1 "
        "harness bank: 0.79 vs 0.39. The same programme from all of Arc's public H1 cells: 0.92 vs 0.69. (b) Sign agreement "
        "of each per-gene programme with the reference fold change on all 2,200 block genes of the 85 h1-test targets "
        "(187,000 pairs; CLEAN4 covers 156k). (c) CLEAN4 + Arc blends, z-weighted. Same line and lab as the Arc source, so an upper "
        "reference for the cross-line board test (458d). One seeded split per point, no intervals.",
        "experiments/462-broad16-lab-signature/work/{halves462,h1halves462,fidlike462,blend462}.json; experiments/458-arc-lab-programme/work/denoise458d.json",
    )


# ----------------------------------------------------------------------------
# F-12  the zero-shot pipeline re-fitted from the panel alone; the two channels read single-channel changes separately
# ----------------------------------------------------------------------------
BOARD_IDS = {  # label -> captured entry id (data/board-entries/<id>.json); ids from each BOARD-RESULT file
    "446b": "UDgjOKL9gI1udAexq6ue", "446bj2": "YvVru3DAF89Jiyr34olR", "446br2": "wE7usCaR6TAlPp16TR0V",
    "453dr": "q3VuTImh4FRGKXH2waCp", "446bj2z": "be90smZUDcIV7j5MgQ03", "457dr": "1FfirZsGsAZvaXYTd39J",
    "446bd2": "R5uE5Rve3R9nWXGbeSLM",
}


def board_entry(tag: str) -> dict:
    return json.load(open(os.path.join(REPO, "data", "board-entries", f"{BOARD_IDS[tag]}.json")))["entry"]


def channel_split(arm: str, base: str) -> tuple[float, float]:
    """(bulk, DE) contributions to the overall change, in overall units: the scaled members' deltas / 6.
    Pure arithmetic on the published scaled members; bulk = pds + mse, DE = nmae + fid + reach + jac."""
    a, b = board_entry(arm), board_entry(base)
    bulk = sum(a[f"score_{k}"] - b[f"score_{k}"] for k in ("pds", "mse")) / 6
    de = sum(a[f"score_{k}"] - b[f"score_{k}"] for k in ("nmae", "fid", "reach", "jac")) / 6
    return bulk, de


def fig12() -> None:
    """Board entries only. Every number is a difference of two published scaled members."""
    fig, (a, b) = plt.subplots(1, 2, figsize=(12.6, 4.3), gridspec_kw={"width_ratios": [1.0, 1.25]})
    # (a) the October path's decomposition on relabelled A/B/C
    steps = [("refit path at CD4 ×1\n453dr − 446b", "453dr", "446b", (-0.0045, -0.0005)),
             ("CD4 ×1.5 restored\n457dr − 453dr", "457dr", "453dr", (0.0005, 0.0035)),
             ("net: refit at the tuned dose\n457dr − 446b", "457dr", "446b", None)]
    xs = np.arange(len(steps))
    for x, (lab, arm, base, band) in zip(xs, steps):
        bulk, de = channel_split(arm, base)
        a.bar(x - 0.17, bulk, 0.32, color=CB["blue"], label="bulk members (pds, mse)" if x == 0 else None)
        a.bar(x + 0.17, de, 0.32, color=CB["orange"], label="DE members (nmae, fid, reach, jac)" if x == 0 else None)
        tot = bulk + de
        a.plot([x - 0.36, x + 0.36], [tot, tot], color=CB["black"], lw=1.4, label="overall change" if x == 0 else None)
        a.text(x, tot + (0.00028 if tot >= 0 else -0.00055), f"{tot:+.4f}", ha="center", fontsize=7)
        if band:
            a.add_patch(plt.Rectangle((x - 0.42, band[0]), 0.84, band[1] - band[0], color=CB["grey"], alpha=0.16, lw=0))
    a.axhline(0, color=CB["black"], lw=0.8)
    a.set_xticks(xs); a.set_xticklabels([s[0] for s in steps], fontsize=7.2)
    a.set_ylabel("change in the overall score (board units)")
    a.set_title("(a) every panel-dependent component re-fitted from the panel alone:\n"
                "the whole bulk cost is one composition choice (grey: preregistered bands)", fontsize=8.4)
    a.legend(frameon=False, fontsize=6.8, loc="upper right")
    # (b) single-channel changes, read by the two channels separately
    arms = [("446bj2 − 446b: per-target bulk key (LODO)", "446bj2", "446b", "bulk"),
            ("453dr − 446b: the refit path (CD4 ×1)", "453dr", "446b", "bulk"),
            ("457dr − 453dr: CD4 ×1.5 on the refit path", "457dr", "453dr", "bulk"),
            ("446bd2 − 446b: a fifth donor line in the bulk", "446bd2", "446b", "bulk"),
            ("446br2 − 446b: zero-gene repair floor −2", "446br2", "446b", "DE"),
            ("446bj2z − 446bj2: exact-rank magnitude pass", "446bj2z", "446bj2", "DE")]
    ys = np.arange(len(arms))[::-1]
    for y, (lab, arm, base, kind) in zip(ys, arms):
        bulk, de = channel_split(arm, base)
        b.barh(y + 0.17, bulk, 0.32, color=CB["blue"])
        b.barh(y - 0.17, de, 0.32, color=CB["orange"])
        for v, dy, col in ((bulk, 0.17, CB["blue"]), (de, -0.17, CB["orange"])):
            if v == 0.0:
                b.text(0.00005, y + dy, "0 (bit-identical)", va="center", fontsize=6.2, color=col)
    b.axvline(0, color=CB["black"], lw=0.8)
    b.set_yticks(ys); b.set_yticklabels([f"{a_[0]}  [{a_[3]} change]" for a_ in arms], fontsize=7)
    b.set_xlabel("contribution to the overall change (board units)")
    b.set_title("(b) a bulk-only change leaves the DE members inside their re-roll noise;\n"
                "an exact-rank DE change leaves the other five members bit-identical", fontsize=8.4)
    save(
        fig,
        "F-12",
        "Unit = scored board entries (A/B/C, 900 perturbations each); every bar is a difference of two published scaled "
        "members divided by 6 (the overall is the unweighted mean of six), split into the bulk members (pds, mse; blue) and "
        "the four DE members (orange). (a) The October day-of path (identity match, CD4 extraction, whitening-operator refit, "
        "bulk-key refit, control-proxy re-level, relabelling) run on A/B/C relabelled as new contexts: it costs −0.0037 "
        "overall at CD4 ×1 (bulk −0.0031, inside its preregistered band), restoring CD4 ×1.5 returns +0.0031 (bulk), and the "
        "net refit at the tuned dose sits −0.0005 from the tuned entry with bulk +0.00001: the DE remainder is reach, inside "
        "the ordering channel's re-roll sd of 0.0045 overall. (b) Six single-channel entries since 446b. Exact board values, "
        "no intervals; the ordering-channel noise between digest-identical arms is sd 0.0045 overall.",
        "data/board-entries/{UDgjOKL9gI1udAexq6ue,YvVru3DAF89Jiyr34olR,wE7usCaR6TAlPp16TR0V,q3VuTImh4FRGKXH2waCp,be90smZUDcIV7j5MgQ03,1FfirZsGsAZvaXYTd39J,R5uE5Rve3R9nWXGbeSLM}.json",
    )


# ----------------------------------------------------------------------------
# F-13  the uncalled-gene hook's bar is lower than the block prior's (517's saved H1 grids)
# ----------------------------------------------------------------------------
def fig13() -> None:
    """Means over the three emission seeds of 517's saved analytic grids (validated against real emitter runs in
    517 §1), plus 517's real-emitter runs as markers. No statistic is recomputed."""
    w = os.path.join(REPO, "experiments", "517-flex-programme", "work")
    grid = read_json(os.path.join(w, "grid517.json"))["rows"]
    hook = read_json(os.path.join(w, "hookgrid517.json"))
    mz = read_json(os.path.join(w, "mzgrid517.json"))
    real = read_json(os.path.join(w, "real517.json"))
    R = [0.0, 0.1, 0.2, 0.3, 0.4]

    def mean(rows, key):
        return float(np.mean([x[key] for x in rows]))

    blk_g = [mean([x for x in grid if x["gated"] and x["q"] == 0.8 and x["w"] == 1.0 and x["tau"] is None and x["r"] == r], "ov") for r in R]
    blk_u = [mean([x for x in grid if (not x["gated"]) and x["q"] == 0.8 and x["w"] == 1.0 and x["tau"] is None and x["r"] == r], "ov") for r in R]
    pd_rows = [[x for x in hook if x["src"] == "arc" and x["gated"] and x["tau"] == 0.05 and x["r"] == r] for r in R]
    pd_ov = [mean(v, "overall_equiv") for v in pd_rows]; pd_nm = [mean(v, "nmae") for v in pd_rows]
    mz_rows = [[x for x in mz if x["arm"] == "hook_g" and x["r"] == r] for r in R]
    mz_ov = [mean(v, "ov") for v in mz_rows]; mz_nm = [mean(v, "nmae") for v in mz_rows]
    c4 = [x for x in hook if x["src"] == "clean4"]
    c4_ov, c4_nm = mean(c4, "overall_equiv"), mean(c4, "nmae")

    fig, (a, b) = plt.subplots(1, 2, figsize=(12.4, 4.3))
    a.plot(R, blk_u, "-o", ms=3.5, color=CB["red"], label="programme as block prior, every target")
    a.plot(R, blk_g, "-o", ms=3.5, color=CB["purple"], label="programme as block prior, design-strong targets (gated)")
    a.plot(R, pd_ov, "-s", ms=3.5, color=CB["sky"], label="uncalled-gene hook, programme P_down (gated)")
    a.plot(R, mz_ov, "-s", ms=3.5, color=CB["blue"], label="uncalled-gene hook, programme mean z, non-essential rows (gated)")
    a.axhline(c4_ov, color=CB["grey"], ls="--", lw=1.0, label="uncalled-gene hook, cross-line consensus (reference)")
    a.axhline(0, color=CB["black"], lw=0.8)
    rm = real["c4-mzg0.05"]; rp = real["c4-arcg0.05"]; rp25 = real["c4-arcg0.05-r0.25"]
    a.plot([0.0], [np.mean([x["ov"] for x in rm])], "D", color=CB["blue"], ms=6, mfc="white", label="real emitter runs (3 seeds)")
    a.plot([0.0, 0.25], [np.mean([x["ov"] for x in rp]), np.mean([x["ov"] for x in rp25])], "D", color=CB["sky"], ms=6, mfc="white")
    a.set_xlabel("r = fraction of the programme's sign information destroyed (random flips)")
    a.set_ylabel("H1 overall-equivalent vs the incumbent (CLEAN4 block, no hook)")
    a.set_title("(a) all six members: the block prior falls below the incumbent at r ≈ 0.17 (every target) … 0.24 (gated);\n"
                "the hooks at r ≈ 0.28 (mean z) … 0.37 (P_down); from r 0.2 the hooks carry a small-gate reach term", fontsize=8.0)
    a.legend(frameon=False, fontsize=6.3, loc="upper right")
    b.plot(R, pd_nm, "-s", ms=3.5, color=CB["sky"], label="hook, programme P_down (gated)")
    b.plot(R, mz_nm, "-s", ms=3.5, color=CB["blue"], label="hook, programme mean z, non-essential rows (gated)")
    b.axhline(c4_nm, color=CB["grey"], ls="--", lw=1.0, label=f"hook, cross-line consensus ({c4_nm:+.4f})")
    b.plot([0.0], [np.mean([x["nmae"] for x in rm])], "D", color=CB["blue"], ms=6, mfc="white", label="real emitter runs (3 seeds)")
    b.plot([0.0, 0.25], [np.mean([x["nmae"] for x in rp]), np.mean([x["nmae"] for x in rp25])], "D", color=CB["sky"], ms=6, mfc="white")
    b.set_xlabel("r = fraction of the programme's sign information destroyed (random flips)")
    b.set_ylabel("H1 Δ nmae, raw (lower is better)")
    b.set_title("(b) nmae, the member the hook moves: it beats the consensus hook\nuntil r ≈ 0.38 (P_down) … 0.42 (mean z)", fontsize=8.2)
    b.legend(frameon=False, fontsize=6.3, loc="lower right")
    save(
        fig,
        "F-13",
        "Unit = h1-test targets (85; nmae is scored on 58), means over three emission seeds of 517's analytic pricer, which "
        "reproduces the production emitter's nmae to 5e-5 and fid to ~10 % (517 §1); open diamonds are real production-emitter "
        "runs (517 §9: mean-z hook nmae −0.0183, P_down hook −0.0149, the P_down hook at r 0.25 −0.0063). The programme is Arc's "
        "H1 CRISPRi screen with h1-test's own targets excluded; r flips a random fraction of its gene signs, emulating a "
        "programme that carries less information about the recipient. As a block prior (red, purple) it replaces the "
        "cross-line consensus on the high-expression block and must beat that consensus; as the uncalled-gene hook (blues) it "
        "replaces the design's own signs on genes the design does not call (0.532 accurate on H1's uncalled gate genes) and "
        "only has to beat those. Crossings are read off the plotted means by linear interpolation. From r 0.2 the analytic "
        "hooks carry a reach term of −0.009 raw from a few small-gate targets, which H1 cannot price for sign priors; that is "
        "why the hook's lead over the consensus hook is shorter in (a), all six members, than in (b), nmae alone. 458d's "
        "preregistered fid bins map onto r ≈ ≤ 0.23 / 0.23–0.25 / 0.25–0.27 / ≥ 0.28 (positive / partial / null / harm; 517 "
        "§0.2, coarse). H1 prices only; the board read of the hook is pending (527).",
        "experiments/517-flex-programme/work/{grid517,hookgrid517,mzgrid517,real517}.json",
    )


def main() -> int:
    print("294 make_figs.py (533 revision) — re-plotting saved tables")
    for fn in (fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9, fig10, fig11, fig12, fig13):
        try:
            fn()
        except Exception as exc:  # keep going; report at the end
            print(f"  !! {fn.__name__} FAILED: {exc!r}", file=sys.stderr)
            raise
    with open(os.path.join(FIGS, "MANIFEST.md"), "w") as fh:
        fh.write("# Figure manifest — experiment 294 (revised by 533, 2026-09-30)\n\n")
        fh.write("Every figure is regenerated by `uv run python make_figs.py`.\n\n")
        for name, caption, sources in MANIFEST:
            fh.write(f"## {name}\n\n**Caption.** {caption}\n\n**Sources.** `{sources}`\n\n")
    print(f"wrote {len(MANIFEST)} figures + figs/MANIFEST.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
