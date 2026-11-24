"""296 -- PREREGISTRATION sec 4: the on-target failure rate of each library's constructs, beside 293's
18.7 % / 4.3 % for Arc's own 2025 panel (recomputed here from 293's saved `kd_all`, not quoted).

    PYTHONDONTWRITEBYTECODE=1 uv run python rates296.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib296 as X  # noqa: E402
from lib296 import L  # noqa: E402

P = L.Logger(X.WORK / "rates296.log")
W = X.WORK
E293 = X.ROOT / "experiments/293-challenge-construct-ceiling/work"


def rate(v, label, extra=None) -> dict:
    v = np.abs(np.asarray([x for x in v if np.isfinite(x)], float))
    if v.size == 0:
        return dict(label=label, n=0)
    r = dict(label=label, n=int(v.size), median_abs_lfc=float(np.median(v)),
             frac_lt_1=float((v < 1.0).mean()), frac_lt_05=float((v < 0.5).mean()))
    if extra:
        r.update(extra)
    P(f"  {label:52s} n {r['n']:6d}   median |lfc| {r['median_abs_lfc']:.3f}   "
      f"|lfc| < 1.0 {100 * r['frac_lt_1']:5.1f} %   < 0.5 {100 * r['frac_lt_05']:5.1f} %")
    return r


def main() -> None:
    t0 = time.time()
    P(f"296 rates296.py -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    panel = [l.split(",")[0] for l in
             (X.ROOT / "data/challenge/pert_counts.csv").read_text().strip().split("\n")[1:]]
    panel_set = set(panel)
    out: dict = {"n_panel_2026": len(panel), "rows": {}}

    # ---------------------------------------------------------------- Arc's 2025 panel (293's own arrays)
    kd25: dict[str, float] = {}
    for sp in ("train", "validation", "test"):
        f = E293 / f"run293_{sp}.json"
        if f.exists():
            kd25.update(json.load(open(f))["kd_all"])
    if kd25:
        out["rows"]["vcc2025_panel"] = rate(list(kd25.values()),
                                            "Arc VCC-2025 panel (293's kd_all, 3 splits, pooled cells)")

    # ---------------------------------------------------------------- K562-GWPS
    jk = json.load(open(W / "k562_296.json"))
    ks = jk["kd_single"]
    pan = [v["lfc"] for t, v in ks.items() if t in panel_set]
    non = [v["lfc"] for t, v in ks.items() if t not in panel_set]
    out["rows"]["k562_panel"] = rate(pan, "K562-GWPS single-construct, 2026-PANEL targets")
    out["rows"]["k562_nonpanel"] = rate(non, "K562-GWPS single-construct, non-panel targets")
    out["rows"]["k562_multi_constructs"] = rate([v["lfc"] for v in jk["kd_cons"].values()],
                                                "K562-GWPS: the individual two-construct constructs")
    cov = sorted(set(ks) & panel_set)
    out["k562_panel_coverage"] = dict(n_panel_on_axis=len(cov), of=len(panel))
    P(f"  [K562-GWPS covers {len(cov)}/{len(panel)} of the 2026 panel as single-construct perturbations "
      f"with a measurable on-target column]")

    # ---------------------------------------------------------------- Nadig
    for wh in ("jurkat", "hepg2"):
        f = W / f"nadig_{wh}_296.json"
        if not f.exists():
            continue
        j = json.load(open(f))
        v = [x["lfc"] for x in j["kd_target"].values() if x["n_constructs"] == 1]
        out["rows"][f"nadig_{wh}_single"] = rate(v, f"Nadig-{wh} single-construct targets")
        out["rows"][f"nadig_{wh}_multi"] = rate([x["lfc"] for x in j["kd_cons"].values()],
                                                f"Nadig-{wh}: the individual two-construct constructs")

    # ---------------------------------------------------------------- Broad-16 (Cas9)
    jb = json.load(open(W / "broad_296.json"))
    full = [v["lfc"] for k, v in jb["kd"].items() if k.endswith("|full")]
    cons = [v["lfc"] for k, v in jb["kd"].items() if k.endswith("|v1") or k.endswith("|v2")]
    out["rows"]["broad_full"] = rate(full, "Broad-16 (Cas9) per (line, target), all cells")
    out["rows"]["broad_constructs"] = rate(cons, "Broad-16 (Cas9) per (line, target, construct)")
    per_line: dict[str, dict] = {}
    for ln in jb["lines"]:
        v = [x["lfc"] for k, x in jb["kd"].items() if k.startswith(f"{ln}|") and k.endswith("|full")]
        v = np.abs(np.asarray([x for x in v if np.isfinite(x)], float))
        if v.size:
            per_line[ln] = dict(n=int(v.size), frac_lt_1=float((v < 1.0).mean()),
                                frac_lt_05=float((v < 0.5).mean()), median=float(np.median(v)))
    out["broad_per_line"] = per_line

    # ---------------------------------------------------------------- CD4
    jc = json.load(open(W / "cd4_296.json"))
    out["rows"]["cd4_guides"] = rate([v["lfc"] for v in jc["kd"].values()],
                                     "CD4 GWCD4i per (condition, guide), my own on-target lfc")
    kdt = pd.read_csv(X.ROOT / "data/external/cd4_cell2026/suppl_tables_2026-08-24/"
                              "guide_kd_efficiency.suppl_table.csv", index_col=0)
    sig = kdt.signif_knockdown.astype(bool)
    out["cd4_authors"] = dict(n=int(len(kdt)), frac_not_signif=float(1 - sig.mean()),
                              by_cond={c: float(1 - g.signif_knockdown.astype(bool).mean())
                                       for c, g in kdt.groupby("culture_condition")},
                              n_high_conf_no_effect=int(kdt.high_confidence_no_effect_guides.sum()))
    P(f"  CD4, the authors' own call: {100 * out['cd4_authors']['frac_not_signif']:.1f} % of "
      f"{len(kdt):,} (guide, condition) rows are NOT a significant knockdown "
      f"({out['cd4_authors']['by_cond']}); {out['cd4_authors']['n_high_conf_no_effect']} rows flagged "
      f"'high confidence no effect'")

    # ---------------------------------------------------------------- indicative: the 2026 panel in K562
    pv = np.abs(np.asarray([v["lfc"] for t, v in ks.items() if t in panel_set
                            and np.isfinite(v["lfc"])], float))
    out["indicative_2026_panel_in_k562"] = dict(
        n=int(pv.size), frac_lt_1=float((pv < 1.0).mean()), frac_lt_05=float((pv < 0.5).mean()),
        note="INDICATIVE ONLY: a different vector, a different line and a different chemistry from the "
             "2026 challenge library; it is a base rate for CRISPRi constructs against these genes, "
             "not a prediction about the challenge's own truths.")
    P(f"  INDICATIVE: of the {pv.size} 2026-panel targets measurable in K562-GWPS, "
      f"{100 * (pv < 1.0).mean():.1f} % have |lfc| < 1.0 and {100 * (pv < 0.5).mean():.1f} % < 0.5")

    # ---------------------------------------------------------------- descriptive tail (POST-HOC)
    desc: dict = {}
    for tag, name in (("k562_296", "K562-GWPS"), ("nadig_jurkat_296", "Nadig-jurkat"),
                      ("nadig_hepg2_296", "Nadig-hepg2")):
        f = W / f"{tag}.json"
        if not f.exists():
            continue
        j = json.load(open(f))
        v = sorted(((c, x["lfc"]) for c, x in j["kd_cons"].items() if np.isfinite(x["lfc"])),
                   key=lambda z: z[1])
        desc[name] = dict(most_efficient=[[c, round(l, 3)] for c, l in v[:20]],
                          least_efficient=[[c, round(l, 3)] for c, l in v[-20:]])
    jb2 = {k: v["lfc"] for k, v in jb["kd"].items() if not k.endswith("|full") and np.isfinite(v["lfc"])}
    vb = sorted(jb2.items(), key=lambda z: z[1])
    desc["Broad-16"] = dict(most_efficient=[[c, round(l, 3)] for c, l in vb[:20]],
                            least_efficient=[[c, round(l, 3)] for c, l in vb[-20:]])
    vc = sorted(((k, v["lfc"]) for k, v in jc["kd"].items() if np.isfinite(v["lfc"])), key=lambda z: z[1])
    desc["CD4-GWCD4i"] = dict(most_efficient=[[c, round(l, 3)] for c, l in vc[:20]],
                              least_efficient=[[c, round(l, 3)] for c, l in vc[-20:]])
    out["descriptive_top20"] = desc

    # promoter / transcript-level construct stems (POST-HOC): the 2025 library's _P1 / _P2 / _P1P2 pattern
    stems: dict = {}
    for tag, name in (("k562_296", "K562-GWPS"), ("nadig_jurkat_296", "Nadig-jurkat")):
        f = W / f"{tag}.json"
        if not f.exists():
            continue
        j = json.load(open(f))
        ex = list(j["kd_cons"])[:8]
        stems[name] = dict(example_construct_ids=ex,
                           has_promoter_suffix=bool(any("_P" in c for c in j["kd_cons"])))
    out["construct_id_shape"] = stems

    out["seconds"] = round(time.time() - t0, 1)
    L.Logger.dump(P, out, W / "rates296.json")
    P(f"\ndone {out['seconds']}s")
    P.save()


if __name__ == "__main__":
    main()
