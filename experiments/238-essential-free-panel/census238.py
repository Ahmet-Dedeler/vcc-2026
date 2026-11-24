"""238 Phase 0 — the FACT, before any scoring: what the essential list is, and how much of
the October asset lives on it.

THE PROBLEM THIS FILE ESTABLISHES.  233 and 236 priced `nadig_hepg2`, `nadig_jurkat`, `rpe1`
and `k562_essential` at +0.0066..+0.0104 six-member composite on `DEF-RUNBOOK` §5's two
synthetic panels.  Both experiments are internally sound (236 reproduces 233 at max |d| 0.0e+00
on 157/157 rows).  But those panels were drawn from DONOR COVERAGE, and donor libraries screen
common-essential genes, so the draws carry 34 / 30 essential genes.  **Arc's one observed panel
carries ZERO** (`DEF-RUNBOOK` §5), and there is an assay reason for that: an essential knockdown
kills the cell, so there is no transcriptional response to score.

Three questions, all answerable from `perts.txt` files and `pert_counts.csv` alone:

  0.1  reproduce `DEF-RUNBOOK` §5's `panel n essential list` row -- **0 / 34 / 30** -- with
       independent code, and STATE WHICH LIST and where it comes from.
  0.3  for each of the four sources, what FRACTION of its carrying rows on each panel is an
       essential gene?  If that is ~100%, the conclusion follows before anything is scored.
  ---  and the eligible-pool arithmetic for Phase 1's essential-free draw.

PROVENANCE OF THE LIST (asserted at runtime, not assumed).  The list is
`data/donors/k562_essential/perts.txt`: the perturbation set of
`data/external/scperturb/ReplogleWeissman2022_K562_essential.h5ad`, i.e. Replogle 2022's K562
**essential-gene** Perturb-seq screen, built by `scripts/build_donor_table.py`.  It is the SAME
list `scripts/build_union_whitener.py --basis-prune essential` uses to drop basis rows and the
same one `def_census.py` reports the panel overlap against, so 0.1 below is a reproduction of
the published row and not a different statistic with the same name.

    uv run python census238.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE))

import panels238 as P  # noqa: E402

# the four sources 233/236 priced, and where their perturbation lists live
OCTOBER_ONLY = {
    "nadig_hepg2": "experiments/166-donor-acquisition/donors/nadig_hepg2",
    "nadig_jurkat": "data/donors/nadig_jurkat",
    "rpe1": "data/donors/rpe1",
    "k562_essential": "data/donors/k562_essential",
}
PRODUCTION = {
    "k562_gwps": "data/donors/k562_gwps",
    "orion_hct116": "data/donors/orion_hct116",
    "orion_hek293t": "data/donors/orion_hek293t",
}
CD4_ROOT = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"
CD4_CONDITIONS = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")

# DEF-RUNBOOK §5's published row, and §5.-2 / 225 §8.3's coverage table -- the known answers
KA_ESSENTIAL_OVERLAP = {"abc": 0, "real": 34, "pess": 30}
KA_COVERAGE = {
    ("abc", "k562_gwps"): 269, ("abc", "orion_hct116"): 290, ("abc", "orion_hek293t"): 294,
    ("real", "k562_gwps"): 300, ("real", "orion_hct116"): 300, ("real", "orion_hek293t"): 300,
    ("pess", "k562_gwps"): 150, ("pess", "orion_hct116"): 291, ("pess", "orion_hek293t"): 293,
    ("abc", "nadig_hepg2"): 0, ("abc", "nadig_jurkat"): 0, ("abc", "rpe1"): 0,
    ("abc", "k562_essential"): 0,
    ("pess", "nadig_hepg2"): 33, ("pess", "nadig_jurkat"): 33, ("pess", "rpe1"): 32,
    ("pess", "k562_essential"): 30,
}


def read_list(rel: str) -> set[str]:
    return {ln.strip() for ln in (ROOT / rel / "perts.txt").read_text().splitlines()
            if ln.strip()}


def sha256_text(rel: str) -> str:
    return hashlib.sha256((ROOT / rel / "perts.txt").read_bytes()).hexdigest()


def cd4_complete() -> set[str]:
    """Production's CD4 rule: a target counts only if all three conditions carry it."""
    per = [{x.strip() for x in (CD4_ROOT / c / "perts.txt").read_text().splitlines()
            if x.strip()} for c in CD4_CONDITIONS]
    return set.intersection(*per)


def main() -> None:
    ess_rel = "data/donors/k562_essential"
    essential = read_list(ess_rel)
    meta = json.loads((ROOT / ess_rel / "meta.json").read_text())
    assert meta["source"].endswith("ReplogleWeissman2022_K562_essential.h5ad"), meta["source"]

    out: dict = {"essential_list": dict(
        path=f"{ess_rel}/perts.txt",
        sha256=sha256_text(ess_rel),
        n=len(essential),
        upstream=meta["source"],
        upstream_n_perturbations=meta["n_perturbations"],
        used_by=["scripts/build_union_whitener.py --basis-prune essential",
                 "scripts/def_census.py (the panel n essential row)",
                 "scripts/make_submission_pm_cd4.py (the pruned operator's provenance)"],
    )}
    print(f"ESSENTIAL LIST  {ess_rel}/perts.txt   n={len(essential)}   "
          f"sha256 {out['essential_list']['sha256'][:16]}...")
    print(f"  upstream: {meta['source']}  ({meta['n_perturbations']} perturbations, "
          f"{meta['n_genes']} genes)\n")

    # ------------------------------------------------------------------ 0.1  the 0/34/30 row
    panels = {k: set(P.read_panel(k)) for k in P.BUNDLES if P.BUNDLES[k].exists()}
    rows = {}
    print(f"{'panel':<8}{'n':>5}{'∩ essential':>13}{'%':>7}   {'known answer':>13}  verdict")
    for k, panel in panels.items():
        n = len(panel & essential)
        ka = KA_ESSENTIAL_OVERLAP.get(k)
        verdict = "n/a (new)" if ka is None else ("PASS" if n == ka else "**FAIL**")
        rows[k] = dict(n_panel=len(panel), n_essential=n,
                       pct=round(100 * n / len(panel), 2), known_answer=ka, verdict=verdict)
        print(f"{k:<8}{len(panel):>5}{n:>13}{100 * n / len(panel):>7.1f}"
              f"{str(ka):>16}  {verdict}")
        if verdict == "**FAIL**":
            raise SystemExit(f"KA-1 FAILED for panel {k}: {n} != {ka}")
    out["ka1_panel_x_essential"] = rows

    # disjointness, asserted
    if "abc" in panels:
        out["disjoint"] = {k: len(panels["abc"] & v) for k, v in panels.items() if k != "abc"}
        for k, v in out["disjoint"].items():
            assert v == 0 or k == "abc", f"panel {k} is not disjoint from abc: {v} shared"

    # ------------------------------------------- 0.3  coverage and ESSENTIAL FRACTION per source
    print("\nCOVERAGE and ESSENTIAL FRACTION of each source's carrying rows")
    print(f"{'panel':<8}{'source':<17}{'carried':>9}{'of which essential':>20}{'frac':>8}"
          f"{'KA':>6}")
    cov: dict = {}
    lists = {**OCTOBER_ONLY, **PRODUCTION}
    for k, panel in panels.items():
        for name, rel in lists.items():
            have = read_list(rel)
            carried = panel & have
            ess = carried & essential
            ka = KA_COVERAGE.get((k, name))
            ok = "" if ka is None else ("PASS" if len(carried) == ka else "FAIL")
            if ok == "FAIL":
                raise SystemExit(f"KA coverage FAILED {k}/{name}: {len(carried)} != {ka}")
            cov[f"{k}|{name}"] = dict(
                carried=len(carried), essential=len(ess),
                frac=(len(ess) / len(carried)) if carried else None,
                known_answer=ka, verdict=ok)
            frac = f"{100 * len(ess) / len(carried):7.1f}%" if carried else "      —"
            print(f"{k:<8}{name:<17}{len(carried):>9}{len(ess):>20}{frac:>8}{ok:>6}")
        c = panel & cd4_complete()
        cov[f"{k}|cd4"] = dict(carried=len(c), essential=len(c & essential),
                               frac=(len(c & essential) / len(c)) if c else None,
                               known_answer=None, verdict="")
        print(f"{k:<8}{'cd4 (complete)':<17}{len(c):>9}{len(c & essential):>20}"
              f"{100 * len(c & essential) / max(len(c), 1):>7.1f}%")
    out["coverage_and_essential_fraction"] = cov

    # ------- the list-free version of the same statement, and the robustness check that matters
    # The exclusion list above is ONE operationalisation of "essential".  The statement that
    # does not depend on it: the four October sources ARE common-essential libraries (Nadig
    # 2025 screened DepMap common essentials in two lines; Replogle's `K562_essential` and
    # `rpe1` are the essential-gene libraries), so the union of their own perturbation lists is
    # the ENTIRE set of targets the October asset can ever carry.  On Arc's one observed panel
    # that union covers ZERO of 300 -- measured, with no list and no assumption.
    union_oct: set[str] = set()
    for rel in OCTOBER_ONLY.values():
        union_oct |= read_list(rel)
    print(f"\nOCTOBER-SOURCE UNION (the asset's whole carrying capacity): {len(union_oct)} "
          f"targets, {len(union_oct & essential)} of them on the essential list "
          f"({100 * len(union_oct & essential) / len(union_oct):.1f}%)")
    uo = {}
    for k, panel in panels.items():
        uo[k] = dict(carried=len(panel & union_oct),
                     essential=len(panel & union_oct & essential))
        print(f"  |{k:5} ∩ union| = {uo[k]['carried']:>3} of 300"
              f"   ({uo[k]['essential']} of those are on the essential list)")
    out["october_union"] = dict(n=len(union_oct),
                                n_essential=len(union_oct & essential), per_panel=uo)

    # ------------------------------------------------------- the eligible pools for Phase 1
    genes = set(P.read_genes("abc"))
    current = set(P.read_panel("abc"))
    pools = {}
    union_cov: set[str] = set()
    for rel in PRODUCTION.values():
        union_cov |= read_list(rel)
    inter_cov: set[str] | None = None
    for rel in PRODUCTION.values():
        s = read_list(rel)
        inter_cov = s if inter_cov is None else (inter_cov & s)
    for label, cover in (("pess_pool(any donor)", union_cov),
                         ("real_pool(all 3 donors)", inter_cov)):
        elig = (genes & cover) - current
        pools[label] = dict(eligible=len(elig), essential_in_pool=len(elig & essential),
                            eligible_essential_free=len(elig - essential))
        print(f"\n{label:<26} eligible {len(elig):>6}   essential {len(elig & essential):>5}"
              f"   essential-free {len(elig - essential):>6}")
    out["eligible_pools"] = pools

    (WORK / "census238.json").write_text(json.dumps(out, indent=2))
    print(f"\n-> {WORK / 'census238.json'}")


if __name__ == "__main__":
    main()
