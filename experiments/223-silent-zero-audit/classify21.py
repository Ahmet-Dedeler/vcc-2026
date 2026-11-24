"""223 Phase 1c — the 21-row classification, assembled from the computed coverage, the
mechanical path scan, and each experiment's OWN recorded result artifact.

Every `reported` number below is quoted from the experiment's own RESULTS.md / result JSON
and the file it came from is named, so the table is checkable rather than remembered.

CLASSES
  VOID     the arm went through a MEDIAN-OVER-ALL-ROWS sizing rule on a <=50%-coverage
           source, so its dose was exactly 0 and its endpoints were bit-identical to its
           baseline.  Its NO-GO tested nothing.
  VALID    a different sizing path, or coverage > 50%.  The NO-GO stands as measured.
  NOT-RUN  no value was ever read (acquisition/metadata only, or a fail-closed gate that
           stopped before truth).  There is no verdict to void.
  UNKNOWN  cannot be determined from what is on disk.

A SEPARATE, weaker flag is carried for arms whose reported delta was EXACTLY ZERO for a
reason other than the sizing bug -- a selector that chose the identity arm.  Those are
`SELECTOR-ZERO`: a real optimiser decision with live machinery, verifiable because the
NON-identity arms in the same table are nonzero.

    uv run python classify21.py
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"

COV = json.loads((WORK / "coverage21.json").read_text())
SCAN = {r["exp"]: r for r in json.loads((WORK / "scan_paths.json").read_text())["mandate"]}

R = "RESULTS.md"

ROWS = [
    dict(exp="055-kolf-direct", source="kolf2",
         path="direct donor BLEND at an explicit fraction 0.10; sizing = median row norm "
              "over the 240 SHARED targets, every row of which carries",
         path_class="explicit-fraction blend (P3-lineage rescale, non-vacuous)",
         reported="held HCT116 -0.001151, held HEK293T +0.002127, frozen K562 +0.003992, "
                  "three-fold macro +0.001656 raw PDS",
         src=R, exact_zero=False, verdict="weak positive, advanced to 057",
         cls="VALID"),
    dict(exp="057-kolf-union-whiten", source="kolf2",
         path="same blend under union whitening, gamma swept 0.00-0.85",
         path_class="explicit-fraction blend",
         reported="gamma 0.65: -0.006869 raw PDS on the frozen K562 holdout; negative at "
                  "every nonzero gamma; unwhitened checkpoint reproduces 055 bit-for-bit",
         src=R, exact_zero=False, verdict="NO-GO (reject KOLF under whitening)",
         cls="VALID"),
    dict(exp="059-feng-ipsc", source="feng_ipsc",
         path="blend at fraction 0.02 on self-significant targets; the scale is the ratio "
              "of median row norms taken over `active_rows` and `feature_mask` -- i.e. "
              "ALREADY the carrying-row rule -- and it RAISES if the median is zero "
              "(`evaluate_lodo.py:377-384`, ValueError 'Feng median effect norm is zero')",
         path_class="carrying-row rule, with a loud zero guard",
         reported="three-fold macro +0.000069837 raw PDS (HCT116 +0.0000150, HEK293T "
                  "+0.0000748, frozen K562 +0.000120)",
         src=R, exact_zero=False,
         verdict="NO-GO (a few rank steps; below resolution)", cls="VALID"),
    dict(exp="079-kolf-local-centering", source="kolf2 is the HELD TRUTH, not the source",
         path="075's reciprocal local centering applied to the v9 effect; KOLF supplies the "
              "truth and the 273-target gallery",
         path_class="not a residual arm at all",
         reported="raw PDS 0.6064291101 -> 0.6067253825, delta +0.0002962724; bootstrap "
                  "[-0.0022893773, +0.0028549881]",
         src=R, exact_zero=False, verdict="NO-GO (fails lower-bound and +0.0015 gates)",
         cls="VALID"),
    dict(exp="085-kolf-pm-gate", source="kolf2 (never opened)",
         path="a PerturbMatch ROUTING rule trained on 064's HCT116/HEK293T outcomes; the "
              "mandatory `always_D` identity fallback was selected, so KOLF was never read",
         path_class="router with an identity fallback",
         reported="always_D 0/0/0 by construction; best non-identity `high:q=0.20` "
                  "+0.000432 / -0.001255 / macro -0.000411; permutation p=1.0",
         src=R, exact_zero=True,
         verdict="NO-GO at the source-only training gate", cls="VALID",
         zero_flag="SELECTOR-ZERO (identity chosen; the non-identity arms are nonzero, so "
                   "the machinery was live)"),
    dict(exp="105-telohaec-residual", source="telohaec",
         path="TeloHAEC as an independent SOURCE in a 50:50 equal blend, on a frozen "
              "300-target NON-panel gallery shared by all four sources; the median row-norm "
              "scale in `source_gate.py:164` is taken over the source's own 2,167 rows, all "
              "of which carry -- it is not vacuous",
         path_class="equal-weight source blend on a non-panel gallery",
         reported="TeloHAEC alone -0.14105 (HCT116) / -0.12654 (HEK293T) raw PDS vs K562; "
                  "equal blend -0.01825 / -0.01731; macro -0.01778; all four gates failed",
         src=R, exact_zero=False,
         verdict="NO-GO; and it explicitly closed the PREREGISTERED residual adapter "
                 "WITHOUT running it", cls="VALID",
         note="THE MANDATE'S FIRST-CHECK ROW. The name says `residual` but no residual was "
              "ever sized: 105 measured direct transfer and then closed the residual arm by "
              "inference. 223 Phase 2 measures the residual arm it never ran."),
    dict(exp="107-jurkat-high-moi", source="GSE249595 high-MOI Jurkat, jointly deconvolved",
         path="source effect added at an explicit blend fraction swept over 0/0.05/0.1/0.2, "
              "selected by nested held-gallery PDS",
         path_class="explicit-fraction blend",
         reported="fraction 0.1 selected on 6/6 galleries; nested macro raw-PDS gain "
                  "+0.004246, bootstrap [+0.001273, +0.006011], above all 99 permutations; "
                  "held PDS +0.001681 (h1) / +0.000116 (HEK293T); directional credit "
                  "-0.000165 / -0.000103",
         src="work/strict-evaluation-results-v2/evaluation-summary.json", exact_zero=False,
         verdict="NO-GO on two of twelve checks: both directional-credit deltas <0 and one "
                 "of six galleries negative (-0.000357)", cls="VALID",
         note="TARGET coverage is 299/300 so the median rule would be live; the binding "
              "limitation is the 349-gene (1.9%) modelled RESPONSE axis, not target "
              "coverage. Its PDS gain was positive and significant on every stated test."),
    dict(exp="115-gse295214-h1", source="gse295214_h1 (Pan 2026 H1)",
         path="acquisition + compact-donor extraction only; no model, no gallery, no truth",
         path_class="source build only", reported="no value endpoint was computed",
         src="README.md", exact_zero=False, verdict="no verdict (build)", cls="NOT-RUN",
         note="the table it built is exactly the source 218 measured returning "
              "`fraction 0.00000` through the median rule, and is half of 222's union."),
    dict(exp="116-dld1-lowmoi", source="dld1_lowmoi",
         path="source-only reliability gate; no held context was ever opened",
         path_class="source-only audit",
         reported="Low1/Low2 split retrieval on the 67 VCC targets: macro PDS 0.947083, "
                  "bootstrap [0.915197, 0.971732], permutation p=0.0001; median target "
                  "cosine 0.0696 (gate 0.2) and broad-axis reliable-row fraction 0 "
                  "(gate 0.5) FAILED",
         src=R, exact_zero=False,
         verdict="NO-GO -- but on a GENE-axis coverage gate (2,536 common measured "
                 "non-panel genes = 13.91% of the axis, against a preregistered 80%), not "
                 "on any measured transfer", cls="VALID",
         note="116's own text: 'scientifically promising for a future, separately frozen "
              "method explicitly designed around sparse-response retrieval'. That is what "
              "the carrying-row residual is. 223 Phase 2 runs it."),
    dict(exp="119-altair-a549", source="altair_a549",
         path="acquisition + compact-source build only",
         path_class="source build only",
         reported="'This artifact is a source representation only ... has not been tested "
                  "for predictive value'", src=R, exact_zero=False,
         verdict="no verdict (build)", cls="NOT-RUN"),
    dict(exp="120-myllia-a549", source="myllia_a549",
         path="acquisition + compact-donor build only", path_class="source build only",
         reported="no value endpoint was computed", src="README.md", exact_zero=False,
         verdict="no verdict (build)", cls="NOT-RUN"),
    dict(exp="121-a549-cross-source-concordance", source="altair_a549 x myllia_a549",
         path="source-vs-source retrieval on a 20-target common gallery, 7,685-coordinate "
              "axis; no VCC truth", path_class="source-only concordance",
         reported="macro PDS 0.452632 (chance 0.5, delta -0.047368); derangement null p "
                  "0.875; bootstrap [0.339, 0.570]; all seven GO criteria failed",
         src=R, exact_zero=False,
         verdict="NO-GO: 'does_not_justify_a_held_context_donor_experiment'", cls="VALID",
         note="the two A549 screens do not agree with EACH OTHER on their 20 common "
              "targets. That is a source-quality result and it survives any sizing rule."),
    dict(exp="126-hcasmc-panel-confirmation", source="IGVF HCASMC",
         path="opaque acquisition; the effect asset must remain an unopened compressed "
              "binary until a new hash-bound authorization exists",
         path_class="sealed acquisition", reported="no value was read, by design",
         src="README.md + VALUE_ACCESS_BOUNDARY.v1.json", exact_zero=False,
         verdict="no verdict (sealed)", cls="NOT-RUN",
         note="223 did NOT open it. 27/300 panel overlap, from its own committed overlap "
              "list."),
    dict(exp="129-gse295214-transfer", source="gse295214_h1 as NEIGHBOURHOOD GEOMETRY only",
         path="Pan effects used to retrieve nearest training targets; 'the query target's "
              "Pan effect was never added to its prediction'. Three 5-degree rotations vs "
              "identity, selected on source-only folds",
         path_class="retrieval geometry, not a residual",
         reported="identity 0/0/0; Pan k=5 macro -0.011010, k=15 -0.016717, k=30 -0.020859 "
                  "raw PDS; all three IMPROVED signed credit (+0.0056 to +0.0136)",
         src=R, exact_zero=True,
         verdict="NO-GO at source-only selection; the 55-target test was never opened",
         cls="VALID",
         zero_flag="SELECTOR-ZERO (identity is the reference arm; the rotations are all "
                   "nonzero and large)"),
    dict(exp="130-telohaec-panel-confirmation", source="telohaec",
         path="metadata preflight only; the value gate never opened",
         path_class="locked", reported="no TeloHAEC expression value was opened",
         src=R, exact_zero=False, verdict="LOCKED, never run", cls="NOT-RUN"),
    dict(exp="131-hues8-hepatocyte-metadata", source="HUES8 ESC / DE / WTC11 hepatocyte",
         path="metadata-only feasibility audit; every object returns HTTP 403/404 to "
              "anonymous access", path_class="never acquired",
         reported="91/300 panel targets in all three lineages, >=100 post-QC cells each; "
                  "no effect matrix exists locally",
         src=R, exact_zero=False,
         verdict="feasibility GO, acquisition NO-GO (access, not coverage)", cls="NOT-RUN"),
    dict(exp="135-telohaec-atlas-confirmation", source="TeloHAEC is the TRUTH; 122's atlas "
                                                        "is the candidate",
         path="122's `equal_mean_agreement` panel atlas scored once against 65 TeloHAEC "
              "panel targets", path_class="not a residual arm at all",
         reported="signed credit +0.004066 (gate +0.005), top-20% sign precision -0.005245 "
                  "(gate +0.010); bootstrap [-0.000456, +0.008541]",
         src=R, exact_zero=False, verdict="terminal NO-GO", cls="VALID"),
    dict(exp="137-jurkat-tangent-cone", source="107's Jurkat source, as a tangent-cone cap",
         path="angular cap on the incumbent, swept over a cap grid, nested selection",
         path_class="rotation with an identity (0 degrees) member",
         reported="nested macro raw-PDS delta -0.000921659; per-fold +0.000256, -0.000256, "
                  "-0.005632, +0.001024, 0.0; signed credit +0.000429; bootstrap "
                  "[-0.004301, +0.005940]",
         src=R, exact_zero=False, verdict="terminal NO-GO", cls="VALID",
         zero_flag="one of five folds selected the 0-degree cap, hence a 0.0 in that fold"),
    dict(exp="138-dld1-replicate-sign", source="dld1_lowmoi (replicate-consistent signs)",
         path="sparse sign gate at k=200, beta=0.05; unsupported targets fall back "
              "bit-exactly to the incumbent",
         path_class="sign gate, explicit k/beta",
         reported="nested macro raw-PDS +0.0000995142, signed credit +0.0000602495; both "
                  "bootstrap intervals strictly positive and above both permutation p95s",
         src=R, exact_zero=False,
         verdict="terminal NO-GO -- real but ~30x below the +0.003 threshold", cls="VALID"),
    dict(exp="139-jurkat-donor-simplex", source="107's Jurkat landmarks as a ROUTER",
         path="landmark-routed convex weights over 132's five source contexts, constrained "
              "to their coordinatewise convex hull; identity is a mandatory member",
         path_class="router with an identity member",
         reported="identity selected in every development fold and on all 613 targets; all "
                  "deltas exactly 0.0",
         src=R, exact_zero=True, verdict="terminal NO-GO", cls="VALID",
         zero_flag="SELECTOR-ZERO. 139's zeros are the only ALL-zero table in the 21. They "
                   "are a router selecting identity, not a sizing bug: 609/613 targets were "
                   "transferable and the routing machinery ran. But this is the one row "
                   "where the zero cannot be cross-checked against a nonzero sibling arm in "
                   "the same table, so it is the weakest VALID in the set."),
    dict(exp="176-cmap-donors", source="LINCS/CMap 2020 (shRNA CGS and CRISPR-KO)",
         path="sign-accuracy weight blend against the h1-test reference DE pool; explicit "
              "weights 0/0.25/0.5/1/2", path_class="explicit-weight sign blend",
         reported="9-line mean landmark sign accuracy 0.534 vs K562's 0.562; trio + CMap "
                  "peaks at +0.005 on 867 landmark genes and nothing on 7,101 inferred; "
                  "<= +0.001 overall on the ~1/3 of targets CMap covers",
         src=R, exact_zero=False,
         verdict="NO-GO for both modalities", cls="VALID",
         note="176 gated its own delivered-space measurement behind gate 3 and gate 3 came "
              "out marginal, so `pds`/cos through the residual path was NEVER measured. "
              "223 Phase 2 measures it."),
]


def main() -> None:
    for r in ROWS:
        c = COV.get(r["exp"], {})
        r["panel_targets"] = c.get("panel_targets")
        r["panel_pct"] = c.get("panel_pct")
        r["gene_axis"] = c.get("n_genes")
        r["median_rule_would_be_vacuous"] = c.get("median_rule_vacuous")
        s = SCAN.get(r["exp"], {})
        r["calls_solve_cd4_fraction"] = "solve_cd4_fraction" in s.get("hits", {})
        r["calls_plolo_solve_fraction"] = "solve_fraction" in s.get("hits", {})
    counts: dict[str, int] = {}
    for r in ROWS:
        counts[r["cls"]] = counts.get(r["cls"], 0) + 1
    out = dict(rows=ROWS, counts=counts,
               void_count=counts.get("VOID", 0),
               headline=("NONE of the 21 experiments was measured at a dose of exactly "
                         "zero. Zero VOID." if counts.get("VOID", 0) == 0 else
                         f"{counts['VOID']} of 21 were VOID."))
    (WORK / "classify21.json").write_text(json.dumps(out, indent=2))

    hdr = (f"{'experiment':38s}{'panel':>6s}{'pct':>6s}{'vac?':>6s}"
           f"{'class':>10s}  path class")
    print(hdr)
    print("-" * 118)
    for r in ROWS:
        print(f"{r['exp']:38s}{str(r['panel_targets']):>6s}"
              f"{str(r['panel_pct']):>6s}"
              f"{('YES' if r['median_rule_would_be_vacuous'] else 'no'):>6s}"
              f"{r['cls']:>10s}  {r['path_class']}")
    print("\ncounts:", counts)
    print(out["headline"])
    print(f"-> {WORK / 'classify21.json'}")


if __name__ == "__main__":
    main()
