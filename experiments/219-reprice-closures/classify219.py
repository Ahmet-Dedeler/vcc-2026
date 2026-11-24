"""219 Phase 1 — classify every terminal verdict into A / B / C.

The mechanical partition first, because it is checkable and it does most of the work:

  8.4's Error 1 is an error in a PRICING CONVENTION -- "price the arm on `mse` through 096's
  `claim` route and carry the DE members at zero".  That convention did not exist before
  experiment 162 (rule 13 was fitted on 2026-09-02, 162 §; AGENTS.md 11.1 rule 13).  An
  experiment that never computed a rule-13 composite has no mse-only price to correct.  So
  the correction can only reach a closure that BOTH quotes a rule-13 floor (+0.0019 /
  +0.0089) or the 096 `claim` route, AND is a direction/profile arm.

Everything else is then hand-adjudicated from the file, with the reason recorded.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
ROWS = json.loads((WORK / "enumeration219.json").read_text())
TERMINAL = {"NO-GO", "VETO", "CLOSED", "NEGATIVE"}

# ---------------------------------------------------------------------------------------
# Hand adjudication.  key -> (class, family, reason).  Only files I read line by line.
# ---------------------------------------------------------------------------------------
HAND = {
 "097-direction": ("B", "direction", "its own recommended arm SHIPPED as a13 and the board scored "
                   "it -0.0203 vs a12 (BOARD-LEDGER). A board point outranks any re-price."),
 "099-donors": ("C", "data", "structural: 3 of 7 prepared donors carry 0/300 panel targets; "
                "kolf2/CD4/feng are 218's territory"),
 "125-hierarchical-panel-posterior": ("B", "direction", "adjudicated by 203; see 203"),
 "150-dual-moment-score-level-hct116": ("B", "profile", "already re-priced by audit 198 (R1); "
                                        "not re-done here"),
 "151-relaxation-consensus-selector": ("B", "profile", "already re-priced by audit 198 (R2); "
                                       "not re-done here"),
 "154-confidence-key-screen": ("ORD", "ordering", "overturned by 186/the board; ordering -> 220"),
 "157-reach-order": ("ORD", "ordering", "ordering channel -> 220"),
 "158-target-confidence": ("B", "direction", "the DE members were measured: nmae -0.0023/-0.0080 "
                           "survives the ORACLE and re-levelling, reach EXACTLY invariant on "
                           "132/132 targets, pds exactly 0. And its cosine input is Orion-pair, "
                           "which 165 measured overstating a profile-geometry gain 12x/48x against "
                           "the board's a23 point."),
 "159-transfer-ceiling": ("B", "direction", "a ceiling measurement, not a priced arm"),
 "161-crossline-residual": ("B", "direction", "three-leg mechanism; context features annihilated "
                            "by the mandatory centring"),
 "164-callset-artefact": ("ORD", "ordering", "call-set/precision -> 220; 181 confirmed"),
 "165-head-protected-massfix": ("B", "profile", "all six priced on h1-test; the arm is measurably "
                                "WORSE than a23 and its premise refuted deterministically. Its "
                                "eps_tail sibling (+0.0013) is a FLOOR question, not a direction "
                                "one, and 165 records its value swings 7.8x with an unpinned "
                                "calibration."),
 "166-donor-acquisition": ("C", "data", "no unexploited panel-covering line exists in public data"),
 "168-sharpening": ("B", "profile", "superseded by 169 (all six) and 188 (the nmae surface)"),
 "169-sharpening-six-member": ("B", "profile", "all six priced; 188 closed the nmae surface "
                               "question AGAINST the arm"),
 "170-gamma-mse": ("B", "operator", "ratio argument; 212 re-closed it on the cosine endpoint"),
 "171-reach-prefix": ("ORD", "ordering", "-> 220"),
 "172-sign-confidence-key": ("ORD", "ordering", "learned ordering key -> 220"),
 "173-sign-repair": ("B", "direction", "all six priced; 178 re-read its pds on the board-validated "
                     "PLOLO frame at -0.0245 raw and the arm is -0.0095 overall"),
 "175-donor-pathway-match": ("B", "direction", "181 re-stated it; 178 measured the sign gain "
                             "consumed by the operator (+0.0043 -> -0.0012)"),
 "176-cmap-donors": ("C", "data", "coverage cap; 181 measured the cap GENEROUS by 40%"),
 "177-resample-emitter": ("ORD", "ordering", "emission-law ordering; 181 CONTRADICTED it on a21"),
 "178-composed-riders": ("B", "direction", "-0.005 to -0.030 on all six; the costs compose "
                         "additively and the gains sub-additively"),
 "179-hy-kt8000-build": ("ORD", "ordering", "181 corrected it to -0.009..-0.004"),
 "180-def-readiness": ("B", "amplitude", "the per-context knee; 211 re-closed it"),
 "182-pretrained-state": ("C", "data", "structurally undeployable"),
 "183-unique-boost": ("B", "direction", "it LOWERS the transfer cosine on 3 of 4 surfaces with "
                      "truth; the multiplier applies to cosine GAINS"),
 "184-tahoe-transfer": ("C", "data", "the operator is sign-preserving and absent in CRISPRi"),
 "185-rank-controlled-emission": ("ORD", "ordering", "rc_hybrid SHIPPED and won +0.0093; "
                                  "rc_hist_a is a live ordering arm -> 220"),
 "188-magnitude-calibration": ("B", "profile", "falsified on a measurement (our |h| carries "
                               "essentially no information about true |r|)"),
 "189-per-gene-reliability": ("ORD", "ordering", "reach +0.00408 +- 0.00424 -> 220; and 214's 77x "
                              "sem tightening under rc_callset is a separate reason to re-open"),
 "190-regulon-ordering": ("DEAD", "ordering", "FALSIFIER 3: a shuffled regulon selected the same "
                          "way reads HIGHER than the real one (0.23699 vs 0.23670)"),
 "192-predicted-callset": ("ORD", "ordering", "call set -> 220 / shipped as rc_callset"),
 "193-whitener-shrinkage": ("B", "operator", "all six priced (+0.000512 best of 55 maps); and an "
                            "operator-geometry arm's multiplier is 0.59-1.23, not 2.51"),
 "194-donor-reestimation": ("B", "direction", "the recoverable part is a MAGNITUDE correction and "
                            "pds_cosine is exactly scale-invariant per row; +0.0003 vs +0.0019"),
 "195-retrieval-transform": ("B", "operator", "-0.00554 pds with a 95% CI excluding zero, and the "
                             "mechanism (dies at the cell-line boundary) is measured"),
 "196-emitted-depth": ("B", "render", "its own abort condition fired (198 §6.2)"),
 "200-production-ablation": ("A", "audit", "3 of 10 rows priced on pds+mse only; and EVERY row "
                             "inherits 214's renderer mis-scaling"),
 "203-reopen-125": ("B", "direction", "cosine +1.6% with a CI including zero, and the marginal "
                    "refund rate e = -0.193 says the arm's extra L1 lands on WRONG-signed genes -- "
                    "so the DE-four propagation cannot be positive for it"),
 "204-repression-weighted-mixing": ("DEAD", "direction", "FALSIFIER 3: equal weighting is the "
                                    "measured MAXIMUM of the transfer cosine; the model's own "
                                    "predicted gain is excluded by the measurement"),
 "205-xcell-pisces": ("C", "data", "the artifacts do not exist (3 files, used_storage 0). Its "
                      "VALUE, however, re-prices 2.3x -> 218"),
 "209-orion-subspace": ("DEAD", "direction", "FALSIFIER 3: the maximiser over all 36 arms on the "
                        "cosine channel IS a permutation control"),
 "208-free-riders": ("A", "render", "rider A's floor turns on the renderer 214 later corrected "
                     "(1.84x) and on whether the emission-law toll applies to a 0.03% profile move"),
 "211-amplitude-reopened": ("B", "amplitude", "all six priced; amplitude does not move the cosine "
                            "(pds and cos are both scale-invariant), so no direction correction "
                            "applies. Its +0.00245 vs +0.0089 is a FLOOR question."),
 "212-gamma-cosine": ("B", "operator", "the ratio argument: the operator trades cosine for "
                      "retrieval at 96.2% of the board's own rate. Instrument-independent."),
 "213-transport-law": ("B", "direction", "a transport law is a per-gene sign-PRESERVING gain "
                       "(0.0000% flips on 7.0M donor entries), so it structurally cannot move the "
                       "sign-accuracy members 8.4's propagation runs through"),
 "214-cheap-positives": ("A", "render", "Arm A prices all six at +0.0054..+0.0088 and is then "
                         "vetoed by a floor defined as the RMS of the UNPRICED channels"),
 "215-leader-decomposition": ("B", "intel", "an analysis of the field, not a priced arm"),
 "187-joint-amplitude": ("B", "amplitude", "amplitude is closed; not a cosine arm"),
 "197-cosine-inversion": ("B", "instrument", "the instrument this audit uses"),
 "198-closure-reaudit": ("B", "audit", "the prior audit"),
 "181-audit-0902": ("B", "audit", "the prior audit"),
 "199-def-runbook": ("B", "ops", "not an arm"),
 "191-rc-hybrid-build": ("ORD", "ordering", "shipped"),
 "202-callset-rider": ("ORD", "ordering", "shipped as rc_callset"),
 "206-cluster-method": ("B", "intel", "not an arm"),
 "207-cd4-dose": ("B", "direction", "SHIPPED as cd4f034 and won +0.0131; the dose ladder is "
                  "216/217's territory"),
 "186-intel-0903": ("B", "intel", "not an arm"),
 "153-de-bench": ("ORD", "ordering", "instrument -> 220"),
 "155-control-covariance": ("B", "direction", "control-cell covariance; closed pre-rule-13"),
 "160-control-subspace": ("B", "direction", "'buy cosine and sell pds' at the measured board "
                          "exchange rate -- a ratio argument"),
 "162-retrodicting-validator": ("B", "instrument", "the instrument"),
 "167-target-expression": ("C", "data", "range restriction by panel construction (AGENTS 2.1)"),
 "010-lodo": ("PRE", "pre-convention", "2026-08-21; its '-0.0089' is a table delta, not rule 13's "
              "floor (hand-checked: zero real rule-13/096 hits)"),
}


def klass(r) -> tuple[str, str, str]:
    e = r["exp"]
    if e in HAND:
        return HAND[e]
    n = int(e.split("-")[0])
    # the strict test: does the file quote a rule-13 FLOOR?  (the 096 keyword alone matches
    # the ordinary English word "claim" and over-counts by ~35, verified by hand on the
    # seven files it flagged: 010 012 014 063 098 116 121 all have ZERO real hits)
    priced_r13 = r["cites_floor_0019"] or r["cites_floor_0089"]
    if n < 162 and not priced_r13:
        return ("PRE", "pre-convention",
                "closed before rule 13 existed (fitted 2026-09-02, exp 162) and never computed an "
                "mse-only composite -- there is no mse-only price for 8.4 to correct")
    if not r["is_direction_flavoured"]:
        return ("B", "non-direction", "not a direction/cosine arm")
    return ("REVIEW", "?", "not hand-adjudicated")


term = [r for r in ROWS if r["verdict"] in TERMINAL]
res = []
for r in term:
    c, fam, why = klass(r)
    res.append(dict(exp=r["exp"], verdict=r["verdict"], klass=c, family=fam, reason=why,
                    priced_members=r["priced_members"], instruments=r["instruments"],
                    cites_r13=bool(r["cites_floor_0019"] or r["cites_floor_0089"]),
                    uses_096=r["uses_096_claim"]))

print(f"terminal verdicts classified: {len(res)}")
print("by class:", dict(Counter(x['klass'] for x in res)))
print()
for c in ("A", "ORD", "DEAD", "C", "B", "PRE", "REVIEW"):
    xs = [x for x in res if x["klass"] == c]
    if not xs:
        continue
    print(f"--- class {c}  ({len(xs)}) ---")
    if c in ("A", "ORD", "DEAD", "C"):
        for x in xs:
            print(f"  {x['exp']:38s} {x['verdict']:8s} {x['reason'][:130]}")
    else:
        print("  " + ", ".join(x["exp"].split("-")[0] for x in xs))
    print()

# how many of the terminal verdicts could the convention error even reach?
reach = [x for x in res if x["cites_r13"]]
print(f"terminal verdicts that quote a rule-13 floor (+0.0019 / +0.0089): {len(reach)}")
(WORK / "classification219.json").write_text(json.dumps(res, indent=1))
print("wrote work/classification219.json")
