#!/usr/bin/env python3
"""Read a scored entry against the shipped incumbent and apply ONE preregistration's reading table.

    uv run python scripts/read_board_entry.py <entry-id | captured.json> [--incumbent <id | captured.json>]
                                              [--prereg 237 | 249 | none]

Sources, in order: a path to a data/board-entries/<id>.json capture (scripts/capture_board_entry.py),
else `vcc --json status <id>`, else data/board-entries/<id>.json. The incumbent defaults to h1union
(iUziqE1gewfcPKzYUDaB) the same way, with its row frozen here as the last resort. The API resolves only
the CURRENT entry, so once the second 2026-09-05 slot publishes, both the 237 entry and h1union resolve
only through their captures.

--prereg 237  (ARM6, experiments/237-gate-validation): the three pre-committed overall bins, verbatim,
              PLUS the `nmae` leg read DIRECTLY off `score_nmae` against the identity's -0.004783, PLUS
              243's raw-`mse` class bins.  Before 2026-09-05 this script inferred the sign of the `nmae`
              prediction from the OVERALL delta ("Δ ≥ +0.0015 ⇒ the identity mis-signed nmae"). That does
              not follow: the other five members can outweigh an `nmae` loss (a mocked entry with overall
              +0.0020 and an `nmae` contribution of -0.0048 was declared a mis-sign). The board reports
              `score_nmae`; the identity is judged on it.
--prereg 249  (cd4svd5, experiments/249-cd4-deflate): falsifiers 1-3 (board-readable), 4-5 cited from
              RESULTS.md (package-side), and the DECISION-TREE branches with an explicit PRECEDENCE
              (E > D > A > B > C) and a HOLD branch for every outcome the tree left uncovered.
--prereg none: the six members and the generic gates only. This is the default: a candidate's
              thresholds are never applied to another candidate's entry.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CAPDIR = REPO / "data" / "board-entries"
H1_ID = "iUziqE1gewfcPKzYUDaB"
# h1union's row, frozen from experiments/222-h1-union-build/BOARD-RESULT.md and its capture.
H1 = {"entry_id": H1_ID, "status": "published", "model_name": "The Nucleus Declined To Comment", "rank": 3,
      "score_avg": 0.21183096391231482, "pds_cosine": 0.8357227796358231, "expr_mse_unbiased_capped_norm": 0.8172630319018631,
      "de_wilcoxon_lfc_nmae": 0.93389, "de_wilcoxon_direction_fidelity_yield_raw": 0.50748,
      "de_wilcoxon_direction_reach_raw": 0.27655, "de_wilcoxon_sig_jaccard": 0.04096,
      "score_pds": 0.74537, "score_mse": 0.17957, "score_nmae": 0.11074, "score_fid": -0.01648,
      "score_reach": 0.22294, "score_jac": 0.02884}
RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
       "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
       "jac": "de_wilcoxon_sig_jaccard"}


CANDIDATES = {  # what each --prereg table may be applied to; anything else is refused, not warned about
    "237": {"sha256": "7a3c2d3e414f8d9dfba0bd1c42f2b93994dedfda73cf207cea9b7aabdc78291c",
            "name": "I'll just keep moving forward.", "package": "submissions/229-stack6-g035.vcc"},
    "249": {"sha256": "6e8b41c6244b32dbe2a9765c7d2f2293b2032c49e59b455982b95bc6bb074533",
            "name": "You can do it too.", "package": "submissions/249-cd4svd5-f034.vcc"},
}


def load(ref: str, *, allow_frozen: bool = False) -> tuple[dict, dict, str]:
    """(entry, wrapper, source). A path wins; then the live API; then the capture dir; then the frozen row."""
    p = Path(ref)
    if p.exists():
        obj = json.loads(p.read_text())
        return (obj.get("entry", obj), obj if "entry" in obj else {}, f"file {p}")
    try:
        r = subprocess.run(["vcc", "--json", "status", ref], capture_output=True, text=True, timeout=120)
        if r.stdout.strip():
            obj = json.loads(r.stdout)
            if obj.get("status"):
                return obj, {}, "live API"
    except Exception:
        pass
    cap = CAPDIR / f"{ref}.json"
    if cap.exists():
        obj = json.loads(cap.read_text())
        return (obj.get("entry", obj), obj if "entry" in obj else {}, f"capture {cap.relative_to(REPO)}")
    if allow_frozen and ref == H1_ID:
        return H1, {}, "frozen row in this script"
    sys.exit(f"cannot resolve {ref}: not a file, not on the API (only the current entry resolves), no capture at {cap}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("entry")
    ap.add_argument("--incumbent", default=H1_ID)
    ap.add_argument("--prereg", choices=["237", "249", "none"], default="none")
    a = ap.parse_args()

    e, ew, esrc = load(a.entry)
    b, _bw, bsrc = load(a.incumbent, allow_frozen=True)
    print(f"entry {e.get('entry_id')} [{esrc}]: status={e.get('status')} rank={e.get('rank')} "
          f"overall={e.get('score_avg')} model={e.get('model_name')!r} panel={e.get('panel_id')} anchors={e.get('anchor_version')}")
    print(f"incumbent {b.get('entry_id')} [{bsrc}]: overall={b.get('score_avg')} panel={b.get('panel_id')}")
    if e.get("score_avg") is None:
        print(f"not scored: status={e.get('status')} error_info={e.get('error_info')}")
        return 1
    if e.get("panel_id") and b.get("panel_id") and e["panel_id"] != b["panel_id"]:
        print(f"WARNING: panels differ ({e['panel_id']} vs {b['panel_id']}); scores are not comparable across panels")

    d = e["score_avg"] - b["score_avg"]
    leg = {}
    print(f"\noverall delta vs incumbent {b['score_avg']:.7f}: {d:+.7f}")
    for k, rk in RAW.items():
        sk = "score_" + k
        leg[k] = (e[sk] - b[sk]) / 6
        print(f"  {k:6} raw {e[rk]:.5f} ({e[rk]-b[rk]:+.5f})   scaled {e[sk]:+.5f} ({e[sk]-b[sk]:+.5f})  -> overall {leg[k]:+.5f}")
    print(f"  sum of member contributions {sum(leg.values()):+.7f} (must equal the overall delta)")
    print("\ngeneric gates (AGENTS.md rule 8.1): shape-only ship floor +0.0019; amplitude-changing +0.0089; "
          "ordering noise sd 0.004502 overall between digest-identical arms; pds/mse deterministic to 16 s.f.")

    if a.prereg in CANDIDATES:
        refusals = guard(a.prereg, e, b, ew)
        if refusals:
            print(f"\nREFUSED: the {a.prereg} reading table does not apply to this entry/baseline pair:")
            for r in refusals:
                print("  -", r)
            return 2
        if a.prereg == "237":
            read_237(e, b, d, leg)
        else:
            read_249(e, b, d, leg)
    else:
        print("\n(no preregistration table applied; pass --prereg 237 or 249 for the armed candidates)")
    return 0


def guard(prereg: str, e: dict, b: dict, wrap: dict) -> list[str]:
    """A candidate's thresholds apply only to that candidate's package, against h1union, on the same panel
    and anchor set. The package is identified by the capture's sha256; without one, by the joke model name."""
    c = CANDIDATES[prereg]
    out: list[str] = []
    sha = (wrap or {}).get("sha256")
    if sha and sha != c["sha256"]:
        out.append(f"capture sha256 {sha[:16]}... is not {prereg}'s package {c['package']} ({c['sha256'][:16]}...)")
    elif not sha and e.get("model_name") != c["name"]:
        out.append(f"no sha256 in the capture and model_name {e.get('model_name')!r} is not {prereg}'s {c['name']!r}")
    elif not sha:
        print(f"note: candidate identified by model_name only (no sha256 in the capture); backfill it with "
              f"scripts/capture_board_entry.py {e.get('entry_id')} --package {c['package']} --sha256 {c['sha256']}")
    if e.get("status") not in ("published", "superseded"):
        out.append(f"entry status {e.get('status')!r} is not a scored state")
    if b.get("entry_id") != H1_ID or abs(float(b.get("score_avg") or 0) - 0.2118310) > 1e-6:
        out.append(f"incumbent is {b.get('entry_id')} at {b.get('score_avg')}, not h1union {H1_ID} at 0.2118310, which the prereg deltas are against")
    for k in ("panel_id", "anchor_version", "partition"):
        if e.get(k) != b.get(k):
            out.append(f"{k} differs: {e.get(k)!r} vs incumbent {b.get(k)!r}; scores are not comparable across panels/anchors")
    return out


def read_237(e: dict, b: dict, d: float, leg: dict) -> None:
    NMAE_PRED = -0.004782935  # 229 six229.json: stack6 overall_nmae_vs_incumbent (overall units = Δscore_nmae/6)
    print("\n=== 237 / ARM6 (stack6 g=0.35 vs h1union) ===")
    print("predictions: 229 six-member composite -0.00219 (= 243's direction legs + the nmae identity leg -0.004783);"
          " 243 class rule +0.0022 ... +0.0038 (point +0.0024), nmae NOT priced")
    bins = [(-9, -0.0015, "<= -0.0015 (P237 s4 row 1)"), (-0.0015, 0.0015, "between (P237 s4 row 2)"), (0.0015, 9, ">= +0.0015 (P237 s4 row 3)")]
    for lo, hi, txt in bins:
        if lo <= d < hi:
            print(f"237 pre-committed overall bin ({d:+.5f}): {txt} -- RECORDED ONLY. Its inference about the nmae identity"
                  " is SUPERSEDED by the direct read below (ADDENDUM-2-READING s1): the overall total does not carry the nmae sign.")
    n = leg["nmae"]
    print(f"\nDIRECT nmae read: realised nmae contribution {n:+.6f} overall (scaled {e['score_nmae']-b['score_nmae']:+.5f}, "
          f"raw {e[RAW['nmae']]-b[RAW['nmae']]:+.5f}) vs identity {NMAE_PRED:+.6f} ({n/NMAE_PRED:+.0%} of it)")
    others = d - n
    print(f"  direction legs (overall minus nmae) realised {others:+.6f} vs 243's +0.0022 ... +0.0038")
    if n <= 0.8 * NMAE_PRED:
        state = "HOLDS"
        if d <= -0.0015:
            final = ("identity VALIDATED at scale on this arm class (sign and magnitude): 228 and ARM6 stay dead on the nmae leg;"
                     " the 48-day programme proceeds on 8.1d/8.1f as written")
            if abs(n - NMAE_PRED) <= 0.2 * abs(NMAE_PRED):
                final += (". P237 s4 falsifier FIRES: leg within 20% of -0.004783, so the two earlier mis-signs were small-magnitude"
                          " artefacts and this slot was mis-spent; record it as such")
        elif d >= 0.0015:
            final = ("identity HOLDS and the total still cleared +0.0015: the gain came from the other five members, so the identity"
                     " is NOT retired; 243's class rule under-priced the direction legs. 228 is re-priced with the identity intact"
                     " (its nmae leg is 2.2x this one)")
        else:
            final = ("identity HOLDS; the total is inconclusive only because the direction legs under-delivered against 243's band."
                     " Nothing is retired; 228 stays dead on the nmae leg")
    elif n < 0:
        state = "PARTIAL"
        final = (f"sign holds, magnitude refuted ({n/NMAE_PRED:.0%} realised): downgrade the identity from a veto to a DISCOUNT;"
                 " re-price 228 with the realised leg scaled by its 2.2x larger raw delta-nmae")
    else:
        state = "MISSIGNED"
        final = ("the identity mis-signed nmae a THIRD time: retire it as a veto, re-derive the leg against the board points,"
                 " re-open 228 on 243's class rule")
        if d <= -0.0015:
            final += (". The overall LOSS came from the direction legs, so what failed here is 243's class rule, not the"
                      " identity's sign; do NOT read P237 s4 row 1 ('validated') from the total")
    print(f"\nFINAL 237 VERDICT [{state}]: {final}")
    m = e[RAW["mse"]]
    mbins = [(0, 0.80768, "above every model (k >= 1.55)"),
             (0.80768, 0.81002, "DOSE class (k 1.25-1.55): the union point was the outlier"),
             (0.81002, 0.81135, "COVERAGE-linear model (k 1.08-1.25)"),
             (0.81135, 0.81392, "NEW-SOURCE class (k 0.75-1.08) -- 243's point prediction 0.81244"),
             (0.81392, 9, "below every model (k < 0.75): the union factor itself was optimistic")]
    for lo, hi, txt in mbins:
        if lo <= m < hi:
            print(f"243 raw-mse class bin ({m:.5f}, Delta raw {m-b[RAW['mse']]:+.5f}): {txt}")
    # Realised k from 243's own ARM6 terms: its price rows are (k_cos -> dcos -> raw mse) on the arm's measured
    # A_c/B_c, and raw mse is linear in dcos at those terms (197's identity), so invert the line through them.
    j = REPO / "experiments/243-transfer-factor/work/arm6_243.json"
    if not (e.get("score_mse", 0) > 0 and b.get("score_mse", 0) > 0):
        print("realised k not computed: score_mse <= 0 on a row (mse is clamped at 0; the identity inverts only rows with score_mse > 0, AGENTS.md)")
        j = None
    try:
        if j is None:
            raise ValueError("skipped")
        a6 = json.loads(j.read_text())
        rows = [(r["k_cos"], r["dcos"], r["raw"]) for r in a6["price"] if all(k in r for k in ("k_cos", "dcos", "raw"))]
        (k1, c1, r1), (k2, c2, r2) = next((x, y) for x in rows for y in rows if abs(x[1] - y[1]) > 1e-6)
        slope = (r2 - r1) / (c2 - c1)            # d raw / d cos, negative
        raw0 = r1 - slope * c1                     # raw mse at zero cosine gain (amplitude already charged)
        per_k = c1 / k1                           # dcos per unit k = rel_plolo_3c * cos_inc
        k_real = (m - raw0) / (slope * per_k)
        print(f"realised k_cos from the measured column: {k_real:.3f} (raw {m:.5f}; 243 terms raw0 {raw0:.5f}, "
              f"d raw/d cos {slope:.3f}, dcos per k {per_k:.6f}); classes: dose 1.31-1.42, coverage-linear 1.08-1.25, new-source 0.94-1.00")
    except ValueError:
        pass
    except Exception as exc:  # the reading above stands without it
        print(f"(realised k not computed: {exc})")
    print("Next: the realised k from the MEASURED cosine (experiments/243-transfer-factor/work/arm6_243.json), never the floor;"
          " local instrument caveat: six229's nmae identity uses fixed top-89 pools and pools 4 folds (KOLF2 in), so this point"
          " tests the instrument on THIS arm class, not every direction arm ever priced.")


def read_249(e: dict, b: dict, d: float, leg: dict) -> None:
    print("\n=== 249 / cd4svd5 (VCC_CD4_DEFLATE_K=5 vs h1union) ===")
    print("predictions: overall +0.0026 ... +0.0040 (upside +0.0086), pds leg +0.0032 / +0.0046 / +0.0091 at 0.7x / 1.0x / 2.0x,"
          " raw pds +0.0086 at 0.7x, raw mse 0.8190-0.8194 (R249 s9; prereg said 0.8183 at h1union's amplitude terms)")
    dp = e[RAW["pds"]] - b[RAW["pds"]]
    dm = e[RAW["mse"]] - b[RAW["mse"]]
    f1 = dp < 0.0040
    f2 = dm > 0.0030
    f3 = not (-0.0010 <= d <= 0.0100)
    print(f"falsifier 1 raw pds < +0.0040: {'FIRES' if f1 else 'no'} (raw pds {dp:+.5f}; pds leg {leg['pds']:+.5f} overall)")
    print(f"falsifier 2 raw mse worse by > +0.0030: {'FIRES' if f2 else 'no'} (raw mse {dm:+.5f}; mse leg {leg['mse']:+.5f} overall)")
    print(f"falsifier 3 overall outside [-0.0010, +0.0100]: {'FIRES' if f3 else 'no'} (overall {d:+.5f})")
    print("falsifier 4 config diff exactly cd4_deflate_k + consequences: PASS off the package (R249 F1, 115/115 pinned)")
    print("falsifier 5 delivered amplitude within 0.5%: PASS off the package, -0.409% (R249 F5) -> branch E does not fire")
    # Precedence: E (package-side, already read) > D > A > B > C > HOLD.
    if f2 or d < 0:
        br = ("D", "overall < 0 or falsifier 2: retire the deflation CLASS on A/B/C; October basis stays UNSET;"
                   " re-derive 243's factor for 're-shape occupied rows' before pricing any re-shaping arm")
    elif dp >= 0.0086 and d >= 0.0026:
        br = ("A", "pds-first class transfers at >= 0.7x. October: ship the PINNED basis (VCC_CD4_DEFLATE_BASIS=V5_abc, 251)."
                   " NEXT SLOT IS NOT AUTOMATIC: k=10 (255) is priced -0.0005 ... -0.0002 over cd4svd5 (1/3 folds) and needs"
                   " its own preregistration on the realised pds factor; +partial (254) is under the gate")
    elif 0.0040 <= dp < 0.0086 and d >= 0.0015:
        br = ("B", "transfer 0.33-0.7x. October: pin the basis. Next slot: hold; re-price k=10 and +partial on the realised"
                   " factor; ship only what clears 1.5x the gate on it")
    elif dp < 0.0040 and d >= 0:
        br = ("C", "falsifier 1: mechanism does not transfer; retire the deflation class; October basis stays UNSET;"
                   " never spend a slot restoring the incumbent")
    else:
        br = ("HOLD", "outcome not covered by A-D (e.g. pds transferred but overall under its bin): hold, reassess with all"
                      " six members, no next-slot escalation and no October decision from this point alone")
    if f3 and br[0] not in ("D", "C"):
        br = ("HOLD", f"falsifier 3 fired (overall {d:+.5f}) on an otherwise-{br[0]} outcome: the prediction band was wrong;"
                      " hold and reassess before any branch action")
    print(f"\nDECISION-TREE branch (precedence E > D > A > B > C > HOLD): {br[0]} -- {br[1]}")
    print("pairing (DECISION-TREE addendum 5 supersedes its 'If 237 lands Delta >= +0.0015' line): 228's per-gene gain becomes a 09-06"
          " candidate only if 237's FINAL verdict is MISSIGNED (or PARTIAL, re-priced with the discount) AND A/B fired here,"
          " re-priced against cd4svd5 as the incumbent.")


if __name__ == "__main__":
    raise SystemExit(main())
