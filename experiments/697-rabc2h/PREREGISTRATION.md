# 697 — R-ABC2H: R-ABC2 + 707's strong-tier hook, on def_day, on A/B/C (preregistered 2026-10-17 ~06:35Z, before R-ABC2's board read)

Asked by the coordinator [2d6f46]: the likely post-R0 replacement entry (OCT22-CHECKLIST §C 6a) is R0 + 707's swap via 708. That
composition has no board read, and the hook has no def_day set-value rehearsal. This entry gives it both, the same standard R0 met.
Slot: 10-19 s2 if K is not adopted (no R-ABC3), else 10-20 s1. Tag rabc2h.

## What is built
R-ABC2's resolved recipe, which must equal 695's `recipe695.json` exactly: R-ABC + θ 0.25. Then `DEF_DAY_DUMP_HOOK=experiments/707-
strong-tier-direction/dumphook707.sh` is added by the driver. docs/october-decisions.json keeps `dump_hook` off, because R0 must not
carry it. Built by def_day.sh on 180's def-bundle with 453's operator 66610be7, exactly as R-ABC2. The hook is env-only (690g), so
def_day.sh is not edited and 690k's proof stands.

## Prediction
pred = R-ABC2's own board read (10-18 s1, entry recorded then) + 707's Δoverall vs ALLSW (**+0.0010330**, sw707 AYbAiHRjuHdil2zCjlsk).
Additivity is the null: the hook acts after the stack dump on the strong tier's direction; θ acts before it, on per-target amplitude.

## Bins on D = R-ABC2H − pred (overall, board; the hook moves the bulk, so pass 1 re-draws: re-roll sd ≈ 0.00025)
| bin | D | consequence |
|---|---|---|
| ABOVE | ≥ +0.0005 | the swap composes with the October recipe; 708's route may use it |
| **COMPOSES** | (−0.0005, +0.0005) | as ABOVE |
| BELOW | ≤ −0.0005 | the 6a route ships WITHOUT the swap (the composition does not hold on A/B/C) |
The read does not change R0 (the hook stays off for R0 either way); it only conditions 708's AMENDMENT 2 swap route.

## Checks (pre-slot)
- The resolved recipe equals R-ABC2's exactly; the only addition is the hook (the driver records it as `added_by_697`).
- The content-aware diff vs R-ABC2 (scripts/def_metadiff.py) is within {bulk_load_sha256, bulk_ctx_scale} and must move
  bulk_load_sha256 (the dump's content). def_day logs "690g: dump hook …dumphook707.sh" and "dump replaced by the hook's".
- CLIFF OK, postbuild PASS, zero-gene 0, novelty rc 0.

## Class and tries
Transfer (707's Arc-panel proxy) composed with scorer mechanics. Tries: 707's own read (1 board arm, WEAK). This entry tests composition
and the day path; it is not a second chance for 707's adoption into R0.
