# LZ access claim: 811's c0_l25 sign check on the sealed reserve (one access)

Written 2026-10-20 09:19Z by the clock, by the lab lead. This is BEFORE any LZ bundle, CD4 bank, operator, dump or score exists, and before
any D/E/F data. Its sha256 goes to the production coordinator for `docs/CONFIRMATION-LEDGER.md` before the first LZ operator fit
(`experiments/800-checker/SEALED.md`, access rule). This is the LZ condition of `docs/LAB-DEF-RUNBOOK.md` §5 (item 5) and
Amendments 1 and 3.

## The candidate (frozen since 10-16 10:49Z)
- 811's **c0_l25**: the many-line consensus, day bank list, λ 0.25, centring off.
- `m811.py` e00192e6…, composite arms by `composite811.py` 6ec25788…, and `hook811.py` c5a5614a… where a dump is involved.
- Its WORKING prices exist, so the access rule's precondition holds:
  - BULK-W +0.0068 (held-out L4–L6);
  - A-only +0.0064, 18/18;
  - ABCp / LABC information reads;
  - and two A/B/C board reads (C1 0.2777649 ABOVE, C0 0.2780624 ABOVE).

## What is priced against what
- **Folds:** L7 and L8 (SEALED.md target shas 97c30e7d… / 32fd27d2…) × the three cancer held lines (K562 GWPS, Orion HCT116, Orion
  HEK293T) = 6 folds.
  - Built on EC2 by 801's panel kit (`experiments/801-lab-folds/kit/kit801.sh`) with A/B/C frames, the L-fold recipe, unchanged.
  - The kit refuses L7/L8 by design. The access adds ONE override keyed to this file's sha256 (its diff recorded in the
    access record), and nothing else changes.
- **Arms:** composite811's B0, B1, C0, C1 on every fold. Two increments:
  - **swap-ON branch:** C1 − B1 (811 on top of R0 + 707's swap);
  - **swap-OFF branch:** C0 − B0 (811 on R0).
  - Both are preregistered here. **The branch that decides is the one production's R0 day path ships** (708's day-of swap rule
    on D/E/F); the other is reported.
- **Statistic:** check800 at its frozen settings: per-fold overall (literal pds + mse via 518's identity), paired target bootstrap
  4,000 draws, seed 800, equal weight per fold, pooled over the 6 folds; pdsl reported. N = 1 on LZ, so there is no best-of-N
  shift.
- **Also reported, not decided on:**
  - each LZ panel's full-library essential-screen count (`814/work/flags/esslib.txt`); arc800 draws hold ~17, so the panels
    are expected to over-read key / direction changes (`docs/LAB-PANEL-REPORT-2026-10-16.md`);
  - the bootstrap CI and folds > 0.

## Reading rule (the only use of this access)
- **PASS (sign):** the deciding branch's pooled point > 0. Runbook §5 item 5 holds.
- **FAIL:** the deciding branch's pooled point ≤ 0. Item 5 fails, so there is no 6a dossier for 811, whatever DEFo reads.
- **Nothing else is read from LZ:** no dose, no λ, no arm choice, no debugging. A failed known answer (engine fold KAs,
  composite811's KAs on the LZ folds) stops the access with no score read, and that is recorded.
- Result to the ledger through the coordinator, whatever it is.

## When and cost
- After DEFo's read (runbook §5) and only if DEFo passes items 1–4. Otherwise LZ stays sealed.
- EC2, ≤ $5. Project=lab, on-demand, terminate-on-shutdown.
