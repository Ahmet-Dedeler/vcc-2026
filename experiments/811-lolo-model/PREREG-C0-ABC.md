# PREREGISTRATION: C0 on A/B/C, the board read of 811's increment WITHOUT 707's swap (the swap-OFF branch)

Written 2026-10-20 02:38Z by the clock. This is BEFORE any C0 number exists: no C0 dump has been written, scored or read. The coordinator
[3c7a7e] offered this read on 10-20 (Amendment 4 of the SIXA decision), to cover SIXA rule 3's swap-OFF branch: if 708 turns the
swap OFF, X2 = R0 + 811 = C0, and C0 − B0 then has fold reads but no board read of its own.

## What is read
- **C0** = R-ABC2's recipe (R0, unswapped stack) plus 811's c0_l25 many-line blend (λ 0.25, day bank list, centring off).
- Built from the lab's UNAPPLIED `experiments/811-lolo-model/dumphook811c0.sh`, which is `hook811.py` (unchanged since C1, sha
  c5a5614a…) applied straight to the stack dump, with no swap707 step. Readiness builds it on def_day as 699, the way 698 built C1.
- Compared against R-ABC2's board read **0.2739849**.

## The prediction (fixed now)
- **Increment source:** 811's information read of C0 − B0 on ABCp, **+0.00269 [+0.00005, +0.00564] 3/3** (composite811 6ec25788…,
  `811/RESULTS.md` §000000). LABC +0.00235 3/3 is reported beside it.
- **ρ (fixed now) = 1.1280**, by the lab's calibration rule (`800-checker/DESIGN.md` §7: a class's ρ is the median of its
  board-pointed known answers), applied to the two direction-class board points on ABCp that now exist:
  - 707: checker +0.00178 ÷ board +0.0010330 = 1.7232;
  - **C1's 811 increment** (10-20 s1): checker C1 − B1 +0.00153 ÷ board +0.0028713 = **0.5329**.
  - **Why the rule, not a single point:** C1 is the only board point for 811 itself, and it says ρ ≈ 0.53. 707 says 1.72. Picking
    either one after seeing C1 would be a choice made with the answer in view. The class median is the lab's standing rule, and
    the two alternatives are reported here so the choice is visible:
    - ρ 0.533 (C1 only) → +0.00505;
    - ρ 1.72 (707 only, the C1 prereg's deflation) → +0.00156.
- **pred (point) = 0.2739849 + 0.00269 ÷ 1.1280 = 0.2739849 + 0.002385 = 0.2763696**
  - pred range over the two class points: 0.2755459 … 0.2790331;
  - the increment's CI ÷ ρ = [+0.00004, +0.00500] (information only).
- **Bins, on D = C0's board overall − pred:**
  - |D| ≤ 0.0005 → AS PREDICTED: C0 in [0.2758696, 0.2768696];
  - D > +0.0005 → ABOVE;
  - **HARM (fixed by the coordinator):** D ≤ −0.0005, i.e. **C0 ≤ 0.2758696**, OR **C0 < 0.2734849** (R-ABC2 − 0.0005).
- **HARM closes X2 on the swap-OFF branch only.** The swap-ON branch (C1, read ABOVE) and the DEFo / DEF reads of the runbook are
  untouched.
- **Noise note:** as for C1, the increment is bulk-only. Through ballast K, the DE members carry the same-law re-roll sd ≈ 0.00025
  overall, inside the ±0.0005 bins.
- **Expectation, recorded:** given C1 (the folds under-read 811 by ~1.9× on top of the swap), C0 may read ABOVE this pred. The
  class-median ρ is deliberately not tuned to C1 alone.

## Read (2026-10-20 08:18Z by the clock; BOARD-LEDGER, entry `4rq3qHIuv8bHAxpTkUeu`, 10-20 s2)
- **C0 board overall: 0.2780624. D = C0 − pred = +0.0016928 → ABOVE.** Not HARM (> 0.2758696 and > 0.2734849).
- C0 − R-ABC2 = **+0.0040775**: bulk +0.003812 (Δmse +0.02314 scaled, Δpds −0.00027), DE +0.000266 (~1 re-roll sd).
- Implied ρ for this increment = ABCp +0.00269 ÷ board +0.0040775 = **0.660**. As with C1 (0.533), the folds under-read 811 on
  the board.
- **The interaction on A/B/C** (4 board points on def_day): swap alone +0.00091; 811 alone +0.00408; both +0.00378. C1 − C0 =
  −0.0003 (about −1.2 re-roll sd): once 811 is in, 707's swap adds nothing on A/B/C. The folds said the same thing in
  direction: the swap absorbed ~45 % of 811's gain on ABCp.
- Both X2 branches (swap ON = C1, swap OFF = C0) now have an A/B/C board read, both ABOVE. 811's 6a decision itself stays the
  DEFo read of the runbook (§5, Amendments 1 and 3).
