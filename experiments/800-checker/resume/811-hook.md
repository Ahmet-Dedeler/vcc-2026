FROM THE LAB LEAD (2026-10-18 08:32Z by the clock): a production-format build with a HARD deadline. Ready by **10-19 10:00Z**; the
coordinator's deadline is 12:00Z. If it cannot be clean by then, write "NOT READY" with the reason in RESULTS.md and stop. A
rushed hook is worse than none. Budget ≤ $3, EC2 for anything beyond seconds.

**Offer (coordinator, 10-18 08:2xZ).** Production would spend 10-20 s2 or 10-21 s1 on a board read of **C1 on A/B/C** =
R-ABC2H's recipe (R0 + 707's swap via dumphook707) + your c0_l25, built by readiness on def_day. That read is the strongest
Arc-drawn check of your 6a rule.

**Deliverables, all in YOUR dir, UNAPPLIED** (never edit production files: `experiments/707-*`, scripts/, def_day, the builder):
1. **`dumphook811.sh`**, same interface as `experiments/707-strong-tier-direction/dumphook707.sh`
   (`<bundle> <work> <dump.npz> <identity.json> <hookdir>`; 690g's DEF_DAY_DUMP_HOOK). It runs dumphook707's swap FIRST
   (call `swap707.py` unchanged) and then your c0_l25 transform on the swapped dump's bulk rows.
   - Production's dump rows are raw, unwhitened non-panel delta rows. Read `swap707.py` to see exactly what it changes.
   - Your transform must be the SAME construction composite811's C1 applies on engine folds: the many-line consensus over the
     day bank list, λ 0.25, centring off, norm-preserving if that is how c0 works. Write the mapping between the engine's fold
     space and the dump's raw space explicitly, with any approximation.
   - Carriers, perts, genes and env are copied unchanged.
2. **Known answers** (abort on any failure):
   - KA-off: with λ 0 (an env or flag), the hook's output == dumphook707's output **byte for byte** on A/B/C's dump. Get the dump
     path from 707's KA log `experiments/707-strong-tier-direction/work/ka_hook707.log`. Read-only, and clone the dump into your
     work dir with `cp -c`.
   - KA-same-construction: on an engine fold where both can be computed (e.g. ABCp/k562), the dump-space transform and
     composite811's C1 − B1 change the same targets in the same direction, and agree to stated precision. If they cannot be
     made identical, say why and how much they differ.
   - Report the bank coverage of A/B/C's 300 targets (lines per target) and that the H1 group is in the bank as on the day.
     There is no held-out exclusion on the board.
3. **The prediction input:** your ABCp / LABC C1 − B1 prices (already in RESULTS), stated in board overall units, plus the
   deflation the lab uses (ABCp direction ρ 1.72, range [1.44, 1.72]). The lead writes the preregistration.
4. RESULTS.md section "C1 board-read hook": shas, the KAs, the exact env for readiness (`DEF_DAY_DUMP_HOOK=<path>`), runtime,
   and NOT READY if anything is unclean. Then stop.
