FROM THE LAB LEAD (2026-10-17 22:35Z by the clock): a build, needed before 10-22. Budget ≤ $2, EC2 only (local only seconds-long).
Deadline 10-21 18:00Z (the release window starts 10-21 20:00Z).

Read `docs/LAB-DEF-RUNBOOK.md` §5 and its **Amendment 1** (just committed). The coordinator's point: 6a allows ONE replacement
attempt, the likely one is R0 + 707's swap, so your c0_l25 must be priced as an INCREMENT on top of what would ship.

**Task.**
1. A new script (say `composite811.py`; m811.py stays byte-identical) that, on any engine801 fold, builds the four arms of
   Amendment 1 by calling the existing constructions unchanged:
   - B0 = inc + 705's θ 0.25 rider;
   - B1 = B0 + 707's swap;
   - C0 = B0 + your c0_l25 blend;
   - C1 = B1 + your c0_l25 blend (swap first, then 811).
   Take 705's and 707's constructions from `experiments/804-known-answers/ka804.py` / `arm702.py` (import; no copy-edit of their
   arithmetic). Your blend comes from m811's functions.
2. **Known answers** (abort on any failure):
   - B0 and B1 reproduce ka804's stored 705 and 707 arm arrays on an L fold and on an ABCp fold (|d| < 1e-12);
   - C0 at λ 0 == B0, and C1 at λ 0 == B1, bit-exact;
   - C0 with the θ rider off == your stored c0_l25 arrays on L1/k562_gwps.
3. **Information read** (counted in wear; decides nothing) of C1 − B1, C1 − B0 and C0 − B0 on ABCp and LABC, through labcloud.sh.
4. Make the kit's 811 job (`experiments/801-lab-folds/kit/run-kit811.sh` and friends) able to run the composite arms on a <TAG>.
   Put the changed kit files in YOUR dir as copies (e.g. `kit/run-kit811c.sh`), not by editing 801's files, and dry-run on LABCk.
5. Report in RESULTS.md (a new section): the script's sha, the KAs, the information read, and the exact 10-22 command. Then stop.
