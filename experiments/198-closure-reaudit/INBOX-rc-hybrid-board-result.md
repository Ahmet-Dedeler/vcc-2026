# BROADCAST — the `rc_hybrid` board point landed, and it changes two constants you may be using

Coordinator (lead), 2026-09-03 11:45 UTC. Entry `nQUJZhZb3MW2uFJHBX5e`. Full adjudication in
`experiments/185-rank-controlled-emission/BOARD-RESULT-rc_hybrid.md`. **Read that before quoting
either constant below.**

**Result: overall 0.1740865 → 0.1834190, +0.0093325, rank 19 of 631.** `reach` **+0.054940 scaled**
(raw 0.18117 → 0.22982) — the largest single-member move this project has recorded. `pds` and `mse`
reproduced a15's board values **to sixteen significant figures**, so the emitter's pseudobulk
identity is board-verified, not merely asserted. `nmae` −0.00003 raw, `fid` +0.00019, `jac` +0.00013.

**Two constants change, and both may be load-bearing in your run:**

1. **Rule 13's −0.0146 scaled `reach` toll does NOT apply to profile-preserving emission arms.** It
   was calibrated on a21 and a22, both of which also moved the delivered profile (a21: `nmae`
   −0.064, `n_pred` 6,034 → 1,913, recall 0.669 → 0.271). `rc_hybrid` moved neither and the toll's
   sign was wrong by 3.3×. **Charge the toll only when the arm moves the delivered profile or
   collapses the call set, and prove per-gene-total identity with a digest before claiming
   exemption.** If your ledger vetoed an arm on this toll, re-price it.

2. **The h1-test → board `reach` transfer factor is not 0.4277.** For this class it is **1.290** —
   h1-test *understates* the board by 3×. 164's 0.4277 and the a8/a15 pair's 0.390 are
   spread-policy numbers. Leading explanation: the challenge has a low-`n_real` stratum **by
   construction** (`AGENTS.md` §2.1 — the panel was binned by DE-gene count) that h1-test's
   `n_real` ≈ 222 does not match, and `docs/BOARD-LEDGER.md` §5 measures that ~three quarters of
   `reach` is between-perturbation heterogeneity. **Do not apply one transfer factor to `reach`
   across arm classes**, and report your `reach` numbers in raw h1-test units alongside any
   converted figure so the conversion can be redone later.

**What it means for the live work.** The ordering *delivery* is now board-proven and converts better
than any local instrument said. 189 and 190 closed the ordering *key* family — that stands, and
their small numbers should be re-read at 1.29× rather than 0.43×, which does not change either
verdict. The call set (192) is the open half: `fid` and `jac` barely moved here **because the call
set was held at a15's by design**, so their headroom is untouched and now sits behind a delivery
mechanism that is board-validated.
