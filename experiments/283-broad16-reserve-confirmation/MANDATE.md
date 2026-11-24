# 283 — The remaining preregistered confirmation: does Broad-16's reproducible line-specific component (F5, 268) hold on the four RESERVED lines? One access, one claim, logged. ≤ 3 hours. Analysis only. No build, no submission.

You are a research subagent. Work in `experiments/283-broad16-reserve-confirmation/`. **Edit nothing outside it** except
appending the claim and its result to `docs/CONFIRMATION-LEDGER.md` as instructed below. No build, no `vcc`, nothing
under `data/`; the reserved lines' single-cell files may be fetched into `work/` (they were deliberately NOT downloaded
by 268 — HTTP range requests, CRC32-verified, as 268 did for the working lines; ≤ 4 GB). `uv run python`. Peak RSS under
~6 GB. **Time box: 3 hours.** `docs/RESEARCH-PROGRAMME-2026-09.md` rule 2 binds you: the reserve is accessed ONCE, for
this claim only, after the claim is written.

## Why
268 measured, on 12 working Broad-16 lines, a line-specific response component that is reproducible across cell halves
and guides (11.5 % of response energy, 4× the null, 79 % surviving across guides), high-rank (27 components for half its
energy) and not basal-state-predictable under leave-one-line-out (F5 in `docs/FINDINGS-2026-09.md`). Four lines —
`JVE127`, `MG63`, `SLR23`, `UMRC3` — were reserved by a fixed rule before any residual was computed
(`experiments/268-broad16-context-residuals/work/reserve.json`, sha256 `1ce639b509d353708feded808914c63233bcee03513ace1c3f47cb93d071b9ce`),
never fitted or selected on, and every reserved line's batch also contains a working line. **The confirmation the
programme owes F5: on the reserved lines, is the line-specific residual reproducible at the size 268 measured, and does
the 12-line shared response predict the reserved lines' responses better than the reserved lines' own residuals
would suggest is possible from basal state?**

## READ FIRST
`docs/RESEARCH-PROGRAMME-2026-09.md` (rule 2), `docs/CONFIRMATION-LEDGER.md` (claim C1's format: claim, data, statistic,
pass rule frozen BEFORE the access, then the result whatever it is), 268's `RESULTS.md` §0–§4 and its scripts
(`load268.py` asserts the reserved columns unread — you will need to lift that assertion INSIDE this experiment's
copy only; `decomp268.py` for the variance decomposition; the split-half unit it used; the LOLO protocol in
`lolo268.py`), `experiments/268-broad16-context-residuals/work/` (the working-line residual map, the shared response,
the basal programmes fitted on the 12 lines).

## Protocol
1. **Write the claim first**, in `PREREGISTRATION.md` (sha256, timestamp) AND append it to `docs/CONFIRMATION-LEDGER.md`
   as **claim C2** in C1's format, BEFORE touching any reserved-line data:
   - C2a (reproducibility): on each reserved line, the split-half cosine of the line-specific residual (response minus
     the 12-line shared response, per target, halves by the same unit 268 used) has a median above 268's OR-control
     null (0.028) with the paired-target bootstrap interval above it, on ≥ 3 of the 4 reserved lines.
   - C2b (the shared response transfers): the 12-line shared response predicts each reserved line's full response
     at a per-target cosine whose median lies within 268's working-line leave-one-line-out range (state that range
     from 268's tables before looking).
   - C2c (basal state does not predict the residual): 268's basal-conditioned model (c), fitted on the 12 working lines
     only, predicts the reserved lines' reproducible residuals no better than the shared-response-only baseline (a) —
     the same null 268 found; state the pass/fail rule for "no better" (difference interval covering 0).
2. **Then the single access**: fetch the four lines' files, build their responses exactly as 268 built the working
   lines' (same pipeline, same gene axis, same centring fitted on the WORKING lines), compute C2a–C2c, and stop
   touching the reserve.
3. **Append the result** to the ledger under C2 (pass / fail per sub-claim, the numbers), and to `RESULTS.md`.

## Deliverable
`PREREGISTRATION.md`, the ledger entries (claim before, result after), `RESULTS.md` with the three sub-claims, their
numbers and intervals, the KA row (the working-line pipeline re-run on one working line reproduces 268's row at 0.0e+00
before the reserve is touched), and one line: **C2a PASS/FAIL, C2b PASS/FAIL, C2c PASS/FAIL**. No model, no
submission, no further use of the reserved lines. Then stop.
