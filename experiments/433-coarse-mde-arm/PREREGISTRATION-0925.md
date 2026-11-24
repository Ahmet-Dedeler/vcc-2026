# 433h / 433k — the 09-25 candidates. Preregistration, 2026-09-23 11:25Z (no 09-24 or 09-25 entry exists)

Which arms go on 09-25 is decided after the 09-24 reads (417e, 433a), by the rule in `docs/NEXT-SLOT-DECISION.md`.
The readings and bins of the two candidates are fixed now.

## 433h (`e7495644…`) = 433a + zero repair −5 + 417h's re-level (A 1.0088 / B 1.05 / C 1.0736)
It is 428h with the coarse MDE. Against 433a it changes two things, each on its own channel:
- **DE channel — the zero repair** (`de_zero_repair`).
  - 433a carries 24 zero-total gated genes, with lfc ≈ −35; 433h has 0.
  - 428's price: Δnmae **0 … −0.008 raw**. It is −0.008 only if every cell sits in a reference gate; on H1 a gated
    zero cost ≈ 0.0013 on the 58-target mean.
  - fid, jac and reach move by the repair's own count noise only (|Δ| ≤ 0.001 / 0.001 / 0.01).
- **Bulk channel — the re-level** (`bulk_ctx_scale`).
  - Predicted +0.0024 overall (417h's price), almost all of it mse. pds moves by less than 0.001.
  - The cliff model behind it is validated by 417a's mse read (k_sh 0.65).
- Paired reading vs 433a:
  - nmae, fid, reach and jac read the repair;
  - pds and mse read the re-level.
  - The two channels are pinned against each other by ballast. If a DE member moves beyond these bins while mse
    reads the re-level cleanly, the repair is doing something unmodelled.
- **Bins:**
  - **nmae:** Δ ≤ −0.002 means the repair mattered (zeros were in the gates); −0.002 < Δ ≤ +0.001 means neutral (it
    stays: it cannot hurt); Δ > +0.001 means stop and diagnose.
  - **mse:** Δ ≤ −0.006 raw means the re-level transfers; Δ > −0.002 means no transfer (re-read the cliff model).
  - **Overall:** +0.001 … +0.005.

## 433k (`f1dc369a…`) = 433h with q 1.25 → 1.5 and 417k's own re-level (A 1.0464 / B 1.092 / C 1.114)
It is 428k with the coarse MDE. Against 433h the DE channel is identical in configuration; only the bulk moves.
- 420's price for the levelled q step: about **+0.004 overall at k_sh 0.55**, mostly pds, plus mse.
  417a's measured transfer puts k_sh between 0.42 and 0.65.
- **Bins (433k − 433h):**
  - **pds:** Δ ≥ +0.004 raw means the dose step transfers; Δ < +0.001 means flat.
  - **mse:** Δ ≤ −0.010 raw means transfer; Δ > 0 means the dose overshoots A's cliff.
  - **DE members:** |Δ| within the ballast's realisation noise (nmae ≤ 0.001, reach ≤ 0.01).
  - **Overall:** +0.002 … +0.006.
- If 433k runs on 09-25 without 433h, read it against 433a. It then carries three changes: the repair, q, and the
  re-level. Say so in the result file.
