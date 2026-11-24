FOLLOW-UP from the lab lead (2026-10-16 02:55Z). Your context is intact. Excellent work on the first pass; the dump708
crash you found is reported to the coordinator. The lab lead now prices the known answers himself (experiments/804-known-
answers/ka804.py runs your engine on EC2; read it so you know what it needs). One addition is needed from you, for the
prospective row 710 (H1 dose: `VCC_H1_UNION_G` 0.70 = H1 ×2.5 vs ALLSW's ×2, board read 10-18 s1; the prediction must be
committed before its submission ~10-18 00:01Z):

1. **Compositions per fold**, on EC2 (operators from S3; charter caps; budget $15):
   - L1–L6 × 3 cancer folds, D/E/F rules: `h2.5c1` beside the incumbent `h2c1` (dump801's `--designs h2.5c1`, i.e. dump708's
     own extra-design path; a full re-dump or a light designs-only pass, whichever is exact).
   - A/B/C × 4 folds, A/B/C rules (ALLSW's stack is CD4 ×1.5): `h2.5c1.5` beside `h2c1.5`, through 438's `abc` binding
     (`FC.final`). KA: the `h2c1.5` rows equal the fold's `fin` (518's `fST`) used by engine801's AbcFold (|d| < 1e-12); if
     518's fST is not h2c1.5, say exactly what it is and stop.
   Deliver as `work/designs/<fold with / as __>.npz` with keys `targets`, `fin_<design>` (rows on the official axis, panel
   columns zero), plus a json with sha256s and the KA. Stage them to the S3 mirror (input/lab/repo/...) the way you staged
   the dumps.
2. Rules update (binding): locally you may only read files and run sub-second < 100 MB aggregations; anything else through
   `experiments/602-council-gate/gate602.sh` or on EC2; and never while `bash experiments/800-checker/window800.sh` prints
   WINDOW (daily 23:30–01:45Z; from the D/E/F release until R0 is submitted). Write `experiments/801-lab-folds/.gitignore`
   covering work/ artefacts (npz, h5ad, logs > 1 MB) if you have not.
3. Write what you did into RESULTS.md and finish with a short message.
