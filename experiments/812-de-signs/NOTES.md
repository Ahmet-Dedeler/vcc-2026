# 812 working notes (state for a resume; RESULTS.md is the record)

## Done (2026-10-16)
- `surf812.py` loader (h1test = 435 zs stats; h1tv = 618 g1 zs stats; key/head/block as lib521). KA: fid from stats ==
  recorded raw fid to 1e-12 (h1test s0 0.5089291277653674, h1tv s0 0.5052857989716864), all 6 seed rows pass.
- `a0class812.py` -> work/a0class812.json (step 1: a0 by class, oracle buys). `pool812.py` -> work/pool812.json
  (CLEAN4 accuracy by CPM x confidence on all non-sig genes; confidence-composed block emulation: fid up, jac down, net ~0).
- `sigcpm812.py` -> work/sigcpm812.json (real-sig genes by CPM: no Wilcoxon-sparsity UP lean worth anything).
- H-NTCB (new mechanism): NTC-guide-specific effects bias the scorer's control mean -> shared sign lean in every
  target's truth. `ntcbias812.py` (H1, guide ids from VCC2025 obs -> work/stage/h1ctrl_guides.tsv sha 51102f53...),
  `ntcvar812.py` (A/B/C between-guide overdispersion). EC2 run 20261016T074407Z-de812-b1 (r7i.xlarge, 40 min cap, $0.23).

- b1 (ntcvar A/B/C: guide R ~ null in B/C, slight excess in A), b2/b3 (ntcbias H1: 5 flagged genes, M0 acc 0.69/0.66
  vs CLEAN4 0.43/0.44, fid +0.0003 raw = +0.00017 overall: real, tiny), b4 (cis812: reach +0.00014 h1test / +0.0021 h1tv
  raw, <= +0.0004 overall), b5 (rates812 exchange rates, head812 head-size composition) running.
- 440 already closed the empty-union jac lever (A/B/C has no null targets; board control copies read jac 0).

## Next
- read b1; if H-NTCB shows accuracy > 0.5 rising with |z| on H1 -> PROPOSAL.md with prediction + falsifier; else record
  the negative.

## Final (08:15Z)
All six EC2 runs terminated and costed ($0.40). RESULTS.md and PROPOSAL.md written: NO-GO on H-NTCB, H-CBLOCK,
H-CISFRONT and H-HEAD; the exchange rates (rates812) are handed to lane 811 as its DE objective. Nothing is left running.

## Follow-up (08:33Z)
skill812 built + KA PASS; step-3 sources run (b7, $0.03). PROPOSAL-u1hook.md (u1 hook from the public per-target mean; worker prior NO-GO, stage 0 held-out read asked). Nothing running.
