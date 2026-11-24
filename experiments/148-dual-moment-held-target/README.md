# Experiment 148 — held-target confirmation of the dual-moment emitter

Status: **preregistered; no confirmation truth opened**

Date: **2026-08-31**

## Question

Does Experiment 147's fixed lambda-1.00 dual-moment emitter transfer from the
40 H1-fast development targets to the other 45 H1-test targets?

This experiment does not tune any emitter parameter. It carries forward:

- lambda `1.00`;
- depth contrast `0.50`;
- within-group quasi spread `0.25`;
- base seed `14700`;
- incumbent pooled columns from gamma `0.65`;
- arithmetic-CPM direction target from gamma `0.45`.

The candidate is built for all 85 H1-test targets because PDS is listwise, but
the primary confirmation cohort is the 45 targets absent from H1-fast.

## Frozen held-target cohort

`GULP1, CD276, TOR1A, PGAM1, ANKZF1, PCCB, ASF1A, MORN2, COASY,
NDUFC2, MED23, RRM2, MLLT1, NFATC3, ZBTB39, PHF12, ARID3A, CSNK1E,
CYCS, APOE, FIBP, ZNF107, CTBP2, DCK, SAFB2, PLAGL2, HDGF, UGP2,
SOX4, MFSD10, FGFR1, NDUFB9, TRAFD1, ZNF92, PIAS1, JMJD1C, STUB1,
SBDS, ZNF664, RBPJ, ACTN1, TNFRSF12A, PBX3, HAX1, NEK2`

The cohort was derived from `h1-test obs - h1-fast obs` before candidate
scoring. No H1-test expression value was inspected to choose it.

## Confirmation gate

On the 45 held targets:

- every per-target raw PDS value must be exactly equal;
- official-anchor projected six-member delta must be at least `+0.004`;
- raw NMAE delta must be at most `-0.004`;
- raw reach delta must be at least `+0.010`;
- raw fidelity delta must be at least `-0.003`;
- raw Jaccard delta must be at least `-0.010`;
- raw normalized MSE delta must be at most `+0.010`.

The full 85-target official-anchor projected average must improve by at least
`+0.005`.

Passing this experiment authorizes a separate cross-context confirmation. It
does not authorize a challenge build or submission.
