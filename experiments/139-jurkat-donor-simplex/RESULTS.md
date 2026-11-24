# Experiment 139 result — terminal NO-GO

Date: **2026-08-30**

## Decision

The high-MOI Jurkat landmark router selected exact identity in every
development fold and on the complete 613-target development set. The
untouched 410-target test remains sealed and must not be opened.

## Frozen development result

- final arm: `identity`;
- nested macro raw-PDS delta: `0.0`;
- nested macro signed-credit delta: `0.0`;
- every held fold selected `identity`;
- final all-development PDS and signed-credit deltas: `0.0`;
- bootstrap intervals and null comparisons are consequently zero;
- promotion: `false`.

This rejects the frozen hypothesis that 138 high-MOI Jurkat landmark
coordinates can usefully route convex weights among Experiment 132's five
allowed source contexts at `rho={0.25,0.50}` and
`temperature={0.10,0.25}`. The result does not test direct use of the Jurkat
effect vector, which remained prohibited.

## Integrity

- 12 source-safe unit tests passed;
- Ruff and diff checks passed;
- 609/613 development targets were transferable and four used exact fallback;
- inner selection used inner-only PDS galleries;
- all candidates remained within the coordinatewise convex hull of the five
  allowed source vectors;
- no test authorization was created and no 410-target value was read.

## Artifact bindings

- prevalue authorization SHA-256:
  `665621539061cb3600ad0b205ec487b32813a117d889be17c57d54d8b608f969`;
- development result SHA-256:
  `c72ca68b3fd52792971e872959faa02f97c23c9f846235fae0f8f1a44c595049`;
- development statistics SHA-256:
  `e9dee51616a2c1dd10a5de8e5a0288daea6a3de7480318f1e36b6bf1fb034325`;
- development artifact manifest SHA-256:
  `73025af13f32f3991b433c125515e7912e1fe9ae7bf37478ba07ef00cc6a079c`;
- terminal candidate-freeze SHA-256:
  `5e41dab49d8868265f691e4f98c73039f962c05442336d16fe01a2bf5ba0217a`.

No count build, `cell-eval2` run, challenge submission, commit, or push was
created by Experiment 139.
