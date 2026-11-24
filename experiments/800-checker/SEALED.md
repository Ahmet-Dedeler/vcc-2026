# Lab reserve LZ (sealed): declaration and access rule

Declared 2026-10-15 23:06Z by the lab lead ([c61d4c]); the draw record was committed at e4cca4b2 (23:03Z) before any
statistic involving these panels existed. To be appended to `docs/CONFIRMATION-LEDGER.md` by the production coordinator
(charter Amendment 1: the lab appends to that ledger only through the coordinator).

## What
Two Arc-like 300-target panels drawn by the lab's rule "arc800" (`experiments/800-checker/DESIGN.md` §3.1;
`experiments/801-lab-folds/draw801.py`, sha256 c4b8507e8f4c6a0f02069ec5766b2ee0b62452494a2df16b2a5246fc7f2a7ead), in the same
run and after the six WORKING panels L1–L6, excluding every panel file on disk at draw time (A/B/C, real, efp, efr, R6 S1–S4,
R7a T1–T4) and each earlier lab draw:

| panel | seed | target-list sha256 (sorted, newline-joined) | file |
|---|---|---|---|
| L7 | 800107 | 97c30e7d11dddc0157219cd62a4312001553a7931b83bce8b3f087ebfddc9619 | `experiments/801-lab-folds/work/sealed/L7/pert_counts.csv`, file sha256 33bb3256ce92f058d7d71c3b7b5e1a3bd15799d6b1e12d4c2af5e0341d8723dc |
| L8 | 800108 | 32fd27d28425b12e102f7c15f16d01b95b10e32da9b37f9312d4df52c50694f1 | `experiments/801-lab-folds/work/sealed/L8/pert_counts.csv`, file sha256 c8d9f1a335a66e9ec9252bd0d5983bc104525b88b128be758d3c0bef77e3ab65 |

Held-out folds when opened: K562-GWPS, Orion HCT116, Orion HEK293T (PLOLO, the bulk engine of 801), i.e. 6 folds.
Nothing else exists for them: no bundle, CD4 bank, operator, dump or score (801's mandate forbids building any).

## Access rule
- One access in total, by the lab lead, for one preregistered final candidate (a dossier's configuration), after its
  WORKING price exists. The claim file (what is priced against what, the statistic, the reading rule) is hashed and the
  hash sent to the coordinator for the ledger BEFORE the first operator fit.
- The access builds L7/L8's bundles, CD4 banks, held-out operators and dumps on EC2, prices the claim with check800 at its
  frozen settings, and appends the result to the ledger through the coordinator, whatever it is.
- No fitting, selection, dose search or debugging may read an LZ score. Known answers run on a WORKING panel first.
