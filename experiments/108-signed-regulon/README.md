# 108 — signed-regulon transfer pre-gate

Status: **completed; direct signed-regulon transform NO-GO**

Date: **2026-08-30**

## Question

Does a signed transcription-factor regulon contain transferable downstream
direction information that unsigned STRING proximity missed in experiments
100 and 104?

This is a mechanism gate, not a submission candidate. Directly applying a
literature network to revise challenge output would sit too close to the
competition's prohibition on revising predictions from literature values. A
positive result here would authorize only a second-stage model trained on
public Perturb-seq labels with regulon-derived features.

## Frozen network

The network is the official OmniPath DoRothEA interaction export:

```text
https://omnipathdb.org/interactions
  ?datasets=dorothea
  &genesymbols=1
  &format=tsv
  &dorothea_levels=A,B,C
  &fields=sources,references,dorothea_level
```

Local artifact:

```text
data/external/omnipath-dorothea/dorothea-abc.tsv
SHA-256 37f4f6959e9faf960837e8ad5a456514bcd081d12dd11161205dad1f3bc8f671
32,628 interactions plus header
```

Only edges with exactly one consensus sign are retained. Confidence weights
are fixed at A=1, B=0.5, C=0.25. Duplicate TF-target edges are summed; zero
net edges are discarded. Each TF regulon row is L2-normalized and must retain
at least five response genes.

## Frozen evaluation

- Source: K562-GWPS only.
- Held contexts: HCT116 and HEK293T.
- Gallery: 300 shared non-panel perturbations chosen by SHA-256 ordering.
- Response axis: official genes shared by K562/HCT116/HEK293T, excluding the
  VCC-300 targets and all 300 gallery targets.
- K562 effects are centered and median-norm-scaled on non-panel source rows.
- Held truth is raw, uncentered log-pseudobulk effect.

The regulon reconstruction uses fixed ridge `lambda=0.01`:

```text
A = X R^T (R R^T + lambda I)^-1
X_reg = A R
```

Two fixed candidates are evaluated:

1. `reweight`: sign-preserving `beta=0.10` reweighting of K562 coordinates by
   agreement with `X_reg`, followed by exact row-norm restoration.
2. `blend10`: 90% K562 plus 10% row-norm-matched `X_reg`, followed by exact
   row-norm restoration.

The primary candidate is `reweight`. Twenty deterministic nulls permute the
reconstructed response coordinates before applying the same reweighting.

Promotion requires reweighting to:

- improve PDS in both held contexts;
- improve macro PDS by at least `+0.003`;
- have a positive target-bootstrap lower bound in both contexts;
- exceed every frozen permutation null in both contexts.

Failure closes this direct signed-regulon transform. It does not rule out
regulon features inside a supervised model.

No count generation, `cell-eval2`, KOLF access, challenge build, submission,
or production edit is authorized.
