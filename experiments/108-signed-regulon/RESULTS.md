# 108 — signed-regulon transfer pre-gate

**Verdict: NO-GO.**

The fixed DoRothEA A/B/C network was reduced to 7,631 consensus-signed input
rows, 89 TF programs, 2,588 aligned signed edges, and 1,309 supported response
genes on a frozen 7,091-gene axis. K562 predictions, the 300-target non-panel
gallery, and 20 response-coordinate permutation nulls were sealed before
opening HCT116 or HEK293T effects.

| held context | K562 PDS | regulon-only | sign-preserving reweight | reweight delta | 10% direction blend delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| HCT116 | 0.68489 | 0.58488 | 0.68424 | **-0.00066** | -0.00196 |
| HEK293T | 0.64746 | 0.58441 | 0.64742 | **-0.00003** | +0.00054 |

Macro sign-preserving reweight delta was **-0.00035**. Its target-bootstrap
interval crossed zero on both contexts, and it did not exceed the frozen null
maximum on either context. All four promotion gates failed.

The network subspace itself is strongly lossy for this retrieval problem.
Adding a small amount of its direction is inconsistent across contexts, while
using it only to reweight matching signs is indistinguishable from noise.

Do not tune the fixed ridge, confidence weights, reweighting strength, or
blend fraction against these outcomes. A later supervised model may still use
regulon-derived features, but this result gives no evidence that a direct
literature-network transform should enter production.

Artifacts:

- `work/manifest.json`
- `work/frozen.npz`
- `work/results.json`
