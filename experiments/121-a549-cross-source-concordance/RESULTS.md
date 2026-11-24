# Experiment 121 results

Status: **complete — NO-GO**

The preregistration, authorization manifest, 10,000 target-label
derangements, and 20,000 paired whole-target bootstrap samples were frozen and
verified before either source expression array was opened. No threshold was
adjusted.

## Primary result

The complete 20-target gallery used one global 7,685-coordinate response axis.
The mask required complete measured support in both sources and excluded every
VCC panel coordinate plus every common target-gene coordinate.

| Direction | PDS | Delta from chance | Own-target cosine |
| --- | ---: | ---: | ---: |
| Altair → Myllia | 0.455263158 | -0.044736842 | -0.029932911 |
| Myllia → Altair | 0.450000000 | -0.050000000 | -0.029932911 |
| Equal macro | **0.452631579** | **-0.047368421** | — |

The derangement-null macro PDS had mean `0.503233026`, p95 `0.575000000`,
maximum `0.653947368`, and one-sided p-value `0.874912509`.

The paired whole-target bootstrap 95% interval for macro PDS was
`[0.339473684, 0.569736842]`.

## Directional signs

| Direction | Top-20% sign precision | Directional sign credit |
| --- | ---: | ---: |
| Altair → Myllia | 0.409141184 | -0.143673756 |
| Myllia → Altair | 0.485230969 | -0.031691401 |

All seven preregistered GO criteria failed. The exact conclusion is:

```text
does_not_justify_a_held_context_donor_experiment
```

This closes direct use of Altair/Myllia A549 cross-source agreement as a
promotion basis. The result does not license post-outcome transformations,
target filtering, coordinate filtering, or threshold changes.

## Frozen hashes

```text
preregistered_protocol.json
7363e1a678a3f13c893886c6ab43ad5641f5b9ad4f1dcb8070b892fc53c781ee

AUTHORIZATION_MANIFEST.json
900a7630411685b149754fd5842d2633d6ca6500f02bfb2b729ccb638c4fba59

work/preregistered-freeze-v1/freeze-manifest.json
042a4cbe3681551e6205fd7fe3f6104fb8c8d488040bdeb75da2080f8478e41e

work/preregistered-freeze-v1/target-derangements.npy
acf951b797ba77b0f92a4e15275285a157f8b0ba421f8b25b37a2b81ae9cef62

work/preregistered-freeze-v1/whole-target-bootstrap.npy
7dd464074c1b50a89809899ce0849c902dc7fab7bf7b52a62b6718fa6e89f331

work/source-only-evaluation-v1/artifact-manifest.json
67a4952e5ec73e0340c57463166e1db9fab0e2ee0bc170e77fb0cbb46470a5b8

work/source-only-evaluation-v1/evaluation-receipt.json
c9514fd61d0141edc49c928200c18d6258e2bec7b0fb70bb331cb92bc77c52db

work/source-only-evaluation-v1/result.json
0fa222defd828a0226a687eb4637f4c34db23924ace6a5f9abf33462c9065f57
```

## Safety

No VCC held perturbation truth, other experiment outcome directory,
`cell-eval2`, count-generation code, challenge-build path, submission path,
commit, or push was accessed.
