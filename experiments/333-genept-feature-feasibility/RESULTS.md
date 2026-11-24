# 333 — feature acquisition succeeds; vector identities remain quarantined

September 9, 2026. **QUARANTINED_IDENTITY_UNRESOLVED** for biological use.
No response model, response-outcome inspection, paid embedding call or
submission. Schema/coverage feasibility is a different result from identity
integrity.

## Acquisition and disclosed recoveries

The official `yiqunchen/GenePT` repository links Zenodo record 10833191,
which declares CC-BY-4.0. The downloaded `GenePT_emebdding_v2.zip` is
574,395,233 bytes, published MD5 `3f6ce4317e3a0091978ae5cb8fbf05a3`
and verified SHA-256
`6193575dcbd7bbf214c8ca3eb518cb3d13272443f98ecd6c26402411a7ca745e`.

Original freeze **63e8a747**. The filename-only recovery **f78057ce**
handles the release's undocumented trailing period on its pickle member.
The original uppercase-name gate failed and is preserved in
`NAME-GATE-FAILURE.md`: 133,736 keys include 11,294 mixed-case keys,
22,635 whitespace-bearing keys and 10,970 uppercasing collisions.
The exact-canonical-only recovery **3a2cf643** performs no uppercasing,
stripping, alias inference or merging. Eight pure tests pass.

Conversion produces 18,187 finite nonzero float32 vectors of dimension
3,072, covering 18,187/18,533 canonical genes and all 300 challenge
targets. Its original schema/coverage status is **FEATURE_FEASIBLE_ONLY**.
Conversion takes 8.40 seconds (9.04 seconds outer), peak RSS
4,772,544,512 bytes.

## Identity audit

Scope **7730234e**, implementation/input freeze **3403336b**.

- 734 exact-vector groups involve 1,545 distinct canonical symbols.
- 34 challenge targets occur in 33 target-linked groups.
- All 33 target-linked groups have identical original float64 vectors
  and share a Python array object after restricted unpickling.
  Float32 conversion is not their cause.
- 716/734 groups have distinct nonempty primary descriptions in the
  released text dictionaries. This includes all 34 affected targets.
- Examples include PRKACA/PRKCA, PDPK1/PDK1, POR/VDAC2, BRD1/BRPF1
  and EZH2/EZH1.

The two source-text dictionaries each have 33,703 primary keys, nonempty
text for all 300 targets and 17,987 of the selected 18,187 canonical
genes. Coverage over the full 18,533-gene axis was not measured here.
The audit does not establish which vector is wrong and does not repair
any vector. Runtime 3.67 seconds, peak RSS 4,583,997,440 bytes.

A bounded read-only source follow-up at 14:45 UTC inspected the official
`input_data/gene_info_table.csv` (Git blob
`03d65efa77c6550c99889007f246294fc660000e`). Its columns are Ensembl ID,
gene name and gene type; it is not an alias-to-vector provenance map.
It assigns different Ensembl IDs to the example gene pairs. The public
text-request and illustrative embedding notebooks do not supply a
release-vector assignment trace that resolves these collisions.
Their code was inspected, not executed. No correction is justified by
this check.

## Artifacts and disposition

- Canonical vector artifact SHA:
  `4c4fa90d90bb57aa6b341087cbc236a66f3bd51924a9671deef0cd89367318a1`
- Canonical result SHA:
  `f1e2b145a26592b854ebad68fc0ed0d7c8ab3dcfe506cee55377d51d222c01a3`
- Canonical input manifest SHA:
  `8cfdc216f11e4f0f237b55401a3cd2211ece69e31376e27339c3ac0070465265`
- Identity audit SHA:
  `be9ee1951ecdf9e6c1b2e048e71e28811f6f780a7adb41c602b47573672cda98`
- Identity input manifest SHA:
  `b03bb353102f873bf7ee379a16418de32172ac1b42b533340ebfcfef53211da9`

Do not use the pretrained vectors for the next biological fit. A direct
text-feature construction would be a different artifact and requires its
own frozen primary-key policy and response-model protocol. This result
does not reject text-derived biological information as a category.
