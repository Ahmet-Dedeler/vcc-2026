# Coordinator review of 304 — accepted corrections, bounded conclusions

304 completed successfully. Its report is an audit, not fresh reserve
confirmation. R4 stays spent and C6 stays VOID.

Independent reconstruction: `scripts/audit_round5.py`, saved output
`work/round5-coordinator-review.json`. Neither reserve HDF5 file is opened.
302's 175 saved TF cosines reproduce to the CSV's rounding bound
(maximum difference 4.8528e-7); its primary median reproduces within 1e-12.

## Corrections accepted

- **No equivalence or latent-mirror bound.** Median observed cosine
  -0.000675, interval [-0.0047, +0.0016], is not evidence of biological
  orthogonality. The KD reliability was borrowed from another bank, and a
  half-sample quantity was mixed with an OE full-sample quantity. Dividing a
  median and its interval by a plug-in median reliability does not establish
  a latent-cosine confidence bound.
- **No established ORF depth law.** The original c ≈ 9,600 was inferred across
  different TFs, not fitted to within-TF depth ladders. The saved-object ladder
  reproduces: split-half r = .04222, .03373, .02641, .02461 at 1, 2, 4, 8
  shards per half. Withdraw the generalization and the 3,100-cell prescription.
- **No established target-specific geometry effect.** The observed
  rho(gen_pc20, |S|) = .382314 reproduces. A separate 2,000-draw uniform
  permutation check gives mean null .28942, central 95% [.15247, .42663],
  plus-one p = .09645. This agrees with 304's 200-derangement diagnostic.
  KD geometry and generic alignability account for much of the observed
  association; the matching-specific effect is not demonstrated.
- **Different isoforms are different interventions.** Their agreement is not
  a noise ceiling for exchangeable repetitions of one construct.
- **Line is not state.** The primary Joung text saved in
  `work/joung_pmc.txt` identifies H1 for the atlas and seven days in STEMdiff
  APEL. Correct DOI: 10.1016/j.cell.2022.11.026. The exact KD exposure/state
  has not been established here.

## Audit wording not adopted without qualification

“The same data refute the depth law” is too broad. The reconstructed ladder
is conditional on the original count normalization, cohort centering, and
sorted nested-shard construction. Six zero-cell shard groups occur among
the full 175 TFs. Log transformation and depth-dependent bias also require
controls before a biological learning-curve claim. The defensible result is
**the current evidence does not establish this ORF depth law**.

“A different assay state is a more likely explanation” is a hypothesis, not
a measured causal comparison. Likewise an audit's NEW/KNOWN tags are scoped
literature judgments, not proof that no prior work exists.

The 58-TF overlap in 303's predecessor is not a numerical error: all 175
scored 302 rows have S_g; 58 overlap R5, of which 45 belong to 302's primary
148-TF subset. Report both populations.

No production setting or submission decision changes. Local findings and
manuscript wording are corrected; nothing is republished externally.
