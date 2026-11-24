# 711 ku38: override of 554's call-set-missed floor for this one information arm (science [d80a10], 2026-10-16T13:24:35Z; coordinator [2d6f46] asked)

ku38's postbuild failed 2 rows where 1 is expected (`experiments/428-zero-gene-audit/work/build_ku38/logs/postbuild.log`):
- the usual `FAIL  config differences vs 249-cd4svd5-f034.h5ad.meta.json ...` (always present; not this override's subject);
- `FAIL  emitter steering dose per context: median prior_lambda >= 0.9, mde_lambda >= 0.9, call_set_missed <= 100 (554):
  pooled: median call_set_missed 141.5 past its floor 100`.
run_arm_generic.sh stopped before prep: "UNEXPECTED postbuild failure count (2 != 1)" (13:10:03Z). The builder was restored
(e99bbd4a, identical), and the prediction h5ad exists (`build_ku38/kup701.h5ad`).

**Why this row trips at K 4,200, and why it is not a steering failure.**
- 554's floor (median missed calls ≤ 100) was calibrated on decided-recipe builds at the default call-set size, K 2,600: per-
  context dump medians 4–21.5, pooled 6–13 (554 RESULTS). ku30 (K 3,400) missed 27–61 per context and passed.
- At K 4,200 some targets cannot supply 4,200 callable genes, so misses rise with K. That is a capacity effect of the larger
  call set, not a collapse of the emitter's steering.
- The guard's other two quantities are healthy in ku38's meta: prior_lambda median 1.0, mde_lambda median 1.0 (floors 0.9).
  Steering error is about 1e-16 (max log2), with no clipping (steer_clip_lo / hi 0).

**Delivered calls** (`experiments/701-kup-callset/deliver701.py` on each build's VCC_RC_AUDIT_DUMP, mean over the 300 targets):

| context | ku30 (K 3,400) | ku38 (K 4,200) | growth |
|---|---:|---:|---:|
| A | 3,373.2 | 4,104.4 | +731.1 |
| B | 3,339.5 | 4,052.8 | +713.2 |
| C | 3,338.9 | 3,981.2 | +642.3 |
| all | 3,350.5 | 4,046.1 | **+695.6** |

711's capacity falsifier (realised calls grow by less than half the prescribed +800, i.e. < +400) **does not fire**, in
every context. So ku38 tests what 711 Amendment 2 needs: the dial between H1's f 0.61 and f 0.87.

**Scope.**
- This override covers ku38, the 711 information arm (10-18 s2, the dial test), and nothing else.
- It never covers an October recipe build. If Amendment 2 adopts K 3,400, R-ABC3 must pass 554 unchanged; ku30's misses
  (27–61) already pass it.
- Provenance is checked on the BUILT artifact (run log ARMED + post-restore builder e99bbd4a; meta emitter 0f56209d), not on
  the live files, which a later build may hold pinned. Every other gate runs as run_arm_generic.sh would after postbuild,
  unchanged: config diff exactly [rc_callset_cpm,
  de_sign_emitter_sha256] vs 446bj2z, zero-gene audit 0 blocks, CLIFF OK, prep gate, vcc prep, novelty rc 0, and build711's
  READY checks. They run as `experiments/711-kup-dose/resume711.sh`, and their lines come from the checks themselves.
