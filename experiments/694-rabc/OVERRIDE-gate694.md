# 694 R-ABC: coordinator override of gate694's config-diff row ([2d6f46], 2026-10-16T03:55:19Z)

gate694 failed one row: config diff vs 524dr2 == [bulk_load_sha256, de_dnmag, de_uprn, sameline], want [de_dnmag, de_uprn]. Checked by hand:
- **bulk_load_sha256:** the two dumps (`694-rabc/work/run/dump_h1x2_cd4x1.npz` vs `658…/dr2/run/dump_h1x2_cd4x1.npz`) have np.array_equal
  bulk, carriers, perts and genes; the only difference is the `env` string's VCC_BULK_DUMP path (each build's own work dir).
- **sameline:** every array of the two sameline.npz artifacts is equal (NaN-aware); only the `provenance` text differs (identity_source path).
  covered_targets 0, rows_changed 0, programme_down_frac identical.
- **F1 bulk identity (bulkid658):** 0 of 900 blocks differ, abs_sum 0. bulk_ctx_scale == 524dr2's (D 1.06 / E 1.1113 / F 1.126).
- Every other gate694 row ok; CLIFF OK; novelty rc 0; package 125c1db94c3843f6d4f422a777f421180d030e065c75604c30d4a6bbd51a48d8 (== the run log's).
So the prediction differs from 524dr2 only by the two intended DE riders; the failing row compared byte hashes of artifacts that embed
their own paths. **Override: R-ABC is READY for 10-16 s1.** The gate gets fixed to compare content, not bytes (asked of [4592ee]); the same
fix is needed for OCT22-CHECKLIST §C 6a's "config diff vs R0" before Oct 22.
