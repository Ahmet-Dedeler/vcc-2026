# 041 — CCLE matching with cross-platform sanity checks

The challenge controls are UMI pseudobulks while CCLE is bulk RNA-seq TPM.
Any apparent nearest neighbor is therefore provisional unless the same metric
recovers known K562 and HCT116 controls across platforms.

This experiment evaluates several platform-robust feature/rank variants, records
the held-out rank of those two known lines, and applies only anchor-validated
methods to contexts A/B/C. Context identities are not used to tune the matcher.
