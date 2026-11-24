# 674 — W2 (W1 + DEF_DAY_CD4_STACK=x1.5): NOT BUILT (dropped by the coordinator, 10-20 ~11:33Z)

- readiness armed w2_674.sh (fd61ce14) to build after W1 ended READY. W1's own driver ended NOT READY at 11:15:45Z (gate673's spec
  gap on bulk_sharpen_pt, see 673 RESULTS), so W2's driver stopped at 11:16:27Z: "NOT READY W1 did not end READY … W2 has no
  reference run". Nothing was built.
- The coordinator then re-gated W1 (READY w1 1337cc9b) and **dropped W2** (dcca185b, "SIXA input record: W1 re-gated, W2 dropped
  (Amendment 3 inactive)"). W2 is not re-armed.
- If W2 is ever revived: gate674 would hit the same bulk_sharpen_pt scale shift (a CD4 dose also rewrites the dump). Allow it only on
  regate673's rule (behavioural fields identical), and point the W1 check at regate673.ready.log, not w1_673.log.
