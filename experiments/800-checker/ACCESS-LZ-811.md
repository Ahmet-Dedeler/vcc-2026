# LZ access record (one access; claim experiments/800-checker/LZ-CLAIM-811.md sha256 3a8f725216075e793a9968b53f2ff0f28f845593ba49e48d3c659f009bd258c2)

- 2026-10-22T07:02:15Z access opened by the lab lead. Override: this wrapper (`lzaccess811.sh` sha256 e990d562b26c987787fbf7d035e996a26815ac52c314fa1db3d85fa5d1f60030);
  tags LZ7 = L7 (file 33bb3256ce92f058d7d71c3b7b5e1a3bd15799d6b1e12d4c2af5e0341d8723dc), LZ8 = L8 (file c8d9f1a335a66e9ec9252bd0d5983bc104525b88b128be758d3c0bef77e3ab65); A/B/C frames; kit801 1ad7545c8d020f18 unchanged.
- 2026-10-22T07:10:00Z kits launched: LZ7=20261022T070218Z-lab801-kit-lz7 LZ8=20261022T070608Z-lab801-kit-lz8 
2026-10-22T07:26:20Z LZ7 kit collected: rc all 0 = 0 failures
2026-10-22T07:26:20Z LZ8 kit collected: rc all 0 = 0 failures
- 2026-10-22T07:32:23Z composite runs launched: 20261022T072630Z-lab811c-kit-lz7 20261022T072928Z-lab811c-kit-lz8  (KA lines read first; ka-lz7/ka-lz8 arrays on disk, unread)
- 2026-10-22T07:46:13Z composite KAs (read before any score): both runs DONE 0, every fold rc 0; KA-705/707 PASS on LZ7 x3 and LZ8 x3 (worst |d| 0.0), ABCp (<= 1.1e-13), L1 (0.0); KA-811 PASS (L1, 0.0); KA-L0 10/10 per run; L1 CD4 compact == 801's.
- 2026-10-22T07:46:13Z pooled run dir run-comp811-LZ (the 6 LZ folds only); aggregation agglz811.py (327b7d02, committed bb0c5cc0 before any LZ score), pairs C0:B0 (deciding) + C1:B1 (reported).
- 2026-10-22T07:48:36Z aggregation 20261022T074617Z-agglz811 (rc 0, $0.008), json `experiments/811-lolo-model/work/run-comp811-LZ/agg811-811-LZ-6048cf1a.json`.
- **RESULT: PASS (sign).** **C0 − B0 (deciding): +0.00999 [+0.00659, +0.01273], 6/6 folds > 0, LT3 +0.00571** (LZ7 K562 / HCT116 / HEK293T, LZ8 K562 / HCT116 / HEK293T: +0.00929 / +0.00912 / +0.01651 / +0.00611 / +0.00811 / +0.01080; pdsl +0.01044). C1 − B1 (reported): +0.00645 [+0.00373, +0.00882], 6/6.
- Essential-screen counts (full library `814/work/flags/esslib.txt`, 24fffbfe…): L7 9, L8 20. These arc800 panels are expected to over-read
  key / direction changes: here 2.2× DEFo's C0 − B0 (+0.00446).
- Nothing else was read from LZ. The kit's ka804 arrays (work/ka-lz7, ka-lz8) and the composite logs' arm lines stayed unread, and there was no
  dose, λ or arm choice and no debugging. **LZ is SPENT** (one access). The result goes to the ledger through the coordinator.
