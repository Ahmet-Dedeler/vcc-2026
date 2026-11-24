# BOARD RESULT — `ARM6` (stack6, six union sources, g = 0.35): **0.2081907, rank 4. Δ = −0.0036403 vs `h1union`. The `nmae` identity is VALIDATED at scale — sign right, magnitude conservative.**

Entry `ydAemVdWmsTx9xhDWcou`, `I'll just keep moving forward.`, description empty. Package
`submissions/229-stack6-g035.vcc` `7a3c2d3e…`. `vcc submit --wait` started 2026-09-05 00:02:21Z (pid 41699,
slot reset + 2 min), upload 2.6 GiB in 12m43s, published 00:28:02Z. Baseline `h1union` = **0.2118310**
(`data/board-entries/iUziqE1gewfcPKzYUDaB.json`). Complete record: `data/board-entries/ydAemVdWmsTx9xhDWcou.json`
(captured by both the second-slot handoff and the coordinator's watcher at 00:28:13–14Z; the API will stop
resolving this id once the 249 entry publishes). Reader transcript:
`reading-ydAemVdWmsTx9xhDWcou.txt` (`scripts/read_board_entry.py … --prereg 237`).

## 1. The six members

| member | raw `h1union` → `ARM6` | Δ raw | scaled Δ | overall contribution |
|---|---:|---:|---:|---:|
| `pds` | 0.83572 → 0.84238 | **+0.00666** | +0.01490 | **+0.00248** |
| `mse` | 0.81726 → 0.81227 | −0.00499 (better) | +0.00525 | **+0.00088** |
| `nmae` | 0.93389 → 0.96016 | +0.02627 (worse) | **−0.04296** | **−0.00716** |
| `fid` | 0.50748 → 0.50757 | +0.00009 | +0.00029 | +0.00005 |
| `reach` | 0.27655 → 0.27698 | +0.00044 | +0.00058 | +0.00010 |
| `jac` | 0.04096 → 0.04099 | +0.00003 | +0.00009 | +0.00002 |
| | | | **sum** | **−0.0036403** |

Same panel and anchor set as the incumbent (`vcc2026-val-1`, `valA-r4+valB-r4+valC-r4`).

## 2. Reading — pre-committed (P237 §4, ADDENDUM-243, ADDENDUM-2-READING), applied by the reader

- **Pre-committed overall bin: row 1, Δ ≤ −0.0015** (recorded; its inference about `nmae` is superseded by the direct read).
- **Direct `nmae` read: −0.00716 overall against the identity's −0.004783 — 150 % of the prediction, same sign.**
  FINAL verdict **HOLDS → VALIDATED at scale on this arm class.** The instrument that was "wrong in sign both times
  it was checked" (P237 §2) was right in sign here and **under-predicted** the loss by a third. 228's per-gene gain
  and `ARM6` stay dead on the `nmae` leg; the 48-day programme proceeds on §8.1d/8.1f as written.
- **P237 §4 falsifier ("Δ ≤ −0.0015 AND the leg within 20 % of −0.004783 ⇒ the slot was mis-spent") does NOT fire**:
  the leg is 50 % off the prediction. The slot bought a genuine calibration point: the identity's magnitude is a
  lower bound on the `nmae` loss for a six-source union, not an over-estimate.
- **Direction legs (overall minus `nmae`): +0.003519**, inside 243's class-rule band **+0.0022 … +0.0038** (point
  +0.0024), at its top. 229's six-member composite −0.00219 had the sign right and was 0.0015 optimistic, entirely on
  the `nmae` leg. Both instruments are confirmed in the direction they were used; neither is retired.
- **Raw `mse` 0.81227 → 243's NEW-SOURCE bin (0.81135–0.81392), point prediction 0.81244**: 0.00017 off.
  **Realised k_cos = 0.962** (new-source band 0.94–1.00). That is the fourth board point for the class model, and it
  lands in the class bin, not the coverage bin (`k = 0.78 + 0.63·cov` would have read 1.08–1.25).

## 3. What the point says beyond the gate

The six-source union at g = 0.35 **bought `pds` +0.00248 and `mse` +0.00088 overall (+0.00336, above the +0.0019
shape-only floor on its own) and paid −0.00716 on `nmae`**. Raw `pds` +0.00666 is almost twice `h1union`'s own
+0.00344 step: the union's *direction* is worth having; its *magnitude* on the Wilcoxon-significant genes is what
the board refuses (raw `nmae` 0.934 → 0.960). Against 243's member-level prices (o_pds +0.00070, o_mse +0.00085):
`mse` landed on its price; `pds` delivered 3.5× the class rule's k_pds 0.705.

**Open hypothesis for the next slot program (not a candidate until preregistered and priced on the now-validated
identity):** an `nmae`-neutral rendering of the same six-source directions — the union's LFC magnitude re-calibrated
on the significant set (188's map-fit family) so the identity's leg reads ≈ 0 — would keep the +0.0034 and drop the
−0.0072. The identity is the right instrument to price it: this point shows it is conservative, not optimistic.

## 4. Board position

Displayed rank 4 (the board shows the latest entry; cosmetic until Oct 22). `h1union` is no longer resolvable on the
API; both records are on disk. The second slot (`cd4svd5`, 249) fired after this capture, at the 00:30Z floor.
