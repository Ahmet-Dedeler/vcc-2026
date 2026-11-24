# Pre-decided: what happens after `h1union`'s board result

Written 2026-09-04 02:25 UTC, **before the build finished and before the result exists.**
Incumbent `cd4f034` = **0.2049196**, rank 7. Leader `cqawesome` = **0.22082**.
Predicted `h1union` = **+0.0055 to +0.0105**, board **0.2104 – 0.2154**.

| branch | condition | next move |
|---|---|---|
| **A** | overall in range **and** raw `pds` **≥ +0.0010** | Rule **and** sources validated. Next slot = **223's ranked stack as riders, re-priced against `h1union`, not against `cd4f034`** — 208's rider rule: the H1 rows are now occupied, so any source sharing them has a smaller marginal. The `g` ladder is the alternative only if 223 returns nothing above the corrected shape-only bar (**+2.5–2.9%** cosine). |
| **B** | overall up, `pds` inside **±0.0010** | The coverage-aware rule is validated **as a code fix**; the 4th-estimate reading is unproven. Next slot = **`g = 0.17`**, not 0.45 — read the `pds` slope off the board rather than off PLOLO. |
| **C** | raw `nmae` rises **> 0.002** | Amplitude overshoot on the carried rows: the union is acting as a **dose increase on 77 rows**. **Halve `g`; add no sources** until the overshoot is understood. This is the failure that killed `cd4f080`. |
| **D** | overall falls | Revert to `cd4f034` as shipped. **Do not spend a slot restoring rank** (§6). **The coverage-aware sizing rule still stands as a code fix** — a solved dose of *exactly zero* is never the right answer — only the sources fall. |
| **E** | byte-identity or the CD4-fraction constraint fails in 222 | **No submission, slot held.** 222's mandate already says this. |

**Unchanged on every branch:** `VCC_ALPHA` **0.982**, `VCC_CD4_TARGET_F` **0.34**, emitter
`rc_callset`, operator sha `e2f65ecf`. One variable per slot.

**And a standing caution about the cushion.** Against `cqawesome`, `pds` + `mse` + `nmae` is
**−0.0317** and we are **+0.0157** ahead on `fid`/`jac`/`reach`. That cushion is an **emission-side**
lever — `rc_hybrid` + `rc_callset` bought +0.0178 of it in two slots without touching the profile —
and the leader sits at `fid` rank 390 and `jac` rank 249 while holding `pds` and `mse` rank 1.
**Any leader who notices the call-set lever takes most of that cushion in one slot. Plan as if the
true gap were 0.03, not 0.016.**
