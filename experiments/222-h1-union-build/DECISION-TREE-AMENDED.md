# AMENDMENT to `DECISION-TREE.md` — reissued band and branch thresholds

**Written 2026-09-04 04:07 UTC. The entry (`iUziqE1gewfcPKzYUDaB`) is `launching`; no score exists,
no member is visible, and nothing about the result is known.** This amendment exists because the
original tree's thresholds were derived from a statistic I have since shown to be inflated, and a
preregistration must be corrected **before** the result, not after it.

## 1. Why the original band is withdrawn

`DECISION-TREE.md` predicts **+0.0055 to +0.0105** (board 0.2104–0.2154). That band was built on the
**4-fold** Δ`pds` of **+0.009166**, which includes `kolf2`. `kolf2` is an iPSC fold, the added source
is **H1, an embryonic stem line**, and `kolf2` gains **+0.021628** — 2.3× the next fold and 3–5× more
favourable, because *both* the source and that recipient are pluripotent. **The three challenge
contexts are Jurkat, HeLa and CAL33 — all cancer lines, none stem-like** (041). 222 reached the same
conclusion independently and excluded `kolf2` from every headline it reported.

The transferable statistic is the **3-cancer** Δ`pds` = **+0.004887**, with relative MACRO cosine
**+7.6803%**.

## 2. The reissued estimate

On 216's **realised** board transfer factors, measured against `cd4f034`'s own board point
(`k_cos` **1.384**, `k_pds` **1.987**) — not on PLOLO rates:

| route | arithmetic | estimate |
|---|---|---:|
| cosine | 7.6803% × 1.384 = 10.63% relative; × board cos 0.09018 = +0.009586; × 0.365 | **+0.00350** |
| `pds` only | +0.004887 × 1.987 = +0.009710 raw; × 2.2138 / 6 | **+0.00358** |

The two agree, but they are **not independent** — 197's inversion makes `mse` a cosine readout, so
the cosine route already contains most of what the `pds` route counts separately. **Take +0.0035 as
the optimistic point, not the sum.**

**Now apply 222's §9 concentration finding**, which the original tree could not: Δ`pds` retains only
**14.6%** under leave-top-3-out while the cosine gain retains **72.7%**. The `pds` half of this arm
rests on three targets per fold and is the half most likely not to transfer.

> ### Reissued: **+0.0010 to +0.0045, central +0.0030.** Board **0.2059 – 0.2094**, central **0.2079**.

## 3. Reissued branch thresholds

| branch | reissued condition | next move |
|---|---|---|
| **A** | overall **≥ +0.0020** and raw `pds` **≥ +0.0010** | Unchanged from A. Rule **and** sources validated. Next slot = 223's stack as riders, re-priced against `h1union` (208's rider rule), each now also owing a leave-top-3-out retention test. |
| **B** | overall **+0.0005 … +0.0020**, or `pds` inside **±0.0010** | Code fix validated, 4th-estimate reading unproven. Next slot = **`g = 0.17`**; read the `pds` slope off the board. |
| **C** | raw `nmae` rises **> 0.002** | Unchanged. Halve `g`, add no sources. (222 measured **+0.00039** on the exact identity, so this is unlikely to fire.) |
| **D** | overall falls | Unchanged. Revert to `cd4f034`; **spend no slot restoring rank**; the sizing rule stands as a code fix. |
| **E** | gate failure in 222 | **Moot — 222 returned BUILDABLE**, all 3 known-answer rows and all 5 falsifiers passed. |
| **F** *(new)* | overall up **and** the realised `pds` gain is a **much smaller fraction of its prediction than the cosine gain is of its own** | **The concentration signature is confirmed on the board.** This is the most informative outcome available and it is a *finding*, not a disappointment: it would establish that leave-top-3-out retention predicts board transfer, and hand **225** a validated ranking rule for every source it prices. Record it in `AGENTS.md` §8 and re-rank 223's stack by retention before any rider flies. |

**Unchanged on every branch:** `VCC_ALPHA` 0.982, `VCC_CD4_TARGET_F` 0.34, emitter `rc_callset`,
operator sha `e2f65ecf`. One variable per slot. And the original tree's closing caution stands
verbatim: **plan as if the true gap to `cqawesome` were 0.03, not 0.019** — our +0.0157 ordering
cushion is an emission-side lever the leader has not yet noticed, and `rc_hybrid` + `rc_callset`
bought it in two slots.
