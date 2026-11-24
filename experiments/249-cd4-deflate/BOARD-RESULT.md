# BOARD RESULT — `cd4svd5` (`VCC_CD4_DEFLATE_K=5`): **0.2160944, rank 3 of 688. Δ = +0.0042634 vs `h1union`. The mechanism transfers (`pds` at 1.2× PLOLO); the cosine cost was under-read 3.8× (falsifier 2 fires). New incumbent.**

Entry `Xv0ayvuoHpsKYWv2iu42`, `You can do it too.`, description empty. Package `submissions/249-cd4svd5-f034.vcc`
`6e8b41c6…`. Fired 2026-09-05 00:30:21Z (pid 77532, after the 237 entry's record was on disk), created 00:30:23Z,
published 00:54:57Z. Baseline `h1union` **0.2118310** (`data/board-entries/iUziqE1gewfcPKzYUDaB.json`). Complete
record `data/board-entries/Xv0ayvuoHpsKYWv2iu42.json` (captured by the handoff at 00:54:59Z and re-validated by
the coordinator's watcher). Reader: `scripts/read_board_entry.py data/board-entries/Xv0ayvuoHpsKYWv2iu42.json --prereg 249`.

## 1. The six members against the preregistration (P249 §3)

| member | raw `h1union` → `cd4svd5` | Δ raw | overall contribution | P249 predicted | verdict |
|---|---:|---:|---:|---:|---|
| `pds` | 0.83572 → **0.85064** | **+0.01492** | **+0.00545** | +0.0032 / +0.0046 / +0.0091 at 0.7 / 1.0 / 2.0× | 1.18× the 1.0× row; raw +0.0149 = **1.22× PLOLO's +0.01226** |
| `mse` | 0.81726 → 0.82181 | **+0.00454** (worse) | **−0.00080** | −0.0003 (raw 0.8190–0.8194, R249 §9) | **2.7× worse than priced** |
| `nmae` | 0.93389 → 0.93397 | +0.00008 | −0.00005 | −0.0001 … +0.0003 | inside |
| `fid` | 0.50748 → 0.50739 | −0.00008 | −0.00006 | −0.0001 | inside |
| `reach` | 0.27655 → 0.27506 | −0.00148 | −0.00026 | −0.0002 | inside |
| `jac` | 0.04096 → 0.04090 | −0.00006 | −0.00003 | ≈ 0 | inside |
| **overall** | 0.2118310 → **0.2160944** | | **+0.0042634** | +0.0026 … +0.0040, upside +0.0086 | above the band's top, under the upside |

Five of six members inside their priced bands; `mse` is the one miss. Same panel and anchors as the incumbent.

## 2. Falsifiers (P249 §4), pre-committed

1. Raw `pds` < +0.0040 — **no** (+0.01492; 3.7× the threshold). The mechanism transfers, at **k_pds 1.22** — above
   `h1union`'s 0.70× and above 243's dose-class band (0.95–1.15).
2. Raw `mse` worse by > +0.0030 — **FIRES** (+0.00454; the line was 0.82026, R249 §9). On 197's identity at the arm's
   own terms (raw at zero cosine change 0.81809, d raw/d cos −1.996): board Δcos −0.001862 = **−1.735 % of the
   measured 0.107332**, against PLOLO's aggregate **−0.4617 %** (macro −0.8087 %). **Realised k_cos = 3.76 on `cos_agg`
   (2.15 on macro)** — a third class, "re-shape occupied rows", far above re-dose (1.36) and new-source (0.94).
3. Overall outside [−0.0010, +0.0100] — **no** (+0.00426).
4. Config diff exactly `cd4_deflate_k` + consequences — **PASS** off the package (R249 F1, 115/115 pinned).
5. Delivered amplitude within 0.5 % — **PASS** off the package (−0.409 %, R249 F5). Branch E does not fire.

## 3. Branch — recorded by the letter, then read on the evidence

**By the letter of `DECISION-TREE.md` and its 23:32Z addendum (precedence E > D > A > B > C): branch D**, because
falsifier 2 fired. Recorded as such.

**What D's text actually asserts, and what the point says.** D was written as "same as C, plus …": *C's* premise is
"the mechanism does not transfer; the arm scored on its cosine legs by luck" — and this point is the exact opposite:
raw `pds` +0.0149 (falsifier 1 not fired, by 3.7×) carried the whole gain **while the cosine legs LOST** (`mse`
−0.0008). D's "retire the deflation CLASS on A/B/C; October basis stays UNSET" inherits C's premise and is
contradicted by the members. D's *own* clause — "the cosine cost was under-read; re-derive 243's class factor for a
re-shape-occupied-rows class from this point before pricing any future re-shaping arm" — is exactly right and is
adopted in full (§2, item 2). Both of A's conditions hold (raw `pds` ≥ +0.0086, overall ≥ +0.0026).

**The October decision does not rest on this tree alone.** `docs/DEF-RUNBOOK.md` §3.1 (250, written 2026-09-04
before this point): *"unset until the 2026-09-05 `cd4svd5` board point is read. If board raw `pds` rises against
0.83572, ship the deflation on D/E/F as the PINNED A/B/C basis, k = 5 (`VCC_CD4_DEFLATE_BASIS=V5_abc`) — never a
day-of re-fit."* Raw `pds` rose 0.83572 → 0.85064. **That trigger fires.** It was pre-committed on the member that
tests the mechanism, independently of the tree, and it agrees with A's October action.

**Decision (lead, 2026-09-05 01:10Z):**
- **October ships the PINNED basis** (`VCC_CD4_DEFLATE_BASIS=V5_abc`, 251: on A/B/C it equals the re-fit to
  2.4e−15; never `VCC_CD4_DEFLATE_K` on D/E/F). Wiring per R251 §7: `DEF_ENV` gains the basis, `SHIPPED_ENV` gains
  `VCC_CD4_DEFLATE_K=5` for A/B/C byte-identity with this package, `ARM_BOARD_DELTA` gains **+0.0042634**,
  `assert_build_carries_arms` gains the row; 240's row 1 is re-run against this package's hash.
- **Every further re-shaping arm is re-priced on k_cos 3.76 (agg) before a slot.** On that factor **k = 10 (255) is
  dead**: its extra −1.08 % aggregate cosine over `cd4svd5` costs ≈ −0.0080 raw `mse` (−0.0014 overall) against
  `pds` +0.00054 × 1.22 ≈ +0.0002 overall. `+partial` (254, cosine +1.16 % agg, `pds` +0.0018) is a different
  operation (rows added, not re-shaped) and needs its own transfer reading before it is priced.
- **The class is not retired and the basis is not left unset.** Recording D by the letter and acting on C's
  premise would be the tree's wording deciding against its own evidence; the preregistration's falsifier 2 says
  what it says — the PLOLO cosine read for this class is refuted — and nothing more.

## 4. Board position

`cd4svd5` **0.2160944 is the incumbent** and the displayed row, rank **3 of 688**. Leader `cqawesome` 0.2245476
(gap **+0.00845**, was +0.01272); 2nd `Aginglab.com` 0.2193913 (moved 18:08Z; gap +0.00330). Two slots on
2026-09-06 (00:00Z); every candidate is re-priced against this entry. `h1union` and the 237 entry no longer resolve
on the API; all three records are in `data/board-entries/`.

## 5. What the two slots bought together

- The `nmae` identity is validated at scale (237) and `nmae` stayed flat here as priced — the instrument that now
  gates direction arms was right twice tonight.
- The board's `pds` transfer for a re-shape of occupied rows is **≥ 1.2×**; its cosine transfer is **3.8×**. The
  first is the largest `pds` step we have ever bought (+0.00545 overall); the second is why the arm netted
  +0.0043 rather than the +0.0086 upside.
- Open hypotheses, each needing a preregistration: an `nmae`-neutral rendering of ARM6's six-source directions
  (237 `BOARD-RESULT.md` §3); a deflation that spends less cosine per unit `pds` (the projector removes 20 % of the
  CD4 energy; a softer shrinkage of the same directions is the obvious knob, priced on 3.76×/1.22×).
