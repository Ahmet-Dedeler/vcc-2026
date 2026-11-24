# 810 — Lane B1: what misdirects the strong tier on held-out lines (decomposition; PROPOSAL-2 awaiting the judge)

Worker: Claude Opus 5.5 under `MANDATE.md` (lab lead [c61d4c]), 2026-10-16 07:06–08:12Z (clock times; this file's first version guessed
08:20Z, which was wrong). Writes only in this directory.
No git, no production code, nothing under `data/`. Every heavy step ran on EC2 through `labcloud.sh`. Local work was
sub-second aggregations of ≤ 5 MB of per-target records; one local gate602 job was queued and then withdrawn unrun.

## 0. Answer
1. **The strong tier's misdirection is the WHITENING of the donor base, not H1/CD4 dilution.** On the top z_src tertile
   (BULK-W, 18 fresh folds), the whitened base points at the held line's truth with median cos 0.019. The unwhitened,
   normrel'd donor mean it is built from (`pre`) points with 0.033. On BULK-ABC cancer it is 0.033 vs 0.046.
   - **H1 and CD4 help the strong tier.** Removing either part loses tier-X on 18/18 BULK-W folds and 3/3 BULK-ABC
     folds (H1 −0.0059 / −0.0060; CD4 −0.0069 / −0.0093).
   - **This refutes 706's guessed mechanism** ("the strong-tier stack rows carry H1 ×2 and CD4 components that do not help
     on cancer truths").
2. **The whitening loss is predictable from public features.** Within the strong tier, `cos(pre,T) − cos(base,T)` rises
   with source strength and donor agreement:
   - BULK-W: Spearman +0.22 (count), +0.16 (agreement), +0.20 (z);
   - BULK-ABC cancer: +0.16, +0.14, +0.23.
   - It is positive for 65 % / 64 % of strong targets.
   - Knockdown depth and per-target noise barely matter (|ρ| ≤ 0.15).
   - The whitening loss is positive on EVERY tier (BULK-W mean per target +0.0096 / +0.0103 / +0.0193, weak / mid /
     strong), and the candidate swap (base → pre, residuals kept) gains at composition level on every tier too (BULK-W
     tier-X +0.0055 / +0.0074 / +0.0203; BULK-ABC +0.0066 / +0.0090 / +0.0092). The "≈ 0 on weak/mid" and the
     strong-agree vs strong-disagree split (+0.0055 vs +0.0005 fold-X, BULK-W) belong to a DIFFERENT construction, the
     pre-only swap (the whole row replaced by `pre`), not to the candidate; the first PROPOSAL.md mixed them up
     (VERDICT-1 item B).
3. **707's U blend works because it partly undoes the whitening, but it also thins H1/CD4.** Re-directing only the base
   toward `pre` (residual parts kept, row norm kept) carries 2–3.5× 707's composition-level cosine gain:
   - strong tier, BULK-W: +33 % vs +11 % fold r;
   - BULK-ABC cancer: +8.6 % vs +2.5 %;
   - kolf2: +3.1 %, where 707's form reads −3.9 %.
4. **Why gating could escape the global γ knife-edge (212, 097, 099 §5.3).** mse reads the energy-weighted aggregate
   cosine; pds is count-weighted. The strong tier holds 59 % (A/B/C) to 73 % (BULK-W) of the realised energy. Unwhitening
   only that third keeps 47–80 % of the global unwhitening's cosine gain at a third of its pds exposure. That is the
   gate's ONLY justification: the cosine gain exists on every tier, so the gate is a bet on pds exposure, tested in
   PROPOSAL-2 against the closed global form `sgu1all` in the same run.
5. **PROPOSAL.md** (P1 β 0.5 / 1, P2 agreement gate) got **REVISE** (`VERDICT-1.md`, judge run 819). **PROPOSAL-2.md**
   (frozen) prices ONE arm, sgu1 (top tertile, β 1), charged N = 6, beside three comparators: inc, 707 and the closed
   global form sgu1all as a known-answer row with a STOP rule. The decision rule is code (`decide810.py`).
   - **Status (08:53Z): decide810 (frozen) returns DE-H1 FIRST** (§9). sgu1's |dK_rel| > 5 % on 12 of 22 folds, all of them
     BULK-W. The bulk legs already fail C1 and C4 on BULK-ABC. sgu1all reproduces the γ closure (the STOP rule did not
     trigger). No DOSSIER.
6. **Known answers: all pass** (§1).
   - New finding: 706's stored `dirswap706-abc.json` is NOT reproduced by its own committed script. The committed path
     gives +1.73 / +2.67 / +1.43 / −2.33 % against the json's +1.64 / +2.56 / +1.79 / −2.30 %. The json was written 46 s
     before the script's last edit.
   - The sign and size of 707's A/B/C evidence are unchanged; the ABC/hek293t point is +1.43 %, not +1.79 %.
7. **707's board read has not landed.** The ledger has no 707 entry at 08:07Z; it is expected 10-17/18. The checker's
   prospective read (804's grade) is BULK-ABC +0.00044 [−0.00092, +0.00178] WEAK and BULK-W +0.00485 SHIP; the lab log
   records the 1b1f4a3c call as NULL ≈ +0.0002. That read is this class's first board point and sets P1's board projection.
8. **Cloud: $0.73 of the lane's $50** (4 launches, below). No 810 instance is running.

## 1. Instruments and known answers
- **Decomposition (`decomp810.py compose`, EC2).** lib801's own calls (bind → build_fold438 → Donors438 → residuals →
  pre_normrel → op.apply → compose) re-run per fold. The composition is split into its additive parts:
  `fin = s·W + s·f_cd4·C + s·f_h1·H + s·f_p·P` (base, cd4, h1, part).
  - L folds use the October rules: h2c1, the fold's held-out operator (801's registry sha256), the panel's CD4 banks.
  - ABC folds use the A/B/C rules: h2c1.5 = 518's fST.
  - Per target it stores the Gram matrix and the truth inner products on the held line's scored genes and on every
    non-panel gene, over base, cd4, h1, part, U (707's raw unit mean of the fold's sources), pre (unwhitened normrel'd
    donor mean), the processed source rows and the raw source rows. It also stores 522's robust strength / count / noise
    per bank, each bank's knockdown of the target (raw delta at the target's own gene), z_src and its tertile
    (fold522.zsrc's rule over the fold's targets).
- **Realised (`decomp810.py realise`).** engine801 (its own KAs on load) prices the October incumbent (pct+P2g) and stores
  x, a, d, pds, pdsl, K per frame and target.

| known answer | folds | worst \|d\| |
|---|---|---|
| parts sum == compose(final) | 22 | 1.8e-15 |
| parts sum == the registry dump's fin (L: dump708 `fin`; ABC: 518 `fST`) | 22 | 3.6e-14 (EC2 vs the Mac-made ABC dumps; L dumps EC2-made: ≤ 2.1e-14) |
| recomposed h1c1 == the dump's fPP | 22 | 3.9e-14 |
| my robust strengths == strength522.npz (abc: 3 banks × value/count/noise, NaN positions) | 4 | 0 |
| engine801 load KAs (PT0, 540c, arms522 …) | 22 | 4.4e-16 |
| KA-706: 706's top tier β 0.25 through dirswap706's committed code path (fold522.Abc / score518) == through the engine with Arm(fin=…): U, w, dX, dpds | 4 ABC | 0, 0, 0, 0 |
| dirswap706-abc.json reproduced by the committed script | 4 ABC | **NO**: \|ΔdX\| 0.022–0.354 %-points (`dev/ka706dbg.py`, run 20261016T073629Z-s810-dbg) |

Thresholds 1e-12, enforced by asserts. The first launch's KA-706 asserted against the json and failed: run a, rc 1 on
4/4 ABC folds. That failure is what exposed the json discrepancy. The assertion now compares the two implementations,
and the json is reported beside them.

## 2. The decomposition (full tables: `work/table810.md`; per-fold: `work/analyse810.json`)
Realised quantities use the incumbent's realised rows. "Composition" means the incumbent's composition rows (before key,
sharpening and render). Tier X is the instrument's aggregate cosine restricted to the tier. ΔXc is a direction-only
replacement at the realised norms. (k/N) counts folds > 0.

**Strong tier (top z_src tertile):**

| quantity | BULK-W (18) | BULK-ABC cancer (3) | ABC/kolf2 |
|---|---|---|---|
| share of realised energy a / truth energy d | 0.73 / 0.35 | 0.59 / 0.33 | 0.62 / 0.33 |
| X_real (tier) | 0.0644 | 0.0754 | 0.0456 |
| median cos realised / composition | 0.0366 / 0.0218 | 0.0733 / 0.0389 | 0.0078 / 0.0084 |
| additive contribution to cos(fin,T): base / cd4 / h1 / part | .0259 / .0078 / .0048 / .0001 | .0329 / .0158 / .0083 / .0008 | .0080 / .0062 / .0138 / .0003 |
| energy share base / cd4 / h1 / part | 0.80 / 0.13 / 0.04 / 0.005 | 0.65 / 0.22 / 0.06 / 0.009 | 0.66 / 0.23 / 0.06 / 0.005 |
| H1 carried (share of targets) | 0.21 | 0.33 | 0.32 |
| median cos to T: base / cd4 / h1 / **pre** / U / src0 / src1 | .019 / .011 / .033 / **.033** / .032 / .023 / .025 | .033 / .020 / .028 / **.046** / .041 / .032 / .033 | .005 / .009 / .036 / .005 / .005 / .006 / .008 |
| donor agreement cos(src0, src1) median | 0.020 | 0.032 | 0.025 |
| ΔXc drop base / cd4 / h1 / part | −.0141 (0/18) / −.0069 (0/18) / −.0059 (0/18) / −.0001 | −.0136 (0/3) / −.0093 (0/3) / −.0060 (0/3) / −.0001 | +.0039 / −.0018 / −.0227 / +.0010 |
| ΔXc pre only / U only | +.0118 (15/18) / +.0098 (16/18) | −.0030 (1/3) / −.0049 (0/3) | −.0236 / −.0233 |
| ΔXc 707-like U blend 0.25 / 0.5 / 1 | +.0068 (18/18) / +.0112 / +.0098 | +.0026 (3/3) / +.0027 / −.0049 | −.0029 / −.0081 / −.0233 |
| **ΔXc base → pre, β 0.5 / 1 (residuals kept)** | **+.0129 (18/18) / +.0203 (18/18)** | **+.0056 (3/3) / +.0092 (3/3)** | **+.0015 / +.0023** |
| in-sample LS oracle in span(parts) / span(all vectors) | +.0123 / +.0752 | +.0084 / +.0600 | +.0212 / +.0645 |
| knockdown median: sources / held line | −0.456 / −0.480 | −0.512 / −0.557 | −0.505 / — |
| literal pds (incumbent) | 0.863 | 0.944 | 0.724 |

**The same quantities on the weak and mid tiers** (`work/table810.md`):
- Removing H1 or CD4 also loses there.
- base → pre gains on composition there too (BULK-W +0.0055 / +0.0074; BULK-ABC +0.0066 / +0.0090), as 212 found
  globally. But the tiers' energy shares are small (0.11 / 0.16 on BULK-W), so in fold units the strong tier carries most
  of it.

**Fold-level, direction-only composition r** (`foldlevel810.py`; realised r is ≈ 0.71× (W) / 0.81× (ABC) of this, by 707's
own pair):

| construction | BULK-W | BULK-ABC cancer | kolf2 |
|---|---|---|---|
| 707 (U blend 0.25, strong) | +11.3 % (18/18) | +2.45 % (3/3) | −3.9 % |
| base → pre 0.5, strong | +21.0 % (18/18) | +5.3 % (3/3) | +2.1 % |
| base → pre 1, strong | +33.3 % (18/18) | +8.6 % (3/3) | +3.1 % |
| base → pre 1, strong ∩ agreement > median | +29.6 % (18/18) | +7.6 % (3/3) | +2.8 % |
| base → pre 1, ALL targets (= global unwhitening, 212's closed lever) | +41.9 % | +18.4 % | +3.9 % |

**Correlates within the strong tier** (Spearman, pooled over folds; `work/explore810-{L,ABC}.txt`):

| outcome | agreement | source strength (count) | held-line strength | key kv | z_src | knockdown (sources) | source noise |
|---|---|---|---|---|---|---|---|
| whitening loss cos(pre,T) − cos(base,T), W / ABC cancer | +0.16 / +0.14 | +0.22 / +0.16 | +0.19 / +0.12 | +0.19 / +0.18 | +0.20 / +0.23 | −0.09 / −0.15 | +0.11 / +0.04 |
| cos(fin,T), W | +0.45 | +0.43 | +0.58 | +0.63 | +0.40 | −0.14 | +0.07 |

**What it says.**
- The strong tier is already the best-directed tier: its realised cos is 3–6× the weak tier's, and it follows the held
  line's own strength (ρ 0.58).
- Its remaining direction error has three parts:
  1. a part the donors carry but the whitener removes (base → pre: +34 % / +15 % of tier-X at the composition level on
     BULK-W / BULK-ABC);
  2. a large part no fixed rule in the span of our vectors reaches. The per-target LS oracle over all vectors gives +0.075
     tier-X, in-sample: which source line resembles the held line for which target. Like 706's per-target amplitude
     oracle, nothing public predicts it yet;
  3. a line-specific remainder (the two sources agree with each other at only cos 0.02–0.03).
- Knockdown quality is not where the error is (|ρ| ≤ 0.15; source and held knockdown medians are within 0.03–0.05 of each
  other on every tier).

## 3. Proposals and verdicts
| file | sha256 | arms | verdict |
|---|---|---|---|
| `PROPOSAL.md` (v1; edited once at ~08:01Z while the judge read it, the judged text was `a421765b…`) | 888576e340975f3d3743aae60a259713c7af29244f4432da0733b1424d060b27 | P1 sgu05 / sgu1, P2 sga1 | **REVISE** (`VERDICT-1.md` b5d87daf607f19603e7d2956561fcbe44fc7f67c4ce794995250cbac538e425b) |
| `PROPOSAL-2.md` (frozen; sha256 in `work/SHA256SUMS.txt`) | — | sgu1 only (+ inc, 707, sgu1all comparators) | **pending (judge)** |

What changed for PROPOSAL-2 (VERDICT-1 items A–E):
- A: `sgu1all` (base → pre, β 1, every target with a donor row) runs in the same run as a known-answer comparator with
  its expected answer from 212 / 097 and a STOP rule; falsifier 2 is now E(sgu1) vs E(sgu1all) in the checker's units.
- B: the mechanism paragraph is corrected (§0 item 2 above); the gate is justified by pds exposure only.
- C: one arm (sgu1), charged N = 6 (e = 1.28 sd) because β 1 and the tertile gate were picked after composition reads of
  six forms on these folds; P2 and β 0.5 are withdrawn. The r projections are labelled as not blind.
- D: the decision rule (STOP / DE-H1 FIRST / CLOSE / HOLD-707 / DOSSIER) is written in PROPOSAL-2 §6 and as code.
- E: this file and the sha256 list below.
- Pricing code: `price810.py` (arms inc, 707, sgu0 KA, sgu1, sgu1all), `cloud/job810p.sh`, `decide810.py`; none has run on
  a candidate. decide810's own plumbing was checked on 804's EXISTING arrays (relabelled 705 / 710 rows in a scratch
  directory, no new pricing, no wear): its points reproduce 804's check800 points exactly (707 W +0.00485, 705 W +0.00092,
  705 ABC +0.00017).
- No DOSSIER.md: nothing is priced.

## 4. Tries and wear
- No check800 evaluation was made, so the wear ledger is untouched.
- Composition-level reads on these 22 folds that shaped the candidate (direction-only, not the instrument's price):
  `work/table810.md`'s 15 ablation columns × 3 tiers (drop-part ×4, U / pre / src0 / src1 only, U blends ×3, base → pre ×2,
  oracles ×2) and `foldlevel810`'s six fold-level forms (U 0.25 strong; base → pre 0.5 / 1 strong; base → pre 1
  strong ∩ agree; base → pre 1 / 0.5 all). sgu1 is therefore charged **N = 6** (PROPOSAL-2 §3).

## 5. Cloud (labcloud.sh; on-demand; Project=lab)
| run | type | what | outcome | actual |
|---|---|---|---|---|
| 20261016T071527Z-s810-decomp-a | r7i.8xlarge | verify + prep L1–L6 + 22 folds | prep rc 7 (801's manifests not shipped: needs `--code-extra experiments/801-lab-folds/work/manifest`); ABC compose PASS, realise KA-706 vs the json FAILED | $0.299 |
| 20261016T072436Z-s810-decomp-b | r7i.8xlarge | verify + prep L1–L6 + 18 L folds | DONE 0 (all KAs) | $0.336 |
| 20261016T073629Z-s810-dbg | r7i.xlarge | dev/ka706dbg.py on ABC k562 / hek | the committed dirswap706 path == the engine exactly; ≠ the json | $0.030 |
| 20261016T074418Z-s810-decomp-c | r7i.2xlarge | 4 ABC folds (KA-706 = two implementations) | DONE 0 | $0.067 |
| **total** | | | | **$0.732** |

- Run b used decomp810.py before the KA-706 edit; the L-fold code path is identical. Run c used sha256 `afe28579…`.
- The full part rows (base, cd4, h1, part, U, pre; 22 folds, 1.34–1.50 × 10^8 bytes each) are kept on S3 under
  `s3://<your-s3-bucket>/output/lab/810parts/<fold>-parts.npz`; keys, sizes and ETags in
  `work/s3parts810.tsv`. Each part set is pinned by the comp json's `fin_sha256` (the parts' sum) and its KA.
- `labcloud.sh status` at finish: no 810 instance running.

## 6. What is left
1. The judge's verdict on PROPOSAL-2.md. On GO: one EC2 run (`cloud/jobs-p.txt`, `cloud/sync810p.txt`, ≈ $1.5), then
   `decide810.py` and `check800.py aggregate` (sgu1 vs inc, family 810-sgu, N = 6), then DOSSIER.md only if PROPOSAL-2 §6
   says so.
2. 707's board read (10-17/18): resolves a HOLD-707 outcome by PROPOSAL-2 §6.
3. For the record: correct 706's ABC numbers to the committed script's +1.73 / +2.67 / +1.43 / −2.33 % (a lead's edit; this
   lane does not touch 706).

## 7. Files (sha256; full list incl. per-fold records in `work/SHA256SUMS.txt`)
| file | sha256 |
|---|---|
| decomp810.py (EC2: compose + realise) | afe285793b81685e38f1b0b96b0b35029d52d5f168a6f683e94e800af0b31359 |
| analyse810.py | 0b510e666e0997bb38d57e4de70060f126ed2ea1a9f4a571032218dd3f4a84b1 |
| table810.py | 21637a6c28845faa80454e90d06791060b5a626779345f3df4d735d07d290f71 |
| foldlevel810.py | 65b4d723209f94dae5d847387f10c5ec4df124c6444795f4865f33c3f583254d |
| explore810.py | cb9fe543bb9c3b5fb578c9812223ddc7dbbd4d5f04b98c8e7db40e1ec704091c |
| price810.py (not run) | da001560468d703fc9d08c004a2b24d1db60b6bc6f19fc49ff087d575f290f77 |
| decide810.py (not run on a candidate) | 724e805907814f7e5b6469c234aca574a47a42476207dc403fc15c641729ec14 |
| dev/ka706dbg.py | b6349476e34867c5f96f99437d4bd8fc68eeaa5aa27ee9c3a4ce32a7ebea090e |
| cloud/job810.sh | 583f1b35c06ecb469f9f7ba2fbde0ebf4bd5a87d7304751c706540913328d62a |
| cloud/job810p.sh | ea9e743fe6a4e4020a63abc10c29047ac0acc4918dbc47134cc435d5cdf7cd91 |
| cloud/jobs-p.txt | ee4fe6836785081e65d925d1b5428046cf1f50069b44793430838b03a5c04d14 |
| cloud/sync810.txt | 56527ed0bd9d6bf6dba331f09890914c4ca45cddff9e354d9c6e9f44e37b64b1 |
| cloud/sync810p.txt | c09135b5289bce2ac128727335a9d397ceb3def604af7fedb3a1e154f16eb41a |
| work/analyse810.json | 1d506cb211605ebaf21dcdc7da2e6d59050af86cf37e3e13d93d2e530dd55b35 |
| work/foldlevel810.json | 4c8d7f90dae18b2e8b5fc34284f7b3ac62ee86b5802c1e2113cbd71149171f4e |
| work/foldlevel810.txt | 82a5e4be5fdb010a650c8df3921bc5b92de2d67f3beb19feb02d5f9af95da864 |
| work/table810.md | 0305e1ba0284f98cdb6935f69a39d011ca6b60e35e109d5b255a1a59be468178 |
| work/explore810-L.txt | e0d9d2131c0cbf094d85a6a3b785938ab9ba4a2b1b60cd0d27a634cb3a666c07 |
| work/explore810-ABC.txt (A/B/C cancer folds, kolf2 out) | 2f3949ae73f5e87161371602173f651a9c64a38f43e3075f977c7034b834c21c |
| work/s3parts810.tsv (S3 parts: keys, sizes, ETags) | 211b868980725795b7f5f45a864e6bd294a41eee43bdfea34d8b01f5f88c432c |
| work/decomp/<fold>-{comp,real}.{npz,json} (22 folds, 4.7 MB) | in work/SHA256SUMS.txt |

## 8. The GO run (VERDICT-2, G1–G4): KA gate FAILED on 1 of 22 folds, so there is NO READ (PROPOSAL-2 §6 rule 0)
Written 08:41Z by the clock.

**Launches:**
- `20261016T081954Z-s810-price-a` (r7i.8xlarge, $0.191). 22 of 22 jobs died before pricing: the instance role may not
  read `s3://…/output/lab/810parts/` (403 on HeadObject).
- `20261016T082954Z-s810-price-b` (r7i.8xlarge, $0.272). The parts were first copied server-side to the input mirror
  (`input/lab/repo/experiments/810-strong-direction/work/parts/`, the same 22 keys and sizes as `work/s3parts810.tsv`).
  The run used `cloud/job810p2.sh`, which is `job810p.sh` with one change: it reads the parts from that synced path.
  Pricing code: frozen `price810.py` da001560….

**G3 count:**
- 21 folds carry all 14 required KA keys: inc == 804 ×4, 707 == 804 ×4, parts sum == F.fin, sgu0 == inc ×5. That is 294
  keys, none missing, worst 6.3e-13.
- The 22nd fold, ABC/orion_hek293t, aborted on `sgu0 == inc: a` = 1.364e-12, above the frozen threshold of 1e-12.
- The likely cause is numerical. sgu0 rebuilds each row as base + cd4 + h1 + part, and on this fold the parts' sum
  differs from F.fin by 3.6e-14 (the largest of the 22; decomp810's KA, RESULTS §1). `a` is an unnormalised ⟨P,P⟩ sum,
  so the gap shows there as an absolute difference just above 1e-12.
- The rule is binding as written: one failed KA means no read. Neither decide810 nor check800 aggregate was run, so the
  wear ledger is untouched. The 21 completed folds' arrays are kept unaggregated in `work/price-b-unread/`.

**Disclosure:**
- Each job's log prints per-arm summary lines before its asserts run.
- I saw ABC/orion_hek293t's lines while diagnosing the failure: sgu1 r +2.81 %, dpds −0.0263, dK +0.022; sgu1all r +15.4 %,
  dpds −0.0533, dK +0.085.
- I have read no other fold's numbers and aggregated nothing.

**Needed from the lead or judge (not done here; the code is frozen):** a ruling on the KA tolerance. Two options:
- re-run unchanged and accept a relative tolerance for the β-0 identity (|d| / max|ref| < 1e-12);
- or have sgu0 return F.fin's own rows at β 0.

Either needs one more launch of about $0.3. ABCp (G1) is not triggered, because no verdict exists.

**Files:**
| file | sha256 |
|---|---|
| `candidate-sgu1.json` (N = 6 declared; G4) | 2de495012c63312e20793d0b0fbaba94bf8d32d16d2c0528874bad9afacb6b30 |
| `cloud/job810p2.sh` | b03329fa2f92989d070f08cb2d2c7b12cae1d109ea82ef306d51e704b11fb96f |
| `cloud/jobs-p2.txt` | 3380173d239712fbed9d09d2d295ebc3b8422cea3b6692479262477eb80c67ac |
| `cloud/sync810p2.txt` | 2ffa814e86f57f0d5c296c45be11808a54f552af228b7671b2fa42c0fd895ff4 |

**G4 notes:**
- Closed list: 063 (`experiments/063-pds-mse-pareto/RESULTS.md:70`: residual de-whitening, a global β blend of the raw
  consensus, pds monotone down) stands beside 097, 099 §5.3, 170 and 212. sgu1all is the same-run stand-in for all of
  them.
- decide810's DOSSIER message says "D1-D4"; it means D1–D3 (PROPOSAL-2 §6). The code is left frozen.
- The "unexpected" flag was not applicable: there is no read.

**Cloud:** lane total $1.20 of $50 (§5 $0.732 + $0.191 + $0.272). No 810 instance is running.

## 9. The GO run, completed under RULING-KA.md, and decide810's verdict
Written 08:53Z by the clock.

**Steps:**
- `price810.py` was patched per `RULING-KA.md` (6c38cf53…): the β-0 identity on x and a is now read relative to max|ref|,
  and the summary lines print only after the KA gate.
  - Old sha256: da001560468d703fc9d08c004a2b24d1db60b6bc6f19fc49ff087d575f290f77.
  - New sha256: 9c33b24fae3bb9cfee135a7d217b2dd495ac238908e878f961914a6d4488fe36.
  - The arithmetic is unchanged; a byte copy of the old file is kept in the scratchpad.
- Only ABC/orion_hek293t was re-run: `20261016T084311Z-s810-price-c`, r7i.xlarge, DONE 0, $0.030.
- The 21 folds from run b were kept as they were. All 22 folds were assembled in `work/price-final/`.

**KA gate (G3):** 22 folds × 14 keys = **308 keys, 0 missing, all < 1e-12**.
- The 14 keys per fold are inc == 804 ×4, 707 == 804 ×4, parts sum == F.fin, and sgu0 == inc on x, a, pds, pdsl and K.
- Worst value: 6.3e-13 (absolute).
- The 21 earlier folds' x and a values were stored as absolute differences ≤ 6.3e-13, so they also pass the relative
  rule.
- ABC/orion_hek293t now reads sgu0 == inc: x 5.7e-15 (relative), a 3.6e-15 (relative), pds / pdsl / K exactly 0.

**decide810.py, frozen (724e8059…), on `work/price-final`** (`work/price-final/decide810.json`, sha256 6a8e2579…):

| set | arm | overall (literal pds) [95 %] | pdsl net | mse leg | pds leg | r | E = −pds/mse | folds + |
|---|---|---|---|---|---|---|---|---|
| BULK-ABC cancer | 707 | +0.00044 [−0.00087, +0.00176] | +0.00064 | +0.00142 | −0.00097 | +1.98 % | 0.69 | 2/3 |
| BULK-ABC cancer | **sgu1** | **−0.00200 [−0.00529, +0.00113]** | −0.00152 | +0.00345 | −0.00545 | +4.83 % | **1.58** | 2/3 |
| BULK-ABC cancer | sgu1all | −0.00987 [−0.01543, −0.00438] | −0.00958 | +0.00747 | −0.01733 | +10.44 % | 2.32 | 0/3 |
| BULK-W | 707 | +0.00485 [+0.00355, +0.00634] | +0.00505 | +0.00570 | −0.00085 | +7.97 % | 0.15 | 18/18 |
| BULK-W | **sgu1** | **+0.00894 [+0.00542, +0.01295]** | +0.00956 | +0.01468 | −0.00574 | +20.5 % | **0.39** | 18/18 |
| BULK-W | sgu1all | +0.00247 [−0.00186, +0.00733] | +0.00465 | +0.01708 | −0.01461 | +23.9 % | 0.86 | 9/18 |
| ABC/kolf2 | sgu1 | +0.00139 [−0.00377, +0.00665] | +0.00187 | +0.00195 | −0.00055 | +2.72 % | 0.28 | 1/1 |

Paired differences (overall):

| set | sgu1 − sgu1all [95 %] | P(E_sgu1 < E_sgu1all) | sgu1 − 707 [95 %] |
|---|---|---|---|
| BULK-ABC cancer | +0.00787 [+0.00312, +0.01235] | 0.86 | **−0.00244 [−0.00479, −0.00022]** |
| BULK-W (G1) | +0.00647 [+0.00338, +0.00924] | 1.000 | +0.00410 [+0.00166, +0.00672] |
| kolf2 | +0.00547 [−0.00345, +0.01418] | 0.62 | +0.00350 [−0.00024, +0.00725] |

**The ladder, in order:**
- **Rule 0, KA gate:** PASS.
- **Rule 1, STOP:** not triggered. sgu1all's BULK-ABC lower bound is −0.0154. The closed global form reads as 212 / 097
  predicted:
  - net −0.0099 against the expected ≤ +0.0005;
  - E 2.32 against the expected ≥ 0.95;
  - so the bulk checker SEES the γ knife-edge. That is a known-answer hit for the instrument on the class every direction
    proposal leans on.
  - The "unexpected" flag (G4) does not apply.
- **Rule 2, DE-H1 FIRST: triggered.** This is the frozen verdict.
  - sgu1: |dK_rel| > 5 % on 12/22 folds, all on BULK-W (worst 11.6 %, L3/orion_hek293t; median 6.0 %).
  - On the four ABC folds sgu1's dK is −0.6 % … +2.2 %.
  - sgu1all: > 5 % on 19/22 (worst 35 %). 707: none (worst 3.0 %).
  - By PROPOSAL-2 §6, the bulk engine cannot price the arm. The ladder stops here.
- **Rule 3, for the record only** (the frozen ladder did not reach it):
  - **C1 fails:** sgu1's BULK-ABC point is −0.00200.
  - **C4 fails:** sgu1 − 707 paired is −0.00244, and its CI excludes 0.
  - C2 holds (E 1.58 < 2.32), C3 holds (+0.0079), C5 holds (BULK-W +0.0089, 18/18).
- **G2 / HOLD-707:** not reached. 707's board read is still not in the ledger at 08:53Z.
- **G1:** not triggered. No DOSSIER.md.

**What it says:**
- **The gating hypothesis is confirmed as a mechanism.** On every surface, gating the unwhitened base to the strong
  tertile improves the pds-for-mse exchange over the closed global form:
  - BULK-ABC: E 1.58 vs 2.32, paired +0.0079;
  - BULK-W: E 0.39 vs 0.86, paired +0.0065.
- **It is not enough on Arc's own panel.** On BULK-ABC, pds still eats 158 % of sgu1's mse gain, and sgu1 loses to 707
  head to head. On the fresh panels it wins clearly: BULK-W +0.0089, and +0.0041 over 707.
  - This repeats the panel split seen for 707 itself (W SHIP, ABC WEAK) and for 438's CD4 dose.
  - The board's panel is the one where the strong targets' unwhitened directions crowd each other's truths: 246 §3.2's
    target-shared directions.
- **The DE-H1 FIRST outcome rests on BULK-W's ballast shift**, where the candidate already wins. A DE-H1 read could change
  the BULK-ABC decision only if the DE members added ≥ +0.0024 to sgu1 relative to 707 through ballast K, where sgu1's K
  moves ≤ 2.2 %.
  - **Recommendation:** do not fund DE-H1 for this family. Record it as **closed for the A/B/C-panel route** (C1 and C4 on
    the deciding surface).
  - Its BULK-W result (+0.0089, E 0.39) is a lead only if D/E/F's panel behaves like the L panels, which is 708's
    day-of question and not this lane's.
  - The lead decides.

**check800 aggregate** (ledger entry appended; `work/price-final/report-sgu1-vs-inc.json`, sha256 077ba2a9…):
- Pooled over all 22 folds: +0.00711 [+0.00427, +0.01048], 21/22 folds positive, leave-top-3-out +0.00004.
  - Its pooled bin is SHIP at N = 1. This pools BULK-W with BULK-ABC and is not the decision surface: decide810's
    per-set read rules.
- The ledger records N = 1. `candidate-sgu1.json` declares N = 6 (G4) and carries both price810 sha256s.

**Cloud:** price-c $0.030. **Lane total $1.23 of $50.** No 810 instance is running (`labcloud.sh status` 08:53Z).

**Files:**
| file | sha256 |
|---|---|
| price810.py (patched per RULING-KA) | 9c33b24fae3bb9cfee135a7d217b2dd495ac238908e878f961914a6d4488fe36 |
| candidate-sgu1.json | 1f388403dc6689b5608c36c7eebe2824a485e13170183e18c98b45eaa97ee8bd |
| cloud/jobs-p3.txt | ea9158111bc1599e6d75fbe1150800f3c7ba5a7555cc12d2f54debb7531739fe |
| work/price-final/decide810.json | 6a8e257938cb495004c8e3105df061160886bd8ff2964d7a962da5e85f0c2a15 |
| work/price-final/report-sgu1-vs-inc.json | 077ba2a9981a0cd2da1001d91e20c9581cb186baefa6b21290261bba7ec785d5 |
| RULING-KA.md / VERDICT-2.md | 6c38cf53… / 63899f9f… |

## 10. ABCp information read (LEAD-DECISION-1; a read, not a reopening)
Written 09:34Z by the clock.

**No re-tune.** Before any number was read, the sha256 of both scripts were recorded and match the GO run's final code:
- price810.py 9c33b24fae3bb9cfee135a7d217b2dd495ac238908e878f961914a6d4488fe36;
- decide810.py 724e805907814f7e5b6469c234aca574a47a42476207dc403fc15c641729ec14.

**Folds and incumbent.**
- The folds are 801's registered `ABCp/k562_gwps`, `ABCp/orion_hct116` and `ABCp/orion_hek293t` (status abc-dayops; day
  operators 48cb88b6…, a0e27947…, 5d982462…). There is no ABCp kolf2.
- `F.panel == "ABC"`, so price810's A/B/C branch applies. The incumbent used is **cap6+P2g(abc)** (ALLSW's bulk), the same as
  on 062's folds.

**Run.** `20261016T092043Z-s810-abcp-a` (r7i.2xlarge, $0.086). Per fold, `cloud/job810abcp.sh` runs:
1. `decompabcp810.py <held> 062`: 801's dump518p route (lib420), parts split by compose_arm's algebra;
2. `decompabcp810.py <held> day`: the same with the day operator;
3. `price810.py ABCp/<held>`.

**Known answers (all < 1e-12; the read was allowed):**
- 062 route:
  - parts sum == compose_arm fST ≤ 2.2e-16;
  - fST == 518's stored dump ≤ 3.6e-14;
  - base, cd4, h1, part, U and pre == the GO run's ABC parts **bit for bit (0.0)** on all three folds. lib420's route and
    decomp810's lib801 route agree.
- Day route: parts sum == fST ≤ 4.4e-16; fST == 801's ABCp dump fST ≤ 1.8e-15.
- `E.load("ABCp/<held>", check=True)`: its load KA (PT0) asserted and passed (rc 0).
- price810 KAs: parts sum == F.fin, and sgu0 == inc on x and a (relative), pds, pdsl and K.
- 804 has no ABCp arrays, so price810's built-in 804 KA had nothing to compare against. In its place, inc and 707 are
  compared with 801's ka-abcp arrays (ka804 unmodified, `../801-lab-folds/work/ka-abcp/bulk/`) on x, a, pds and pdsl.
- **3 folds × 14 keys = 42 keys, 0 missing, worst 2.8e-13** (`work/abcp/ka-abcp810.json`).

**Results** (overall board units, literal pds deciding; E = pds share of the mse gain = −pds leg / mse leg):

| arm | BULK-ABCp (day operator): point [95 %] | folds > 0 | E | pdsl net | BULK-ABC (062 operator): point [95 %] | folds > 0 | E |
|---|---|---|---|---|---|---|---|
| 707 | **+0.00178 [+0.00050, +0.00316]** | 3/3 | 0.21 | +0.00208 | +0.00044 [−0.00087, +0.00176] | 2/3 | 0.69 |
| sgu1 | **+0.00012 [−0.00317, +0.00346]** | 2/3 | 0.97 | +0.00085 | −0.00200 [−0.00529, +0.00113] | 2/3 | 1.58 |
| sgu1all | −0.00488 [−0.01058, +0.00090] | 1/3 | 1.52 | −0.00368 | −0.00987 [−0.01543, −0.00438] | 0/3 | 2.32 |
| sgu1 − sgu1all (paired) | +0.00500 [−0.000002, +0.00981]; P(E_sgu1 < E_sgu1all) 0.92 | | | | +0.00787 [+0.00312, +0.01235]; P 0.86 | | |
| sgu1 − 707 (paired) | −0.00166 [−0.00402, +0.00071] | | | | −0.00244 [−0.00479, −0.00022] | | |

- mse legs on ABCp: 707 +0.00226 (r +3.17 %); sgu1 +0.00463 (r +6.47 %); sgu1all +0.00935 (r +13.07 %).
- dK_rel on ABCp:
  - sgu1: +0.6 % (k562), +4.9 % (hct116), **+7.4 % (hek293t)**;
  - 707: ≤ 2.1 %;
  - sgu1all: up to 24.8 %.

**decide810 applied as frozen to ABCp.**
- Setup: decide810 hard-codes the A/B/C fold names. Its run directory `work/price-abcp-decide/` therefore holds the ABCp
  folds under the names `ABC__<held>`, plus the 18 L folds and 062's kolf2 (reported only), all unchanged copies.
- Verdict (`decide810.json` b89348d0…): **DE-H1 FIRST**. The trigger is BULK-W's dK, as before, and now also ABCp's
  hek293t fold (7.4 %).
- STOP did not fire: sgu1all's lower bound is −0.0106, and the closure still reads (E 1.52 ≥ 0.95).
- Past that gate, the ladder would **CLOSE on C4**: sgu1 − 707 is −0.00166 on ABCp too.
  - C1 (+0.00012 > 0), C2 (0.97 < 1.52), C3 (+0.0050) and C5 hold.
  - D1 fails regardless: the point after e(6) is −0.0020.
- G1 is not met, so there is **no DOSSIER**. The family stays closed for the A/B/C-panel route.

**The answer to the question this read buys:**
- **Yes, partly.** With A/B/C's operator replaced by the day operator, sgu1's BULK-ABC read rises from −0.0020 to +0.0001,
  and its pds share falls from 1.58 to 0.97. The W/ABC split shrinks from 0.0109 to 0.0088, so about a fifth of it was 062's
  whitener.
- **The rest is the panel.** The day operator lifts 707 by as much (+0.0004 → +0.0018), so sgu1 still loses to 707 there.

**Wear rows** (check800 aggregate on `work/abcp/`, one per arm; reports `work/abcp/report-<arm>-vs-inc.json`):

| arm | family | checker bin on ABCp (its own N) | N declared (BULK-ABC reads of this arm) |
|---|---|---|---|
| sgu1 | 810-sgu (`candidate-sgu1-abcp.json` 532553cc…) | HARM after best-of-6 (ledger N 6) | **7** = the N = 6 selection branch, charged with the 062 read, + this ABCp read |
| 707 | 707-swap | SHIP (N 1) | **2** = 804's 062 read + this ABCp read. 801's ka-abcp holds the same 707 arrays to 2.8e-13: not a separate configuration |
| sgu1all | 810-sgu1all-ka | HARM (N 1) | **2** = the GO run's 062 read + this ABCp read (a known-answer arm, not a try of 810-sgu) |

**One by-product for the lead:** 707 itself reads SHIP on ABCp: +0.00178 [+0.0005, +0.0032], 3/3 folds, E 0.21, against
WEAK on 062. Its board read is pending (10-17/18). This is the day-recipe instrument's prospective call for 707, made after the
fact here and labelled as such.

**Cloud:** $0.086. **Lane total $1.32 of $50.** No 810 instance is running (`labcloud.sh status` 09:34Z).

**Files:**
| file | sha256 |
|---|---|
| decompabcp810.py | 245240b3de51f371f103d2c483732cefea9b33d23a2776968874bca2c2038655 |
| cloud/job810abcp.sh | 8ab5c3140e66ab75a5b2f73d54d4f44bdbabb1f1d8c5d51edb37584ed53f0554 |
| cloud/jobs-abcp.txt | 2a615df6e4571d0aea34de7ccb124e492e30aa36c4d96098555ca953b205c202 |
| cloud/sync810abcp.txt | 44080a89b0b3ddf3d9185b049e107524f5fd9031a8cd36d296ccf5416eee4089 |
| candidate-sgu1-abcp.json | 532553cc9f36a58569014505d980d5430b97e1f76bad167f3c18c6ee1b2e5cc3 |
| work/abcp/ka-abcp810.json | a310dc7186a94c671ec4a5fc89c661c0d3a7f10284cb587b5909a69d8a09e8d2 |
| work/abcp/report-{sgu1,707,sgu1all}-vs-inc.json | b16432d0… / 54a76ee5… / c0ae7111… |
| work/price-abcp-decide/decide810.json | b89348d078f79aa0eb64819ae7e52d5d48680b284e970349096867ac1ca833d4 |
| work/abcp/parts/ABCp__{k562_gwps,orion_hct116,orion_hek293t}-parts.json (fin_sha256 231decac… / 7f469f16… / 5b1ead19…) | b0866c7a… / c87c7295… / 240f39b2… |

The part rows themselves (npz) stayed on the instance and were not uploaded. They are rebuilt by `decompabcp810.py` in about
30 s per fold.

## 11. LABC instrument read (lab lead, 10-16 ~10:55Z; an instrument read, not a reopening)
Written 11:10Z by the clock.

**What was read.** Three arms, sgu1, sgu1all and 707, each against inc, on 801's `LABC/<held>` folds:
- the folds are A/B/C's 300 targets built with the fresh-fold (L) recipe and D/E/F rules (registry status `labc`);
- inc = pct+P2g, because `F.panel` is "LABC", so price810's non-ABC branch applies;
- the code is price810.py **9c33b24fae3bb9cfee135a7d217b2dd495ac238908e878f961914a6d4488fe36**, the same as the GO run and the
  ABCp read. No re-tune.

**Run.** `20261016T105743Z-s810-labc-a` (r7i.2xlarge, DONE 0, **$0.080**).
- `cloud/job810labc.sh`, run on EC2 only, because the Mac was under the coordinator's local-job hold:
  - It arranges the LABC inputs in the EC2 checkout: data/challenge's bundle files; A/B/C's production CD4 banks, which are
    byte-identical to 801's LABC banks (801 §14); and ABCp's day operators, the ones the LABC dumps used.
  - Then, per fold, it runs `decomp810.py compose` (afe28579…) and price810.
  - Then labc810.py's KA gate, then `check800.py aggregate --no-wear` for each arm.

**Known answers: all pass. The read was allowed.**
- decomp810: the parts sum to the LABC dump's `fin` at ≤ 8.9e-16. compose(final) == the dump's fin and h1c1 == fPP, both
  exactly 0.0, which proves the input arrangement.
- engine801's load KA (PT0) passed.
- labc810.py gate: **3 folds × 14 keys = 42 keys, 0 missing, worst 1.05e-15** (`work/labc/ka-labc810.json`). The keys per
  fold:
  - parts sum == F.fin, and sgu0 == inc on x, a, pds, pdsl, K;
  - inc and 707 == 801's ka-labc arrays (ka804 unmodified) on x, a, pds, pdsl. 804 has no LABC arrays.

**Results** (check800, paired target bootstrap 4,000; overall = mse leg + literal pds leg, in board units):

| arm vs inc | LABC point [95 %] | folds > 0 | mse leg | pds leg | pds share of the mse gain | per fold (k562 / hct116 / hek293t) |
|---|---|---|---|---|---|---|
| sgu1 | **−0.00015 [−0.00367, +0.00342]** | 2/3 | +0.00494 | −0.00508 | 1.03 | +0.00096 / +0.00337 / −0.00477 |
| sgu1all | −0.00605 [−0.01183, +0.00015] | 1/3 | +0.00982 | −0.01587 | 1.62 | −0.00521 / −0.01361 / +0.00068 |
| 707 | +0.00115 [−0.00012, +0.00252] | 3/3 | +0.00187 | −0.00073 | 0.39 | +0.00104 / +0.00231 / +0.00008 |

- check800 labels this block "BULK-W" because it treats every non-ABC fold that way. These are the three LABC folds only.
- The arrays for panelrank800 / overread800 are at `work/labc/bulk/LABC__<held>/{inc,sgu1,sgu1all,707}.npz`. They use
  price810's standard-format keys, the same as `work/price-final`.

**Beside the other A/B/C-panel reads:**

| arm | BULK-ABC (062) | BULK-ABCp | LABC | BULK-W |
|---|---|---|---|---|
| sgu1 | −0.00200 | +0.00012 | −0.00015 | +0.00894 |
| sgu1all | −0.00987 | −0.00488 | −0.00605 | +0.00247 |
| 707 | +0.00044 | +0.00178 | +0.00115 | +0.00485 |

LABC lands with ABCp, not with the L folds. With the panel held fixed, the L recipe and rules leave sgu1 at about 0 and 707
small and positive. The gap to BULK-W is carried by the panel's targets, as 801 §14 found for the key arm.

**Wear.**
- check800 ran on EC2 with `--no-wear`, because the local hold forbade running it on the Mac.
- The three wear rows are in `work/wear810-labc.jsonl` (sha256 91bebb16…), in the ledger's field format, for the lead to
  append to 800's ledger. N, counted as A/B/C-panel reads of the arm including this one:

| arm | family | N | made up of |
|---|---|---|---|
| sgu1 | 810-sgu | 8 | the 6-try selection branch + 062 + ABCp + LABC |
| 707 | 707-swap | 3 | 062 + ABCp + LABC |
| sgu1all | 810-sgu1all-ka | 3 | 062 + ABCp + LABC |

**Cloud:** $0.080 (the launch cap was $0.46; the actual cost is under the $0.3 grant). **Lane total $1.40 of $50.** No 810 instance is
running.

**Files:**
| file | sha256 |
|---|---|
| labc810.py | d8decc952f6a3fa296e75ca342e2cb4cd02c9bdcfadc2e18c6dde569978b1b47 |
| cloud/job810labc.sh | bb191d9328132df5172c7862c48aa7f438b6ad3e92bc7160a368f260661a8e95 |
| cloud/jobs-labc.txt | 3313e5cfb735edc98d38d3f7db9897838a293ef18e22b1a9541acfdefa5ee37e |
| cloud/sync810labc.txt | e2ba431b26e4fd187d3a074e058730237f85213c555acbd63fda92149ce0b1b9 |
| work/labc/ka-labc810.json | a26d1da61aa64784034125500d76067987bf13d6a6d8296868053b379030672f |
| work/labc/report-{sgu1,sgu1all,707}-vs-inc.json | f1be0bb7… / bd8e8f40… / cfb97500… |
| work/wear810-labc.jsonl | 91bebb168fe733ed1ec46ee7c8da63867449861d8220922a71131d1baa16b8dc |
| work/labc/bulk/*/*.npz | in work/SHA256SUMS.txt |

## 12. A-only instrument read (lab lead, 10-16 12:44Z; an instrument read, not a reopening)
Written 13:00Z by the clock.

**What was read.**
- **Arms:** sgu1 and 707 vs inc on 801's 18 `A<k>/<held>` folds. These are the A-only arc800 draws A1–A6 with the
  essential screens removed, on the L recipe, D/E/F rules, inc = pct+P2g (registry status `a-only`).
- **Code:** price810.py **9c33b24fae3bb9cfee135a7d217b2dd495ac238908e878f961914a6d4488fe36** (unchanged; no re-tune) and
  decomp810.py afe28579….
- **Run:** `20261016T124635Z-s810-aonly-a`, r7i.8xlarge, 24/24 jobs rc 0, **$0.337**. The launch cap was $1.97; the actual
  cost is within the $1.5 grant.
- **Job script:** `cloud/job810aonly.sh`. Per panel, `prep` arranges 801's L-recipe inputs in the EC2 checkout:
  - 199's def-bundle controls, sha256-checked;
  - `cd4_801.py` (it fell back to the relaxed gates as 801's own extraction did; coverage 283–284);
  - frames801.py;
  - operators from a server-side S3 copy of `output/lab/ops801/A<k>` into this lab's mirror
    (`input/lab/repo/experiments/810-strong-direction/work/aops/`), each checked against 801's registry sha256
    (`work/aops/registry-aonly.json`).

  801's S3 write probe was not used. Per fold, `fold` runs decomp810 compose, then price810.
- **Aggregation:** check800, run locally (seconds) with `--no-wear`.

**Known answers: all pass. The read was allowed** (`aonly810.py`, `work/aonly/ka-aonly810.json`).
- 18 folds × 16 keys = **288 keys, 0 missing, worst 4.2e-15**. The keys per fold:
  - decomp810's parts sum == the A dump's fin, and h1c1 == fPP. These prove the EC2-side CD4 extraction and operators equal
    what 801's dumps used.
  - price810's parts sum == F.fin, and sgu0 == inc on x, a, pds, pdsl and K.
  - inc and 707 == 801's ka-aonly arrays (ka804 unmodified) on x, a, pds and pdsl.
- engine801's load KA (PT0) passed on every fold.

**Results** (check800, paired target bootstrap 4,000; overall = mse leg + literal pds leg, in board units; 18 folds pooled
with equal weight):

| arm vs inc | A-only pooled point [95 %] | folds > 0 | mse leg | pds leg | pds share of the mse gain | r | leave-top-3-out | per held line (k562 / hct116 / hek293t; folds > 0) |
|---|---|---|---|---|---|---|---|---|
| sgu1 | **+0.00869 [+0.00443, +0.01297]** | 17/18 | +0.01556 | −0.00687 | 0.44 | +21.8 % | −0.00039 | +0.0018 (5/6) / +0.0059 (6/6) / +0.0184 (6/6) |
| 707 | **+0.00472 [+0.00324, +0.00624]** | 17/18 | +0.00598 | −0.00126 | 0.21 | +8.4 % | +0.00168 | +0.0023 (5/6) / +0.0034 (6/6) / +0.0085 (6/6) |
| sgu1 vs 707 (paired, extra) | +0.00301 [+0.00015, +0.00544] | 13/18 | | | | | −0.00269 | |

- check800 labels the block "BULK-W" because it treats every non-ABC fold that way. These are the 18 A folds only.
- Arrays: `work/aonly/bulk/A<k>__<held>/{inc,sgu1,707}.npz`, using price810's standard-format keys, as in `work/price-final`.
  Reports: `work/aonly/report-{sgu1,707}-vs-inc.json` and `report-sgu1-vs-707.json`.

**Beside the other panels** (pooled points):

| arm | BULK-ABC (062) | BULK-ABCp | LABC | A-only (18) | BULK-W (L1–L6, 18) |
|---|---|---|---|---|---|
| sgu1 | −0.00200 | +0.00012 | −0.00015 | +0.00869 | +0.00894 |
| 707 | +0.00044 | +0.00178 | +0.00115 | +0.00472 | +0.00485 |

Removing the essential screens from the arc800 draw does not move the fresh panels toward A/B/C on these direction arms.
A-only reads like BULK-W (sgu1 +0.0087 vs +0.0089; 707 +0.0047 vs +0.0049). sgu1's gain is concentrated, though: its
leave-top-3-out mean per fold is −0.0004, against +0.0017 for 707.

**Wear.** Two rows in `work/wear810-aonly.jsonl`, in the ledger's field format, for the lead to merge:
- sgu1: family 810-sgu, N = 1 on the A-only surfaces. The family's N = 6 selection charge was made on BULK-W / BULK-ABC.
- 707: family 707-swap, N = 1. 801's ka-aonly holds the same arrays, so it is not a separate configuration.

**Cloud:** $0.337. **Lane total $1.74 of $50.** No 810 instance is running.

**Files:**
| file | sha256 |
|---|---|
| aonly810.py | fc0a3778b306377ea0e0cfd0742df02f776e379ebc2a90e74337e4e21755b23b |
| cloud/job810aonly.sh | 936f79e8d0cbc1794a7dc55890c70a82cbd129d2e28f9ef559315f6944b0261b |
| cloud/jobs-aonly.txt | bea99a3b27a6d5c7ba4579472d60cec62cf5a8c35fe23767cba875f9a297721c |
| cloud/sync810aonly.txt | 5c48228a42fbd95f1b05ea1a800d1b5ffee6126137ddb33c51bf620e271812a7 |
| work/aops/registry-aonly.json | 2cc7f23d4a5fdd4a73021fe8107c756803a9bffa2b107fde4d93800bbb947b8a |
| work/aonly/ka-aonly810.json | 5f10758bb6fddc2c14a8f56027d90519b84a924d894d332936e5b1cfe822bfeb |
| work/aonly/report-sgu1-vs-inc.json / report-707-vs-inc.json / report-sgu1-vs-707.json | 2aa478f1… / 3f780ad0… / 342927aa… |
| work/wear810-aonly.jsonl, work/aonly/bulk/*/*.npz | in work/SHA256SUMS.txt |
