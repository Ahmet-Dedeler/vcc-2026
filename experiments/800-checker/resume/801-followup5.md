FROM THE LAB LEAD (2026-10-16 11:49Z by the clock): follow-up 5. This is the corrected default draw rule ("A-only"), built under the
DESIGN §3.1 amendment of today. Budget ≤ $5, EC2 only (the local-job HOLD may still be on; window800 tells you).
Deadline: 10-17 18:00Z.

**Why.** arc800 removed only K562's essential list. FINDINGS §8.1g (Arc's stated criterion) also requires removing the RPE1,
HepG2 and Jurkat essential screens. Your L1–L6 carry 14–24 such targets each, and 814 shows they carry most of the
direction over-read (`experiments/814-panel-rule/RESULTS.md`; DESIGN §3.1 amendment).

**Task.**
1. **Draw six panels "A1–A6"** with draw801's arc800 procedure plus ONE change: the pool also removes
   `data/donors/rpe1/perts.txt` ∪ `data/donors/nadig_jurkat/perts.txt`. If 814 names a different HepG2 list before you
   draw, add it; otherwise note that it is absent.
   - Independent draws from the same full pool (overlap with each other and with L1–L6 is allowed and recorded).
   - Seeds 801301–801306. Never L7/L8's seeds, and never their targets.
   - Asserts at draw time: 0 K562-essential, 0 essential-screen, every quota met.
   - Record each panel's sha256, its overlap with L1–L6 and with A/B/C, and the strength-bin shares, BEFORE any fold is built.
2. **Build the 18 folds** (A1–A6 × K562 GWPS / Orion HCT116 / Orion HEK293T) with your L-fold recipe unchanged: CD4 banks, a
   fresh day-recipe 3-donor operator fit per panel and held line, dump801, the engine KAs. Register them as `A<k>/<held>`,
   status `a-only`.
3. **Run ka804's L arm set** on all 18 (unmodified, `KA804_OUT` outside the repo). Write standard arrays to
   `work/ka-aonly/{bulk,meta}/` and mirror them to S3. I then ask 810 (sgu1) and 811 (c0_l25) for their arms.

**Preregistered reading (the lead's bar, written now, before any number):**
- panelrank800 with A1–A6 in place of L1–L6. **PASS if ABCp's mean percentile is ≥ 0.35 and its bottom count ≤ 6/21**; LABC
  is reported beside it.
- Per-class over-read factors are re-measured with overread800.
- A class whose ABCp percentile stays < 0.25 keeps a factor > 1.

Report in RESULTS.md §15 with runs and $, then stop.
