# 461 — the October path with Arc's denoised programme (`DEF_DAY_PROGRAMME_SOURCE=h1_2025_full`, flags `pt`)

Rehearsal 461dr started 11:58:54Z 09-26 (`chain461.sh`), after 458d's build.
- def_sameline: `--programme-source h1_2025_full` gave D/E/F the programme (275 rows, coverage 0.975); mode ON for
  3 contexts with `auto:pt`.
- The dump's bulk and carrier arrays are identical to 453dr's. The files differ only in the recorded env string: the
  dump path, and q 1.25 in 453's record because it ran before 455 moved the driver default to 1.75. The dump is taken
  before sharpening, so the re-level should retrace 453's.

## Result (461dr, 11:58Z → 13:20Z): PASS
- The re-level retraced 453dr's **exactly**:
  - lvl0 ρ 0.8407 / 0.7703 / 0.7346;
  - lvl1 at 1.0793 / 1.1275 / 1.1546;
  - CONVERGED at lvl2, **D 1.0859 / E 1.1383 / F 1.1662** (453's factors), ρ 0.9793 / 0.9794 / 0.9795.
- So the programme switch leaves the bulk untouched, bit for bit, as designed.
- Zero-gene audit 0 in 0 blocks (lfc < −5: 17; > +3: 3, max +3.95). `vcc prep` verified 360,000 cells. The builder was
  restored (`4feae27d`).
- Package `submissions/461dr-october-path-arcprog-abc-q1.75.vcc`, sha256 `f2191696…`, novel.
- `gate461.py`: **ALL GATES PASS**; the config diff vs 453dr is exactly `sameline`. The first run failed on
  `bulk_load_sha256`: the dump file records its env (run path; 453's recorded q 1.25, before 455). The gate now accepts
  that only when the dump arrays are identical, which they are. `work/gates461.log` is written.
- It is a 10-01 candidate only if 458d reads positive (`464-slot-1001/chain1001.sh`).
