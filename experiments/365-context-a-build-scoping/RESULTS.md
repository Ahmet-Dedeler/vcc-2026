# 365 — A-only arm is scopable to one insertion point + one flag; `replace` vacuous on A/B/C

Status: **COMPLETE** (read-only scoping; no code, no builds, no slots).

## Known-answer reproduction (own code, this session)

| check | record | remeasured | match |
|---|---|---|---|
| nadig_jurkat median row norm | `meta.json` 12.371221542358398 | 12.371221542358398 (numpy, mmap) | **exact** |
| nadig_jurkat ∩ A/B/C panel | 0/300 (225 §8.3, 235) | 0/300 | **exact** |
| nadig_jurkat genes ∩ official axis | 8,284/18,533 (362 PLAN cites 8,284 mapped) | 8,284 | **exact** |
| bank shape/dtype/finiteness | 2270×8881 float32 | same, all finite | **exact** |
| 323 raw members | `metrics-cap15.json` | pds 0.88964 / mse 0.97495 / jac 0.07476 / **nmae 0.95084** / fid 0.49755 / reach 0.25435 | read OK |
| 343 headline | NMAE 0.950836, sign 0.588/0.600 | file re-read, values as cited | read OK |

(Note: `wc -l` reports 8880/2269 on the bank's txt files — missing trailing
newline, not missing rows; Python split gives the meta's 8881/2270.)

## What was measured
Read-only reads: builder composition flow (`build_v9_pre_whiten_effects` :1410 →
PM :1598 → whiten :1864 → CD4 merged add + preserve :2700–2860 →
`build_production_effects` :2943 → `run()` :4792 → writer :4370, per-context branch
:4502–4509, knee precedent :5315–5427); `def_pipeline.py` (`day_env`, arms guard
:171–248, postbuild `--allow-config-diff whitener_id` :614); `def_postbuild.py`
(`ALLOWED_CONFIG_DIFFS` :81, incumbent = 249 meta :120); `run319p.sh` (319p = cap
−1.5 + partial 0.34); `patch254.py` + 254 RESULTS §1 (inert-flag pattern);
DEF-RUNBOOK §§1–4; 346/362/363 PLANs (the three constructions). Light probes only
(< 1 GiB RSS, numpy/pandas mmap).

## Results
Delivered: `PLAN.md` (mandate, inputs, falsifiers) and `DESIGN.md` (§0 bank facts +
constraints C1/C2; §1 insertion point `:5325` with before/after normalisation
enumeration; §2 `VCC_CTX_A_PRIOR` flag, artifact format, loader validation, three
applicators; §3 invariant table; §4 test plan; §5 80–95 min estimate).

| design question | answer |
|---|---|
| Insertion point | `run()` :5325, `context_effects={"A":…, "B":effects, "C":effects}` → writer :5433; knee precedent |
| Normalisations acting on correction | preserve NO; A's CPM gate YES; render+rc_callset YES; panel asserts upstream (untouched) |
| Inertness proof | 254 KA-1 pattern (staged + package sha vs 319p meta) |
| Config-diff shape | new `ctx_a_prior*` keys only (+ consequences); `--allow-config-diff` extended; arms guard extended |
| `replace` on A/B/C | **vacuous (0/300 carry) — builder must raise, never silent-zero** |
| Time | ≈ 80–95 min decision→upload; fits iff gate opens immediately + patch pre-staged |

## Verdict against the PLAN.md kill criterion
Not killed: the insertion point is statable as one file:line with all
normalisations before/after enumerated, the arm confines to one config key plus
consequences, and nothing touches shared whitener/preserve/emitter semantics. The
C1 falsifier fired as designed — against `replace`-on-A/B/C specifically, not the
doc: it is recorded as a raise-requirement, and `sign_block`/`add_vector` are
unaffected (target-independent).

## Limits
- H1-panel reasoning only where cited; no A/B/C forecast is made (per mandate §5).
- No retrieval/cosine/proxy gain is claimed as a score gain; modes are "engineering
  feasible", gated by 346/362/363 — if none clears its gate, build nothing.
- Line numbers are HEAD-as-read; the implementer re-verifies anchors (patch script
  must use exact-anchor `replace_once`, 254 pattern) and hash-gates builder bytes.
- "100 of A" in the mandate's artifact sketch is a placeholder: on A/B/C the carry
  is 0, on October it is re-derived from that panel — never hardcoded.

## Exact commands to reproduce
```bash
memory_pressure | tail -1
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 python3 -c "
import numpy as np, pandas as pd, json
panel = set(pd.read_csv('data/challenge/pert_counts.csv')['target_gene'].astype(str))
print('overlap:', len(set(open('data/donors/nadig_jurkat/perts.txt').read().split()) & panel), '/300')
d = np.load('data/donors/nadig_jurkat/delta.npy', mmap_mode='r')
print('shape:', d.shape, 'median norm:', float(np.median(np.linalg.norm(d, axis=1))))
print('finite:', bool(np.isfinite(d).all()))"
python3 -c "import json; print(json.load(open('experiments/323-current-family-validation/work/metrics-cap15.json'))['raw'])"
grep -n "context_effects = None" scripts/make_submission_pm_cd4.py
grep -n "ALLOWED_CONFIG_DIFFS = " scripts/def_postbuild.py
```

## Peak memory and wall time
Peak RSS < 1 GiB (mmap + pandas CSV reads only). Wall time ≈ 35 min (reads +
probes + three docs). No scratch files retained; nothing written outside
`experiments/365-context-a-build-scoping/` (contains `PLAN.md`, `DESIGN.md`,
`RESULTS.md`). No subprocess over 4 GiB was started (gate check: 72% free at start).
