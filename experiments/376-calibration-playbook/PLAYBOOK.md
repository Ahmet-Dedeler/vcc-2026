# PLAYBOOK — same-line calibration artifact for a newly revealed context (Oct 22)

Scope: context X is revealed as cell line L, and we hold a public genome-wide
CRISPRi bank for L. Build a 373-`blend_rows` calibration artifact for X in minutes,
hand it to the PI for builder wiring. This playbook builds the artifact only —
it never runs the production builder (rule 1).

Tool: `experiments/376-calibration-playbook/build_calibration_artifact.py`
(estimator primitives imported verbatim from 341/334; validated bit-identical on
the Jurkat known answer, 376 RESULTS.md K1).

## 0. Identify the bank (5 min, human)
Map L to a bank directory with `delta.npy / perts.txt / genes.txt / b_ctrl.npy`:
`data/donors/` (k562_gwps, orion_hct116, orion_hek293t, nadig_jurkat, rpe1,
k562_essential) or `experiments/166-donor-acquisition/donors/` (nadig_hepg2).
No bank for L → stop; this playbook does not apply (HeLa and CAL33 have none —
that is why A got one and B/C did not). Check the per-line worth in 376
RESULTS.md §"Results" before spending anything: HepG2 +0.057 and Jurkat +0.042
( clean nulls) are the good cases; HCT116 +0.007 is marginal; RPE1's +0.100 is
~half shared-component (shuffled +0.055 — read the caveat, do not quote +0.100).

## 1. Preview the axes, no fit (~5 s)
```
memory_pressure | tail -1   # require >=50% free
uv run python - <<'PY'
import sys; sys.path.insert(0,'experiments/376-calibration-playbook')
from pathlib import Path
import build_calibration_artifact as B
ax = B.build_axes(Path('data/donors/<L>').resolve(),
                  [Path(d).resolve() for d in
                   ('data/donors/k562_gwps','data/donors/orion_hct116',
                    'data/donors/orion_hek293t')],
                  Path('data/challenge/pert_counts.csv'),
                  Path('data/challenge/gene_names.csv'))
print({k: len(v) if isinstance(v,list) else v for k,v in ax.items()
       if k not in ('canonical','r_perts','r_genes','panel')})
print('recipient panel overlap:', ax['recip_panel_overlap'])
PY
```
Read: `training` (need ≥ a few hundred; Jurkat had 2,006), `covered` (fraction of
canonical that the artifact can touch; Jurkat 4,951 = 26.7%), `test` (panel
targets with full donor consensus; Jurkat 259/300), `recip_panel_overlap`
(excluded from training automatically — K562 carries 269, that is normal, not an
error). If the recipient IS one of the three donors it is dropped automatically
(2 donors remain; K562/HCT116/HEK293T rows of the value table are exactly this
case). Expect `covered` to be small when the recipient's gene universe is small
or training is huge (K562-as-recipient: 682 genes — the map still builds, it just
touches 3.7% of the axis).

## 2. Build the artifact (~20 s–2 min, ≤12 GiB)
```
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python \
  experiments/376-calibration-playbook/build_calibration_artifact.py \
  --recipient data/donors/<L> \
  --donors data/donors/k562_gwps data/donors/orion_hct116 data/donors/orion_hek293t \
  --panel data/challenge/pert_counts.csv --axis data/challenge/gene_names.csv \
  --beta 0.35 --context X --out experiments/376-calibration-playbook/work/ctx_X_b035.npz \
  --heldout-seed <any int, e.g. 2210910>
```
Measured: Jurkat full build (fit + poison + held-out) 19.6 s, peak 8.4 GiB;
slowest bank (k562_essential) ~12 s + held-out, peak 11.2 GiB across all seven
(Table §"Results"). Serial only while a production build runs.

## 3. Validate (the PI's checklist — all must pass)
1. **Poison PASS** in stdout (`poison PASS byte-identical`) — recipient panel rows
   never read. Abort otherwise.
2. **Panel-gene exclusion**: the tool asserts zero panel symbols on the covered
   axis (it refuses to write otherwise — the 373 builder refuses panel support
   too). Confirm in the sidecar: `panel_genes_dropped` lists what was removed.
3. **Key set**: `targets genes matrix beta mode context norm_rule`, uncompressed
   `np.savez` (~10 MB for Jurkat-scale). `matrix` must be the UNBLENDED ridge
   prediction; `beta` the blending weight the builder applies; `mode=blend_rows`;
   `context=X`. Loader rejects anything else (373 `patch373.py` E5b/E5c).
4. **Held-out row** (sidecar `heldout_within_bank`, printed too): same−diag gain
   positive on BOTH splits with a ≤0 shuffled control. If shuffled is positive
   (RPE1 precedent), discount the gain by the shuffled floor before quoting.
5. **Config-key gate** at build time (375 `run319pA.sh` pattern): the finished
   package's meta must differ from its parent in EXACTLY the prior keys
   (`ctx_a_prior_mode/sha256/context/beta/targets` for A; new-context keys per
   the PI's builder wiring) — diff the artifact meta, not the command line.
6. **β choice per line**: 0.35 is retrieval-safe for Jurkat (+0.0051, resolved)
   and flat for HCT116/K562, but costs retrieval for HepG2 (−0.0100), HEK293T
   (−0.0118) and RPE1 (−0.0309) on 341's frame (376 value table). For those
   lines a lower β (or no ship) is the default — re-read the row before wiring.

## 4. Hand to the PI
Report: artifact path + sha256, sidecar path, the four numbers (training /
covered / query / held-out gain + shuffled), retrieval column, wall time and
peak RSS. The PI wires the builder flag (375's `VCC_CTX_A_PRIOR` pattern) and
owns the builder change, the five-key diff, and the slot decision. Do not run
`def_pipeline.py`, `make_submission*.py`, `vcc prep` (dry-run or real), or
`vcc submit` from this playbook.
