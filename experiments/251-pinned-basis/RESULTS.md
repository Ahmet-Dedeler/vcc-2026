# 251 — The PINNED deflation basis is in the builder: `VCC_CD4_DEFLATE_BASIS=V5_abc` applies A/B/C's gene-space `V_5` (250's October reading) on any panel, never a day-of re-fit. Flag unset is byte-inert; the four known-answer rows below. **PINNED PATH READY — verdict line at the bottom (§10).**

**Date 2026-09-04.** Deterministic. **No `vcc submit`, no `cell-eval2`, no model name, no
description, no slot proposal.** Files changed outside this directory: `scripts/make_submission_pm_cd4.py`
(+293 / −0) and `scripts/def_postbuild.py` (+66 / −0, the mandate's item 3, §1.3), both committed as
**`c321f76a0add0b74f71732c1de04204a86fa85f3`** with the known-answer hashes in the message
(`work/commit-message.txt`); plus the pinned artifact under `data/donors/cd4_deflation_basis/`
(force-added past `.gitignore`'s `*.npz`, 1.74 MB, so the pin is reproducible from the repository).
`scripts/def_pipeline.py` untouched (sha `ba0d6533…` before and after). Every ≥ 15 GiB stage ran
alone under the mandate's gate (`memory_pressure` ≥ 40 % free AND `vm_stat` free+inactive ≥ 30 GiB AND
zero other builders AND zero `vcc prep`, the working-tree builder hashed to this experiment's bytes
before each staging row); peak RSS per stage in §8.

| row | required | measured | verdict |
|---|---|---|---|
| **artifact** | `V5_abc.npz`: gene-space `V_5` from 249's staged A/B/C residual, a path and a sha | `data/donors/cd4_deflation_basis/V5_abc.npz` **sha256 `ba7ade3fc4f59de5c2de5d8bf270f243ab8e8ddfd0728be9e2c8c71c9b18122b`**, 1,743,704 B; `V` (18,533 × 5) array sha `dabf88f7…` **== 250's `basis250-production.npz[V5_abc]` bit for bit**; `U_5` sha **`f06e0262…` == the shipped `cd4svd5` audit**; source `experiments/249-cd4-deflate/work/stage/cd4-off.npz` sha `3a99dd3b…` | **PASS** |
| **KA 1** flags unset, staging | `effects_array_sha256` `abeb0749…` | **`abeb0749f8ecc9410fcec0d9d0cebd2569c97ca456d73f9e6b3bfd359164c315`** through the changed builder (19.95 GiB); **and** through the pre-251 byte copy with `VCC_CD4_DEFLATE_BASIS=V5_abc` set (unknown to it) | **PASS** |
| **KA 1** flags unset, package | `.h5ad` `8df4dadd…` on 12,648,807,800 B | **`8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6` on 12,648,807,800 B — MATCH** through the changed builder, both flags unset (`run251ka.sh`, 1329.8 s / 17.30 GiB; `def_postbuild` 32/32, `moved=[]`; no `cd4_deflate_*` key, no `deflation` block) | **PASS** |
| **re-fit untouched** | `VCC_CD4_DEFLATE_K=5` still reproduces 249's staged arm | **`a28d5be2e9cb3ed17d11a4c5f40f470893f0c968a1ae87b0b3b358321ead5d11`** bit for bit (20.11 GiB) | **PASS** |
| **KA 2** on A/B/C the pinned basis IS the re-fit basis | `VCC_CD4_DEFLATE_BASIS=V5_abc` reproduces `VCC_CD4_DEFLATE_K=5`'s staged effects bit for bit + the same energy numbers | staged **`aad5936597ce636b6832179dcec85cbd9dceafefa08d904ac7c88c30a9a34f25`** vs `a28d5be2…`: **max \|Δ\| 2.36e−15** on entries up to 1.70 (relative 1.4e−15), `cd4.fraction` **0.34104860341941284 identical**, `energy_removed_frac` **0.19866046506100588 identical**, captured 0.19823595471029617 vs top-k 0.1982359547102955 (6.7e−16), union fraction identical; **NUMERIC identity, NOT byte identity** — one projector in two floating-point associations (`U(UᵀR_c)` vs `(R_cV)Vᵀ`), the same shape as 247 §3.4's carrying-row row (4.4e−16 there) | **PASS at numeric identity; byte identity NO (§3)** |
| **KA 3** 233's `real` and 238's `efr` through the builder's staging path | 250's pinned rows to the digit: `pds` +0.0153 / +0.0107 at −0.1 … −2.0 % cosine; diagnostic 0.10–0.14 | **`real` Δ`pds` +0.015323919, agg +0.06 %, macro +0.48 %; `efr` +0.010649077, −1.99 %, −1.36 %** — every per-fold `cos_agg` / macro / `pds` / fraction at **0.0e+00** vs 250's json (24/24 rows), effect matrices == 250's `compose` bit for bit; diagnostic **0.11713736981577881 / 0.12493656042166369** (250's to the last digit) | **PASS** |
| **KA 4** the rehearsal on 238's `ef-bundle`, `VCC_CD4_DEFLATE_BASIS=V5_abc` exported and nothing else | meta carries the pinned block, `def_postbuild` passes, the diagnostic prints, `vcc prep --dry-run --contexts D,E,F` clean | `scripts/def_pipeline.py` preflight → … → postbuild → prep, every stage rc 0; 244's operator accepted by the skip branch on its axis; **meta: `config.cd4_deflate_basis = V5_abc`, `effects.cd4.deflation.mode = pinned`, 186 carried rows, energy removed 0.13556624003063233, captured 0.1333778923231684 (= 250's `efp` reading to the last digit)**; the builder printed the diagnostic in `repin-panel-support` and `build`; `[arms]` 6 OK; **`def_postbuild` 40/40 through the driver** (the 8 new rows + the echo line); **`✓ dry run`**, 360,000 cells, D/E/F 120,000 each, targets verified; builder byte-restored to `10d0177e…`, `scripts/` clean | **PASS** |
| unit / env / suites | — | `unit251` **47/47**; `envcheck251` **10/10**; the builder's 57-test suite **56/56 under both builders, 0 regressions**; 235's hardening tests **9/9**; `def_postbuild` on the 222 / 249 / 248 metas **31/31 · 32/32 · 32/32 under both scripts** | **PASS** |

---

## 0. Why this exists, in four lines

250 measured that 246/249's deflation transfers to a disjoint 300-target panel in its `pds` mechanism but
that a **day-of re-fit** of the projector on the new panel's rows removes panel-specific signal too (cosine
cost 4–14× A/B/C's, price ≤ 0 on both low-coverage draws), while the **pinned A/B/C basis** beats the re-fit
on 4/4 disjoint panels (`pds` +0.0153 / +0.0037 / +0.0073 / +0.0107, price +0.0061 / +0.0009 / +0.0014 /
+0.0030, never negative). Pinning needed a builder path distinct from `VCC_CD4_DEFLATE_K` (rule 10), a
249-style known answer on A/B/C, and a rehearsal. This is that path. Whether it ships is the lead's call
after the 2026-09-05 `cd4svd5` board point.

---

## 1. The change

### 1.1 The artifact — `data/donors/cd4_deflation_basis/V5_abc.npz` (`basis251.py`, `work/basis251.log`)

From 249's staged A/B/C CD4 residual before deflation (`experiments/249-cd4-deflate/work/stage/cd4-off.npz`,
sha256 `3a99dd3b37bf356ab2b99cb15c3b809b551ae7d0c315c24727f5dd84a24b894c`; 300 rows, 291 carried, 276
complete; the object whose k = 5 deflation the shipped `cd4svd5` meta audits), 250's `lib250.gene_basis`
verbatim: rows = `target_condition_count > 0`, `R` masked to the gene support, `R_c` = per-gene centred over
the carried rows and re-zeroed outside the support, `G = R_c R_cᵀ`, `eigh`, `U_5` = top-5 eigenvectors,
`V_5 = R_cᵀ U_5 / s_5`. Stored: `V` (18,533 × 5, float64, orthonormal to 4.0e−15), `U` (291 × 5), `s`, the top-15
eigenvalues, `k`, the official gene axis (18,533 symbols, `gene_names.csv` sha `25bfa667…`), the 291 carried
targets, the A/B/C per-gene mean (provenance only — the pinned path centres on the DAY's rows), and a
provenance json (source path + three array shas, construction, application, date, threads, numpy). Sidecar
`V5_abc.json` repeats every sha.

| known answer (abort on failure) | measured |
|---|---|
| the builder's own `deflate_cd4_residual(res, 5)` on the source reproduces the shipped audit | n_carried 291, energy_removed_frac 0.19866046506100588, top_k 0.1982359547102955, `U_5` sha **`f06e0262cec2affb6c0594d2cebf8285888fd36ed73b518d514c7f0129b74e4c`** — 4/4 identical |
| this script's `U_5` hashes to `f06e0262…`; `V_5` == 250's `basis250-production.npz[V5_abc]` | **bit for bit** (under the builder's 4-thread BLAS default; 250 fitted it at 6 threads — same bytes) |
| gene axis == `data/challenge/gene_names.csv` (18,533, sha `25bfa667…`); source target order == `pert_counts.csv` | PASS |
| pinned vs re-fit arithmetic on the source's own rows | max \|Δ\| **7.8e−15** on entries up to 2.82, **not bitwise** — the pre-declared reason KA 2 is numeric |

Registry entry in the builder (rule 11, every value a plain string so `def_postbuild.py` reads it as a literal):

```
"V5_abc": {"path": "data/donors/cd4_deflation_basis/V5_abc.npz",
           "sha256": "ba7ade3fc4f59de5c2de5d8bf270f243ab8e8ddfd0728be9e2c8c71c9b18122b",
           "v_sha256": "dabf88f7f95ffc4dac109de640ad29b943773ec4e45c5ddb2a4231b0003cce61",
           "u_k_sha256": "f06e0262cec2affb6c0594d2cebf8285888fd36ed73b518d514c7f0129b74e4c", "k": "5"}
```

### 1.2 `scripts/make_submission_pm_cd4.py` (+293 / −0 on 248's `2632b92` bytes `f7a12a70…` → `10d0177e…`; `patch251.py`, five exact-anchor edits, `--check` reproduces the working tree from the byte copy)

- **`CD4_DEFLATION_BASIS_REGISTRY`** + **`VCC_CD4_DEFLATE_BASIS`** (import time, after `VCC_CD4_DEFLATE_K`):
  unset → `None`, inert. Refused at import: both flags set (*"two variables (rule 10); set exactly one"*), an
  unregistered name (case-sensitive), `K` invalid alongside. `envcheck251` 10/10 (§1.4).
- **`load_cd4_deflation_basis(name, *, official_genes=None)`** — fail-closed on every pin, in order: unregistered
  name; missing file; file sha256 ≠ the pin; `k` ≠ the pin; shape ≠ (n_genes, k) / non-finite / not orthonormal
  (max |VᵀV − I| > 1e−9); `V` array sha256 ≠ the pin; provenance `U_k` sha ≠ the pin; and when `official_genes`
  is given (`run()` passes the bundle's axis BEFORE the three-donor model is loaded), a gene axis that is not
  this bundle's, symbol for symbol. Cached per name; **every call that names an axis verifies it** (the first
  draft returned a cache verified against one bundle for another — caught by `unit251`'s "wrong gene ORDER"
  row and fixed before commit).
- **`deflate_cd4_residual_pinned(cd4, basis, *, basis_info)`** — 250's `lib250.deflate_gene_basis` verbatim
  (`unit251`: bit-identical): on the carried rows, per-gene centring on THIS panel's rows, re-zero outside the
  support, `P = R_c V`, `R_k = R_c − P Vᵀ`, re-zero outside the support. No fit, no eigh, no seed. Returns a
  `CD4Prepared` with the SAME mode / votes / counts / details (so `effects.cd4.mode`, `supported_entries`,
  `fully_supported_targets` cannot move) plus the audit: `mode: pinned`, `basis`, `basis_path`, `basis_sha256`,
  `basis_v_sha256`, `basis_u_k_sha256`, `k`, `n_genes`, `orthonormality_err`, `gene_axis_verified`,
  `n_carried`, `energy_removed_frac`, **`captured_frac_of_centred`** (250's census diagnostic),
  `own_top_k_energy_frac_of_centred` (the day's own top-k share, recorded because 250 §6.2 measured it NOT
  to discriminate), `centred_energy_frac_of_raw`, `provenance`.
- **`apply_frozen_cd4_residual`** — one `elif CD4_DEFLATE_BASIS is not None:` branch after 249's `if
  CD4_DEFLATE_K is not None:` (mutually exclusive at import, so exactly one can run; unset, neither does):
  load → deflate → assert the panel coordinates are still zero → print the diagnostic → then
  `solve_cd4_dose` (247's router) sizes the deflated object, 249's *resolved* construction. The union
  residuals are never touched.
- **`run()`** — with the flag on, the basis is hashed and its axis verified against `spec.genes` right after
  the frozen-input hashing, before the 17 GiB model load (`[cd4-basis] …` line).
- **Meta** — `effects.cd4.deflation` (the audit above) and `config.cd4_deflate_basis`, written only when on
  (the `h1_union_g` convention): a default build's config block and audit are byte-identical to every prior
  build's. 249's re-fit audit keys are untouched (§3: `k5` reproduces `a28d5be2…`).

### 1.3 `scripts/def_postbuild.py` (+66 / −0; the mandate's item 3 — exactly what was touched)

Without this, 247's rule-10 row *"config differences vs the incumbent are only […]"* refuses a pinned build:
`moved` gains `cd4_deflate_basis`, and the driver's `postbuild` stage passes only `--allow-config-diff
whitener_id`. The shipped-arms rows (`alpha`, `cd4_target_f`, `h1_union_g`, `h1_union_sources`, `emitter`)
were never affected. Four exact-anchor edits (`patch251.py`, `--check` reproduces):

1. `builder_literal_or_empty()` — `builder_literal` returning `{}` when the source predates the table.
2. A `251:` block before the 240 arms rows, run **only when `config.cd4_deflate_basis` is present**: eight rows —
   the name is in the builder's `CD4_DEFLATION_BASIS_REGISTRY` (parsed as a literal, never imported);
   `config.cd4_deflate_k` ABSENT; `effects.cd4.deflation.mode == "pinned"` and `.basis` == the config name;
   `.basis_sha256` == the pin; `.basis_v_sha256` and `.k` == the pin's; **sha256 of the file at the pin's path
   == the pin** (read the pin, hash the file, as the operator rows do); `gene_axis_verified` is `true`;
   `captured_frac_of_centred` is a number in (0, 1). Plus one row that fires when the audit says `pinned`
   but the config carries no name. The key is added to the allowed diff set only inside this block.
3. `allowed |= extra_allowed` in the rule-10 section.
4. An echo line in the summary: `cd4  deflation PINNED 'V5_abc' sha … k=5: N carried rows, energy removed …,
   captured … of the day's centred CD4 energy (250: 0.10–0.14 disjoint, 0.198 A/B/C), own top-k share …
   (recorded, not gated)`.

A build without the key sees none of these rows: 222's meta **31/31**, 249's arm meta **32/32**
(`--allow-config-diff cd4_deflate_k`), 248's KA meta **32/32** — identical counts under the pre-251 byte copy
(`postbuild251.py`, `work/postbuild251-pre.log`). The firing rows are demonstrated on the rehearsal meta in §5.

### 1.4 Unit level, import rules, suites (`work/unit251.log`, `envcheck251.log`, `suite251.log`, `tests235.log`)

`unit251.py` (two builders in one process, the working tree vs the pre-251 byte copy): **47/47**.
OFF: output array, CD4 audit and union audit bit-identical to pre-251 with and without an H1-style extra,
no `deflation` key. ON (synthetic residual with complete / partial / untouched rows): pinned(own gene basis)
== the K re-fit to **2.7e−15** (audit energies to 1e−12), == `lib250.deflate_gene_basis` **bit for bit**;
mode / votes / counts untouched; zero outside the support and on the partial rows; a foreign basis on
full-support rows obeys Pythagoras to 1e−12 (a projector, not a norm cut) with `R_k V` = 0 to 1e−10; a wrong
n_genes / k = 0 / 1-D basis and an empty residual are refused. Through `apply_frozen_cd4_residual`: the
fraction re-solves, rows CD4 does not carry equal OFF's × one scalar (4.4e−16), carried rows move in
direction, the union's frame is untouched, a basis with mass on a panel gene leaves the panel coordinates at
whitened × scale, pinned(own basis) == the K re-fit through the whole entry point to 1.1e−15, and the re-fit
audit is 249's. LOADER: `V5_abc` loads and == 250's `V` bit for bit; a wrong gene order, an unregistered
name, a wrong file sha, a tampered `V` (re-pinned file sha, not orthonormal), a wrong `k`, a wrong `V` sha,
a wrong provenance `U_k` sha and a missing file are each refused. REAL (249's `cd4-off.npz` on a one-hot
stub base whose median row norm is exactly the shipped 4.229295604454299): the re-fit path re-solves the
shipped fraction **0.34104860341941284 bit for bit** and the pinned path gives the same fraction (|Δ| 0.0)
and the same output to **2.5e−15**.

`envcheck251.py` (ten subprocess imports): unset / `V5_abc` / `K=5` / blank → ok; both set → refused ("two
variables"); `nope` and `V5_ABC` → refused (unregistered); `K=0`+basis → refused on `K`; the pre-251 byte
copy imports with the flag set and has no attribute. **10/10.**

Regression suites: the builder's 57-test suite **56/56 under both builders, 0 regressions, 0 fixed** (the one
failure predates 222, 247 §5.1); 235's hardening tests **9/9** (1.06 GiB).

---

## 2. Known-answer row 1 — the flag unset is byte-inert

### 2.1 Staging (`stage251.py` = 249's `stage249.py` + the `run()`-style basis verification; `stage251.sh`, each row alone and gated; `work/stage/stage251.log`, verbatim)

```
2026-09-04T15:52:55Z gate: memory_pressure free 82%  free+inactive 31.4 GiB  other builders 0  vcc preps 0  builder intended
2026-09-04T15:53:12Z END   pre  rc=0  peak RSS 13.70 GiB  array_sha256 abeb0749f8ecc9410fcec0d9d0cebd2569c97ca456d73f9e6b3bfd359164c315  cd4.fraction 0.32207201067578145
2026-09-04T15:53:30Z END   off  rc=0  peak RSS 19.95 GiB  array_sha256 abeb0749f8ecc9410fcec0d9d0cebd2569c97ca456d73f9e6b3bfd359164c315  cd4.fraction 0.32207201067578145
2026-09-04T15:53:48Z END   k5  rc=0  peak RSS 20.11 GiB  array_sha256 a28d5be2e9cb3ed17d11a4c5f40f470893f0c968a1ae87b0b3b358321ead5d11  cd4.fraction 0.34104860341941284
2026-09-04T15:54:05Z END   pinned  rc=0  peak RSS 19.74 GiB  array_sha256 aad5936597ce636b6832179dcec85cbd9dceafefa08d904ac7c88c30a9a34f25  cd4.fraction 0.34104860341941284
2026-09-04T15:54:10Z KA1(staging) pre MATCH abeb0749f8ecc9410fcec0d9d0cebd2569c97ca456d73f9e6b3bfd359164c315
2026-09-04T15:54:10Z KA1(staging) off MATCH abeb0749f8ecc9410fcec0d9d0cebd2569c97ca456d73f9e6b3bfd359164c315
2026-09-04T15:54:10Z KA2(re-fit untouched) k5 MATCH a28d5be2e9cb3ed17d11a4c5f40f470893f0c968a1ae87b0b3b358321ead5d11
```

| row | builder | env | `effects_array_sha256` | `effects.cd4.fraction` | peak RSS |
|---|---|---|---|---|---:|
| `pre` | pre-251 byte copy (HEAD `89cb22a`, `f7a12a70…`) | shipped arms + `VCC_CD4_DEFLATE_BASIS=V5_abc` | **`abeb0749…` MATCH** (the flag does not exist there) | 0.32207201067578145 | 13.70 GiB |
| `off` | working tree `10d0177e…` | shipped arms, flags unset | **`abeb0749…` MATCH** | 0.32207201067578145 | 19.95 GiB |
| `k5` | working tree | + `VCC_CD4_DEFLATE_K=5` | **`a28d5be2…` = 249's staged arm bit for bit** (the re-fit path is 249's) | 0.34104860341941284 | 20.11 GiB |
| `pinned` | working tree | + `VCC_CD4_DEFLATE_BASIS=V5_abc` | `aad59365…` — §3 | 0.34104860341941284 | 19.74 GiB |

(The `pre` row's 13.7 GiB is the same computation with the donor pages already cached; 249 measured 19.9.)

### 2.2 The full A/B/C package (`run251ka.sh` → `build251.sh` = 222's `build222.sh` + the two passthroughs, both unset, gated, alone; `work/build/run251ka.log`, verbatim)

```
2026-09-04T15:58:42Z gate: memory_pressure free 77%  free+inactive 30.3 GiB  other builders 0  vcc preps 0  builder intended
2026-09-04T15:58:42Z KA BUILD start (both flags unset)
2026-09-04T16:20:52Z KA BUILD rc=0 wall 1329.79 s peak RSS 17.30 GiB
2026-09-04T16:21:15Z KA h5ad sha256 8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6 bytes 12648807800  (want 8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6 12648807800)  MATCH
2026-09-04T16:21:16Z KA postbuild rc=0 32/32 checks pass 0 FAIL lines
2026-09-04T16:21:16Z DELETED .../work/build/251-ka-h1union-g028.h5ad (byte duplicate of submissions/222-h1union-g028.h5ad; meta kept in .../work/build)
```

| what | measured | required |
|---|---|---|
| `work/build/251-ka-h1union-g028.h5ad` sha256 / bytes | **`8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6` / 12,648,807,800 — MATCH** | the shipped `h1union` (240 row 1, 247 §3.4, 249 §2.2, 248 §2.2) |
| `output.nnz` / `effects.final_sha256` / `effects.cd4.fraction` / union fraction | 1,577,411,397 / `abeb0749…` / 0.32207201067578145 / 0.19207161964444322 | the shipped meta's |
| `config.cd4_deflate_basis` / `config.cd4_deflate_k` / `effects.cd4.deflation` | **absent / absent / absent** | absent (the `h1_union_g` convention) |
| `def_postbuild.py` against the 222 incumbent meta | **32/32 checks pass, `moved=[]`** | — |
| `hashes.implementation` for the builder | `10d0177e…` (commit `c321f76`'s bytes) | changed, as 235's / 247's / 249's / 248's did; the h5ad did not |
| cost | **1329.8 s**, peak RSS **17.30 GiB** | 249: 1345 s / 20.00 GiB; 248: 1478 s / 16.72 GiB |

The h5ad was re-hashed and deleted as a verified byte duplicate of `submissions/222-h1union-g028.h5ad`; the
meta is kept in `work/build/`. The `.vcc` was not re-prepped: it is a deterministic transform of the h5ad
(240 row 1, 235 row 2, 249 §2.2 all took `8df4dadd…` to `f13ce40b…`), so the byte identity rests on the h5ad
row above, as 247's and 248's did.

---

## 3. Known-answer row 2 — on A/B/C the pinned basis IS the re-fit basis (`compare251.py`, `work/compare251.json`)

```
  pinned  array sha256 aad5936597ce636b6832179dcec85cbd9dceafefa08d904ac7c88c30a9a34f25
  re-fit  array sha256 a28d5be2e9cb3ed17d11a4c5f40f470893f0c968a1ae87b0b3b358321ead5d11  (249 stage-svd5.json a28d5be2e9cb3ed1...)
  k5 (this builder)    a28d5be2e9cb3ed17d11a4c5f40f470893f0c968a1ae87b0b3b358321ead5d11  == 249's bitwise: True
  bitwise equal: False
  max |d| 2.359e-15 on entries up to 1.705 (relative 1.384e-15); 5,244,114 of 5,559,900 entries differ, on 300 rows (291 CD4-carried, 9 not carried; max |d| on non-carried rows 4.441e-16)
  pinned vs OFF (abeb0749) max |d| 5.341e-01 (the flag ACTS)
  cd4.fraction pinned 0.34104860341941284 vs re-fit 0.34104860341941284 (|d| 0.0e+00)
  preserve 0.9115882187626341 vs 0.9115882187626339; union fraction 0.19207161964444322 vs 0.19207161964444322
  deflation pinned: n_carried 291, energy_removed 0.19866046506100588, captured 0.19823595471029617, own top-k 0.19823595471029556, basis V5_abc file ba7ade3fc4f59de5 V dabf88f7f95ffc4d U_k f06e0262cec2affb axis verified True
  deflation re-fit: n_carried 291, energy_removed 0.19866046506100588, top-k 0.1982359547102955, U_k f06e0262cec2affb
  |d energy_removed| 0.0e+00; |captured - top-k| 6.7e-16
ROW 2 byte identity: NO (expected: two associations of one projector)   ROW 2 numeric identity (<= 1e-12 on the matrix, 1e-14 on the fraction): PASS
```

**What "bit for bit" turned out to mean, measured.** The two paths apply the same projector to the same
291 rows: the re-fit computes `R_c − U_5 (U_5ᵀ R_c)` (row space, `U_5` from `eigh`), the pinned path computes
`R_c − (R_c V_5) V_5ᵀ` (gene space, `V_5 = R_cᵀ U_5 / s_5`). IEEE arithmetic associates the two products
differently, so the deflated residuals agree to **7.8e−15** (`basis251.py`, on the source rows) and the staged
effect matrices to **2.4e−15** (relative 1.4e−15), the fraction, the energy numbers and the union's fraction
are identical to the last digit, the preserve scalar moves by 2e−16, and the non-carried rows move by that
scalar (4.4e−16). No construction that is a genuine gene-space projector on a foreign panel can reproduce
the row-space bytes on A/B/C — the only way to get `a28d5be2…` from the pinned path would be a branch that
detects A/B/C's own rows and runs 249's arithmetic, which would make this row test nothing. The pre-declared
reading in `compare251.py` (written before the row ran) is numeric identity at ≤ 1e−12, and it is met with
three orders of magnitude to spare; the same reading 247 §3.4 used for the carrying-row rule on A/B/C
(4.4e−16, "numeric identity at one ULP"). **The unit-level preview (§1.4 REAL) gave 2.5e−15 on a stub base
with the same conclusion before any 20 GiB stage ran.**

Consequence for the calibration argument (250 §4.5 iii): the `cd4svd5` board point on 2026-09-05 calibrates
the pinned object to floating-point order — the h5ad a pinned A/B/C build would emit is not byte-identical to
`249-cd4svd5-f034.h5ad`, but its effect matrix differs by 1e−15 relative, far below anything the emitter's
integer draw or the board can resolve.

---

## 4. Known-answer row 3 — 250's pinned rows through the builder on the disjoint panels (`repro251.py`, one (panel, fold) per process; `work/repro251-*.{json,log}`, `repro251-summary.log`)

162's PLOLO as 250 ran it (`lib250` byte copy — `bind250` / `lib250` / `gate250` sha-identical to 250's, caches
under this experiment's `work/plolo/`; 233's per-panel re-fit operator, axis-asserted; 233's / 238's re-extracted
CD4 banks; 247's carrying-row sizing at `g_cd4` = 0.3443655278649641), with the WHOLE of the builder's
`apply_frozen_cd4_residual` (pinned deflate → carrying-row size → merge with the H1 union → one preserve)
substituted for 250's `deflate_gene_basis` + `solve_cd4` + `compose`. The builder imported with October's env
(`def_pipeline.DEF_ENV` on top of the shipped arms) plus `VCC_CD4_DEFLATE_BASIS=V5_abc`.

| panel / fold | KA-A `BASE` vs 250 | KA-B builder pinned deflation vs `lib250.deflate_gene_basis` | KA-D diagnostic vs `basis250.json` | KA-C arm vs 250 `pinABC5` | fraction base → arm (= 250's) | arm: macro / agg / Δ`pds` |
|---|---|---|---|---|---|---|
| `real` / k562_gwps | 0.0e+00 | bit-identical; audit identical | 0.11713736981577881 (|Δ| 0.0) | **0.0e+00** | 0.382842 → 0.408034 | −2.79 % / −2.32 % / +0.009275 |
| `real` / orion_hct116 | 0.0e+00 | bit-identical | 0.0 | **0.0e+00** | 0.363300 → 0.387206 | +1.55 % / +0.41 % / +0.017514 |
| `real` / orion_hek293t | 0.0e+00 | bit-identical | 0.0 | **0.0e+00** | 0.363814 → 0.387754 | +2.67 % / +2.08 % / +0.019182 |
| `efr` / k562_gwps | 0.0e+00 | bit-identical | 0.12493656042166369 (|Δ| 0.0) | **0.0e+00** | 0.349800 → 0.367736 | +3.04 % / +0.58 % / +0.013103 |
| `efr` / orion_hct116 | 0.0e+00 | bit-identical | 0.0 | **0.0e+00** | 0.352992 → 0.371091 | −4.04 % / −3.60 % / +0.012014 |
| `efr` / orion_hek293t | 0.0e+00 | bit-identical | 0.0 | **0.0e+00** | 0.352441 → 0.370511 | −3.08 % / −2.96 % / +0.006830 |

Three cancer folds, against 250's own json rows recomputed the same way (`work/repro251-summary.log`):

| | `real`: this run | 250 | \|Δ\| | 250 §4.3 | `efr`: this run | 250 | \|Δ\| | 250 §4.3 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| rel MACRO cosine | **+0.004781411** | +0.004781411 | 0.0e+00 | +0.48 % | **−0.013599386** | −0.013599386 | 0.0e+00 | −1.36 % |
| rel AGGREGATE cosine | **+0.000576161** | +0.000576161 | 0.0e+00 | +0.06 % | **−0.019933889** | −0.019933889 | 0.0e+00 | −1.99 % |
| Δ`pds` | **+0.015323919** | +0.015323919 | 0.0e+00 | +0.01532 | **+0.010649077** | +0.010649077 | 0.0e+00 | +0.01065 |

`pds` positive on 3/3 folds on both panels; Δq_w +0.00000 / −0.00033; the union's solved fraction untouched
on every cell (|Δ| 0.0); effect matrices vs 250's `compose` max |Δ| **0.0** on all six; 24/24 known-answer
rows at 0.0e+00; peak RSS 4.94–5.10 GiB, 7–13 s per cell. **The builder's October path IS 250's pinned
construction to the bit**, and the mandate's digits (+0.0153 / +0.0107 at −0.1 … −2.0 % — 250 quotes the
aggregate as +0.06 % on `real`, i.e. inside "−0.1 …", and −1.99 % on `efr`) reproduce.

---

## 5. Known-answer row 4 — the rehearsal (`rehearsal251.sh` → `work/rehearsal/`; `rehearsal251.log`, `driver.log`, `driver-prep.log`, `logs/<stage>.log`, `timings.json`)

The driver on 238's `ef-bundle` (D/E/F, 300 targets, 0 in A/B/C, 0 essential, K562 146/300, CD4 union
186/300 — 244's panel) with **`VCC_CD4_DEFLATE_BASIS=V5_abc` exported and nothing else** (the driver's
defaults supplied the arms: `VCC_ALPHA 0.982 · VCC_CD4_TARGET_F 0.34 · VCC_H1_UNION_G 0.28 · VCC_CD4_SIZING
carrying · VCC_CD4_TARGET_G 0.3443655278649641`, `VCC_UNION_SOURCES` unset, emitter `rc_callset`, the shipped
flags), 244's fitted operator handed to `--whitener-out`, a throwaway frozen table, `--cd4-thresholds relaxed`
(the same bundle 244 §2.2 confirmed the mapping for), `--restore-builder-after`, gate 30 GiB; committed
builder `c321f76` (`10d0177e…`), `git status --short scripts/` empty before the pin (244 §8 edit 3).

```
2026-09-04T16:21:16Z ARMED pid=26080  builder before: 10d0177e...  git: c321f76
2026-09-04T16:21:16Z gate: memory_pressure free 79%  free+inactive 31.5 GiB  other builders 0  vcc preps 0
=== [whitener] exists and its axis matches this bundle (metadata:whitener_union_pruned_def.json); skipping the fit: .../244-october-rehearsal/work/run1/whitener_union_pruned_def.npz
[cd4-basis] V5_abc: data/donors/cd4_deflation_basis/V5_abc.npz sha256 ba7ade3fc4f59de5c2de5d8bf270f243ab8e8ddfd0728be9e2c8c71c9b18122b (V dabf88f7f95ffc4d, U_k f06e0262cec2affb), k=5, gene axis verified against the bundle (18533 genes)
[cd4] deflation pinned V5_abc (file sha256 ba7ade3fc4f59de5, k=5): 186 carried rows, energy removed 0.135566, captured 0.133378 of the day's centred CD4 energy (250: 0.10-0.14 on disjoint panels, 0.198 on A/B/C; < 0.05 or > 0.19 is unmeasured), own top-k share 0.231608
[cd4] sizing=carrying g=0.344366: solved 0.3600962339 on 182 carrying rows (60.7%); the all-rows rule at f=0.34 would solve 0.4149939363281411 (over-dose x1.16725)
--- [build] 1333.8s, peak RSS 20.05 GiB, rc=0
[arms] the built config against the configuration asked for:
  OK  config.alpha: asked 0.982, the artifact carries 0.982
  OK  config.cd4_target_f: asked 0.34, the artifact carries 0.34
  OK  config.h1_union_g: asked 0.28, the artifact carries 0.28
  OK  config.h1_union_sources: asked ['gse295214_h1', 'h1_2025'], the artifact carries ['gse295214_h1', 'h1_2025']
  OK  config.cd4_sizing: asked 'carrying', the artifact carries 'carrying'
  OK  config.cd4_target_g: asked 0.3443655278649641, the artifact carries 0.3443655278649641
cd4       deflation PINNED 'V5_abc' sha256 ba7ade3fc4f59de5… k=5: 186 carried rows, energy removed 0.13556624003063233, captured 0.1333778923231684 of the day's centred CD4 energy (250: 0.10-0.14 disjoint, 0.198 A/B/C), own top-k share 0.23160821107385518 (recorded, not gated)
  PASS  config differences vs 222-h1union-g028.h5ad.meta.json are only ['cd4_deflate_basis', 'cd4_sizing', 'cd4_target_g', 'contexts', 'expect_contexts', 'targets', 'whitener_id']: moved=['cd4_deflate_basis', 'cd4_sizing', 'cd4_target_g', 'contexts', 'expect_contexts', 'whitener_id']
40/40 checks pass
--- [postbuild] 0.8s, peak RSS 0.13 GiB, rc=0
✓ dry run — no .vcc written   cells: 360000   contexts: column 'context' → D=120000, E=120000, F=120000   targets: verified against the official list
restored scripts/make_submission_pm_cd4.py from pin-artifacts.builder-before.py (sha256 10d0177e33e53e9f…): feaa1509f5c1e163… -> 10d0177e33e53e9f… (identical)
2026-09-04T16:46:55Z builder after restore: 10d0177e...  (intended 10d0177e...)  git:
```

| # | stage | wall s | peak RSS GiB | rc | observed |
|---|---|---:|---:|---:|---|
| 1 | preflight | 2.9 | 0.32 | 0 | 66/66 |
| 2 | identity | 9.7 | 2.74 | 0 | D JURKAT, E HELA, F CAL33 |
| 3 | cd4-extract (relaxed) | 1.6 | 0.32 | 0 | 183 / 185 / 184, union 186 |
| 4 | census | 0.3 | 0.10 | 0 | gate pass |
| 5 | whitener | skipped | — | — | 244's `373df206…` accepted on axis `f060e926…` (metadata json), 247 B#18 |
| 6 | repin-frozen | 1.3 | 0.23 | 0 | throwaway table |
| 7 | pin-artifacts | 0.6 | 0.05 | 0 | `union_pruned_def` + 15 CD4 hashes; builder `10d0177e…` → `feaa1509…` |
| 8 | repin-panel-support | 12.5 | 19.96 | 0 | **the pinned diagnostic printed**; `panel_support` {133, 182} |
| 9 | build | **1333.8** | **20.05** | 0 | `[cd4-basis]` verified axis; `[cd4] deflation pinned …`; `[arms]` 6 OK |
| 10 | postbuild | 0.8 | 0.13 | 0 | **40/40** (32 of 244's + the 8 pinned-basis rows), `moved` = the five 244 saw + `cd4_deflate_basis` |
| 11 | prep (`--dry-run --contexts D,E,F`) | 17.5 | 16.70 | 0 | `✓ dry run`, 360,000 cells, D/E/F 120,000 each, targets verified |
| — | restore-builder | 0.1 | 0.03 | 0 | byte-identical to the pre-pin copy = commit `c321f76` |

**The meta** (`work/rehearsal/lantern-under-the-stairs.h5ad.meta.json`; the h5ad `76a389c7b0f476f7bf8829445718c0bd96df0a1ba475589b54e9b78817621084`,
12,623,641,976 B, nnz 1,574,625,808, 360,000 rows, kept, gitignored): `config.cd4_deflate_basis = "V5_abc"`,
`config.cd4_deflate_k` absent, `cd4_sizing carrying` at `g` 0.3443655278649641, `whitener_id union_pruned_def`,
operator `373df206…` axis `f060e926…`; `effects.cd4.deflation` = `{mode pinned, basis V5_abc, basis_path
data/donors/cd4_deflation_basis/V5_abc.npz, basis_sha256 ba7ade3f…, basis_v_sha256 dabf88f7…, basis_u_k_sha256
f06e0262…, k 5, n_genes 18533, orthonormality_err 4.0e−15, gene_axis_verified true, n_carried 186,
energy_removed_frac 0.13556624003063233, captured_frac_of_centred 0.1333778923231684,
own_top_k_energy_frac_of_centred 0.23160821107385518, centred_energy_frac_of_raw 0.99976, provenance {…}}`;
`effects.cd4.fraction` **0.3600962338656195** on 182 carrying rows (244's un-deflated 0.3417130399758081, ×1.054 —
the resolved sizing on the deflated object, as on A/B/C's ×1.059), over-dose factor of the all-rows rule 1.16725
(244's un-deflated object 1.17016); the H1 union untouched (`carrying_rows` 28, `solved_fraction`
0.20198547321112115 = 244's to the digit); merged preserve 0.9531004116394588 (244: 0.9553649891807288).

**The diagnostic through the builder on the day's own re-extracted banks equals 250's `efp|by_abc|5`
(0.1333778923231684) to the last digit and the day's own top-5 share equals 250's 0.23160821107385526 to
8e−17** — the builder's CD4 object on this bundle IS 250's production-shaped object (250 KA-P1's finding on
A/B/C, now on a disjoint panel too).

`postbuild251.py --rehearsal` (`work/postbuild251.log`, **16/16**): the rehearsal meta **40/40** under the changed
script; under the pre-251 byte copy **31/32, FAIL** on rule 10 (`moved` carries `cd4_deflate_basis`) — the row the
mandate's item 3 is about; against a pre-251 builder SOURCE (`--builder`, no registry) **33/40, FAIL** on the
registry row; tampered copies each fail the named row — basis sha moved, `mode` not pinned, `cd4_deflate_k`
added, name unregistered, audit pinned without the config key, `gene_axis_verified` false, `V` sha moved.

Process note, on the record: after the build, the prep gate read free+inactive 24.6–25.6 GiB for two
minutes (the fresh 12.6 GB h5ad's clean page cache, 247 §3.4's symptom); rather than wait out the script's
2-hour poll I dropped the clean cache with `sudo purge` (24.7 → 34.9 GiB, `memory_pressure` 79 % free) and the
gate opened on the next poll. No process was touched; the two armed submit scripts (pids 26557 / 80453) slept
through everything.

---

## 6. The captured-energy diagnostic on every panel touched (`effects.cd4.deflation.captured_frac_of_centred`)

| panel | CD4 carried rows | captured by the pinned `V_5` | own top-5 share (recorded, not gated) | energy removed by the pinned deflation | where |
|---|---:|---:|---:|---:|---|
| A/B/C (production, 291 rows) | 291 / 300 | **0.19823595471029617** | 0.19823595471029556 | 0.19866046506100588 | §3, the staged `pinned` row; == 250's `abc|by_abc|5` and the shipped top-k 0.1982359547102955 to 6.7e−16 |
| `real` (199's realistic draw, 300 rows) | 300 / 300 | **0.11713736981577881** | 0.20248526907148018 | 0.12004522940489226 | §4, all three folds; == 250's `basis250.json` to the last digit |
| `efr` (238's essential-free realistic draw) | 300 / 300 | **0.12493656042166369** | 0.1984391612775082 | 0.1284920356342738 | §4, all three folds; == 250's |
| `ef-bundle` = 238's `efp` (essential-free pessimistic draw, the rehearsal) | 186 / 300 | **0.1333778923231684** | 0.23160821107385518 | 0.13556624003063233 | §5, the rehearsal meta through the driver on the day's own re-extracted banks; == 250's `efp|by_abc|5` to the last digit |

250's reading (§6.2): 0.10–0.14 on the four disjoint panels, 0.198 in-sample, 0.0005 random; < 0.05 means the
hub is absent from the day's residual, > 0.19 means a panel that behaves like A/B/C's own rows; and the own
top-5 share (0.198–0.232 everywhere) does not discriminate and must not gate.

---

## 7. What the lead has to decide, and what the driver does NOT yet do

- **Nothing in `scripts/def_pipeline.py` changed** (the mandate's item 3). With `VCC_CD4_DEFLATE_BASIS=V5_abc`
  exported, the driver passes it through to every builder subprocess (`{**os.environ, **resolved}`) but does
  **not print it** in its "configuration this run will build" block and `assert_build_carries_arms` does
  **not assert it** — the same visibility gap `VCC_CD4_DEFLATE_K`, `VCC_UNION2_*` and `VCC_MASS_FIX` have
  today. The build log's `[cd4-basis] …` / `[cd4] deflation pinned …` lines and `def_postbuild`'s echo are the
  places it shows. If the lead ships it, the three-line change is: `DEF_ENV` gains
  `"VCC_CD4_DEFLATE_BASIS": "V5_abc"`, `ARM_BOARD_DELTA` gains its board delta once the 2026-09-05 point is
  read, and `assert_build_carries_arms` gains a `cd4_deflate_basis` row keyed on `resolved.get("VCC_CD4_DEFLATE_BASIS")`.
- `def_postbuild.py` accepts the key only with every pin verified (§1.3); it does not require it. Requiring it
  on D/E/F is a one-line addition once the runbook row says SHIPS.

---

## 7b. Proposed `DEF-RUNBOOK` §3 decision-table row and the census trigger (text only — the runbook is not edited)

| flag / choice | SHIPS | validated by | if D/E/F look different |
|---|---|---|---|
| **`VCC_CD4_DEFLATE_BASIS=V5_abc`** (251) — the PINNED A/B/C basis, k = 5, in place of any `VCC_CD4_DEFLATE_K` | **Unset until the 2026-09-05 `cd4svd5` board point is read.** If it reads board raw `pds` up against the incumbent's 0.83572 (246 predicted +0.009 … +0.025 at 0.7 … 2.0×), export **`VCC_CD4_DEFLATE_BASIS=V5_abc`** on D/E/F (add it to `def_pipeline.DEF_ENV` and to `assert_build_carries_arms`, §7) and **never `VCC_CD4_DEFLATE_K`** (the two are refused together at import). Artifact `data/donors/cd4_deflation_basis/V5_abc.npz` sha `ba7ade3f…`, `V` sha `dabf88f7…`, `U_5` sha `f06e0262…` (= the object the board point scores, to 1e−15); the builder refuses a file, a `V`, a `k`, a provenance sha or a gene axis that is not the pin's. Composes with 247's carrying-row sizing: the deflation happens before `solve_cd4_dose` and the deflated object is sized at `g_cd4` (251 KA 3: the fraction rises ×1.066 / ×1.051 on `real` / `efr`). | 250 (pinned A/B/C basis on four disjoint panels: `pds` +0.0153 / +0.0037 / +0.0073 / +0.0107 at −0.1 … −2.0 % cosine, price +0.0061 / +0.0009 / +0.0014 / +0.0030, never negative, better than the re-fit on 4/4); **251**: the builder path == 250's construction to 0.0e+00 on `real` / `efr` (24/24 rows), byte-inert unset (`abeb0749…`, `8df4dadd…`), numeric identity with `cd4svd5` on A/B/C (2.4e−15), the rehearsal on 238's ef-bundle through the driver (§5) | **Census trigger — CD4 carried rows on the D/E/F panel (`effects.cd4.fully_supported_targets` + the partial rows = `effects.cd4.deflation.n_carried`): < 250/300 → a day-of re-fit reads ≤ 0 (250: 174 / 186 carried → −0.0002 / −0.0026) and must not be used; the pinned path stays ≥ 0 there (+0.0009 / +0.0014) and is the only deflation that may ship on such a draw.** Read `effects.cd4.deflation.captured_frac_of_centred` off the meta (the builder prints it and `def_postbuild` echoes it): 0.10–0.14 measured on the four disjoint draws (0.117 / 0.103 / 0.133 / 0.125), 0.198 on A/B/C; **< 0.05 or > 0.19 is unmeasured territory — stop and read 250 §6.2 before shipping**. `own_top_k_energy_frac_of_centred` (0.198–0.232 on every panel, transferring or not) is recorded and must NOT gate. Never change k on the day. |

The trigger in one line: **`n_carried < 250` ⇒ re-fit forbidden (reads ≤ 0), pinned allowed (reads ≥ 0);
`captured_frac_of_centred ∉ [0.05, 0.19]` ⇒ unmeasured, do not ship without a reading.**

## 7c. Honest limits

- **KA 2 is numeric, not byte, identity** (§3): 2.4e−15 on the staged effect matrix, the fraction and the
  energy numbers identical. A gene-space projector cannot reproduce a row-space product's bytes; the
  pre-declared criterion was ≤ 1e−12 and the measurement is 2.4e−15. Anyone who needs the literal
  `a28d5be2…` bytes has them from `VCC_CD4_DEFLATE_K=5`, which is untouched (`k5` row).
- **No board member was measured and no `cell-eval2` was run.** KA 3 reproduces 250's PLOLO rows through
  the builder; it is not new evidence that the `pds` gain transfers. The 2026-09-05 `cd4svd5` point will
  calibrate the object.
- **The rehearsal reuses 244's operator** through the driver's skip branch (accepted on its `axis_sha256`,
  247 B#18); the whitener stage's 361 s were not re-run, so the "29.7 min" of 244 is not re-measured here.
  Everything downstream of the operator (repin-panel-support, build, postbuild, prep) ran through the driver.
- **The driver does not print or assert the new variable** (§7); until `DEF_ENV` carries it, an operator
  who forgets to export it gets the un-deflated stack silently — the same class as 240's finding, and the
  reason the runbook row above says to add it to the driver on the day it ships.
- **`pess` / `efp` were not run through the builder's PLOLO path** (the mandate named `real` and `efr`);
  the ef-bundle rehearsal (= 238's `efp`, 186 carried rows) covers the low-coverage diagnostic through the
  builder (0.133378 = 250's `efp|by_abc|5` reading), not its price.
- **The artifact's npz bytes carry a zip timestamp** (`np.savez`), so re-running `basis251.py` reproduces
  `V` bit for bit but not the file's sha256; the registry pins the file as written (`--force` refuses to
  overwrite without a re-pin), and the `V` array sha `dabf88f7…` is the content pin.

## 8. Cost and memory (every ≥ 15 GiB stage alone; the gate logged before each)

| stage | wall | peak RSS |
|---|---:|---:|
| `basis251.py` | 1.2 s | 0.69 GiB |
| `unit251.py` / `envcheck251.py` / `suite251.py` / 235's tests | 1.3 s / 9 s / 56 s / ~2 min | 0.71 / 0.23 / 0.34 / 1.06 GiB |
| `repro251.py` × 6 cells | 7–13 s each | 4.94–5.10 GiB |
| staging `pre` / `off` / `k5` / `pinned` (each alone) | 17 / 12 / 12 / 12 s incl. inputs | 13.70 / 19.95 / 20.11 / 19.74 GiB |
| `def_postbuild` × 6 (KA-1..3, both scripts) | ~3 s each | 0.12 GiB |
| KA package build (`run251ka.sh`, both flags unset) + `def_postbuild` | **1329.8 s** + 1 s | **17.30 GiB** (the meta's own `peak_rss_gib` 17.30) |
| rehearsal (`rehearsal251.sh`): repin-panel-support / build / postbuild / prep | 12.5 / **1333.8** / 0.8 / 17.5 s; launch → prepped artifact **25.6 min** (16:21:16 → 16:46:55Z, incl. the 2-minute gate hold) | 19.96 / **20.05** / 0.13 / 16.70 GiB |
| `postbuild251.py --rehearsal` (16 `def_postbuild` runs) | ~50 s | 0.13 GiB |

Disk: `work/` holds 12 GB — the rehearsal h5ad (12.6 GB, kept for re-verification, gitignored), the four
staged 44 MB effect matrices and their CD4 / union residual npz, the KA build's meta (its h5ad deleted on
match). `submissions/` was not written to.

---

## 9. Files

```
basis251.py        the artifact writer + its three known-answer rows                -> data/donors/cd4_deflation_basis/V5_abc.{npz,json}, work/basis251.{json,log}
patch251.py        the builder change (5 edits) and the def_postbuild change (4) as exact-anchor replacements; --check reproduces the working tree
envcheck251.py     the import-time rules, ten subprocess imports incl. the pre-251 byte copy   -> work/envcheck251.{json,log}
unit251.py         47 unit rows: OFF bit-identical, ON == lib250 / == the re-fit, loader refusals, the real residual on a stub base  -> work/unit251.{json,log}
suite251.py        the builder's 57 tests under both builders                          -> work/suite251.{json,log}; 235's tests -> work/tests235.log
stage251.py/.sh    KA rows 1-2 at staging level (pre / off / k5 / pinned), gated, alone -> work/stage/
compare251.py      row 2: pinned vs 249's re-fit, both readings                        -> work/compare251.json
repro251.py        row 3: 250's pinABC5 rows through the builder on real / efr        -> work/repro251-<panel>-<fold>.{json,log}, repro251-summary.{json,log}
postbuild251.py    def_postbuild's known answers (both scripts; the rehearsal meta and its tampered copies)  -> work/postbuild251{-pre,}.{json,log}
build251.sh        222's build flags + the two passthroughs;  run251ka.sh  the gated package-level row -> work/build/
rehearsal251.sh    KA row 4, the driver on 238's ef-bundle, two gated invocations       -> work/rehearsal/
bind250.py, lib250.py, gate250.py   byte copies of 250's (sha-identical), so the instrument's caches land here
work/make_submission_pm_cd4.PRE.py, work/def_postbuild.PRE.py   the git controls (HEAD 89cb22a) taken before the first edit
work/builder.intended.sha256        the gate's reference (10d0177e...)   work/builder251.diff   the change as git saw it
work/commit-message.txt             the builder commit's message with the known-answer hashes
```

```bash
cd ~/Code/virtual-cell-challenge-2026
uv run python experiments/251-pinned-basis/basis251.py
uv run python experiments/251-pinned-basis/patch251.py --check experiments/251-pinned-basis/work/make_submission_pm_cd4.PRE.py experiments/251-pinned-basis/work/def_postbuild.PRE.py
uv run python experiments/251-pinned-basis/envcheck251.py; uv run python experiments/251-pinned-basis/unit251.py
uv run python experiments/251-pinned-basis/suite251.py; uv run python experiments/235-october-hardening/tests/test_october_hardening.py
for p in real efr; do for f in k562_gwps orion_hct116 orion_hek293t; do uv run python experiments/251-pinned-basis/repro251.py --panel $p --fold $f; done; done
uv run python experiments/251-pinned-basis/repro251.py --summary
bash experiments/251-pinned-basis/stage251.sh                 # pre / off / k5 / pinned, gated; ends with compare251.py
uv run python experiments/251-pinned-basis/postbuild251.py    # then, after the rehearsal: --rehearsal work/rehearsal/<name>.h5ad.meta.json
bash experiments/251-pinned-basis/run251ka.sh; bash experiments/251-pinned-basis/rehearsal251.sh
```

---

## 10. Verdict

**PINNED PATH READY.** `VCC_CD4_DEFLATE_BASIS=V5_abc` is in the builder (commit `c321f76`, builder
`10d0177e…`), pinned by path and three shas to the artifact `data/donors/cd4_deflation_basis/V5_abc.npz`
(`ba7ade3f…`; `V` `dabf88f7…` == 250's bit for bit; `U_5` `f06e0262…` == the shipped `cd4svd5` audit), refused
together with `VCC_CD4_DEFLATE_K`, refused on any pin or gene-axis mismatch, byte-inert unset at staging
(`abeb0749…`) and package (`8df4dadd…`) level with the re-fit path untouched (`a28d5be2…`); on A/B/C it
reproduces `cd4svd5`'s effects at **numeric, not byte, identity** (max |Δ| 2.4e−15, fraction and energies
identical — one projector in two floating-point associations, §3); on 233's `real` and 238's `efr` it
reproduces 250's pinned rows through 162's instrument at **0.0e+00** on every cell (`pds` +0.015324 / +0.010649,
aggregate cosine +0.06 % / −1.99 %, diagnostic 0.1171 / 0.1249); and the Oct-22 driver ran it end to end on
238's ef-bundle with the flag exported and nothing else (`def_postbuild` 40/40 through the driver, the pinned
block and the diagnostic 0.13338 in the meta, `vcc prep --dry-run --contexts D,E,F` clean, the builder
restored). Whether it ships is the 2026-09-05 `cd4svd5` board point's call; the runbook row and the census
trigger to adopt are in §7b, and the three-line driver change that makes the flag printed and asserted is in
§7. **Not submitted; nothing under `submissions/` written.**
