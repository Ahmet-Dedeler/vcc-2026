# 429: put 423's exact-U down-gene pass (`dn_k250_zs`) into production and test it stacked with 426's MDE and 428's zero repair. Mandate from the PI, 2026-09-23 06:00Z

**Read first.**
- `experiments/423-rank-preserving-magnitude/RESULTS.md`, especially §6 (the rule), §8 (the builder spec, not applied)
  and §5 (float32 safety);
- `experiments/426-mde-magnitude/RESULTS.md` (the MDE: an A-channel magnitude add, now `VCC_DE_MDE` in production);
- `experiments/428-zero-gene-audit/RESULTS.md` (`VCC_DE_ZERO_REPAIR`);
- `docs/FINDINGS.md` §0–§1.

**Why.**
- 423's `dn_k250_zs` is a zero-sum carrier↔carrier move on down-regulated "mass-holder" genes, kept inside each cell's
  control interval. It gave nmae −0.0026 on H1 with the other five members bit-identical: the cleanest DE-magnitude
  arm we have.
- 426's MDE works earlier, on the class split before the |z| design, and gives nmae −0.0096. Nobody has measured
  whether the two stack.

**Territory.** Only `experiments/429-zs-production/`.
- No git, no `vcc`, no `cell-eval2`. Do not edit `scripts/`, `src/` or any other experiment's files. The PI applies
  production patches.
- H1 runs are strictly serial, ≤ 8 GB RSS each, behind the memory gate (`memory_pressure` ≥ 40 % free and free+inactive
  ≥ 20 GiB). A production build (~20 GB) may be running at the same time.
- Nobody else works on the DE channel today.

**Tasks.**
- **(a) Stacked H1 test.** Copy 426's harness (`run426.py`, `k_emit426.py`, `drive426.sh`) into your directory. Add the
  zs post-pass from `rp423.py` (vendored or imported read-only) after the ballast block. Report 3 seeds (0, 1, 2) with
  `--bulk ps` for:
  - base (must reproduce 416's `base_s{S}_bps` stats bit-for-bit: the known answer);
  - base + zs (must reproduce 423's §7 production-bulk numbers for `dn_k250_zs`: the second known answer);
  - MDE (`--mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3`; must reproduce 426's numbers: the third known answer);
  - MDE + zs;
  - MDE + zs + `--zero-repair -5`. On the base carriers this should equal MDE + zs: there are no zeros on H1 base.

  Give all six members, with sd over seeds and overall at board slopes (pds 2.2138, mse −1.051, nmae −1.6532,
  fid 3.3645, reach 1.1292, jac 2.7263; overall = Σ slope·Δ/6). Say whether zs keeps U exactly on top of the MDE's
  class splits (423's verify).
- **(b) Production patch (prepared, not applied).**
  - `src_rank_mag.py`: the zs path vendored per 423 §8, with a module docstring. The PI copies it to
    `src/vcc2026/rank_mag.py`.
  - `patch_zs429.py`: adds `VCC_DE_MAG=dn_k250_zs` to `scripts/make_submission_pm_cd4.py`. Anchor it on exact text
    and refuse on missing or duplicate anchors (see `experiments/428-zero-gene-audit/patch_zr428.py` for the house
    style). Include 423 §8's scorer-order control block. Run order after ballast: zero repair, then zs; say why. The
    config key is written only when set.
  - `test_zs429.py`, with a `TEST_BUILDER` / `TEST_MODULE` dry-copy mode (copy the builder into `scripts/.dry429_*.py`
    for the parse test and delete the copy after; the builder resolves ROOT from its own path):
    - parse and refusals;
    - the module's invariants on real H1 blocks;
    - a known-answer comparison of the vendored module against `rp423.py` on the same blocks.
  - Apply the patch only to copies. Report the shas.
- **(c) Priced recommendation.**
  - The stacked arm's H1 table.
  - Its predicted board range, with the ×1 … ×2.3 pure-magnitude reasoning from 423 §6. Correct that reasoning if the
    09-24 board read of 426a has landed by then (`docs/BOARD-LEDGER.md`).
  - Whether zs should ride the next arm with the MDE, or wait.

**Deliverable.** `experiments/429-zs-production/RESULTS.md` covering:
- the known answers;
- the table;
- the patch/test shas;
- the recommendation;
- a compliance note (memory, serial runs, files written).
