"""428: production patch for the zero-gene repair (`VCC_DE_ZERO_REPAIR=<floor_log2>`, e.g. -5).

- Installs src/vcc2026/zero_repair.py: `zero_repair` copied verbatim from experiments/428-zero-gene-audit/repair428.py,
  whose H1 test is in RESULTS.md.
- Builder: parses the knob (requires `VCC_PERT_SCALE_BALLAST`, the ballast path, where it runs right after
  `ballast_rows`), applies the repair to every block, records `de_zero_repair` in the config only when set, and adds the
  per-block `zero_repaired` / `zero_counts_moved` to the rc diagnostics.
- Unset is a no-op: no code runs and no key is written.

usage: uv run python patch_zr428.py [--builder PATH] [--module PATH]   (defaults: the production files, in place)
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def rep(s: str, old: str, new: str, what: str) -> str:
    n = s.count(old)
    if n != 1:
        raise SystemExit(f"{what}: anchor found {n} times (need exactly 1)")
    return s.replace(old, new)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--builder", default=str(ROOT / "scripts/make_submission_pm_cd4.py"))
    ap.add_argument("--module", default=str(ROOT / "src/vcc2026/zero_repair.py"))
    a = ap.parse_args()
    mod = Path(a.module)
    if mod.exists():
        raise SystemExit(f"{mod} exists already")
    src = (HERE / "repair428.py").read_text()
    mod.write_text(src.replace('"""428: zero-gene repair.', '"""Zero-gene repair (experiment 428; copied verbatim from '
                                                               'experiments/428-zero-gene-audit/repair428.py).', 1))
    b = Path(a.builder)
    s = b.read_text()
    if "VCC_DE_ZERO_REPAIR" in s:
        raise SystemExit(f"{b} already carries the 428 patch")
    anchor = "DE_MDE_ENV = os.environ.get(\"VCC_DE_MDE\", \"\").strip()\n"
    knob = '''# `VCC_DE_ZERO_REPAIR=<floor_log2>` -- arm 428 (2026-09-23): a gated gene (control mean CPM > 5) that a block's 400 cells
# never express has a delivered arithmetic-mean lfc of log2(1e-9 / ctrl) ~ -35. The carrier design's knockdown was so
# deep that pass 1 emitted nothing, and ratio-mode ballast copies the carriers' composition. 426a's package had 24 such
# (context, target, gene) cells on A/B/C; on H1, 61 of them cost nmae +0.081, and this repair brought it back to +0.0005
# (experiments/428-zero-gene-audit). The repair gives the carriers n = max(1, round(ctrl_cpm * 2**floor * sum(carrier
# libs) / 1e6)) single counts in distinct random carrier cells. Each cell pays with one count from its own largest
# sub-gate gene, so every library and every other tested gene is unchanged. The gene stays called down. It runs after
# ballast_rows. Requires VCC_PERT_SCALE_BALLAST. Unset is a no-op: no code runs and no key is written.
DE_ZERO_REPAIR_ENV = os.environ.get("VCC_DE_ZERO_REPAIR", "").strip()
DE_ZERO_REPAIR = None
if DE_ZERO_REPAIR_ENV:
    try:
        DE_ZERO_REPAIR = float(DE_ZERO_REPAIR_ENV)
    except ValueError:
        raise SystemExit(f"VCC_DE_ZERO_REPAIR must be a float log2 floor: {DE_ZERO_REPAIR_ENV!r}")
    if not -10.0 <= DE_ZERO_REPAIR <= -1.0:
        raise SystemExit(f"VCC_DE_ZERO_REPAIR must be in [-10, -1]: {DE_ZERO_REPAIR_ENV!r}")
_ZR_GENE_POS: dict[str, int] = {}
'''
    s = rep(s, anchor, knob + anchor, "builder knob")
    # the ballast-path requirement is checked where PERT_SCALE_BALLAST is known (after its own parsing)
    g_old = "BALLAST_LB = 995_000\n"
    g_new = ("if DE_ZERO_REPAIR is not None and PERT_SCALE_BALLAST is None:\n"
             "    raise SystemExit(\"VCC_DE_ZERO_REPAIR needs VCC_PERT_SCALE_BALLAST (it runs after ballast_rows)\")\n" + g_old)
    s = rep(s, g_old, g_new, "builder guard")
    blk_old = '''                            if block.shape[0] != int(spec.cells_per_pert):
                                raise AssertionError(f"{context}/{target}: ballast block has {block.shape[0]} rows")
'''
    blk_new = blk_old + '''                            if DE_ZERO_REPAIR is not None:
                                # 428: repair gated genes with a zero block total (delivered lfc ~ -35)
                                from vcc2026 import zero_repair as _zr  # noqa: PLC0415
                                if not _ZR_GENE_POS:
                                    _ZR_GENE_POS.update({str(g): i for i, g in enumerate(spec.genes)})
                                block, _zr_diag = _zr.zero_repair(
                                    np.asarray(block), int(spec.cells_per_pert) - int(ballast_k),
                                    np.asarray(rc_control.keep, dtype=bool), np.asarray(rc_control.ref_lin, dtype=np.float64),
                                    _ZR_GENE_POS.get(str(target), -1),
                                    np.random.default_rng(RC_JITTER_SEED_BASE + 7919 + row // spec.cells_per_pert),
                                    floor_log2=DE_ZERO_REPAIR)
                                rc_diag = {**rc_diag, **_zr_diag}
'''
    s = rep(s, blk_old, blk_new, "builder block")
    cfg_old = '''            **({} if CARRIER_GAMMA is None else {"carrier_gamma": CARRIER_GAMMA}),
'''
    cfg_new = cfg_old + '''            **({} if DE_ZERO_REPAIR is None else {"de_zero_repair": {"floor_log2": DE_ZERO_REPAIR,
                                                                    "where": "after ballast_rows; carriers only; paid by each cell's largest sub-gate gene"}}),
'''
    s = rep(s, cfg_old, cfg_new, "builder config")
    b.write_text(s)
    print(f"module {mod}: sha {hashlib.sha256(mod.read_bytes()).hexdigest()}")
    print(f"builder {b}: patched, sha {hashlib.sha256(b.read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
