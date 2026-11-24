"""429: production patch for 423's exact-U down-gene pass (`VCC_DE_MAG=dn_k250_zs`).

- Installs src/vcc2026/rank_mag.py: a byte-for-byte copy of experiments/429-zs-production/src_rank_mag.py, whose sha256
  must equal MODULE_SHA256 below (the file the H1 tests ran). The builder pins the same sha.
- Builder (scripts/make_submission_pm_cd4.py), five anchored edits; each anchor must occur exactly once:
  1. knob: parses `VCC_DE_MAG` (`<sel>_k<k*100>_zs`, zs mode only), needs VCC_DE_SIGN_PRIOR; pins the module sha;
  2. guard: needs VCC_PERT_SCALE_BALLAST (the pass runs after ballast_rows);
  3. 423 section 8's scorer-order control block per context (count / (lib / 1e6), the evaluator's order; +1.6 GB per
     context, built right after rc_control while the count matrix is still alive; released with rc_control); the
     module's sha is checked here, before any block is emitted;
  4. the pass on every block, after ballast_rows and AFTER the zero repair (428), before to_csr; its per-block
     diagnostics (`zs_*`, flat numerics) join rc_diag;
  5. config: `de_mag` is written only when the knob is set.
- Requires the 428 patch (the zero-repair anchors). Unset is a no-op: no code runs, no control block is built, no key
  is written, and the pass consumes no RNG, so a knob-on build differs from knob-off only in the moved carrier counts.

Run order after ballast: zero repair, then zs, because
  (a) the repair gets first claim on each carrier's sub-gate counts (it pays one count per repaired cell from the cell's
      largest sub-gate gene; zs then budgets from what is left and can never overdraw it; the reverse order could leave
      a carrier with nothing to pay with, 423 saw ~7 budget-exhausted carriers per block);
  (b) the pass's in-block asserts (per-gene totals, libraries, every gene's U, ballast rows, no sign flip) then certify
      the final delivered carriers; nothing edits the block after them;
  (c) zs moves only genes with a finite nonzero delivered lfc, so it sees the repaired genes as what they are after
      repair (finite, down) rather than as -35 zeros; in practice their single counts are frozen (a 1-count cell cannot
      drop to 0 without tying the control zeros), so the repair's delivered floor is kept.

usage: uv run python patch_zs429.py [--builder PATH] [--module PATH]   (defaults: the production files, in place)
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
MODULE_SHA256 = "ed6728327a15f039d89b05bfc8f65281f138242ea914544116d54ec99ecde638"


def rep(s: str, old: str, new: str, what: str) -> str:
    n = s.count(old)
    if n != 1:
        raise SystemExit(f"{what}: anchor found {n} times (need exactly 1)")
    return s.replace(old, new)


KNOB = '''# `VCC_DE_MAG=dn_k250_zs` -- arms 423/429 (2026-09-23): 423's exact-U magnitude pass on the DE channel's arithmetic
# means. The K-class carriers deliver a strongly down-designed, low-expressed gene by parking its total on a few
# "mass-holder" carriers above the control max (every other carrier a zero tied with the control zeros): the Wilcoxon sees
# a strong down gene, the arithmetic mean reports only the bulk total (H1: median |l|/|r| 0.46, sign accuracy 0.65-0.74 on
# these genes). For every down gene the pass moves counts from low-library carriers to high-library carriers, each cell
# strictly inside its own interval between adjacent control values (evaluator order, count / (lib / 1e6)), each library
# held by the cell's own sub-gate genes. Per-gene totals (pds, mse's profile half), every library, the ballast rows and
# every gene's rank sum U (calls, signs) are unchanged; only the arithmetic mean moves (nmae). H1: 423 section 7, nmae
# -0.0023 with the other five members bit-identical; 429, stacked on VCC_DE_MDE and VCC_DE_ZERO_REPAIR. Runs after
# ballast_rows and AFTER the zero repair (experiments/429-zs-production/RESULTS.md: the repair gets first claim on each
# carrier's sub-gate counts, and the pass's asserts certify the final carriers). Module src/vcc2026/rank_mag.py, sha pinned
# here and checked before the first block. Requires VCC_PERT_SCALE_BALLAST and VCC_DE_SIGN_PRIOR (the K-class carriers).
# Only dn_k250_zs is priced. Unset is a no-op: no code runs, no control block is built, no key is written.
DE_MAG_ENV = os.environ.get("VCC_DE_MAG", "").strip()
DE_MAG = None
DE_MAG_MODULE_PATH = ROOT / "src/vcc2026/rank_mag.py"
DE_MAG_MODULE_SHA256 = "__MODULE_SHA256__"
if DE_MAG_ENV:
    import re as _re  # noqa: PLC0415
    _dm = _re.fullmatch(r"(all|q50|q80|up|dn|q80up|q80dn)_k(\\d{3,4})_zs", DE_MAG_ENV)
    if _dm is None:
        raise SystemExit(f"VCC_DE_MAG must look like '<sel>_k<k*100>_zs', e.g. dn_k250_zs (zs mode only): {DE_MAG_ENV!r}")
    if not 100 <= int(_dm.group(2)) <= 1000:
        raise SystemExit(f"VCC_DE_MAG k must be in k100..k1000: {DE_MAG_ENV!r}")
    if DE_SIGN_PRIOR is None:
        raise SystemExit("VCC_DE_MAG needs VCC_DE_SIGN_PRIOR (the K-class carriers it moves between library classes)")
    DE_MAG = DE_MAG_ENV
'''

CTRL = '''                if DE_MAG is not None:
                    # 423/429: the zs pass keeps every carrier inside its interval between the controls' values in the
                    # EVALUATOR's arithmetic (count / (lib / 1e6)); rc_control's count * (1e6 / lib) differs by ulps.
                    if not DE_MAG_MODULE_PATH.exists():
                        raise SystemExit(f"VCC_DE_MAG needs {DE_MAG_MODULE_PATH} (experiments/429-zs-production/src_rank_mag.py)")
                    _dm_sha = hashlib.sha256(DE_MAG_MODULE_PATH.read_bytes()).hexdigest()
                    if _dm_sha != DE_MAG_MODULE_SHA256:
                        raise SystemExit(f"rank_mag.py sha256 {_dm_sha} != pinned {DE_MAG_MODULE_SHA256}")
                    from vcc2026 import rank_mag as _rank_mag  # noqa: PLC0415
                    rc_scorer_control = _rank_mag.scorer_control_block(control_matrix, min_cpm=RC_SCORER_MIN_CPM)
                    print(f"[de_mag] {context}: {DE_MAG}; scorer-order control block, {int(rc_scorer_control.keep.sum())} genes "
                          f"(differs from the emitter gate on {int((rc_scorer_control.keep != rc_control.keep).sum())})", flush=True)
                    memory_checkpoint(f"de_mag scorer control block {context}", max_rss_gib)
'''

PASS = '''                            if DE_MAG is not None:
                                # 423/429: exact-U zero-sum carrier<->carrier magnitude pass on down genes, after the zero
                                # repair (see the knob's comment); asserts totals, libraries, U per gene and the ballast rows
                                from vcc2026 import rank_mag as _rank_mag  # noqa: PLC0415
                                block, _dm_diag = _rank_mag.apply(np.asarray(block), int(spec.cells_per_pert) - int(ballast_k),
                                                                  rc_scorer_control, DE_MAG)
                                rc_diag = {**rc_diag, **_dm_diag}
'''

CFG = '''            **({} if DE_MAG is None else {"de_mag": {"variant": DE_MAG, "module_sha256": DE_MAG_MODULE_SHA256,
                                                    "where": "after ballast_rows and the zero repair; carriers only; zero-sum "
                                                             "between library classes inside each cell's control interval "
                                                             "(evaluator-order control block)"}}),
'''


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--builder", default=str(ROOT / "scripts/make_submission_pm_cd4.py"))
    ap.add_argument("--module", default=str(ROOT / "src/vcc2026/rank_mag.py"))
    a = ap.parse_args()
    src = HERE / "src_rank_mag.py"
    got = hashlib.sha256(src.read_bytes()).hexdigest()
    if got != MODULE_SHA256:
        raise SystemExit(f"{src} sha256 {got} != the tested {MODULE_SHA256}")
    mod = Path(a.module)
    b = Path(a.builder)
    s = b.read_text()
    if "VCC_DE_MAG" in s:
        raise SystemExit(f"{b} already carries the 429 patch")
    if "VCC_DE_ZERO_REPAIR" not in s:
        raise SystemExit(f"{b} lacks the 428 patch (VCC_DE_ZERO_REPAIR); the pass runs after the zero repair")
    if mod.exists():
        raise SystemExit(f"{mod} exists already")
    # 1. knob, after the zero repair's knob
    a1 = "_ZR_GENE_POS: dict[str, int] = {}\n"
    s = rep(s, a1, a1 + KNOB.replace("__MODULE_SHA256__", MODULE_SHA256), "builder knob")
    # 2. guard, where PERT_SCALE_BALLAST is known
    a2 = '    raise SystemExit("VCC_DE_ZERO_REPAIR needs VCC_PERT_SCALE_BALLAST (it runs after ballast_rows)")\n'
    s = rep(s, a2, a2 + 'if DE_MAG is not None and PERT_SCALE_BALLAST is None:\n'
                        '    raise SystemExit("VCC_DE_MAG needs VCC_PERT_SCALE_BALLAST (it runs after ballast_rows and the zero repair)")\n',
            "builder guard")
    # 3. the scorer-order control block: initialised with rc_control, built after its checkpoint, released with it
    a3 = "            rc_control = None\n            if emitter in RC_EMITTERS:\n"
    s = rep(s, a3, "            rc_control = None\n            rc_scorer_control = None   # 429: VCC_DE_MAG's evaluator-order "
                   "control block\n            if emitter in RC_EMITTERS:\n", "builder control init")
    a3b = '                memory_checkpoint(f"rc control block {context}", max_rss_gib)\n'
    s = rep(s, a3b, a3b + CTRL, "builder control block")
    a3c = "            del rc_control\n            rc_control = None\n            gc.collect()\n"
    s = rep(s, a3c, "            del rc_control\n            rc_control = None\n            rc_scorer_control = None\n"
                    "            gc.collect()\n", "builder control release")
    # 4. the pass, after the zero repair inside the ballast path
    a4 = "                                rc_diag = {**rc_diag, **_zr_diag}\n"
    s = rep(s, a4, a4 + PASS, "builder block")
    # 5. config key, only when set
    a5 = ('''            **({} if DE_ZERO_REPAIR is None else {"de_zero_repair": {"floor_log2": DE_ZERO_REPAIR,\n'''
          '''                                                                    "where": "after ballast_rows; carriers only; paid by each cell's largest sub-gate gene"}}),\n''')
    s = rep(s, a5, a5 + CFG, "builder config")
    mod.write_bytes(src.read_bytes())
    b.write_text(s)
    print(f"module {mod}: sha {hashlib.sha256(mod.read_bytes()).hexdigest()}")
    print(f"builder {b}: patched, sha {hashlib.sha256(b.read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
