"""426: production patch for the MDE-structured DE magnitude (`VCC_DE_MDE="c,cap,min_count"`).

Emitter (src/vcc2026/rc_emit_k.py, prior mode): an optional `mde_add` array in the prior spec. It is copied verbatim
from experiments/426-mde-magnitude/k_emit426.py, where it was measured, without the debug print. Absent `mde_add`,
every line that runs is the old one.
Builder (scripts/make_submission_pm_cd4.py):
- parses the knob and requires VCC_DE_SIGN_PRIOR;
- computes the per-context add from the context's own controls;
- passes it in the lfc_target dict;
- records `de_mde` in the config only when set;
- re-pins DE_SIGN_EMITTER_SHA256 to the patched emitter.

usage: uv run python patch_mde426.py [--emitter PATH] [--builder PATH]   (defaults: the production files, in place)
Refuses when any anchor is missing or not unique, or when a file already carries the patch.
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


def block_from_k426(start: str, end: str) -> str:
    """The exact text between two anchors in k_emit426.py (start included, end excluded)."""
    s = (HERE / "k_emit426.py").read_text()
    i = s.index(start)
    j = s.index(end, i)
    return s[i:j]


def patch_emitter(path: Path) -> str:
    s = path.read_text()
    if "mde_add" in s:
        raise SystemExit(f"{path} already carries the 426 patch")
    move_old = "            move = np.where(blk, _t2(kf) - E0v, 0.0)\n"
    move_new = block_from_k426(move_old, "            w_ = (nk * Lk).astype(np.float64)\n")
    s = rep(s, move_old, move_new, "emitter move")
    dem_old = ("            def _demand(lam):\n"
               "                A_t_ = np.where(blk, A0 * np.power(2.0, lam * move), np.nan)\n")
    dem_new = block_from_k426("            def _demand(lam, lam_m=0.0):\n", "                Mk_ = base_.copy()\n")
    s = rep(s, dem_old, dem_new, "emitter demand")
    dose_old = '            dose_rec = {"dose_nominal": 1.0, "dose_delivered": float(lam)}\n'
    dose_new = block_from_k426("            if move_m is not None and lam > 0.0:\n", dose_old) + dose_old
    s = rep(s, dose_old, dose_new, "emitter dose")
    s = rep(s, '"prior_mass_residual": float(_R(kf))})', '"prior_mass_residual": float(_R(kf)), **mde_rec})', "emitter rec")
    path.write_text(s)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def patch_builder(path: Path, emitter_sha: str) -> None:
    s = path.read_text()
    if "VCC_DE_MDE" in s:
        raise SystemExit(f"{path} already carries the 426 patch")
    old_pin = 'DE_SIGN_EMITTER_SHA256 = "dbbb91210a2d808fd03ae291ae3dc1479a30e58c675c5dc9536da6c266e5c09c"\n'
    s = rep(s, old_pin, f'DE_SIGN_EMITTER_SHA256 = "{emitter_sha}"   # 426: + optional mde_add (dbbb9121... before)\n',
            "builder pin")
    knob_anchor = "DE_SIGN_SHRINK_FLOOR = 0.002\n"
    knob = knob_anchor + '''# `VCC_DE_MDE="c,cap,min_count"` -- arm 426 (2026-09-23): MDE-structured magnitude on the DE channel's arithmetic
# means (A-channel only). nmae reads only the reference's significant genes, whose |lfc| is bounded below by the minimum
# detectable effect, ~ 1/sqrt(control mean count per cell). Adding sign x min(c / sqrt(count), cap) to every steered
# gene's delivered LFC therefore gains (2a - 1) x add per gene (a = sign accuracy), as long as add < |lfc|. Genes below
# min_count counts/cell get nothing. Per-gene totals (the bulk: pds, mse's profile half), the pass-2 |z| design (calls,
# jac, reach order) and every delivered sign are unchanged. The add is mass-balanced by one factor on the heavier side.
# H1 (426, 410b construction): see experiments/426-mde-magnitude/RESULTS.md. Requires VCC_DE_SIGN_PRIOR.
# Unset is a no-op: the lfc_target dict and the config are exactly as before.
DE_MDE_ENV = os.environ.get("VCC_DE_MDE", "").strip()
DE_MDE = None
if DE_MDE_ENV:
    try:
        DE_MDE = tuple(float(x) for x in DE_MDE_ENV.split(","))
    except ValueError:
        raise SystemExit(f"VCC_DE_MDE must be 'c,cap,min_count' floats: {DE_MDE_ENV!r}")
    if len(DE_MDE) != 3 or not (0.0 < DE_MDE[0] <= 0.5 and 0.0 < DE_MDE[1] <= 0.5 and 0.0 <= DE_MDE[2] <= 10.0):
        raise SystemExit(f"VCC_DE_MDE must be 'c,cap,min_count' with 0 < c <= 0.5, 0 < cap <= 0.5, 0 <= min_count <= 10: "
                         f"{DE_MDE_ENV!r}")
    if DE_SIGN_PRIOR is None:
        raise SystemExit("VCC_DE_MDE needs VCC_DE_SIGN_PRIOR (the K-class prior-mode emitter)")
'''
    s = rep(s, knob_anchor, knob, "builder knob")
    ctl_old = '''                memory_checkpoint(f"rc control block {context}", max_rss_gib)
'''
    ctl_new = '''                memory_checkpoint(f"rc control block {context}", max_rss_gib)
                if DE_MDE is not None:
                    # 426: the add from this context's own controls (scorer gate axis), c / sqrt(expected count per cell)
                    _cnt = np.asarray(rc_control.ref_lin, dtype=np.float64) * float(np.median(libraries)) / 1e6
                    de_mde_add = np.minimum(DE_MDE[0] / np.sqrt(np.maximum(_cnt, 1e-3)), DE_MDE[1])
                    de_mde_add = np.ascontiguousarray(np.where(_cnt >= DE_MDE[2], de_mde_add, 0.0))
                    print(f"[de_mde] {context}: c {DE_MDE[0]} cap {DE_MDE[1]} min_count {DE_MDE[2]}; add median "
                          f"{float(np.median(de_mde_add)):.4f}, {int((de_mde_add > 0).sum())}/{de_mde_add.size} genes", flush=True)
'''
    s = rep(s, ctl_old, ctl_new, "builder control block")
    init_old = "            de_sign_classes = None\n            de_sign_full = None\n"
    s = rep(s, init_old, init_old + "            de_mde_add = None\n", "builder init")
    lt_old = '''                                **({} if DE_SIGN_PROTECT_Q is None else {"protect_q": DE_SIGN_PROTECT_Q}),
                            },'''
    lt_new = '''                                **({} if DE_SIGN_PROTECT_Q is None else {"protect_q": DE_SIGN_PROTECT_Q}),
                                **({} if de_mde_add is None else {"mde_add": de_mde_add, "mde_scope": "all"}),
                            },'''
    s = rep(s, lt_old, lt_new, "builder lfc_target")
    cfg_old = '''                    **({} if DE_SIGN_PROTECT_Q is None else {"de_sign_protect_q": DE_SIGN_PROTECT_Q}),
'''
    cfg_new = cfg_old + '''                    **({} if DE_MDE is None else {"de_mde": {"c": DE_MDE[0], "cap": DE_MDE[1], "min_count": DE_MDE[2],
                                                            "scope": "all", "law": "min(c/sqrt(ref_lin*median_lib/1e6), cap)"}}),
'''
    s = rep(s, cfg_old, cfg_new, "builder config")
    path.write_text(s)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--emitter", default=str(ROOT / "src/vcc2026/rc_emit_k.py"))
    ap.add_argument("--builder", default=str(ROOT / "scripts/make_submission_pm_cd4.py"))
    a = ap.parse_args()
    before = hashlib.sha256(Path(a.emitter).read_bytes()).hexdigest()
    if before != "dbbb91210a2d808fd03ae291ae3dc1479a30e58c675c5dc9536da6c266e5c09c":
        raise SystemExit(f"emitter sha {before} is not the pinned dbbb9121...")
    sha = patch_emitter(Path(a.emitter))
    patch_builder(Path(a.builder), sha)
    print(f"emitter {a.emitter}: dbbb9121... -> {sha}")
    print(f"builder {a.builder}: patched, sha {hashlib.sha256(Path(a.builder).read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
