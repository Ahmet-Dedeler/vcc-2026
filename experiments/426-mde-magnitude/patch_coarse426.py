"""426h: production patch for the count-balanced MDE (`VCC_DE_MDE="c,cap,min_count,coarse"`).

Emitter (src/vcc2026/rc_emit_k.py), copied from experiments/426-mde-magnitude/k_emit426.py where it was measured:
- `split_coarse_extremes`, which moves a gene's counts only between the lowest- and highest-library classes, in whole
  counts, with no polish;
- in prior mode, with `mde_coarse` in the spec, the non-block MDE genes take targets relative to their unsteered split's
  delivered A and are split coarse. Block genes keep `split_for_A`.
Without `mde_coarse`, every line that runs is the old one.
Builder: VCC_DE_MDE takes an optional 4th field `coarse`, which passes `mde_coarse` and records `"coarse": true` in the
`de_mde` config; the 3-field form writes exactly the old config. DE_SIGN_EMITTER_SHA256 is re-pinned.

Why: 429a's per-block dump shows context A delivering no MDE (`mde_lambda` median 0). split_for_A's polish drains the
middle class, ~6,400 counts per block on H1, and A has only ~4k of slack under the capacity margin. The coarse path cuts
the per-class shift ~17x at the same H1 nmae (426 RESULTS, 426h).

usage: uv run python patch_coarse426.py [--emitter PATH] [--builder PATH]
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def rep(s, old, new, what):
    n = s.count(old)
    if n != 1:
        raise SystemExit(f"{what}: anchor found {n} times (need exactly 1)")
    return s.replace(old, new)


def k426_block(start: str, end: str) -> str:
    s = (HERE / "k_emit426.py").read_text()
    i = s.index(start)
    return s[i:s.index(end, i)]


def patch_emitter(path: Path) -> str:
    s = path.read_text()
    if "split_coarse_extremes" in s:
        raise SystemExit(f"{path} already carries the coarse patch")
    fn = k426_block("def split_coarse_extremes(M0, A_t, Lk, n):", "def split_for_A(Ms, A_t, Lk, nk, caps, A_tol=1e-12):")
    s = rep(s, "def split_for_A(Ms, A_t, Lk, nk, caps, A_tol=1e-12):",
            fn + "def split_for_A(Ms, A_t, Lk, nk, caps, A_tol=1e-12):", "emitter fn")
    s = rep(s, "            base_ = split_massprop(Ms, w_)\n", '''            base_ = split_massprop(Ms, w_)
            # 426h (`mde_coarse`): a non-block gene's MDE target is taken relative to the arithmetic mean its UNSTEERED
            # (mass-proportional) split already delivers, and its move uses only the extreme-class channel (no polish),
            # so an A-balanced set of MDE moves is also count-balanced per class.
            mde_coarse = move_m is not None and bool(lfc_target.get("mde_coarse", False))
            A_ms = arith_of_split(base_, Lk, n) if mde_coarse else None
''', "emitter base")
    old_t = '''                if lam_m > 0.0:
                    A_t_ = np.where(steer_set, A0 * np.power(2.0, lam * move + lam_m * move_m), np.nan)'''
    new_t = '''                if lam_m > 0.0 and A_ms is not None:
                    A_t_ = np.where(blk, A0 * np.power(2.0, lam * move + lam_m * move_m),
                                    np.where(steer_set, A_ms * np.power(2.0, lam_m * move_m), np.nan))
                elif lam_m > 0.0:
                    A_t_ = np.where(steer_set, A0 * np.power(2.0, lam * move + lam_m * move_m), np.nan)'''
    s = rep(s, old_t, new_t, "emitter targets")
    old_sp = '''                Mk_ = base_.copy()
                st_ = np.isfinite(A_t_)
                if st_.any():
                    Mk_s_, A_ach_, chi_, clo_ = split_for_A(Ms[st_], A_t_[st_], Lk, nk, caps[st_])
                    Mk_[st_] = Mk_s_
                else:
                    A_ach_, chi_, clo_ = np.array([]), np.array([]), np.array([])'''
    new_sp = k426_block('''                Mk_ = base_.copy()
                st_ = np.isfinite(A_t_)
                if st_.any() and lam_m > 0.0 and mde_coarse:''', "                ok_ = bool(np.all(Mk_.sum(0) - (nk * Lk - margin_) <= 0))")
    new_sp = new_sp.rstrip("\n")
    s = rep(s, old_sp, new_sp, "emitter split")
    path.write_text(s)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def patch_builder(path: Path, old_em: str, new_em: str) -> None:
    s = path.read_text()
    if "de_mde_coarse" in s or "coarse\" in DE_MDE" in s:
        raise SystemExit(f"{path} already carries the coarse patch")
    s = rep(s, f'DE_SIGN_EMITTER_SHA256 = "{old_em}"', f'DE_SIGN_EMITTER_SHA256 = "{new_em}"', "builder pin")
    old_parse = '''    try:
        DE_MDE = tuple(float(x) for x in DE_MDE_ENV.split(","))
    except ValueError:
        raise SystemExit(f"VCC_DE_MDE must be 'c,cap,min_count' floats: {DE_MDE_ENV!r}")'''
    new_parse = '''    _mde_parts = DE_MDE_ENV.split(",")
    DE_MDE_COARSE = len(_mde_parts) == 4 and _mde_parts[3].strip() == "coarse"   # 426h: count-balanced non-block moves
    if len(_mde_parts) not in (3, 4) or (len(_mde_parts) == 4 and not DE_MDE_COARSE):
        raise SystemExit(f"VCC_DE_MDE must be 'c,cap,min_count' or 'c,cap,min_count,coarse': {DE_MDE_ENV!r}")
    try:
        DE_MDE = tuple(float(x) for x in _mde_parts[:3])
    except ValueError:
        raise SystemExit(f"VCC_DE_MDE must be 'c,cap,min_count' floats: {DE_MDE_ENV!r}")'''
    s = rep(s, old_parse, new_parse, "builder parse")
    s = rep(s, '''DE_MDE_ENV = os.environ.get("VCC_DE_MDE", "").strip()
DE_MDE = None
''', '''DE_MDE_ENV = os.environ.get("VCC_DE_MDE", "").strip()
DE_MDE = None
DE_MDE_COARSE = False
''', "builder init")
    s = rep(s, '''                                **({} if de_mde_add is None else {"mde_add": de_mde_add, "mde_scope": "all"}),''',
            '''                                **({} if de_mde_add is None else {"mde_add": de_mde_add, "mde_scope": "all"}),
                                **({"mde_coarse": True} if (de_mde_add is not None and DE_MDE_COARSE) else {}),''', "builder lfc_target")
    old_cfg = '''                    **({} if DE_MDE is None else {"de_mde": {"c": DE_MDE[0], "cap": DE_MDE[1], "min_count": DE_MDE[2],
                                                            "scope": "all", "law": "min(c/sqrt(ref_lin*median_lib/1e6), cap)"}}),'''
    new_cfg = '''                    **({} if DE_MDE is None else {"de_mde": {"c": DE_MDE[0], "cap": DE_MDE[1], "min_count": DE_MDE[2],
                                                            "scope": "all", "law": "min(c/sqrt(ref_lin*median_lib/1e6), cap)",
                                                            **({"coarse": True} if DE_MDE_COARSE else {})}}),'''
    s = rep(s, old_cfg, new_cfg, "builder config")
    path.write_text(s)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--emitter", default=str(ROOT / "src/vcc2026/rc_emit_k.py"))
    ap.add_argument("--builder", default=str(ROOT / "scripts/make_submission_pm_cd4.py"))
    a = ap.parse_args()
    old_em = hashlib.sha256(Path(a.emitter).read_bytes()).hexdigest()
    if old_em != "72be7e166d8673b64ecebec81f28234c348b02c16879d5ce6937c18785c0f61d":
        raise SystemExit(f"emitter sha {old_em} is not 72be7e16...")
    new_em = patch_emitter(Path(a.emitter))
    patch_builder(Path(a.builder), old_em, new_em)
    print(f"emitter {a.emitter}: 72be7e16... -> {new_em}")
    print(f"builder {a.builder}: patched, sha {hashlib.sha256(Path(a.builder).read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
