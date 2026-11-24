"""496 (478's route): an optional `mde_add_sign` in the emitter's lfc_target. A sign source for the MDE add on complement
genes (neither head nor block, i.e. the uncalled genes), e.g. the CLEAN4 per-gene prior where |P_down - 0.5| >= tau;
0 = keep the delivered sign. Only the add's direction changes: the gene's base LFC is untouched, and fu/fd re-balance the
mass exactly as before (478 priced this construction analytically: H1 nmae -0.0015, reach +0.002..+0.0035 at tau 0.05).
Absent key = the emitter is unchanged (no code path runs, no diagnostic key is written).
usage: uv run python experiments/496-uncalled-sign-hook/patch_emitter496.py <rc_emit_k.py> [--check]
"""
import sys
from pathlib import Path

path = Path(sys.argv[1]); src = path.read_text()
old1 = """                base_t = np.where(blk, _t2(kf), E0v)
                sgn = np.sign(base_t)
                elig = (Ms > 0) & np.isfinite(E0v) & np.isfinite(madd) & (madd > 0) & (sgn != 0)
"""
new1 = """                base_t = np.where(blk, _t2(kf), E0v)
                sgn = np.sign(base_t)
                sign_rec = {}
                if "mde_add_sign" in lfc_target:
                    # 496 (478's route): the add's direction on complement genes (not head, not block: the uncalled genes)
                    # comes from a sign source, e.g. the CLEAN4 prior where confident; 0 = keep the delivered sign. The base
                    # LFC is untouched; fu/fd below re-balance the mass as for any sign vector.
                    asg = np.asarray(lfc_target["mde_add_sign"], dtype=np.float64)
                    if asg.shape != (n_keep,):
                        raise ValueError("mde_add_sign must match the gated gene axis")
                    use_s = ~head & ~blk & (asg != 0)
                    sign_rec = {"mde_n_sign_cover": int(use_s.sum()), "mde_n_sign_flip": int((use_s & (asg != sgn)).sum())}
                    sgn = np.where(use_s, asg, sgn)
                elig = (Ms > 0) & np.isfinite(E0v) & np.isfinite(madd) & (madd > 0) & (sgn != 0)
"""
old2 = """                           "mde_scope_code": float({"all": 0, "block": 1, "nonblock": 2}[scope])}
"""
new2 = """                           "mde_scope_code": float({"all": 0, "block": 1, "nonblock": 2}[scope]), **sign_rec}
"""
for o, n in ((old1, new1), (old2, new2)):
    if src.count(o) != 1:
        raise SystemExit(f"anchor count {src.count(o)} != 1: {o[:60]!r}")
if "mde_add_sign" in src:
    raise SystemExit("already patched")
src = src.replace(old1, new1).replace(old2, new2)
compile(src, str(path), "exec")
if "--check" in sys.argv:
    print("patch applies cleanly (not written)")
else:
    path.write_text(src)
    print(f"patched {path}")
