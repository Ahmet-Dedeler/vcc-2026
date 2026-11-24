"""317 — is `--effect-scale s` the same knob as `VCC_ALPHA -> 0.982*s`?

`profile_from_effect`'s own docstring says `effect_scale` "composes multiplicatively with the
VCC_ALPHA env var, which scales the same". That is only exact if everything between `alpha` and
the profile is SCALE-EQUIVARIANT. This checks the two steps that could break it, on the real
production functions, with a stub operator (no 1.29 GB load needed: the claim is about scaling,
not about the operator's values).

  step 1  `apply_v9_whitening` — its preserve factor is `median|E_cols| / median|W E_cols|`,
          a ratio of norms of the same input, so it must be INVARIANT under E -> cE, and the
          output must be exactly c times the output at E.
  step 2  `EXP062.add_post_whiten_residual` at a SOLVED FRACTION — the residual is added at a
          fraction of the row's own energy, so it too must be equivariant.
"""
import os, sys, importlib.util
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[2]
os.environ.setdefault("VCC_ALPHA", "0.982")
os.environ.setdefault("VCC_MASS_FIX", "minnorm")
os.environ.setdefault("VCC_WHITEN_GAMMA", "0.65")
os.environ.setdefault("VCC_WHITEN_LAM", "0.0")
os.environ.setdefault("VCC_CD4_TARGET_F", "0.10")
sys.path.insert(0, str(REPO / "src"))
spec = importlib.util.spec_from_file_location("prod317", REPO / "scripts/make_submission_pm_cd4.py")
B = importlib.util.module_from_spec(spec); sys.modules["prod317"] = B
spec.loader.exec_module(B)
print(f"[317] builder ALPHA={B.ALPHA} MASS_FIX={B.MASS_FIX} GAMMA={B.WHITEN_GAMMA}")

rng = np.random.default_rng(0)
G, T, C = 400, 30, 250
E = rng.normal(size=(T, G)) * rng.lognormal(0, 1, size=(1, G))
cols = np.sort(rng.choice(G, C, replace=False))


class StubWhitener:
    """A linear operator with the same call signature; linearity is all this test needs."""
    def __init__(self, n, seed=1):
        r = np.random.default_rng(seed)
        A = r.normal(size=(n, n)) / np.sqrt(n)
        self.M = A @ A.T + np.eye(n)

    def apply(self, X, gamma, lam, renorm=False):
        w, V = np.linalg.eigh(self.M + lam * np.eye(self.M.shape[0]))
        return X @ (V * np.power(w, -gamma / 2.0)) @ V.T


W = StubWhitener(C)
out1, a1 = B.apply_v9_whitening(E, W, cols, gamma=0.65, lam=0.0)
for c in (0.70, 0.85, 1.20, 1.45, 2.00):
    out2, a2 = B.apply_v9_whitening(E * c, W, cols, gamma=0.65, lam=0.0)
    rel_scale = abs(a2["operator_preserve_scale"] / a1["operator_preserve_scale"] - 1.0)
    rel_out = float(np.abs(out2 - c * out1).max() / max(np.abs(c * out1).max(), 1e-300))
    print(f"  c={c:<5} preserve-scale invariant: rel {rel_scale:.2e}   "
          f"output == c * output(E): max rel {rel_out:.2e}   "
          f"{'PASS' if rel_scale < 1e-12 and rel_out < 1e-12 else 'FAIL'}")

# step 2 — the CD4 residual at a solved fraction, on the production object
cd4 = B.EXP062.CD4Prepared(mode="complete_mean3",
                           mean_effect=rng.normal(size=(T, G)) * 0.1,
                           gene_votes=np.ones((T, G)),
                           target_condition_count=np.full(T, 3),
                           condition_details=[])
r1 = B.EXP062.add_post_whiten_residual(out1, cd4, 0.0947, cols)[0]
print("\n  `add_post_whiten_residual` at a FIXED fraction is equivariant in the whitened input:")
for c in (0.70, 1.45, 2.00):
    r2 = B.EXP062.add_post_whiten_residual(c * out1, cd4, 0.0947, cols)[0]
    # the residual bank is NOT scaled with the candidate, so the sum is c*whitened + resid,
    # i.e. equivariance holds only if the bank is scaled too — check both readings
    rel_both = float(np.abs(r2 - c * r1).max() / max(np.abs(c * r1).max(), 1e-300))
    print(f"    c={c:<5} residual added to c*E vs c*(residual added to E): max rel "
          f"{rel_both:.2e}")
