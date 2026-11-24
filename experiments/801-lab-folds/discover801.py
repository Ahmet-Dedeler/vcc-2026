"""801 -- find the data closure of dump801 (and engine801's A/B/C path) WITHOUT the heavy operator: the whitening operator
is replaced by an identity stub (its path and json are still recorded), so the run stays light (< 2 GB) and its NUMBERS
ARE MEANINGLESS. Its only output is the list of repo files every code path opens (guard801's audit hook), which the cloud
runs stage. Known-answer tolerances are disabled here on purpose (TOL = inf); the real runs keep 1e-12.

    GUARD801_CLOSURE=work/closure/<name>.txt uv run python experiments/801-lab-folds/discover801.py dump <panel> <fold>
    GUARD801_CLOSURE=work/closure/abc-<fold>.txt uv run python experiments/801-lab-folds/discover801.py abc <fold>
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
from pathlib import Path  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import guard801  # noqa: E402,F401


class _Eye:
    """An identity matrix that never exists: M @ I == M, I.T == I."""
    __array_ufunc__ = None

    def __rmatmul__(self, M):
        return M

    @property
    def T(self):
        return self


class StubOperator:
    """Identity stand-in for plolo.Operator: records the operator files, never loads the matrix."""

    def __init__(self, path, *, name=None):
        path = Path(path).resolve()
        self.path = path
        guard801.CLOSURE.add(str(path.relative_to(Path(guard801.ROOT))))
        js = path.with_suffix(".json")
        if js.exists():
            js.read_text()
        self.name = f"stub:{path.stem}"

    def _check_axis(self, wcols):
        import numpy as np
        self.evals = np.ones(len(wcols))
        self.evecs = _Eye()
        return None

    def apply(self, M, wcols, gamma, lam, *, check_axis=True):
        return M.copy()


def main():
    mode = sys.argv[1]
    if mode == "dump":
        panel, fold = sys.argv[2], sys.argv[3]
        import dump801 as D
        D.L.plolo.Operator = StubOperator
        D.TOL = float("inf")
        D.OUT = HERE / "work" / "discovery"
        D.OUT.mkdir(parents=True, exist_ok=True)
        D.FOLD_OPS = "held"
        D.main(panel, fold, None, [], None, ka522=True)
    elif mode == "dumpL":                       # a WORKING panel, with dump801's own CLI spec (held ops never read)
        panel, fold = sys.argv[2], sys.argv[3]
        import dump801 as D
        D.L.plolo.Operator = StubOperator
        D.TOL = float("inf")
        D.OUT = HERE / "work" / "discovery"
        D.OUT.mkdir(parents=True, exist_ok=True)
        D.FOLD_OPS = "held"
        W = HERE / "work" / panel
        folds = ["k562_gwps", "orion_hct116", "orion_hek293t"]
        spec = dict(bundle=str((W / "bundle").resolve()), rules="def", cd4_root=str((W / "cd4/cd4_compact").resolve()),
                    full_op=str(Path("none").resolve()), folds=folds, meta=None,
                    synthetic=not (W / "bundle/frames_src.txt").exists(),  # kit --frames: the bundle's own controls
                    fold_ops={h: (W / "ops" / f"held_{h}.npz").resolve() for h in folds})
        D.main(panel, fold, spec, [], None, ka522=False)
    elif mode == "abc":
        sys.path.insert(0, str(Path(guard801.ROOT) / "experiments/522-october-bulk-key"))
        import fold522
        S = fold522.Abc(sys.argv[2], check=True)
        print("abc", S.held, S.n, "targets")
    else:
        raise SystemExit(mode)


if __name__ == "__main__":
    main()
