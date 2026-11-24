#!/usr/bin/env python3
"""690g — validate a def_day dump hook's outputs before any level build reads them (installed as scripts/def_dumphook_check.py).

    uv run python scripts/def_dumphook_check.py <stack_dump.npz> <hook_dump.npz> <bundle> [<bulk_ctx.npz>]

The hook's dump replaces the stack dump as VCC_BULK_LOAD, so it must be the same object with new numbers:
  - the keys bulk, carriers, perts, genes, env all present; perts and genes identical arrays; env byte-identical (598's
    tercilec reads the dump key's carried flags from it);
  - bulk and carriers the same shape and dtype as the stack dump's, every value finite.
The optional per-context npz (672's VCC_BULK_LOAD_CTX) must hold perts and genes identical to the dump's, contexts equal to
the bundle manifest's, and one finite float64 bulk_<ctx> of the dump's bulk shape per context.
Prints one line per check and, last, per context the cosine of the hook's bulk to the stack dump's (mean over targets),
so the log shows how far the hook moved the bulk. Exit 0 iff every check passes.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

NEED = ("bulk", "carriers", "perts", "genes", "env")


def main() -> int:
    if len(sys.argv) not in (4, 5):
        print(__doc__)
        return 2
    base_p, hook_p, bundle = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    ctx_p = Path(sys.argv[4]) if len(sys.argv) == 5 else None
    bad: list[str] = []
    base = np.load(base_p, allow_pickle=False)
    hook = np.load(hook_p, allow_pickle=False)
    miss = [k for k in NEED if k not in hook.files]
    if miss:
        print(f"FAIL hook dump lacks {miss}")
        return 1
    for k in ("perts", "genes"):
        ok = hook[k].shape == base[k].shape and bool(np.all(hook[k] == base[k]))
        print(f"{'PASS' if ok else 'FAIL'} {k} identical to the stack dump's ({hook[k].shape[0]})")
        ok or bad.append(k)
    ok = str(hook["env"]) == str(base["env"])
    print(f"{'PASS' if ok else 'FAIL'} env byte-identical")
    ok or bad.append("env")
    for k in ("bulk", "carriers"):
        ok = hook[k].shape == base[k].shape and hook[k].dtype == base[k].dtype and bool(np.isfinite(hook[k]).all())
        print(f"{'PASS' if ok else 'FAIL'} {k} shape {hook[k].shape} dtype {hook[k].dtype}, finite")
        ok or bad.append(k)
    if ctx_p is not None:
        z = np.load(ctx_p, allow_pickle=False)
        man = json.loads((bundle / "manifest.json").read_text())
        ctxs = [str(c) for c in man["contexts"]]
        okc = "contexts" in z.files and sorted(map(str, z["contexts"])) == sorted(ctxs)
        print(f"{'PASS' if okc else 'FAIL'} bulk_ctx contexts == the bundle's {ctxs}")
        okc or bad.append("ctx contexts")
        for k in ("perts", "genes"):
            ok = k in z.files and z[k].shape == base[k].shape and bool(np.all(z[k] == base[k]))
            print(f"{'PASS' if ok else 'FAIL'} bulk_ctx {k} identical to the dump's")
            ok or bad.append(f"ctx {k}")
        for c in ctxs:
            key = f"bulk_{c}"
            ok = key in z.files and z[key].shape == base["bulk"].shape and z[key].dtype == np.float64 \
                and bool(np.isfinite(z[key]).all())
            print(f"{'PASS' if ok else 'FAIL'} bulk_ctx {key} shape/dtype/finite")
            ok or bad.append(key)
    if not bad:
        b0, b1 = base["bulk"], hook["bulk"]
        n = np.linalg.norm(b0, axis=1) * np.linalg.norm(b1, axis=1)
        cos = np.where(n > 0, (b0 * b1).sum(1) / np.where(n > 0, n, 1), 1.0)
        print(f"moved: shared bulk mean cosine to the stack dump {cos.mean():.4f} (min {cos.min():.4f})")
        if ctx_p is not None:
            for c in ctxs:
                bc = z[f"bulk_{c}"]
                nc = np.linalg.norm(b0, axis=1) * np.linalg.norm(bc, axis=1)
                cc = np.where(nc > 0, (b0 * bc).sum(1) / np.where(nc > 0, nc, 1), 1.0)
                print(f"moved: bulk_{c} mean cosine to the stack dump {cc.mean():.4f}")
    print("dumphook check " + ("PASS" if not bad else f"FAIL {bad}"))
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
