"""417 builder window 2b: bulk-matrix I/O so any BULK design computed elsewhere can ride the ballast.

VCC_BULK_DUMP=<path.npz>  writes the finished bulk matrix (after VCC_PERT_SCALE / VCC_BULK_GENE_GAIN, i.e. exactly what
                          the writer would render for pds/mse), the carrier matrix, perts, genes and the VCC_* env, prints
                          its sha256, and exits 0 before emission (a dump run makes no h5ad).
VCC_BULK_LOAD=<path.npz>:<sha256>  (requires VCC_PERT_SCALE_BALLAST) replaces the bulk matrix with the pinned file's
                          "bulk" after checking sha, perts and genes; the carriers stay this build's. The config records
                          path, sha and the dump's env; the audit records max |loaded - own bulk| (0 for a same-env dump).
Unset is a no-op. Apply only when no build has the builder pinned. Apply AFTER patch_next417.py.
usage: uv run python patch_bulkio417.py [--write]      (PATCH_TARGET=<path> to dry-run on a copy)
"""
import difflib
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
P = Path(os.environ.get("PATCH_TARGET", ROOT / "scripts/make_submission_pm_cd4.py"))
src = P.read_text()
out = src


def sub(old, new, label):
    global out
    n = out.count(old)
    if n != 1:
        raise SystemExit(f"[{label}] expected exactly one match, found {n}")
    out = out.replace(old, new)


sub('''if BULK_CTX_SCALE is not None and PERT_SCALE_BALLAST is None:
    raise SystemExit("VCC_BULK_CTX_SCALE is a bulk-only arm and needs VCC_PERT_SCALE_BALLAST (welded it reaches the DE channel)")
''', '''if BULK_CTX_SCALE is not None and PERT_SCALE_BALLAST is None:
    raise SystemExit("VCC_BULK_CTX_SCALE is a bulk-only arm and needs VCC_PERT_SCALE_BALLAST (welded it reaches the DE channel)")
# `VCC_BULK_DUMP` / `VCC_BULK_LOAD` (417): the bulk matrix as a pinned file, so a bulk design built with a different
# composition (418) rides the ballast while the carriers keep this build's design. DUMP writes and exits before emission;
# LOAD="<path.npz>:<sha256>" swaps the bulk after every effect-level bulk transform. Unset is a no-op.
BULK_DUMP = os.environ.get("VCC_BULK_DUMP", "").strip() or None
BULK_LOAD = os.environ.get("VCC_BULK_LOAD", "").strip() or None
if BULK_LOAD is not None:
    if ":" not in BULK_LOAD:
        raise SystemExit(f"VCC_BULK_LOAD must be '<path.npz>:<sha256>': {BULK_LOAD!r}")
    _bl_path, _bl_sha = BULK_LOAD.rsplit(":", 1)
    BULK_LOAD = (Path(_bl_path), _bl_sha.strip().lower())
    if len(BULK_LOAD[1]) != 64:
        raise SystemExit(f"VCC_BULK_LOAD sha must be a full sha256: {BULK_LOAD[1]!r}")
    if PERT_SCALE_BALLAST is None:
        raise SystemExit("VCC_BULK_LOAD replaces the ballast bulk and needs VCC_PERT_SCALE_BALLAST")
if BULK_DUMP is not None and BULK_LOAD is not None:
    raise SystemExit("VCC_BULK_DUMP and VCC_BULK_LOAD are exclusive")
''', "knobs")

sub('''    del donors, pm, cd4_banks, cd4, h1_union, union2, operator, model
    gc.collect()
    memory_checkpoint("before challenge emission", args.max_rss_gib)
''', '''    # --- `VCC_BULK_DUMP` / `VCC_BULK_LOAD` (417): the bulk matrix as a pinned file.
    if BULK_DUMP is not None:
        carriers = effects_unscaled if (PERT_SCALE is not None and PERT_SCALE_BALLAST is not None) else effects
        np.savez_compressed(BULK_DUMP, bulk=np.asarray(effects, dtype=np.float64), carriers=np.asarray(carriers, dtype=np.float64),
                            perts=np.asarray(list(map(str, spec.perts))), genes=np.asarray(list(map(str, spec.genes))),
                            env=np.asarray(json.dumps({k: v for k, v in sorted(os.environ.items()) if k.startswith("VCC_")})))
        print(f"[bulk-dump] wrote {BULK_DUMP} sha256={file_sha256(Path(BULK_DUMP))} -- exiting before emission", flush=True)
        raise SystemExit(0)
    if BULK_LOAD is not None:
        if PERT_SCALE is None or context_effects is not None:
            raise SystemExit("VCC_BULK_LOAD needs the shared-matrix ballast path (VCC_PERT_SCALE + VCC_PERT_SCALE_BALLAST)")
        got = file_sha256(BULK_LOAD[0])
        if got != BULK_LOAD[1]:
            raise SystemExit(f"VCC_BULK_LOAD sha mismatch: {got} != {BULK_LOAD[1]}")
        with np.load(BULK_LOAD[0], allow_pickle=False) as z:
            if list(map(str, z["perts"])) != list(map(str, spec.perts)) or list(map(str, z["genes"])) != list(map(str, spec.genes)):
                raise SystemExit("VCC_BULK_LOAD perts/genes differ from the bundle's")
            loaded = np.asarray(z["bulk"], dtype=np.float64)
            dump_env = json.loads(str(z["env"]))
        if loaded.shape != np.asarray(effects).shape or not np.isfinite(loaded).all():
            raise SystemExit(f"VCC_BULK_LOAD matrix shape/finiteness: {loaded.shape}")
        effect_audit["bulk_load"] = {"path": str(BULK_LOAD[0]), "sha256": got, "dump_env": dump_env,
                                     "max_abs_diff_vs_own_bulk": float(np.max(np.abs(loaded - np.asarray(effects, dtype=np.float64))))}
        print(f"[bulk-load] {BULK_LOAD[0]} sha={got[:16]} max|loaded-own|={effect_audit['bulk_load']['max_abs_diff_vs_own_bulk']:.3e}", flush=True)
        effects = loaded

    del donors, pm, cd4_banks, cd4, h1_union, union2, operator, model
    gc.collect()
    memory_checkpoint("before challenge emission", args.max_rss_gib)
''', "io")

sub('''            **({} if BULK_CTX_SCALE is None else {"bulk_ctx_scale": BULK_CTX_SCALE}),
''', '''            **({} if BULK_CTX_SCALE is None else {"bulk_ctx_scale": BULK_CTX_SCALE}),
            **({} if BULK_LOAD is None else {"bulk_load": str(BULK_LOAD[0].name), "bulk_load_sha256": BULK_LOAD[1]}),
''', "config")

sys.stdout.writelines(difflib.unified_diff(src.splitlines(True), out.splitlines(True), "a/" + P.name, "b/" + P.name, n=1))
if "--write" in sys.argv:
    P.write_text(out)
    print("\nWROTE", P)
