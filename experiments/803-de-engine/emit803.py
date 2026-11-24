#!/usr/bin/env python3
"""803 — launch ONE step of the real-emitter DE harness, importing the existing drivers verbatim (never edited) and
redirecting every output into an 803 run directory. One process per step, so each step's peak RSS is its own.

  emit <surface> <hook> <callset_cpm> <outdir> <run496 args ...>
      surface h1test: hook none   -> 618/run618e.py mode h1ka (no swaps: 435's harness = run496 + the production emitter)
                      hook c4h    -> 612/run612.py mode c4kgall (the plain CLEAN4 hook, tau 0.05, = 537's cpm+c4h)
                      hook idall3 -> 612/run612.py mode idall3 (the ungated v3 identity hook)
      surface h1tv:   hook none / c4h / idall3 / cp -> 618/run618e.py mode tv / tvc4h / tvidall3 / tvcp (618's faithful
                      tv key, ref_tv, the four CLEAN4 banks; cp = 528's call-propensity block, 618's G3)
      callset_cpm: 2200 = the shipped CallSet(400, 2200, 1.9); any other value patches run400.CALLSET's CPM-block size,
                      the harness twin of the builder's VCC_RC_CALLSET_CPM (523).
  rnup <surface> <outdir> <rnup630 args ...>
      630/rnup630.py (the UP pass, optionally + 631's DN pass, on run496's dumped pre-pass blocks); for h1tv its reference
      dir is ref_tv (as 618/rnup618.py does); outputs into <outdir>.

Every step writes <outdir>/manifest-<step>.json at exit: argv, python / platform, the git blob sha1 of every imported
module under the repo, every file the process opened under the repo (sha256 for files < 256 MB; size + mtime above:
the big inputs are sha-checked once per instance by cloud/setup803.sh), wall seconds and peak RSS.
"""
from __future__ import annotations

import atexit
import hashlib
import json
import os
import platform
import resource
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
T0 = time.time()
ORIG_ARGV = list(sys.argv)
STEP = sys.argv[1]
assert STEP in ("emit", "rnup"), STEP
OPENED: dict[str, str] = {}


def _audit(event, args):
    if event == "open" and args and isinstance(args[0], (str, bytes, os.PathLike)):
        try:
            p = os.path.abspath(os.fsdecode(args[0]))
        except Exception:  # noqa: BLE001
            return
        mode = args[1] if len(args) > 1 and isinstance(args[1], str) else "r"
        if p.startswith(str(ROOT)) and p not in OPENED:
            OPENED[p] = mode


sys.addaudithook(_audit)


def git_blob(path: Path) -> str:
    b = path.read_bytes()
    return hashlib.sha1(b"blob %d\0" % len(b) + b).hexdigest()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(outdir: Path, extra: dict):
    mods = {}
    for name, m in list(sys.modules.items()):
        f = getattr(m, "__file__", None)
        if f and str(Path(f).resolve()).startswith(str(ROOT)) and f.endswith(".py"):
            p = Path(f).resolve()
            mods[str(p.relative_to(ROOT))] = git_blob(p)
    files = {}
    for p, mode in sorted(OPENED.items()):
        pp = Path(p)
        if not pp.is_file() or str(pp).startswith(str(outdir.resolve())) or pp.suffix in (".py", ".pyc"):
            continue
        if any(c in mode for c in "wax+"):
            continue
        st = pp.stat()
        rel = str(pp.relative_to(ROOT))
        files[rel] = ({"sha256": sha256(pp), "size": st.st_size} if st.st_size < (256 << 20)
                      else {"size": st.st_size, "mtime": st.st_mtime, "sha256": "checked-by-setup803"})
    rec = dict(step=STEP, argv=ORIG_ARGV, driver_argv=sys.argv, python=sys.version, platform=platform.platform(), wall_s=round(time.time() - T0, 1),
               peak_rss_kb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, modules=mods, inputs=files,
               env={k: v for k, v in os.environ.items() if k.startswith(("NUMBA", "OMP", "OPENBLAS", "POLARS", "RUN4", "VCC"))},
               **extra)
    (outdir / f"manifest-{STEP}.json").write_text(json.dumps(rec, indent=1, default=str))


def install_bulk_transform(M, spec: dict, surface: str, extra: dict):
    """§8: an arm's `bulk_transform` = {module (repo-relative path), sha256, args}. `module.transform(eff, genes, targets,
    surface, **args) -> (rows, info)` is applied to the surface's composition rows (the cap15 matrix run496 passes to
    ps_bulk), i.e. BEFORE the key scaling, the sharpening and the ballast render, where production applies a stack-dump
    transform. The DE carriers (run496's eff_de_native) are untouched. The rows must keep the shape and be finite."""
    import importlib.util  # noqa: PLC0415
    import numpy as np  # noqa: PLC0415
    path = ROOT / spec["module"]
    got = sha256(path)
    assert got == spec["sha256"], f"bulk_transform module sha {got} != the arm's {spec['sha256']}"
    ms = importlib.util.spec_from_file_location("bt803_" + path.stem, path)
    mod = importlib.util.module_from_spec(ms)
    ms.loader.exec_module(mod)
    rec = {"module": spec["module"], "sha256": got, "git_blob": git_blob(path), "args": spec.get("args", {})}
    extra["bulk_transform"] = rec
    _ps = M.ps_bulk

    def ps_bulk(eff, genes, targets):
        rows, info = mod.transform(eff, genes, targets, surface, **spec.get("args", {}))
        rows = np.asarray(rows)
        assert rows.shape == np.shape(eff) and np.isfinite(rows).all(), "bulk_transform changed the shape or made non-finite rows"
        rec["info"] = info
        rec["identical_to_input"] = bool(np.array_equal(rows, eff))
        return _ps(rows, genes, targets)

    M.ps_bulk = ps_bulk


def install_sharpen(M, spec: str, extra: dict):
    """Follow-up 1 (bulk sharpening on DE-H1, 599's wiring): the ballast's bulk profile is rendered with a per-target q by the
    builder's own calls. q comes from the live builder's `sharpen_q_by_tercile(bulk, carriers, genes, perts, pt)` with
    bulk = run496.ps_bulk's key-scaled rows and carriers = the cap15 design rows (production: the dump's bulk and carrier
    rows), and the render is `profile_from_effect(eff_bk[i], b_ctrl, min_cpm=10, sharpen_q=q_i, protect=panel)` (the
    builder's profile_bulk call, without its effect_scale / gene_scale / shrink_center / ctx-scale terms, which are off in
    the harness). Before replacing it, the profile run496 passed is asserted equal to the unsharpened one recomputed.
    spec: noop (wrapper installed, q None: must equal the plain run) | scalar:<q> | tercile:<a,b,c> | stercile:<a,b,c>."""
    import numpy as np  # noqa: PLC0415
    B = M.H.builder389()
    st: dict = {}
    rec = {"spec": spec, "builder_sha256": sha256(ROOT / "scripts/make_submission_pm_cd4.py"), "per_target": []}
    extra["sharpen"] = rec
    _ps = M.ps_bulk

    def ps_bulk(eff, genes, targets):
        eb = _ps(eff, genes, targets)
        tg = [str(t) for t in targets]
        st.update(eb=eb, genes=genes, targets=tg)
        if spec == "noop":
            q = None
        elif spec.startswith("scalar:"):
            q = np.full(len(tg), float(spec.split(":", 1)[1]))
        else:
            mode, qs = spec.split(":", 1)
            assert mode in ("tercile", "stercile"), spec
            q, audit = B.sharpen_q_by_tercile(eb, eff, [str(g) for g in genes], tg, {"mode": mode, "q": [float(x) for x in qs.split(",")]})
            q = np.asarray(q, dtype=np.float64)
            rec["audit"] = audit
        st["q"] = q
        rec["q"] = None if q is None else q.tolist()
        return eb

    _load = M.E.load_reference

    def load_reference():
        out = _load()
        st["r"] = out[0]
        return out

    _bb = M.ballast_block

    def ballast_block(profile_de, profile_bulk, libs, emit_rng, spread, key, i, emit_c, kmax, k_force=None):
        r = st["r"]
        if "eff_bk" not in st:
            assert [str(t) for t in r["targets"]] == st["targets"], "reference target order != the bulk's"
            st["eff_bk"] = M.R.map_axis(st["eb"], st["genes"], r["genes"])
            panel = set(st["targets"])
            st["panel"] = np.array([str(g) in panel for g in r["genes"]])
        pb0 = B.profile_from_effect(st["eff_bk"][i], r["b_ctrl"], min_cpm=10.)
        assert np.array_equal(pb0, profile_bulk), f"target {i}: run496's bulk profile != the recomputed unsharpened one"
        if st["q"] is not None:
            profile_bulk = B.profile_from_effect(st["eff_bk"][i], r["b_ctrl"], min_cpm=10., sharpen_q=float(st["q"][i]),
                                                 protect=st["panel"])
        block, dg = _bb(profile_de, profile_bulk, libs, emit_rng, spread, key, i, emit_c, kmax, k_force=k_force)
        rec["per_target"].append({"i": int(i), "target": st["targets"][i],
                                  "q": None if st["q"] is None else float(st["q"][i]), "K": int(dg["ballast_K"])})
        return block, dg

    M.ps_bulk = ps_bulk
    M.E.load_reference = load_reference
    M.ballast_block = ballast_block


def main():
    if STEP == "emit":
        surface, hook, cpm, outdir = sys.argv[2], sys.argv[3], int(sys.argv[4]), Path(sys.argv[5]).resolve()
        args = sys.argv[6:]
        outdir.mkdir(parents=True, exist_ok=True)
        extra = dict(surface=surface, hook=hook, callset_cpm=cpm)
        atexit.register(write_manifest, outdir, extra)
        if surface == "h1tv" or hook == "none":
            mode = {"none": "tv", "c4h": "tvc4h", "idall3": "tvidall3", "cp": "tvcp"}[hook] if surface == "h1tv" else "h1ka"
            d618 = ROOT / "experiments/618-h1-trainval-surface"
            # run618e writes into 618/work/<NAME>: NAME = the relative path to our run dir, so nothing lands in 618
            name = os.path.relpath(outdir, d618 / "work")
            sys.path.insert(0, str(d618))
            sys.argv = [str(d618 / "run618e.py"), mode, name] + args
            import run618e as W  # noqa: PLC0415  (read-only)
            M = W.M
            extra.update(driver="618/run618e.py", mode=mode, driver_rec=W.REC)
        else:
            d612 = ROOT / "experiments/612-key-gated-hook"
            mode = {"c4h": "c4kgall", "idall3": "idall3"}[hook]
            sys.path.insert(0, str(d612))
            sys.argv = [str(d612 / "run612.py"), mode] + args
            import run612 as W  # noqa: PLC0415  (read-only; its __main__ block is not run, so 612's dir gets no json)
            M = W.M
            M.WORK = outdir
            M.H.WORK = outdir
            M.HERE = outdir
            extra.update(driver="612/run612.py", mode=mode,
                         driver_rec={k: v for k, v in W.REC.items() if k != "per_target"})
        assert Path(M.WORK).resolve() == outdir, (M.WORK, outdir)
        if cpm != 2200:
            assert tuple(M.H.CALLSET) == (400, 2200, 1.9), M.H.CALLSET
            M.H.CALLSET = (400, cpm, 1.9)
        extra["callset"] = list(M.H.CALLSET)
        if os.environ.get("DE803_BULK_TRANSFORM", ""):
            install_bulk_transform(M, json.loads(os.environ["DE803_BULK_TRANSFORM"]), surface, extra)
        if os.environ.get("DE803_SHARPEN", ""):
            install_sharpen(M, os.environ["DE803_SHARPEN"], extra)
        M.main()
    else:
        surface, outdir = sys.argv[2], Path(sys.argv[3]).resolve()
        args = sys.argv[4:]
        outdir.mkdir(parents=True, exist_ok=True)
        extra = dict(surface=surface)
        atexit.register(write_manifest, outdir, extra)
        sys.path.insert(0, str(ROOT / "experiments/630-rank-neutral-up"))
        import rnup630 as U  # noqa: PLC0415  (read-only)
        if surface == "h1tv":
            U.R.WORK = ROOT / "experiments/618-h1-trainval-surface/work/ref_tv"
        U.WORK = outdir
        extra.update(driver="630/rnup630.py", refdir=str(U.R.WORK))
        sys.argv = [str(ROOT / "experiments/630-rank-neutral-up/rnup630.py")] + args
        U.main()


if __name__ == "__main__":
    main()
