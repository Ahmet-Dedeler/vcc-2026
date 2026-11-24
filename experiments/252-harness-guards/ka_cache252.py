"""252 item #6 -- the fold cache, keyed and validated, with one REAL rebuild.

  A. `fold_cache_status` on this experiment's byte copies of 162's four caches (no key):
     every one must validate as `legacy` with the six recomputed checks listed.
  B. the same caches under changed inputs -- sources, op, PM_BLEND, MIN_CPM -- must report
     `miss` with the disagreement named (no rebuild is attempted here; `probes252.py` shows
     `build_fold` itself refusing them via the rebuild sentinel).
  C. a real rebuild of the k562_gwps fold into an EMPTY cache dir (`work/rebuild/`) through
     the unchanged build path, `limit_gib=10`, then: (i) every array bit-identical to 162's
     cache; (ii) the new file carries `cache_key` + `inputs`; (iii) a second `build_fold` is a
     keyed `hit`; (iv) with the source list changed the keyed cache reports `miss` naming the
     differing input.

    /usr/bin/time -l uv run python ka_cache252.py     -> work/ka_cache.json ; exit 1 on FAIL
"""
from __future__ import annotations

import json
import resource
import shutil
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(HERE))

import plolo  # noqa: E402
from gate252 import wait_for_memory  # noqa: E402

ORIG = ROOT / "experiments/162-retrodicting-validator/work"
ARRAYS = ("pre", "truth", "b_ctrl", "gate", "scored", "wcols", "targets")


def main() -> None:
    t0 = time.time()
    out: dict = dict(A={}, B={}, C={}, failures=[])

    def fail(msg):
        print(f"  *** {msg}", flush=True)
        out["failures"].append(msg)

    # ---- A. legacy validation on the byte copies
    plolo.WORK = WORK / "plolo"
    for held in plolo.FOLDS:
        s = plolo.fold_cache_status(held)
        out["A"][held] = s
        print(f"A {held}: {s['status']} validated={s['validated']} key={s['key'][:12]}", flush=True)
        if s["status"] != "legacy":
            fail(f"A: {held} did not validate as legacy: {s}")

    # ---- B. changed inputs must miss, naming the disagreement
    held = "k562_gwps"
    saved = dict(plolo.FOLDS[held])
    saved_pm, saved_cpm = plolo.PM_BLEND, plolo.MIN_CPM
    cases = {}
    try:
        plolo.FOLDS[held]["sources"] = ("orion_hek293t",)
        cases["sources -> ('orion_hek293t',)"] = plolo.fold_cache_status(held)
        plolo.FOLDS[held].update(saved)
        plolo.FOLDS[held]["op"] = plolo.HELD_OPS / "held_orion_hek293t.npz"
        cases["op -> held_orion_hek293t.npz"] = plolo.fold_cache_status(held)
        plolo.FOLDS[held].update(saved)
        plolo.PM_BLEND = 0.25
        cases["PM_BLEND -> 0.25"] = plolo.fold_cache_status(held)
        plolo.PM_BLEND = saved_pm
        # MIN_CPM 10 -> 5 leaves this fold's gate IDENTICAL (every gene the k562 bank carries is
        # above 5 CPM), so a legacy cache is rightly still valid there; the miss is demonstrated
        # at the smallest threshold that actually changes the gate (20, measured by probes252).
        plolo.MIN_CPM = 5.0
        cases["MIN_CPM -> 5.0 (gate unchanged: still valid)"] = plolo.fold_cache_status(held)
        plolo.MIN_CPM = 20.0
        cases["MIN_CPM -> 20.0"] = plolo.fold_cache_status(held)
        plolo.MIN_CPM = saved_cpm
        cases["restored"] = plolo.fold_cache_status(held)
    finally:
        plolo.FOLDS[held].update(saved)
        plolo.PM_BLEND, plolo.MIN_CPM = saved_pm, saved_cpm
    for k, s in cases.items():
        out["B"][k] = dict(status=s["status"], reason=s["reason"], key=s["key"][:12])
        print(f"B {k}: {s['status']} -- {s['reason']}", flush=True)
        want = "legacy" if (k == "restored" or "still valid" in k) else "miss"
        if s["status"] != want:
            fail(f"B: {k}: expected {want}, got {s['status']}")

    # ---- C. one real rebuild, gated, into an empty dir
    print(f"C gate: {wait_for_memory(tag=' rebuild')}", flush=True)
    rb = WORK / "rebuild"
    if rb.exists():
        shutil.rmtree(rb)
    rb.mkdir()
    for c in "ABC":                                   # bctrl caches: read-only copies
        shutil.copyfile(ORIG / f"bctrl-challenge-{c}.npy", rb / f"bctrl-challenge-{c}.npy")
    plolo.WORK = rb
    t1 = time.time()
    fold = plolo.build_fold(held, limit_gib=10.0)
    out["C"]["rebuild_seconds"] = round(time.time() - t1, 1)
    out["C"]["rebuild_audit_cache"] = fold["audit"]["cache"]
    z_new = np.load(rb / f"fold-{held}.npz", allow_pickle=True)
    z_old = np.load(ORIG / f"fold-{held}.npz", allow_pickle=True)
    ident = {}
    for k in ARRAYS:
        x, y = z_old[k], z_new[k]
        ident[k] = bool(x.shape == y.shape and x.dtype == y.dtype and np.array_equal(x, y))
        if not ident[k]:
            fail(f"C: rebuilt {k} differs from 162's cache")
    out["C"]["arrays_bit_identical_to_162"] = ident
    out["C"]["new_cache_has_key"] = "cache_key" in z_new.files and "inputs" in z_new.files
    a_old, a_new = json.loads(str(z_old["audit"])), json.loads(str(z_new["audit"]))
    out["C"]["audit_fields_identical"] = {k: a_old[k] == a_new.get(k) for k in a_old}
    print(f"C rebuilt in {out['C']['rebuild_seconds']}s; arrays identical {ident}; "
          f"audit identical {out['C']['audit_fields_identical']}", flush=True)
    s = plolo.fold_cache_status(held)
    out["C"]["second_call"] = dict(status=s["status"], key=s["key"][:12])
    if s["status"] != "hit":
        fail(f"C: keyed cache did not hit on the second call: {s}")
    fold2 = plolo.build_fold(held)
    out["C"]["second_build_fold_cache"] = fold2["audit"]["cache"]
    try:
        plolo.FOLDS[held]["sources"] = ("orion_hek293t",)
        s = plolo.fold_cache_status(held)
    finally:
        plolo.FOLDS[held].update(saved)
    out["C"]["keyed_cache_sources_changed"] = dict(status=s["status"], reason=s["reason"])
    print(f"C keyed cache, sources changed: {s['status']} -- {s['reason']}", flush=True)
    if s["status"] != "miss" or "sources" not in str(s["reason"]):
        fail(f"C: keyed cache did not miss on changed sources: {s}")

    out["seconds"] = round(time.time() - t0, 1)
    out["peak_rss_gib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3
    (WORK / "ka_cache.json").write_text(json.dumps(out, indent=1, default=str))
    print(f"KA-cache: {len(out['failures'])} failures, {out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB "
          f"-> {'PASS' if not out['failures'] else 'FAIL'}")
    if out["failures"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
