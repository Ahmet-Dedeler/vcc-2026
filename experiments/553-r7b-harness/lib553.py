#!/usr/bin/env python3
"""553 — the surface library: the H1 DE instrument (323 -> 400 -> 410b -> 416 -> 426 -> 429 -> 435 -> 496 -> 528 -> 537)
generalised from h1-test's frozen inputs to a SURFACE = (control bundle, target list, gene axis, truth path).

Nothing here re-implements the emitter or the scorer. Every statistic is computed by the SAME function the H1 instrument
calls, on the surface's own inputs:

  control side (allowed before any access; RESERVE.md section 6):
    `controlref(...)`     evaluate323.reference()'s control half, verbatim: rc_emit.control_block_from_counts (the emitter's
                          sorted control block), exact_de323.control_block (the scorer's), the control count sums, the
                          control jackknife, libraries, b_ctrl. Any number of controls.
  truth side (the ACCESS only; `access553.py`):
    `truth_reference(...)` evaluate323.reference()'s truth half: E.block_stats per target (sums, jackknife, Wilcoxon p / BH q /
                          lfc against the scorer's control block) and E.de_frame. Any number of truth cells per target.
    `score_stats(...)`    evaluate323.score() with the reference passed in and the prediction's cell counts explicit (400 per
                          target, as a submission is; on h1-test identical to the truth's 400).
    `anchors(...)`        anchors323.run()'s baseline (cell-eval2's generic-response profile, dispersed emission, seed 0) and
                          5-split replicate (cell-eval2's _derive_seeds / _disjoint_halves / _nmae_ref_from_tables), computed
                          group by group so it fits in < 8 GB on h1-test (323's version peaked at 19 GB). Scored with
                          anchors323.fast_score itself (read-only import).
    `scaled(...)`         cell-eval2's own score_one with each member's catalog policy anchored at the measured replicate
                          (score.py's _replicate_entries), i.e. (u - b) / (r - b) with the mse clamp and nmae's -6 floor.

SEALING. The R7b truth file is opened only through `open_truth()`, which refuses it unless the access has been armed by
`access553.py r7b --owner-go <file>` in the same process. Nothing in this module computes a statistic of a perturbed cell
otherwise.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("POLARS_MAX_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
sys.dont_write_bytecode = True

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E323 = ROOT / "experiments/323-current-family-validation"
R7B = ROOT / "experiments/552-r7-reserve/work/r7b"
for _p in (str(ROOT / "src"), str(E323)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

CTRL = "non-targeting"
TS = 50000.0
MEMBERS = ("pds", "mse", "nmae", "fid", "reach", "jac")        # evaluate323.MEMBERS' order
FULL = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
        "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
        "jac": "de_wilcoxon_sig_jaccard"}
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
EMITTER_SHA = "0f56209d57e8385abd08a7a15692ae3241b9373b436b4921546e546d571e5876"
N_PRED_CELLS = 400
# no source of R7b's own line may enter any donor mix, bank, prior, registry or key (mandate item 3)
FORBIDDEN_SOURCES = ("sunshine", "sunshinehein2023", "calu3", "calu-3", "calu_3")

_ACCESS_ARMED = {"r7b": False}


# ---------------------------------------------------------------------------------------------------------------- utils
def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def arr_sha(a) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1, default=_jsonable) + "\n")
    tmp.replace(path)


def _jsonable(x):
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, Path):
        return str(x)
    raise TypeError(type(x))


def rss_gb() -> float:
    import resource
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9     # macOS: bytes


def log(msg: str) -> None:
    print(f"[553 {time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------------------------------------------------------- surfaces
CTRL_HALF_RE = r"^(?P<base>h1r7(?:s\d+)?|r7bd\d+)c(?P<k>\d)(?P<h>[AB])$"


def surface_spec(name: str) -> dict:
    """The surface's inputs (paths only). `sealed` surfaces never read their truth outside the armed access."""
    import re
    m = re.fullmatch(CTRL_HALF_RE, name)
    if m:
        # a CONTROL HALF of a surface (`prep553.py chalf <base> <k>`): half of the base's controls (split k, half A or
        # B), the prediction side rebuilt from them and emitted again, scored against the base's target cells with the
        # same half of the controls -- the control re-draw that the within-surface band cannot see (RESULTS). On a
        # rehearsal base the half file carries the target cells too; on an R7b base it holds controls only and the
        # access keeps the truth's controls whose obs names are in it (the half's control sums are then checked).
        base = m.group("base")
        spec = surface_spec(base)
        d = WORK / name
        spec.update(name=name, ctrl_half=dict(base=base, split=int(m.group("k")), half=m.group("h")),
                    control_h5ad=d / f"{name}.h5ad", controlref=d / "controlref553.npz",
                    controlref_json=d / "controlref553.json")
        if base.startswith("h1r7"):
            spec.update(truth_h5ad=d / f"{name}.h5ad")
        else:
            spec.update(ctrl_names_json=d / f"{name}.json")
        return spec
    if name == "h1":
        return dict(
            name="h1", sealed=False,
            control_h5ad=ROOT / "data/harness/h1-test.h5ad",
            truth_h5ad=ROOT / "data/harness/h1-test.h5ad",
            # 323's arrays: the control half of evaluate323.reference(); controlref() reproduces them bit for bit (KA-C)
            controlref=E323 / "work/reference.npz",
            controlref_json=E323 / "work/reference.json",
            effects=E323 / "work/effects.npz", compose=E323 / "work/compose.json", effects_key="cap15",
            operator=E323 / "work/whitener_h1_pruned.npz", op_columns=E323 / "work/stage.npz",
            key="keys528",
        )
    if name == "h1cb":
        # KA-CB2: h1's inputs through R7b's code path -- the sealed scoring stub and the October key form (528's
        # callblock formula on h1-test's own controls, work/h1cb/callblock_h1cb.npz); `lcap-kept` must equal 528's run
        spec = surface_spec("h1")
        spec.update(name="h1cb", sealed=True, key="callblock528")
        return spec
    if name.startswith("h1r7"):
        # the R7b-SHAPED dress rehearsal on h1-test (not a reserve): a subset of its controls and R7b's cells-per-target
        # distribution (metadata), built by `prep553.py h1r7 build` (h1r7, seed 55307) or `prep553.py h1r7 <name> --seed
        # <s>` (further draws, e.g. h1r7s2); carriers / operator / bulk key are h1's
        spec = surface_spec("h1")
        d = WORK / name
        spec.update(name=name, sealed=False, control_h5ad=d / f"{name}.h5ad", truth_h5ad=d / f"{name}.h5ad",
                    controlref=d / "controlref553.npz", controlref_json=d / "controlref553.json", key="callblock528")
        return spec
    if name.startswith("h1t7"):
        # ablation of the rehearsal: h1's prediction side unchanged (all 18,400 controls, 323's reference) against a truth
        # with a rehearsal draw's target cells (h1t7 <- h1r7, h1t7s2 <- h1r7s2, ...) and ALL of h1-test's controls
        # (`prep553.py h1t7 build` / `prep553.py h1t7 h1r7s2`): h1 -> h1t7 isolates the truth's cell counts; h1t7 -> h1r7
        # the control count (prediction side and the truth's DE controls)
        spec = surface_spec("h1")
        spec.update(name=name, truth_h5ad=WORK / name / f"{name}.h5ad")
        return spec
    if name == "fx":
        # test-only synthetic surface (test553.py builds it): the access script end to end in seconds
        d = WORK / "test553/fx"
        return dict(name="fx", sealed=False, control_h5ad=d / "fx.h5ad", truth_h5ad=d / "fx.h5ad",
                    controlref=d / "controlref553.npz", controlref_json=d / "controlref553.json")
    if name.startswith("r7bd"):
        # R7b with a human-library depth floor (a CANDIDATE amendment, not the reserve as built), PREDICTION SIDE ONLY: the
        # floored controls (work/<name>/context_R_<floor>.h5ad, `prep553.py subset <name>`); carriers, operator and bulk
        # key are R7b's (same panel). The truth is never floored here (truth_rows refuses): if the owner adopts the floor,
        # 552 rebuilds the surface with it and a plain spec points at the new files.
        floor = int(name[len("r7bd"):])
        spec = surface_spec("r7b")
        d = WORK / name
        spec.update(name=name, depth_floor=floor, control_h5ad=d / f"context_R_{floor}.h5ad",
                    controlref=d / "controlref553.npz", controlref_json=d / "controlref553.json")
        # if 552 adopts the floor, it writes the floored surface here (its own build record with sha256s); the access then
        # reads that truth as an ordinary sealed surface (and still checks its controls == the floored controls above)
        fl = ROOT / f"experiments/552-r7-reserve/work/r7b_floor{floor}"
        if (fl / "truth_r7b.h5ad").exists() and (fl / "build_r7b.json").exists():
            spec.update(truth_h5ad=fl / "truth_r7b.h5ad", floored_build_json=fl / "build_r7b.json")
        return spec
    if name == "r7b":
        return dict(
            name="r7b", sealed=True,
            control_h5ad=R7B / "bundle/context_R.h5ad",
            truth_h5ad=R7B / "truth_r7b.h5ad",
            build_json=R7B / "build_r7b.json",
            bundle=R7B / "bundle",
            controlref=WORK / "r7b/controlref553.npz",
            controlref_json=WORK / "r7b/controlref553.json",
            effects=WORK / "r7b/effects553.npz", compose=WORK / "r7b/compose553.json", effects_key="cap15",
            operator=WORK / "r7b/whitener_r7b_pruned.npz", op_columns=WORK / "r7b/stage553.npz",
            cd4=WORK / "r7b/cd4_compact_r7b",
            key="callblock528",
        )
    raise ValueError(f"unknown surface {name!r}")


def assert_clean_sources(paths) -> None:
    """Refuse any input whose path names R7b's own line or study (sunshine / Calu-3)."""
    bad = [str(p) for p in paths if any(f in str(p).lower() for f in FORBIDDEN_SOURCES)]
    # the surface's OWN control bundle lives under 552-r7-reserve/work/r7b (its path does not name the study)
    if bad:
        raise SystemExit(f"a forbidden source entered an input list: {bad}")


def r7b_build_record() -> dict:
    spec = surface_spec("r7b")
    return json.loads(Path(spec["build_json"]).read_text())


def check_r7b_outputs(check_truth_sha: bool = True) -> dict:
    """The R7b bundle's files against build_r7b.json's sha256 record. The truth file is HASHED only (allowed by the mandate),
    never parsed."""
    rec = r7b_build_record()
    out = {}
    for rel, want in rec["outputs"].items():
        if rel == "truth_r7b.h5ad" and not check_truth_sha:
            continue
        got = sha256(R7B / rel)
        out[rel] = dict(sha256=got, ok=got == want)
        if got != want:
            raise SystemExit(f"R7b {rel}: sha256 {got} != build_r7b.json {want}")
    return out


def open_truth(spec: dict):
    """The ONLY way this directory's code opens a truth file. R7b's is refused unless the access was armed in-process."""
    if spec["sealed"] and not _ACCESS_ARMED.get(spec["name"], False):
        raise PermissionError(f"{spec['name']}'s truth is SEALED (RESERVE.md sections 3 and 6): only `access553.py "
                              f"{spec['name']} --owner-go <file>` may open it")
    return Path(spec["truth_h5ad"])


def arm_access(name: str) -> None:
    _ACCESS_ARMED[name] = True


# ---------------------------------------------------------------------------------------------------------------- h5ad reading
def _labels(h, path):
    """evaluate323.labels (categorical / nullable string / plain), verbatim semantics."""
    import h5py
    obj = h[path]
    if isinstance(obj, h5py.Group):
        if "categories" in obj:
            cats = obj["categories"].asstr()[:]
            codes = obj["codes"][:]
            if np.any(codes < 0):
                raise ValueError(f"missing categorical label at {path}")
            return cats[codes].astype(str)
        if "values" in obj:
            if "mask" in obj and np.any(obj["mask"][:]):
                raise ValueError(f"missing nullable-string label at {path}")
            return obj["values"].asstr()[:].astype(str)
        raise ValueError(f"unsupported string encoding at {path}")
    return obj.asstr()[:].astype(str)


def read_labels_genes(path: Path):
    import h5py
    with h5py.File(path, "r") as h:
        lab = _labels(h, "obs/target_gene")
        key = h["var"].attrs.get("_index", "_index")
        genes = _labels(h, f"var/{key}")
    return lab, genes


def read_rows(path: Path, rows: np.ndarray | None, chunk_rows: int = 4096):
    """CSR rows `rows` (strictly ascending) of an h5ad's X, as a scipy CSR with the file's dtypes; rows None = every row.
    Memory: the output is preallocated from indptr and filled through two fixed chunk buffers (h5py read_direct), so only
    the selected rows plus ~2 chunks are ever resident (h1-test's 18,400 controls of 52,400: 0.9 GB, not the 2.7 GB file)."""
    import h5py
    import scipy.sparse as sp
    with h5py.File(path, "r") as h:
        X = h["X"]
        if X.attrs.get("encoding-type") != "csr_matrix":
            raise ValueError(f"{path}: X is not CSR")
        n, g = (int(x) for x in X.attrs["shape"])
        ipd = X["indptr"][:].astype(np.int64)
        if rows is None:
            data = X["data"][:]
            ind = X["indices"][:]
            return sp.csr_matrix((data, ind, ipd), shape=(n, g))
        rows = np.asarray(rows, dtype=np.int64)
        if rows.size and (np.any(np.diff(rows) <= 0) or rows[0] < 0 or rows[-1] >= n):
            raise ValueError("rows must be strictly ascending and in range")
        lens = ipd[rows + 1] - ipd[rows]
        indptr = np.zeros(rows.size + 1, dtype=np.int64)
        indptr[1:] = np.cumsum(lens)
        data = np.empty(int(indptr[-1]), dtype=X["data"].dtype)
        ind = np.empty(int(indptr[-1]), dtype=X["indices"].dtype)
        mx = max((int(ipd[min(s + chunk_rows, n)] - ipd[s]) for s in range(0, n, chunk_rows)), default=0)
        bd = np.empty(mx, dtype=X["data"].dtype)
        bi = np.empty(mx, dtype=X["indices"].dtype)
        k = 0
        for s in range(0, n, chunk_rows):
            e = min(n, s + chunk_rows)
            j0, j1 = np.searchsorted(rows, [s, e])
            if j0 == j1:
                continue
            a, b = int(ipd[s]), int(ipd[e])
            X["data"].read_direct(bd, np.s_[a:b], np.s_[0:b - a])
            X["indices"].read_direct(bi, np.s_[a:b], np.s_[0:b - a])
            for r in rows[j0:j1]:
                lo, hi = int(ipd[r]) - a, int(ipd[r + 1]) - a
                data[indptr[k]:indptr[k + 1]] = bd[lo:hi]
                ind[indptr[k]:indptr[k + 1]] = bi[lo:hi]
                k += 1
        assert k == rows.size
        return sp.csr_matrix((data, ind, indptr), shape=(rows.size, g))


def int_colsums(x, chunk: int = 2048) -> np.ndarray:
    """np.asarray(x.astype(np.float64).sum(0)).ravel().astype(np.int64) without the full float64 copy: the counts are
    integers, so float64 row-chunk partial sums are exact (< 2^53) and the int64 result is identical in any order."""
    out = np.zeros(x.shape[1], dtype=np.float64)
    for lo in range(0, x.shape[0], chunk):
        out += np.asarray(x[lo:lo + chunk].astype(np.float64).sum(0)).ravel()
    if not np.all(out == np.floor(out)) or out.max(initial=0.0) >= 2 ** 53:
        raise ValueError("non-integer or overflowing count sums")
    return out.astype(np.int64)


def control_block_keep(xc, keep, min_cpm: float = 5.0, chunk: int = 1000):
    """exact_de323.control_block with the gene universe FIXED to `keep`: the same normalized_counts / log1p / sort /
    ctrl_tiesums arithmetic on a subset of the controls (a control half), scored on the full reference's gated genes.
    With keep = control_block(xc).keep it is control_block(xc) exactly (test553 T8)."""
    import exact_de323
    from vcc2026 import rc_emit
    n, g = xc.shape
    keep = np.asarray(keep, dtype=bool)
    mean = np.zeros(g, dtype=np.float64)
    for lo in range(0, n, chunk):
        mean += exact_de323.normalized_counts(xc[lo:lo + chunk]).sum(0)
    mean /= n
    cs = np.empty((int(keep.sum()), n), dtype=np.float64)
    for lo in range(0, n, chunk):
        hi = min(lo + chunk, n)
        cs[:, lo:hi] = np.log1p(exact_de323.normalized_counts(xc[lo:hi])[:, keep]).T
    cs.sort(axis=1)
    return rc_emit.RcControlBlock(keep=keep, cs=cs, tie=rc_emit.ctrl_tiesums(cs), ref_lin=mean[keep], n_controls=n,
                                  min_cpm=min_cpm)


# ---------------------------------------------------------------------------------------------------------------- control side
def controlref(control_h5ad: Path, *, emission: bool = True, check=None):
    """evaluate323.reference()'s control half on any control file (any number of controls). Returns a dict with the
    reference.npz key names (keep, cs, tie, ref_lin, emit_cs, emit_tie, emit_ref_lin, libraries, b_ctrl) plus
    ctrl_sums / ctrl_jk / n_controls / genes. The file's rows labelled `non-targeting` are the controls."""
    import exact_de323
    from cell_eval2.moments import jackknife_correction
    from vcc2026 import rc_emit
    lab, genes = read_labels_genes(control_h5ad)
    rows = np.flatnonzero(lab == CTRL)
    xc = read_rows(control_h5ad, rows)
    log(f"controls {xc.shape[0]} x {xc.shape[1]} ({xc.nnz:,} nnz, {xc.dtype}) from {control_h5ad.name}")
    out = {"genes": genes}
    if emission:
        ec = rc_emit.control_block_from_counts(xc)
        out.update(emit_cs=ec.cs, emit_tie=ec.tie, emit_ref_lin=ec.ref_lin, emit_keep=ec.keep)
        del ec
        if check is not None:           # compare, then drop the big block before the next one is built
            check(out, ("emit_cs", "emit_tie", "emit_ref_lin"))
            out.pop("emit_cs")
    c = exact_de323.control_block(xc)
    if emission and not np.array_equal(c.keep, out["emit_keep"]):
        raise AssertionError("emitter and scorer CPM gate differ")
    out.pop("emit_keep", None)
    cp = int_colsums(xc)
    cj = jackknife_correction(xc, np.zeros(xc.shape[0], dtype=np.int64), 1, TS, chunk=32)[0]
    out.update(keep=c.keep, cs=c.cs, tie=c.tie, ref_lin=c.ref_lin,
               libraries=np.asarray(xc.sum(1)).ravel().astype(float),
               b_ctrl=np.log1p(TS * cp / cp.sum()), ctrl_sums=cp, ctrl_jk=np.float64(cj),
               n_controls=np.int64(xc.shape[0]))
    return out


def library_classes(libraries) -> tuple:
    """Production's emitter classes (make_submission_pm_cd4.py `de_sign_classes`): equal thirds at the context's
    control-library q25 / q50 / q75. On h1-test this is run400's CLASSES3 constant exactly (asserted by the harness)."""
    q = np.quantile(np.asarray(libraries, dtype=np.float64), [0.25, 0.5, 0.75])
    return ((1 / 3, float(q[0])), (1 / 3, float(q[1])), (1 / 3, float(q[2])))


H1_CLASS_LEVELS = (16078.75, 20190.5, 25330.5)      # h1-test's control-library q25 / q50 / q75 (run400's CLASSES3)


def flex_classes(libraries) -> tuple:
    """The `flex` class rule: the context's median library x H1's relative class levels (q25/q50, 1, q75/q50). A deviation
    from production's rule for surfaces whose control-library quartiles are far apart; on h1-test it equals production
    (the q75 level differs by one ulp and every class layout is identical)."""
    med = float(np.quantile(np.asarray(libraries, dtype=np.float64), 0.5))
    q25, q50, q75 = H1_CLASS_LEVELS
    return ((1 / 3, med * (q25 / q50)), (1 / 3, med * 1.0), (1 / 3, med * (q75 / q50)))


def load_controlref(spec: dict, *, need_cs: bool = True) -> dict:
    """The control-side arrays of a surface, from 323's reference.npz (H1) or this directory's controlref553.npz (R7b).
    Truth arrays are NEVER read here (on H1, `truth_arrays` reads them separately for the inline KA score)."""
    path = Path(spec["controlref"])
    info = json.loads(Path(spec["controlref_json"]).read_text())
    want = info["statistics_sha256"]
    got = sha256(path)
    if got != want:
        raise ValueError(f"{path}: sha256 {got} != recorded {want}")
    keys = ["genes", "keep", "tie", "ref_lin", "emit_tie", "emit_ref_lin", "libraries", "b_ctrl"]
    if need_cs:
        keys += ["cs", "emit_cs"]
    with np.load(path, allow_pickle=False) as z:
        out = {k: z[k] for k in keys}
        if spec["name"] in ("h1", "h1cb") or spec["name"].startswith("h1t7"):     # 323's reference.npz
            out["ctrl_sums"] = z["sums"][-1]
            out["ctrl_jk"] = z["jk"][-1]
            out["n_controls"] = int(z["counts"][-1])
            out["targets"] = z["targets"]
        else:
            out["ctrl_sums"] = z["ctrl_sums"]
            out["ctrl_jk"] = z["ctrl_jk"]
            out["n_controls"] = int(z["n_controls"])
            out["targets"] = z["targets"]
    return out


def surface_targets(spec: dict) -> np.ndarray:
    """The target list: H1 = the reference's (sorted labels); R7b = the bundle's pert_counts.csv (== build_r7b.json)."""
    if spec["name"] in ("h1", "h1cb", "fx") or spec["name"].startswith(("h1r7", "h1t7")) and not spec.get("ctrl_names_json"):
        with np.load(spec["controlref"], allow_pickle=False) as z:
            return z["targets"].astype(str)
    import pandas as pd
    t = np.asarray(pd.read_csv(Path(spec["bundle"]) / "pert_counts.csv")["target_gene"], dtype=str)
    rec = r7b_build_record()
    if list(t) != list(rec["targets"]) or list(t) != sorted(t):
        raise ValueError("R7b pert_counts.csv differs from build_r7b.json's sorted target list")
    return t


# ---------------------------------------------------------------------------------------------------------------- truth side
def csr_mean0(x) -> np.ndarray:
    """scipy's CSR `.astype(float64).mean(axis=0)` bit for bit (scipy 1.18: every element times 1/n, then the axis-0 sum,
    a CSC mat-vec that adds each row's values into y in row order), streamed row by row with no scaled copy. Checked equal
    to scipy on 3,000 h1-test control rows and a 400-cell target group (RESULTS.md)."""
    n = x.shape[0]
    inv = 1.0 / n
    y = np.zeros(x.shape[1], dtype=np.float64)
    ip, ind, dat = x.indptr, x.indices, x.data
    for r in range(n):
        a, b = ip[r], ip[r + 1]
        y[ind[a:b]] += dat[a:b].astype(np.float64) * inv
    return y


EXCLUDED = "__excluded__"


def truth_rows(spec: dict, *, xt: bool = True, exclude=()) -> dict:
    """The truth held as its control rows (xc) and, if asked, its target rows (xt), both CSR in file order; `pos` maps a
    file row to its local index in xc or xt. xc + xt is the whole matrix (h1-test: 0.94 + 1.8 GB), read by the
    preallocating reader; nothing else is ever materialised whole. Opens the truth only through open_truth()."""
    if spec.get("depth_floor") and not spec.get("floored_build_json"):
        raise PermissionError(f"{spec['name']} is a prediction-side variant: a depth floor on the TRUTH is an amendment for "
                              "552 to build (experiments/552-r7-reserve/work/r7b_floor<floor>/), not something this "
                              "access applies")
    if spec.get("floored_build_json"):
        want = json.loads(Path(spec["floored_build_json"]).read_text())["outputs"]["truth_r7b.h5ad"]
        if sha256(Path(spec["truth_h5ad"])) != want:
            raise SystemExit("the floored truth's sha256 differs from its build record")
    path = open_truth(spec)
    lab, genes = read_labels_genes(path)
    if spec.get("ctrl_names_json"):
        # an R7b-based control half: the truth's controls outside the half are treated as absent (never read)
        import h5py
        keep_names = set(json.loads(Path(spec["ctrl_names_json"]).read_text())["control_obs_names"])
        with h5py.File(path, "r") as h:
            key = h["obs"].attrs.get("_index", "_index")
            names = _labels(h, f"obs/{key}")
        out = (lab == CTRL) & ~np.isin(names, list(keep_names))
        lab = np.where(out, EXCLUDED, lab)
        if int((lab == CTRL).sum()) != len(keep_names):
            raise ValueError("the truth's controls do not contain the half's control cells by name")
    if exclude:
        # excluded targets' cells are treated as absent from the file (their rows are never read)
        lab = np.where(np.isin(lab, list(exclude)), EXCLUDED, lab)
    crow = np.flatnonzero(lab == CTRL)
    trow = np.flatnonzero((lab != CTRL) & (lab != EXCLUDED))
    pos = np.full(lab.size, -1, dtype=np.int64)
    pos[crow] = np.arange(crow.size)
    pos[trow] = np.arange(trow.size)
    return dict(path=path, lab=lab, genes=genes, crow=crow, trow=trow, pos=pos, xc=read_rows(path, crow),
                xt=read_rows(path, trow) if xt else None, exclude=tuple(exclude))


def control_moments(T: dict, *, expect_ctrl_sums=None):
    """evaluate323.reference()'s control row: exact integer column sums, the control count and the control jackknife
    (on the float32 control rows, as evaluate323; the kernel upcasts to float64 itself, so anchors323's float64 call gives
    the same bits). Returns (sums as float64, count, jackknife)."""
    from cell_eval2.moments import jackknife_correction
    xc = T["xc"]
    cp = int_colsums(xc)
    if expect_ctrl_sums is not None and not np.array_equal(cp, expect_ctrl_sums):
        raise AssertionError("the truth's control cells differ from the prediction side's control bundle")
    cj = jackknife_correction(xc, np.zeros(xc.shape[0], dtype=np.int64), 1, TS, chunk=32)[0]
    return cp.astype(np.float64), int(xc.shape[0]), float(cj)


def truth_reference(spec: dict, T: dict, cmom, *, expect_keep=None, expect_targets=None):
    """evaluate323.reference() generalised (any number of controls and of truth cells per target): the scorer's control
    block from the control rows, each target's E.block_stats against it, the control row from `cmom` (control_moments).
    Returns (ref, c, de) with ref holding reference.npz's truth keys (perts, sums [targets + control], counts, jk, p, q,
    lfc, genes, targets, keep)."""
    import evaluate323 as E
    import exact_de323
    t0 = time.monotonic()
    lab, genes, xc = T["lab"], T["genes"], T["xc"]
    targets = np.array(sorted(set(lab) - {CTRL, EXCLUDED}))
    if expect_targets is not None and list(map(str, expect_targets)) != list(map(str, targets)):
        raise ValueError("truth targets differ from the prediction side's target list")
    c = exact_de323.control_block(xc)
    if expect_keep is not None and not np.array_equal(c.keep, expect_keep):
        raise AssertionError("the truth's control gate differs from the prediction side's control reference")
    xt, pos = T["xt"], T["pos"]
    sums, jks, ps, qs, lfcs, counts = [], [], [], [], [], []
    for i, target in enumerate(targets):
        rows = xt[pos[np.flatnonzero(lab == target)]]
        if rows.shape[0] < 2:
            raise ValueError(f"{target}: fewer than 2 truth cells")
        s, j, p, q, l = E.block_stats(rows, c)
        sums.append(s); jks.append(j); ps.append(p); qs.append(q); lfcs.append(l)
        counts.append(rows.shape[0])
        if (i + 1) % 20 == 0:
            log(f"truth reference {i + 1}/{len(targets)} ({time.monotonic() - t0:.1f}s, peak RSS {rss_gb():.2f} GB)")
    cp, n_ctrl, cj = cmom
    ref = dict(genes=genes, targets=targets, perts=np.append(targets, CTRL), sums=np.vstack([*sums, cp.astype(np.int64)]),
               counts=np.array([*counts, n_ctrl]), jk=np.array([*jks, cj]), p=np.asarray(ps), q=np.asarray(qs),
               lfc=np.asarray(lfcs), keep=c.keep)
    de = E.de_frame(targets, genes[c.keep], ps, qs, lfcs)
    return ref, c, de


def score_stats(ref: dict, keep, real_de, sums, jks, p, q, lfc, *, pred_counts=None):
    """evaluate323.score() with the reference passed in: the same cell-eval2 calls in the same order. `pred_counts` is the
    prediction's cells per group (default 400 per target + the real control count); the bulk_lognorm mse reads only the
    jackknife, so on h1-test (400 truth cells) this is evaluate323.score() exactly."""
    import polars as pl
    import evaluate323 as E
    import run323 as R
    from cell_eval2.config import EvalConfig
    from cell_eval2.de import prepare_de
    from cell_eval2.metrics.discrimination import discrimination_score
    from cell_eval2.prep import bulk_lognorm_means
    from cell_eval2.run import dispatch_de_metrics
    T = len(ref["targets"])
    if pred_counts is None:
        pred_counts = np.append(np.full(T, float(N_PRED_CELLS)), float(ref["counts"][-1]))
    pred_sums = np.vstack([sums, ref["sums"][-1]])
    pred_jks = np.append(jks, ref["jk"][-1])
    pm = bulk_lognorm_means(pred_sums.astype(float), TS)
    rm = bulk_lognorm_means(ref["sums"].astype(float), TS)
    pdsv = discrimination_score(
        pred_bulk=(ref["perts"], pm), real_bulk=(ref["perts"], rm),
        genes=ref["genes"], pert_col="target_gene", control=CTRL,
        distance="cosine", rank_denominator="n-1", tie_policy="midrank",
        exclude_target_gene=True, exclusion_scope="panel", control_source="real")
    msekit = R.load_module("msekit553", ROOT / "experiments/096-mse-floor/msekit.py")
    mse, num, den = msekit.score(ref["perts"], pred_sums, np.asarray(pred_counts, dtype=np.float64), pred_jks,
                                 ref["perts"], ref["sums"], np.asarray(ref["counts"], dtype=np.float64), ref["jk"],
                                 ref["genes"])
    pred = E.de_frame(ref["targets"], ref["genes"][keep], p, q, lfc)
    prepared = prepare_de(pred, real_de, control=CTRL, sort_by="abs_log2_fold_change",
                          p_adj_threshold=.05, nan_lfc_policy="mask")
    cfg = EvalConfig.from_preset("vcc2026")
    rows = dispatch_de_metrics([FULL[m] for m in ("nmae", "fid", "reach", "jac")], prepared, cfg)
    rows += [{"perturbation": t, "metric": FULL["pds"], "value": v} for t, v in pdsv.items()]
    raw = {"mse": mse}
    for short in MEMBERS:
        if short != "mse":
            vals = [row["value"] for row in rows if row["metric"] == FULL[short]]
            raw[short] = float(np.nanmean(vals))
    return raw, rows, {"numerator": num, "denominator": den}


# ---------------------------------------------------------------------------------------------------------------- anchors
def _group_stats(get_rows, labels, genes, ctrl_block=None, pre=None):
    """anchors323.stats() one group at a time (no AnnData copies, no float64 copies): `get_rows(label)` returns the
    group's rows IN THE ORDER anchors323 would read them (a boolean mask on the half's / prediction's row order).
    Identical numbers: the control group's sums are exact integer column sums in float64 (anchors323 sums a float64 copy
    of the same integer counts), other groups are summed as anchors323 does; the jackknife kernel and de_from_counts
    upcast to float64 themselves. The DE control block is built from the control
    group's rows exactly as anchors323.stats does, unless it is given (then it must be that same block).
    `pre[CTRL] = (sums, count, jk)` reuses the control group's moments computed on identical rows."""
    import evaluate323 as E
    import exact_de323 as D
    from cell_eval2.moments import jackknife_correction
    perts = np.array(sorted(labels), dtype=str)
    ctrl = D.control_block(get_rows(CTRL)) if ctrl_block is None else ctrl_block
    sums, jk, counts, ps, qs, lfcs, targets = [], [], [], [], [], [], []
    for target in perts:
        if pre is not None and target in pre:
            assert target == CTRL, "only the control group's moments are reused"
            s_, n_, j_ = pre[target]
            sums.append(np.asarray(s_, dtype=np.float64)); counts.append(int(n_)); jk.append(float(j_))
            continue
        rows = get_rows(target)
        # the control group (large, integer truth counts): exact integer column sums, no float64 copy; every other group
        # (small; the baseline's rows are fractional) exactly as anchors323 sums it
        sums.append(int_colsums(rows).astype(np.float64) if target == CTRL
                    else np.asarray(rows.astype(np.float64).sum(0)).ravel())
        counts.append(rows.shape[0])
        jk.append(jackknife_correction(rows, np.zeros(rows.shape[0], dtype=np.int64), 1, TS, chunk=32)[0])
        if target != CTRL:
            _, p, q, lfc = D.de_from_counts(rows, ctrl)
            ps.append(p); qs.append(q); lfcs.append(lfc); targets.append(target)
        del rows
    de = E.de_frame(np.array(targets), genes[ctrl.keep], ps, qs, lfcs)
    return {"perts": perts, "genes": genes, "sums": np.array(sums), "counts": np.array(counts), "jk": np.array(jk),
            "de": de}, ctrl


def _obs_frame(path: Path):
    """The truth's obs as anndata decodes it (categorical target_gene: groupby order = category order)."""
    import h5py
    from anndata.io import read_elem
    with h5py.File(path, "r") as h:
        return read_elem(h["obs"])


def _full_rows(T: dict, groups: dict):
    xc, xt, pos = T["xc"], T["xt"], T["pos"]

    def get(g):
        return xc if g == CTRL else xt[pos[groups[g]]]
    return get


def anchors_base(spec: dict, T: dict, *, ctrl_block=None, ctrl_moments=None):
    """anchors323.run()'s full stats and dispersed seed-0 baseline, one group at a time. Returns (baseline raw, full
    stats, diagnostics).
    Profile: cell-eval2's own generic_response_profile, run on a one-row-per-group AnnData holding each group's exact
    pseudobulk mean (the per-group `X[rows] -> float64 -> mean` of its _grouped_means, reproduced by csr_mean0).
    Emission: build_baseline_prediction's single Generator draws each non-control group's control rows in sorted-label
    order; each group's rows are regenerated on demand (same draws, same float32 scaling, same eliminate_zeros), and the
    prediction's control rows are the truth's own (so its control block and moments are the truth's).
    `ctrl_block` / `ctrl_moments`: truth_reference's control block and control moments (identical by construction)."""
    import anndata as ad
    import pandas as pd
    import anchors323 as A
    from cell_eval2.baseline import _emission_scale, generic_response_profile
    open_truth(spec)
    genes_all = T["genes"]
    t0 = time.monotonic()
    labels = np.asarray(T["lab"]).astype(str)
    uniq = sorted(set(labels.tolist()) - {EXCLUDED})
    groups = {g: np.flatnonzero(labels == g) for g in uniq}
    get = _full_rows(T, groups)
    pre = None if ctrl_moments is None else {CTRL: ctrl_moments}
    full, full_ctrl = _group_stats(get, uniq, genes_all, ctrl_block=ctrl_block, pre=pre)
    log(f"anchors: full stats ({time.monotonic() - t0:.1f}s, peak RSS {rss_gb():.2f} GB)")
    means = np.vstack([csr_mean0(get(g)) for g in uniq])
    pseudo = ad.AnnData(X=means, obs=pd.DataFrame({"target_gene": uniq}, index=[f"g{i}" for i in range(len(uniq))]),
                        var=pd.DataFrame(index=pd.Index(genes_all)))
    prof = generic_response_profile(pseudo, pert_col="target_gene", control=CTRL)
    ctrl_pb = means[uniq.index(CTRL)]
    del pseudo, means
    scale, scale_diag = _emission_scale(prof.values, ctrl_pb, genes_all)
    xc = T["xc"]
    rng = np.random.default_rng(0)
    order_seen = []

    def base_rows(g):
        order_seen.append(g)
        draw = rng.integers(0, xc.shape[0], size=int(groups[g].size))
        out = xc[draw]
        out.data = (out.data * scale[out.indices]).astype(np.float32)
        out.eliminate_zeros()
        return out
    ci = int(np.flatnonzero(full["perts"] == CTRL)[0])
    pre_b = {CTRL: (full["sums"][ci], full["counts"][ci], full["jk"][ci])}
    bstats, _ = _group_stats(base_rows, uniq, genes_all, ctrl_block=full_ctrl, pre=pre_b)
    if order_seen != sorted(set(labels[(labels != CTRL) & (labels != EXCLUDED)].tolist())):
        raise AssertionError("baseline draws were not taken in cell-eval2's sorted non-control order")
    base = A.fast_score(bstats, full)
    del bstats, full_ctrl
    log(f"anchors: baseline {json.dumps(base)} ({time.monotonic() - t0:.1f}s, peak RSS {rss_gb():.2f} GB)")
    return base, full, dict(baseline_emission=scale_diag, seconds_base=time.monotonic() - t0)


def replicate_halves(spec: dict, split_index: int, *, n_splits: int = 5, base_seed: int = 0, exclude=()):
    """The two disjoint halves of anchors323.replicate()'s split `split_index`: _derive_seeds(base_seed, n_splits)[i], then
    ceiling._disjoint_halves' permutations (one Generator over the obs groupby order; floor(n/2) cells per group, the odd
    one dropped). Returns ({label: original rows in the half's order} for A, the same for B)."""
    from cell_eval2.anchor import _derive_seeds
    path = open_truth(spec)
    obs = _obs_frame(path)
    gidx = obs.groupby("target_gene", observed=True).indices
    r_ = np.random.default_rng(_derive_seeds(base_seed, n_splits)[split_index])
    ga, gb = {}, {}
    for pert, idx in gidx.items():
        if str(pert) in exclude:
            continue                     # as if absent from the file: no permutation is drawn for it
        perm = r_.permutation(np.asarray(idx))
        h = perm.size // 2
        if h < 1:
            continue
        ga[str(pert)] = perm[:h]
        gb[str(pert)] = perm[h:2 * h]
    if CTRL not in ga or set(ga[CTRL].tolist()) & set(gb[CTRL].tolist()):
        raise AssertionError("replicate halves share control cells (or the control has < 2 cells)")
    return ga, gb


def replicate_half_stats(spec: dict, groups: dict):
    """anchors323.stats() of one replicate half, reading ONLY that half's rows from the truth file (any size: the half's
    own control block, its own control moments); each group's rows in the half's permuted order."""
    path = open_truth(spec)
    genes_all = read_labels_genes(path)[1]
    rows_all = np.sort(np.concatenate(list(groups.values())))
    x = read_rows(path, rows_all)
    loc = {g: np.searchsorted(rows_all, r) for g, r in groups.items()}
    st, _ = _group_stats(lambda g: x[loc[g]], list(groups), genes_all)
    return st


def replicate_score(sa: dict, sb: dict, full_de):
    """anchors323.replicate()'s per-split score: fast_score(sb, sa, own_controls=True), nmae replaced by
    _nmae_ref_from_tables on the full-reference gate (323's ANCHOR-NOTE.md)."""
    import anchors323 as A
    from cell_eval2.de import resolve_target_genes
    from cell_eval2.lfc_nmae_ref import _nmae_ref_from_tables
    resolution = resolve_target_genes(full_de, sorted(full_de["target"].unique().to_list()))
    raw = A.fast_score(sb, sa, own_controls=True)
    full_gate = _nmae_ref_from_tables(full_de, sa["de"], sb["de"], p_adj_threshold=.05, min_gate_size=10,
                                      target_resolution=resolution)
    raw["nmae"] = float(full_gate["nmae_ref_raw"].mean())
    return raw, int(full_gate.height)


def anchors_rep(spec: dict, T: dict, full, *, n_splits: int = 5, base_seed: int = 0):
    """anchors323.replicate() in one process (fixtures and small surfaces; the access runs each half in its own process):
    per split the halves (replicate_halves), each half's stats against its own controls, replicate_score. Returns
    (replicate raw mean, per-split rows, seconds)."""
    t0 = time.monotonic()
    open_truth(spec)
    xc, xt, pos = T["xc"], T["xt"], T["pos"]
    splits = []
    for i in range(n_splits):
        ga, gb = replicate_halves(spec, i, n_splits=n_splits, base_seed=base_seed, exclude=T.get("exclude", ()))

        def getter(gd):
            return lambda g: (xc if g == CTRL else xt)[pos[gd[g]]]
        sa, _ = _group_stats(getter(ga), list(ga), T["genes"])
        sb, _ = _group_stats(getter(gb), list(gb), T["genes"])
        raw, cohort = replicate_score(sa, sb, full["de"])
        from cell_eval2.anchor import _derive_seeds
        splits.append({"split": i, "seed": int(_derive_seeds(base_seed, n_splits)[i]), "raw": raw,
                       "nmae_full_gate_cohort": cohort})
        log(f"anchors: split {i + 1}/{n_splits} {json.dumps(raw)} ({time.monotonic() - t0:.1f}s, peak RSS {rss_gb():.2f} GB)")
        del sa, sb
    rep = {m: float(np.mean([s["raw"][m] for s in splits])) for m in MEMBERS}
    return rep, splits, time.monotonic() - t0


def anchors(spec: dict, T: dict, *, n_splits: int = 5, base_seed: int = 0, ctrl_block=None, ctrl_moments=None):
    """anchors_base + anchors_rep: (record, full) with record = baseline raw, replicate raw (mean of the splits), the
    per-split rows."""
    base, full, diag = anchors_base(spec, T, ctrl_block=ctrl_block, ctrl_moments=ctrl_moments)
    rep, splits, sec = anchors_rep(spec, T, full, n_splits=n_splits, base_seed=base_seed)
    return dict(baseline=base, replicate=rep, replicate_splits=splits, base_seed=base_seed, n_splits=n_splits,
                seconds=diag["seconds_base"] + sec, baseline_emission=diag["baseline_emission"]), full


def scaled(raw: dict, base: dict, rep: dict) -> dict:
    """cell-eval2's replicate-scale score per member: score_one with the catalog policy anchored at the replicate
    (score._replicate_entries), i.e. (u - b) / (r - b) with the mse [0, 1] clamp and nmae's linear -6 floor."""
    from dataclasses import replace
    from cell_eval2.catalog import CATALOG
    from cell_eval2.scoring import score_one
    out = {}
    for m in MEMBERS:
        spec = CATALOG[FULL[m]]
        policy = replace(spec.scoring, anchor=float(rep[m]), allow_negative_baseline=False)
        out[m] = float(score_one(float(raw[m]), float(base[m]), policy))
    out["overall"] = float(np.mean([out[m] for m in MEMBERS]))
    return out


def scaled_simple(raw: dict, base: dict, rep: dict) -> dict:
    """report323.scaled(): (u - b) / (r - b), mse clipped to [0, 1], nmae floored at -6 (the cross-check of `scaled`)."""
    out = {}
    for m in MEMBERS:
        v = (raw[m] - base[m]) / (rep[m] - base[m])
        if m == "mse":
            v = min(max(v, 0.0), 1.0)
        if m == "nmae":
            v = max(v, -6.0)
        out[m] = float(v)
    out["overall"] = float(np.mean([out[m] for m in MEMBERS]))
    return out
