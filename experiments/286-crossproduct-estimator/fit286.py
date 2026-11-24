"""286 -- ONE DONOR LINE PER PROCESS: the E1 / E2 filters from that line's cross-product moments, applied to the
line's PRODUCTION un-whitened rows for all 300 panel targets.  The fold sweep then only reads small arrays.

    PYTHONDONTWRITEBYTECODE=1 uv run python fit286.py orion_hct116

Per source covariance (the arm and its three controls) two symmetric eigenproblems are solved in float32 with LAPACK
`evd`: eigh(S-hat) (the PSD projection, E2's truncation, the mechanism reading) and eigh(L^-1 S+ L^-T) (E1's exact
Wiener gains).  Sources:
    real       the line's own S-hat, N-hat                                        the ARM
    shufpair   S-hat, N-hat from the shuffled target pairing (control Cs)
    wrong      ANOTHER line's S-hat on this line's columns, own N-hat (control Cw)
    rand       an energy-matched random PSD covariance, own N-hat (control Cr)
Mirror (control Cm, G' = 2I - G_phi) is free from the arm's own factors.

The k562_gwps production row carries the PerturbMatch 0.5 blend on the shared columns exactly as
`lib263.prepared_sources` builds it (asserted against the fold's staged rows by the sweep, KA-C).

Outputs: work/filt/<line>.npz -- one (300, 18533) float32 row block per arm, plus the top-100 eigenvectors of S-hat,
its spectrum, the Wiener gains and the audits; work/filt/<line>.json.
"""
from __future__ import annotations

import gc
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/286-crossproduct-estimator"
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(HERE))
import plolo  # noqa: E402
import est286 as E  # noqa: E402

plolo.WORK = WORK / "plolo"
E062, E064 = plolo.E062, plolo.E064
OFFICIAL, PANEL, PANEL_SET = plolo.OFFICIAL, list(plolo.PANEL), plolo.PANEL_SET
WRONG_SOURCE = {"k562_gwps": "orion_hct116", "orion_hct116": "orion_hek293t", "orion_hek293t": "orion_hct116"}
N_TOPVEC = 100
LOG: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def rss() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def production_rows(line: str) -> tuple[np.ndarray, np.ndarray, dict]:
    """The line's production un-whitened rows for all 300 panel targets, with the PerturbMatch 0.5 blend on
    k562_gwps's shared columns (lib263.prepared_sources, verbatim arithmetic)."""
    z = np.load(WORK / "prod" / f"{line}.npz", allow_pickle=True)
    rows = np.asarray(z["rows"], dtype=np.float64)
    avail = np.asarray(z["avail"], bool)
    aud = dict(raw_median=float(z["raw_median"]), mapped=int(z["mapped"]), pm_rows=0)
    if line == "k562_gwps":
        official_pos = {g: i for i, g in enumerate(OFFICIAL)}
        panel_positions = np.asarray([official_pos[g] for g in PANEL if g in official_pos], dtype=np.int64)
        bank = E062.Bank.load(line, plolo.DONOR_PATHS[line])
        pm_bank = E064.load_perturbmatch(plolo.PM_PATH)
        pm = E062.prepare_main_source(pm_bank, PANEL, OFFICIAL, panel_positions, chunk_rows=256, max_rss_gib=9.0)
        shared = E064.shared_k562_pm_columns(bank, pm_bank, OFFICIAL, PANEL_SET)
        both = avail & np.asarray(pm.target_available, bool)
        r = np.flatnonzero(both)
        if r.size:
            rows[np.ix_(r, shared)] = ((1 - plolo.PM_BLEND) * rows[np.ix_(r, shared)]
                                       + plolo.PM_BLEND * pm.effects[np.ix_(r, shared)])
        aud["pm_rows"] = int(both.sum())
        aud["pm_shared_cols"] = int(len(shared))
        del pm, pm_bank, bank
        gc.collect()
    return rows, avail, aud


def main(line: str) -> None:
    t0 = time.time()
    P(f"[{line}] est286.py sha256 {E.sha256_self()}")
    (WORK / "filt").mkdir(parents=True, exist_ok=True)
    meta = json.loads((WORK / "moments" / f"{line}.json").read_text())
    assert meta["est286_sha256"] == E.sha256_self(), "est286.py changed since staging"
    cols = np.load(WORK / "moments" / f"{line}_cols.npy")
    m, n_fit = int(cols.size), int(meta["n_rows"])
    rows, avail, aud_prod = production_rows(line)
    X = rows[:, cols].copy()
    P(f"  production rows {rows.shape}, {int(avail.sum())}/300 available, PM-blended rows {aud_prod['pm_rows']}; "
      f"m = {m}, n_fit = {n_fit}; rss {rss():.2f} GiB")

    N = np.load(WORK / "moments" / f"{line}_N.npy")                    # float32 (m, m)
    out: dict = {}
    aud: dict = dict(line=line, m=m, n_fit=n_fit, est286_sha256=E.sha256_self(), prod=aud_prod,
                     wrong_source=WRONG_SOURCE[line], sources={})

    def source(tag: str):
        if tag == "real":
            return np.load(WORK / "moments" / f"{line}_S.npy"), N, None
        if tag == "shufpair":
            return (np.load(WORK / "moments" / f"{line}_S_shuf.npy"),
                    np.load(WORK / "moments" / f"{line}_N_shuf.npy"), None)
        if tag == "wrong":
            other = WRONG_SOURCE[line]
            oc = np.load(WORK / "moments" / f"{other}_cols.npy")
            So = np.load(WORK / "moments" / f"{other}_S.npy")
            pos = {int(c): i for i, c in enumerate(oc.tolist())}
            keep = np.asarray([pos[int(c)] for c in cols if int(c) in pos], dtype=np.int64)
            mine = np.asarray([i for i, c in enumerate(cols) if int(c) in pos], dtype=np.int64)
            S = np.zeros((m, m), dtype=np.float32)
            S[np.ix_(mine, mine)] = So[np.ix_(keep, keep)]
            del So
            gc.collect()
            return S, N, dict(other=other, shared_cols=int(mine.size), of=int(m))
        if tag == "rand":
            Sr = np.load(WORK / "moments" / f"{line}_S.npy")
            tr = float(np.trace(Sr.astype(np.float64)))
            del Sr
            gc.collect()
            return E.random_psd(m, n_fit, max(tr, 1e-12)), N, dict(target_trace=tr)
        raise ValueError(tag)

    for tag in ("real", "shufpair", "wrong", "rand"):
        t1 = time.time()
        S, Nj, extra = source(tag)
        lam, V, a_eig = E.psd_eigh(S, log=P)
        del S
        gc.collect()
        if tag == "real":
            out["topvec"] = np.ascontiguousarray(V[:, ::-1][:, :N_TOPVEC], dtype=np.float32)
            out["lam"] = lam.astype(np.float64)
            out["nu"] = E.direction_noise(Nj, V).astype(np.float64)
            lam_top = lam[::-1][:N_TOPVEC].copy()
        Sp = E.s_plus(lam, V)
        del V                            # E2 uses the saved top-100 eigenpairs; the full V is not needed again
        gc.collect()
        L, Q, k, a_w = E.wiener_factors(Sp, Nj, 1.0, log=P)
        del Sp
        gc.collect()
        for phi in E.PHIS:
            if tag != "real" and phi == 0.0:
                continue
            out[f"e1:phi{phi:g}" + ("" if tag == "real" else f":{tag}")] = \
                E.apply_e1(X, L, Q, k, float(phi)).astype(np.float32)
        if tag == "real":
            out["gains"] = (k / (1.0 + k)).astype(np.float64)
            for phi in (0.10, 0.25):
                out[f"e1:phi{phi:g}:mirror"] = E.apply_e1(X, L, Q, k, float(phi), mirror=True).astype(np.float32)
            for r in E.RANKS:
                out[f"e2:r{r}"] = E.apply_e2(X, L, lam_top, out["topvec"], int(r)).astype(np.float32)
        aud["sources"][tag] = dict(eig=a_eig, wiener=a_w, extra=extra, seconds=round(time.time() - t1, 1))
        del L, Q, k, lam
        gc.collect()
        P(f"  [{tag}] done {time.time()-t1:.0f}s, rss {rss():.2f} GiB")

    # scatter the filtered columns back onto the official axis (columns outside `cols` pass through unchanged)
    blocks = {}
    for key, val in out.items():
        if key in ("topvec", "lam", "nu", "gains"):
            blocks[key] = val
            continue
        full = rows.astype(np.float32).copy()
        full[:, cols] = val
        blocks[key] = full
    blocks["E0"] = rows.astype(np.float32)
    blocks["avail"] = avail
    blocks["targets"] = np.asarray(PANEL)
    blocks["cols"] = cols
    np.savez_compressed(WORK / "filt" / f"{line}.npz", **blocks)
    aud["arms"] = sorted(k for k in blocks if k not in ("topvec", "lam", "nu", "gains", "avail", "targets", "cols"))
    aud["seconds"] = round(time.time() - t0, 1)
    aud["peak_rss_gib"] = rss()
    (WORK / "filt" / f"{line}.json").write_text(json.dumps(aud, indent=1, default=str))
    (WORK / "filt" / f"{line}.log").write_text("\n".join(LOG) + "\n")
    P(f"[{line}] done {(time.time()-t0)/60:.1f} min, {len(aud['arms'])} arms, peak rss {rss():.2f} GiB")


if __name__ == "__main__":
    for ln in sys.argv[1:]:
        main(ln)
