#!/usr/bin/env python3
"""def_identity_prior.py: the identity-augmented CLEAN4 block prior (council 600: 639 d4cd5f, 622 7772dc, 604 4d0207).

Biology (639, 622): a gene a line expresses above the other lines' controls tends to fall under ANY perturbation
(perturbation erodes the line's identity programme). Combined with CLEAN4, it beats CLEAN4 leave-one-line-out on the
gate genes' signs (636's instrument: hct116 / hek293t / h1, CIs above 0; Broad-16 12 lines, mean +0.009 at lambda 0.5).

For every context c of a bundle:
  references(c) = the bundle's other contexts + the pinned reference profiles (A/B/C challenge controls, H1 Flex 2025),
                  ONE per distinct identity call (bundle contexts first), and never one whose call equals c's own
                  (622e: a same-line reference kills the signal). Only Flex profiles are used (622: 3' references hurt).
                  SISTER guard (626, 7772dc): drop a reference whose Spearman similarity to the context, computed on
                  log1p-control-CPM DEVIATIONS from the pool mean, exceeds 0.3 (most similar first), never leaving fewer
                  than 2 references. On Arc's {A, B, C, H1} every pairwise similarity is negative, so it never fires
                  there; on the Broad-16 it fires only on the true sister pair, at a cost of -0.0001.
  I_c           = log1p(ctrl CPM_c) - mean over references of log1p(ctrl CPM), on the genes c and every reference measure
                  (a NaN in a profile = not measured; H1's pinned profile carries only its > 5 CPM genes, as 480's gate file).
  score         = zc + lam * z(I_c), on CLEAN4-COVERED genes only (a covered gene with no I_c gets zc alone), where
                  zc = (P_C4 - 1/2) / sd(P_C4 - 1/2) is CENTRED AT P = 1/2 (639's validated form, not mean-centred) and
                  z(I_c) is mean/sd over the genes where I_c is defined; positive -> down.
  P (v3)        = SIDE from sign(score) (639 LOLO / 639d's validated rule); CONFIDENCE |P - 1/2| = CLEAN4's own |P - 1/2|
                  distribution on covered genes, assigned by the rank of |score| (stable; ties broken by CLEAN4's own
                  |P - 1/2|, then gene index, so lam 0 returns CLEAN4's P bit for bit). v2 rank-mapped the SIGNED score onto
                  CLEAN4's P, which kept CLEAN4's down-count and undid part of the identity down-shift (most in A).
                  Identity-ONLY genes (CLEAN4 undefined) get NO prior at any lam: 639d (d4cd5f, 1f794019) measured that
                  extension's signs right only 0.27-0.34 on HEK293T's gate genes and 0.49 on HCT116's uncalled ones, so
                  they keep CLEAN4's pinned psign and no P, exactly as without this module (639's block prereg (ii)).
                  psign = -1 where P > 0.5, +1 where P < 0.5, CLEAN4's pinned psign where P == 0.5 or undefined.
KA (hard, asserted): lam 0 returns CLEAN4's psign on every gene and CLEAN4's P_down on its covered genes.

As a library: `identity_prior(...)`. As a script (rehearsals / inspection):
  uv run python scripts/def_identity_prior.py --bundle <dir> --identity <identity.json> --lam 0.5 --out <npz>
writes a VCC_SAMELINE-format artifact (programme_<c>, pdown_<c>, w 0, no rows) for every context of the bundle.
On the day it runs inside `def_sameline.py --programme-source identity:<lam>` (patch 604), which keeps any same-line rows.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1] if Path(__file__).resolve().parent.name == "scripts" else Path(__file__).resolve().parents[2]
CLEAN4_PATH = ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz"
CLEAN4_SHA256 = "7e532f84a2425176a4582a2e42e458bc3fd25bc3267f890b27a530565df5936e"
REFS_PATH = ROOT / "experiments/604-identity-prior/work/flex_refs604.npz"
SISTER_SIM_MAX = 0.3   # 626 (7772dc): Broad-16's one true sister pair reads 0.343; non-sisters max 0.256
DEFAULT_LAM = 0.5      # the council's ship point (622g: Broad-16 mean +0.009 at 0.5; lambda 1 only as a dose partner)
REFS_SHA256 = "07762ab11be6461e084cd0ed64c1dd9d36029d2fc5b130d664bbfeeea3111598"   # pinned by patch604 (flex_refs604.npz)
# 625's CLEAN4-strength octiles (7772dc; v8 design): `octile` 0 = P_C4 < 0.440 (strongest UP prior) ... 7 = P_C4 > 0.583
# (strongest DOWN); -1 = not covered. NOT the `module` array (that is the Hadamard column order).
MODULES625_PATH = ROOT / "experiments/625-module-readout/work/modules625_v8.npz"
MODULES625_SHA256 = "c24e3304ffd63d77d76b4a74341e09a2ee8396d7e372955e9c6feee167133dc1"


def file_sha256(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def ctrl_cpm(h5ad_path: Path, genes: list[str]) -> np.ndarray:
    """Per-cell CPM mean of a context file's cells (the scorer's per-cell CPM; 639's abc_profiles), float64."""
    import anndata as ad
    import scipy.sparse as sp
    a = ad.read_h5ad(h5ad_path, backed="r")
    if list(map(str, a.var_names)) != list(genes):
        raise SystemExit(f"{h5ad_path}: gene axis differs from gene_names.csv")
    X = sp.csr_matrix(a.X[:])
    lib = np.asarray(X.sum(1)).ravel().astype(np.float64)
    if (lib <= 0).any():
        raise SystemExit(f"{h5ad_path}: a control cell has zero counts")
    return np.asarray((sp.diags(1e6 / lib) @ X).mean(0)).ravel().astype(np.float64)


def _z(v: np.ndarray) -> np.ndarray:
    out = np.full(v.shape, np.nan)
    m = np.isfinite(v)
    if m.sum() > 1 and v[m].std() > 0:
        out[m] = (v[m] - v[m].mean()) / v[m].std()
    return out


def load_clean4(genes) -> tuple[np.ndarray, np.ndarray]:
    if file_sha256(CLEAN4_PATH) != CLEAN4_SHA256:
        raise SystemExit("CLEAN4 prior sha256 moved")
    with np.load(CLEAN4_PATH, allow_pickle=False) as z:
        if list(map(str, z["genes"])) != list(map(str, genes)):
            raise SystemExit("CLEAN4 gene axis differs from the official gene order")
        return np.asarray(z["P_down"], dtype=np.float64), np.asarray(z["psign"], dtype=np.float64)


def load_octiles(genes) -> np.ndarray:
    """625's octile per gene of `genes` (matched BY NAME; -1 = in no octile), from the sha-pinned module file."""
    if file_sha256(MODULES625_PATH) != MODULES625_SHA256:
        raise SystemExit(f"{MODULES625_PATH}: sha256 differs from the pinned {MODULES625_SHA256}")
    with np.load(MODULES625_PATH, allow_pickle=True) as z:   # object gene names; the file is sha-pinned above
        pos = {str(g): i for i, g in enumerate(z["genes"])}
        oc = np.asarray(z["octile"], dtype=np.int64)
    return np.array([oc[pos[g]] if g in pos else -1 for g in map(str, genes)], dtype=np.int64)


def flip_vector(genes, P_c4: np.ndarray, flip_octiles: str | None) -> np.ndarray | None:
    """Covered genes whose CLEAN4 prior is reversed (P -> 1 - P) by an 8-character 0/1 octile mask; None = no flip."""
    if flip_octiles is None:
        return None
    if len(flip_octiles) != 8 or set(flip_octiles) - {"0", "1"}:
        raise SystemExit(f"--flip-octiles must be 8 characters of 0/1 in octile order 0..7 (got {flip_octiles!r})")
    if "1" not in flip_octiles:
        return None
    oc = load_octiles(genes)
    chosen = np.array([flip_octiles[k] == "1" for k in range(8)])
    return np.isfinite(P_c4) & (oc >= 0) & chosen[np.clip(oc, 0, 7)]


def identity_prior(genes, ctx_cpm: dict, ctx_calls: dict, ref_cpm: dict, ref_calls: dict, lam: float,
                   P_c4: np.ndarray, ps_c4: np.ndarray, flip_octiles: str | None = None) -> dict:
    """Per context: {psign, pdown, prov}. ctx_cpm / ref_cpm: name -> CPM vector on `genes` (NaN = not measured).
    flip_octiles (625, 7772dc): an 8-character 0/1 mask over 625's CLEAN4-strength octiles. The arm is the BASE (no-flip)
    artifact with ONLY the masked octiles' SIDES re-decided: their prior becomes P_C4 -> 1 - P_C4 (the pinned psign negated
    where P != 1/2), their side is sign(zc' + lam * zi) with zc keeping the unflipped prior's scale, and they keep the
    BASE's confidence |P - 1/2|. Every non-octile gene's pdown and psign are byte-identical to the base (15:05Z, 7772dc's
    option b: a tau > 0 hook then sees no change outside the octile). At lam 0 this is exactly the flipped prior (KA).
    None or all zeros = no flip, bit-identical (the provenance gains no flip keys)."""
    fv = flip_vector(genes, P_c4, flip_octiles)
    # zc's SCALE always comes from the unflipped prior (14:55Z fix): recomputing sd(P - 1/2) from the flipped prior would
    # rescale zc for EVERY covered gene and, at lam > 0, move non-octile genes' sides too, which contaminates a module
    # readout. With the original scale, only the masked octiles' scores change. At lam 0 the scale cannot matter (side and
    # rank-mapped confidence are scale-free), and without a flip it is the same number as before (bit-identical).
    zc_scale = float((P_c4[np.isfinite(P_c4)] - 0.5).std())
    P_eff, ps_eff = P_c4, ps_c4
    if fv is not None:
        P_eff = np.where(fv, 1.0 - P_c4, P_c4)
        ps_eff = np.where(fv & (P_eff != 0.5), -ps_c4, ps_c4)
    covered = np.isfinite(P_c4)
    cov_idx = np.flatnonzero(covered)
    out = {}
    for c in ctx_cpm:
        own = ctx_calls.get(c)
        refs, seen = [], {own} if own is not None else set()
        # 607 (cloud review): a pinned reference whose NAME equals a bundle context's (a rehearsal reusing ids A/B/C with
        # other calls) must not overwrite that context in the pool dict; it is renamed pin:<name>. With no collision every
        # name, and so every artifact byte, is unchanged.
        for name, prof, call in ([(o, ctx_cpm[o], ctx_calls.get(o)) for o in ctx_cpm if o != c]
                                 + [((r if r not in ctx_cpm else f"pin:{r}"), ref_cpm[r], ref_calls.get(r)) for r in ref_cpm]):
            key = call if call is not None else f"__{name}"
            if key in seen:
                continue
            seen.add(key)
            refs.append((name, prof))
        if not refs:
            raise SystemExit(f"identity prior: context {c} has no reference left after the guard")
        # 626 sister guard: similarity of log1p deviations from the pool mean (pool = the context + its references)
        pmask = np.isfinite(ctx_cpm[c]).copy()
        for _, prof in refs:
            pmask &= np.isfinite(prof)
        Lp = {n: np.log1p(p[pmask]) for n, p in [(c, ctx_cpm[c])] + refs}
        pool_mean = np.mean(list(Lp.values()), axis=0)
        def _rank(v):
            r = np.empty(len(v)); r[np.argsort(v, kind="stable")] = np.arange(len(v)); return r
        dc = _rank(Lp[c] - pool_mean)
        sims = {n: float(np.corrcoef(dc, _rank(Lp[n] - pool_mean))[0, 1]) for n, _ in refs}
        dropped = []
        for n in sorted(sims, key=lambda k: -sims[k]):
            if sims[n] > SISTER_SIM_MAX and len(refs) > 2:
                refs = [(m, p) for m, p in refs if m != n]
                dropped.append(n)
        mask = np.isfinite(ctx_cpm[c]).copy()
        for _, prof in refs:
            mask &= np.isfinite(prof)
        I = np.full(len(genes), np.nan)
        I[mask] = np.log1p(ctx_cpm[c][mask]) - np.mean([np.log1p(p[mask]) for _, p in refs], axis=0)
        zc4 = np.full(len(genes), np.nan)
        zc4[covered] = (P_c4[covered] - 0.5) / zc_scale   # v3: centred at P = 1/2 (639d's zc); scale of the UNFLIPPED prior
        zi = _z(I)
        has_i = np.isfinite(zi)
        score = np.where(covered, zc4, 0.0) + lam * np.where(has_i, zi, 0.0)
        scored = covered   # 639d: identity-only genes get NO prior at any lam (the extension's signs are wrong-leaning)
        P = np.full(len(genes), np.nan)
        # v3 (d4cd5f 13:30Z): side = sign(score); confidence = CLEAN4's own |P - 1/2| values, handed out by rank of |score|.
        conf = np.abs(P_c4[cov_idx] - 0.5)
        src = cov_idx[np.lexsort((cov_idx, conf))]                               # CLEAN4 genes by |P - 1/2| (stable)
        dst = cov_idx[np.lexsort((cov_idx, conf, np.abs(score[cov_idx])))]       # genes by |score|, ties: own |P - 1/2|
        Psrc = P_c4[src]
        side_src, side_dst = np.sign(Psrc - 0.5), np.sign(score[dst])
        # 607: a score of exactly 0 has no side; it gets P = 1/2, so psign keeps that gene's own pinned CLEAN4 sign (never
        # occurs on the A/B/C artifacts; at lam 0 it is exactly the P == 1/2 genes, which keep 1/2)
        P[dst] = np.where(side_dst == 0, 0.5, np.where(side_src == side_dst, Psrc, 1.0 - Psrc))
        ido = scored & ~covered   # empty by construction since 639d; kept so the provenance field stays explicit
        if fv is not None:
            # 625 flip (option b): re-decide ONLY the octile genes' sides from the flipped prior; keep the base's |P - 1/2|
            zcf = zc4.copy()
            zcf[fv] = (P_eff[fv] - 0.5) / zc_scale
            side_f = np.sign(np.where(covered, zcf, 0.0) + lam * np.where(has_i, zi, 0.0))
            P_base = P.copy()
            turn = fv & (side_f != 0) & (np.sign(P_base - 0.5) != side_f)
            P[turn] = 1.0 - P_base[turn]
            P[fv & (side_f == 0)] = 0.5
        psign = ps_eff.copy()
        psign[np.isfinite(P) & (P > 0.5)] = -1.0
        psign[np.isfinite(P) & (P < 0.5)] = 1.0
        if lam == 0.0:
            if not (np.array_equal(P[covered], P_eff[covered]) and np.array_equal(psign, ps_eff)):
                raise SystemExit(f"identity prior KA FAILED for {c}: lam 0 must reproduce CLEAN4 (flipped where masked) exactly")
        prov = {"lam": lam, "own_call": own, "references": [n for n, _ in refs], "identity_genes": int(mask.sum()),
                "sister_similarity": {n: round(v, 4) for n, v in sims.items()}, "sister_dropped": dropped,
                "identity_only_scored": int(ido.sum()), "scored": int(scored.sum()),
                "down_frac": float((psign < 0).mean()), "down_frac_clean4": float((ps_eff < 0).mean()),
                "flips_vs_clean4": int((psign != ps_eff).sum())}
        if fv is not None:   # flips_vs_clean4 / down_frac_clean4 above are then relative to the FLIPPED prior
            prov.update({"flip_octiles": flip_octiles, "flipped_genes": int(fv.sum()), "modules625_sha256": MODULES625_SHA256})
        out[c] = {"psign": psign, "pdown": P, "prov": prov}
    return out


def load_refs(path: Path = REFS_PATH, sha: str | None = REFS_SHA256) -> tuple[list, dict, dict]:
    if sha is not None and file_sha256(path) != sha:
        raise SystemExit(f"{path}: sha256 differs from the pinned {sha}")
    with np.load(path, allow_pickle=False) as z:
        genes = list(map(str, z["genes"]))
        names = list(map(str, z["names"]))
        calls = dict(zip(names, map(str, z["calls"])))
        cpm = {n: np.asarray(z[f"cpm_{n}"], dtype=np.float64) for n in names}
    return genes, cpm, calls


def main() -> None:
    import csv
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--identity", required=True, help="def_identity.py's JSON (its consensus `call` per context)")
    ap.add_argument("--lam", type=float, required=True)
    ap.add_argument("--refs", default=str(REFS_PATH))
    ap.add_argument("--refs-sha", default=REFS_SHA256)
    ap.add_argument("--out", required=True)
    ap.add_argument("--flip-octiles", default=None, help="625: 8 characters of 0/1 (octile order 0..7) reversing those "
                    "CLEAN4 octiles' prior before scoring; default none (bit-identical)")
    a = ap.parse_args()
    bundle = Path(a.bundle)
    with open(bundle / "gene_names.csv") as fh:
        genes = [r["gene_name"] for r in csv.DictReader(fh)]
    with open(bundle / "pert_counts.csv") as fh:
        perts = [r["target_gene"] for r in csv.DictReader(fh)]
    calls = dict(json.loads(Path(a.identity).read_text())["call"])
    rgenes, rcpm, rcalls = load_refs(Path(a.refs), a.refs_sha)
    if rgenes != genes:
        raise SystemExit("reference profiles' gene axis differs from the bundle's")
    ctx_cpm = {c: ctrl_cpm(bundle / f"context_{c}.h5ad", genes) for c in sorted(calls)}
    P_c4, ps_c4 = load_clean4(genes)
    res = identity_prior(genes, ctx_cpm, calls, rcpm, rcalls, a.lam, P_c4, ps_c4, flip_octiles=a.flip_octiles)
    art = {"genes": np.asarray(genes), "perts": np.asarray(perts), "contexts": np.asarray(sorted(res))}
    for c, r in res.items():
        art[f"programme_{c}"] = r["psign"]
        art[f"programme_ok_{c}"] = np.asarray(True)
        art[f"w_{c}"] = np.asarray(0.0)
        art[f"pdown_{c}"] = r["pdown"]
        art[f"targets_{c}"] = np.zeros(0, dtype=np.int64)
        art[f"rows_{c}"] = np.zeros((0, len(genes)), dtype=np.float32)
        art[f"signs_{c}"] = np.zeros((0, len(genes)), dtype=np.int8)
    prov = {"source": "def_identity_prior", "lam": a.lam, "refs": str(a.refs), "refs_sha256": file_sha256(Path(a.refs)),
            "clean4_sha256": CLEAN4_SHA256, "contexts": {c: r["prov"] for c, r in res.items()}}
    art["provenance"] = np.asarray(json.dumps(prov))
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out, **art)
    sha = file_sha256(Path(a.out))
    Path(str(a.out) + ".json").write_text(json.dumps({"sha256": sha, **prov}, indent=1))
    print(f"[identity-prior] wrote {a.out} sha256={sha} " + json.dumps({c: r["prov"] for c, r in res.items()}), flush=True)


if __name__ == "__main__":
    main()
