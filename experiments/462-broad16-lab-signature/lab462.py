"""462 — are perturbation-response programmes lab signatures? A matched-target test on the Broad-16 release (one lab, Cas9
knockout, 12 working cancer lines, the same targets) against other labs' CRISPRi screens of the SAME targets.

450h-l (whole screens, then matched sets on the Weissman/Xaira/Satija banks) found that a bank's per-gene programme agrees
across lines within a lab (r ~0.5) far more than across labs within a line (K562: ~0.1-0.2). 458's board arm (09-30)
bets on it: Arc's own H1 programme as Arc's cancer contexts' block prior. The Broad-16 adds 66 within-lab, cross-line
pairs on one fixed target set, and the same targets are screened by K562-GWPS, Orion HCT116/HEK293T and VIPerturb K562.

programme(bank) = per-gene P_down over the common target set T: the fraction of T's responses < 0 among the non-zero
ones (def_sameline's 396 rule). Agreement = Pearson r of (P_down - 0.5) over a common gene axis; sign agreement on
genes where both |P - 0.5| >= 0.05. Reliability = target split-halves (P_down on two disjoint random halves of T, 50
splits, Spearman-Brown to full length), so a cross-bank r can be disattenuated by each bank's own noise.

Broad responses: 268's `load268.build_responses` (log1p(5e4 p / sum p) minus the same line's 10 olfactory-receptor KO
controls, genes measured in every working line, the 100 Broad target genes removed, mean control CPM >= 5), working lines
only (268's reserve), targets with >= 6 cells in every working line.

usage: uv run python experiments/462-broad16-lab-signature/lab462.py   -> work/lab462.json, work/lab462.md
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "experiments/268-broad16-context-residuals/work"))
sys.path.insert(0, str(ROOT / "scripts"))
import load268 as B  # noqa: E402
import def_sameline as S  # noqa: E402

S.BANK_DIRS.update({"orion_hek293t": ROOT / "data/donors/orion_hek293t", "rpe1": ROOT / "data/donors/rpe1"})
LAB = {"k562_gwps": "Weissman", "k562_essential": "Weissman", "nadig_jurkat": "Weissman", "nadig_hepg2": "Weissman",
       "rpe1": "Weissman", "orion_hct116": "Xaira", "orion_hek293t": "Xaira", "viperturb_k562": "Satija"}
LINE = {"k562_gwps": "K562", "k562_essential": "K562", "nadig_jurkat": "Jurkat", "nadig_hepg2": "HepG2", "rpe1": "RPE1",
        "orion_hct116": "HCT116", "orion_hek293t": "HEK293T", "viperturb_k562": "K562"}
MAIN = ("k562_gwps", "orion_hct116", "orion_hek293t", "viperturb_k562")          # the matched set's other-lab banks
WIDE = ("k562_essential", "nadig_jurkat", "nadig_hepg2", "rpe1")                   # a smaller second set adds these
SEED, NSPLIT = 462, 50
t0 = time.time()
log = lambda *x: print("[462]", *x, f"({time.time() - t0:.0f}s)", flush=True)


def pdown(M):
    """M: (targets, genes) responses, NaN = unmeasured. P_down over the non-zero finite entries per gene."""
    fin = np.isfinite(M) & (M != 0)
    dn = ((M < 0) & fin).sum(0)
    nz = fin.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(nz > 0, dn / np.maximum(nz, 1), np.nan), nz


def rel_split(M, rng):
    """Target split-half reliability of the programme, Spearman-Brown corrected to the full target count."""
    n = M.shape[0]
    rs = []
    for _ in range(NSPLIT):
        p = rng.permutation(n)
        a, _ = pdown(M[p[: n // 2]])
        b, _ = pdown(M[p[n // 2:]])
        ok = np.isfinite(a) & np.isfinite(b)
        rs.append(np.corrcoef(a[ok] - .5, b[ok] - .5)[0, 1])
    r = float(np.mean(rs))
    return r, 2 * r / (1 + r) if r > -1 else np.nan


def pair_stats(pa, pb):
    ok = np.isfinite(pa) & np.isfinite(pb)
    r = float(np.corrcoef(pa[ok] - .5, pb[ok] - .5)[0, 1])
    strong = ok & (np.abs(pa - .5) >= .05) & (np.abs(pb - .5) >= .05)
    sa = float(np.mean(np.sign(pa[strong] - .5) == np.sign(pb[strong] - .5))) if strong.any() else np.nan
    return r, sa, int(ok.sum()), int(strong.sum())


def run(bank_names, tag, H, genes, lines, targets0):
    banks = {b: S.load_bank(b) for b in bank_names}
    T = [t for t in targets0 if all(t in banks[b][1] for b in bank_names)]
    log(f"{tag}: matched targets {len(T)} (Broad usable {len(targets0)}; banks {list(bank_names)})")
    R, ok, avail = B.build_responses(H, genes, lines, T, keys=("full",), min_cpm=5)
    gsym = np.asarray(genes)[ok]
    full = R["full"].astype(np.float64)                                      # (T, L, G)
    full[~avail] = np.nan
    progs, mats = {}, {}
    for j, l in enumerate(lines):
        mats[f"broad:{l}"] = full[:, j, :]
    for b in bank_names:
        delta, bperts, bgenes, _ = banks[b]
        pos = {g: i for i, g in enumerate(bgenes)}
        col = np.array([pos.get(g, -1) for g in gsym])
        rows = np.asarray(delta[[bperts.index(t) for t in T]], dtype=np.float64)
        M = np.where(col >= 0, rows[:, np.maximum(col, 0)], np.nan)
        mats[b] = M
    # a common gene axis: measured (non-zero somewhere) in every bank
    cover = np.ones(len(gsym), bool)
    for k, M in mats.items():
        cover &= (np.isfinite(M) & (M != 0)).any(0)
    rng = np.random.default_rng(SEED)
    rel = {}
    for k, M in mats.items():
        progs[k], _ = pdown(M[:, cover])
        rel[k] = rel_split(M[:, cover], rng)
    names = list(mats)
    lab = lambda k: "Broad" if k.startswith("broad:") else LAB[k]
    line = lambda k: k.split(":")[1] if k.startswith("broad:") else LINE[k]
    pairs = []
    for a, b in itertools.combinations(names, 2):
        r, sa, n, ns = pair_stats(progs[a], progs[b])
        den = np.sqrt(max(rel[a][1], 1e-9) * max(rel[b][1], 1e-9))
        same_lab, same_line = lab(a) == lab(b), line(a) == line(b)
        ptype = ("same line, same lab" if same_lab and same_line else "same line, different lab" if same_line
                 else "different line, same lab" if same_lab else "different line, different lab")
        pairs.append(dict(a=a, b=b, lab_a=lab(a), lab_b=lab(b), type=ptype, r=r, r_disatt=r / den, sign_agree=sa,
                          n_genes=n, n_strong=ns))
    summary = {}
    for ptype in sorted({p["type"] for p in pairs}):
        for labs in sorted({tuple(sorted((p["lab_a"], p["lab_b"]))) for p in pairs if p["type"] == ptype}):
            sel = [p for p in pairs if p["type"] == ptype and tuple(sorted((p["lab_a"], p["lab_b"]))) == labs]
            summary[f"{ptype} | {'×'.join(labs)}"] = dict(
                pairs=len(sel), r_mean=float(np.mean([p["r"] for p in sel])),
                r_median=float(np.median([p["r"] for p in sel])),
                r_disatt_mean=float(np.mean([p["r_disatt"] for p in sel])),
                sign_agree_mean=float(np.nanmean([p["sign_agree"] for p in sel])))
    frac_down = {k: float(np.nanmean(v < .5)) for k, v in progs.items()}
    return dict(tag=tag, targets=T, n_targets=len(T), n_genes=int(cover.sum()), banks=list(bank_names), lines=list(lines),
                reliability={k: dict(split_half_r=v[0], spearman_brown=v[1]) for k, v in rel.items()},
                frac_down=frac_down, summary=summary, pairs=pairs)


def main():
    lines = list(B.WORKING)
    H, genes = B.load_halves(lines)
    targets0 = [t for t in B.usable_targets(H, lines, min_cells=6) if t not in B.CONTROLS]
    log(f"Broad working lines {len(lines)}; usable non-control targets {len(targets0)}")
    out = {"main": run(MAIN, "main", H, genes, lines, targets0),
           "wide": run(MAIN + WIDE, "wide", H, genes, lines, targets0)}
    (HERE / "work").mkdir(exist_ok=True)
    (HERE / "work/lab462.json").write_text(json.dumps(out, indent=1))
    md = []
    for k, res in out.items():
        md.append(f"### {k}: {res['n_targets']} matched targets, {res['n_genes']} genes, banks {res['banks']} + "
                  f"{len(res['lines'])} Broad lines\n")
        md.append("| pair type | labs | pairs | mean r | median r | disattenuated r | sign agreement |\n|---|---|---:|---:|---:|---:|---:|")
        for key, v in res["summary"].items():
            ptype, labs = key.split(" | ")
            md.append(f"| {ptype} | {labs} | {v['pairs']} | {v['r_mean']:+.3f} | {v['r_median']:+.3f} | "
                      f"{v['r_disatt_mean']:+.3f} | {v['sign_agree_mean']:.3f} |")
        md.append("\nreliability (Spearman-Brown, target split-halves): " + ", ".join(
            f"{k.replace('broad:', '')} {v['spearman_brown']:.2f}" for k, v in res["reliability"].items()))
        md.append("\nfraction of genes down-programmed: " + ", ".join(
            f"{k.replace('broad:', '')} {v:.2f}" for k, v in res["frac_down"].items()) + "\n")
    (HERE / "work/lab462.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
