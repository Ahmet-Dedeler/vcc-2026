"""220 — DELIVER each ordering key through the shipped `rc_callset` emitter on h1-test.

WHY THIS EXISTS.  189 and 190 both measured `reach(key)` COUNTERFACTUALLY: they reordered a
frozen DE table and asked what `reach` would have been.  Under a15's `quasi` emitter that was
the only option, because the delivered order was emission noise.  Under `rc_callset` the
delivered order IS the key (`rc_emit.py`:983 -- delivered = prescribed U unsuppressible), so
the key can be SHIPPED and read back.  That is the one thing 214's finding genuinely buys the
ordering channel, and it did not exist when 189/190 ran.

THE ARM CLASS.  Only `key` moves.  `rc_hybrid_block` asserts `block.sum(0) == totals` for
every block, and `totals` comes from `floor(x)` with `x` the artifact's own intended profile,
which no key touches.  So the per-gene panel totals are IDENTICAL by construction and this is
rule 13's `emission-law, profile-PRESERVING` class -- toll 0.  The digest is printed for all
85 blocks rather than asserted in prose.

THE ONE CONFOUND, AND IT IS MEASURED.  `rc_callset` selects the head of its prescribed call
set as the top `n_head = 400` by `|key|` (`rc_emit.py`:917).  So a key change also moves the
CALL SET, which moves `fid` and `jac`.  Every arm reports the head-set Jaccard against the
incumbent's so the ordering half and the call-set half are never conflated.

KNOWN-ANSWER ROW (rule 10, abort): arm `absH` at seed 1 IS `rc_callset` and must reproduce
`202/work/de202.rc_callset.npz` bit for bit on `pred_p`, `pred_padj`, `pred_lfc`, `pred_z`.

Donors: none opened.  This is the frozen a15 h1-test artifact reconstructed by 196's `Recon`,
exactly as 214 did it.  `data/donors/h1_2025` is never touched.

    ARMS=... uv run python emit220.py
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
for p in ("experiments/196-emitted-depth", "experiments/185-rank-controlled-emission",
          "experiments/190-regulon-ordering", "src"):
    sys.path.insert(0, str(ROOT / p))

import harness185 as H                                                    # noqa: E402
import recon196                                                           # noqa: E402
import vcc2026.rc_emit as RC                                              # noqa: E402
from vcc2026.donor_transfer import emit_counts                            # noqa: E402

import common220 as C                                                     # noqa: E402

REF191 = ROOT / "experiments/191-rc-hybrid-build/work"
KA_TABLE = ROOT / "experiments/202-callset-rider/work/de202.rc_callset.npz"

n, PHI, DISPERSION = 400, 4.0, 0.1
N_HEAD, N_CPM, Z_CAP = 400, 2200, 1.90       # rc_callset's prescribed call set
KEEP_TOP, KEEP_SPREAD = 8000, 0.25           # a15's tiered spread policy
SEEDS = (90_000, 7_000_000, 11_111)          # 200/214's matched-seed ladder
MAXT = int(os.environ.get("MAXT", "0"))
TAG = os.environ.get("OUTTAG", "de220")


# ------------------------------------------------------------------ the key library
def build_keys(H0, genes_full, keep, targets):
    """Every delivered key, on the FULL gene axis, carrying `H`'s sign.

    A key is `H * m` with `m > 0` a per-(target, gene) multiplier, so `sign(key) = sign(H)`
    and the emitted SIGNS -- hence sign accuracy, hence `nmae`'s `q_w` -- are untouched.
    Only the ORDER and the top-400 head move.  That is what makes this an ordering arm.
    """
    G_full = genes_full.size
    theta_full = C.theta_on(genes_full)
    keys = {"absH": H0.copy()}
    meta = {"absH": dict(kind="incumbent |H|")}

    def add(name, mult_full, **info):
        keys[name] = H0 * mult_full
        meta[name] = info

    # -- 189: the per-gene cross-line reliability prior, gene-only multiplier ----------
    for lam in (0.25, 0.50, 1.00):
        m = np.exp(np.clip(lam * theta_full, -50, 50))[None, :]
        add(f"k189_l{int(lam*100):03d}", np.broadcast_to(m, H0.shape),
            kind="189 |H| x exp(lam theta)", lam=lam)

    # -- 189 shuffled controls: theta permuted across the WHOLE gene axis (189 F1) -----
    for j, sd in enumerate((11, 22, 33, 44, 55, 66, 77, 88, 99, 111, 122, 133)):
        th = np.random.default_rng(sd).permutation(theta_full)
        m = np.exp(np.clip(0.50 * th, -50, 50))[None, :]
        add(f"shuf189_{j+1}", np.broadcast_to(m, H0.shape),
            kind="189 SHUFFLED theta", lam=0.50, seed=sd)

    # -- 189 section 6.2's conditional form: reorder only inside the top X% of |H| ------
    for frac, lam in ((0.10, 1.0), (0.25, 0.5)):
        m = np.ones(H0.shape)
        gate = np.zeros(H0.shape, bool)
        for i in range(H0.shape[0]):
            a = np.abs(H0[i])
            k = max(1, int(round(frac * G_full)))
            gate[i, np.argpartition(-a, k - 1)[:k]] = True
        mm = np.broadcast_to(np.exp(np.clip(lam * theta_full, -50, 50))[None, :], H0.shape)
        m = np.where(gate, mm, 1.0)
        add(f"k189_top{int(frac*100):02d}", m,
            kind="189 conditional: top-X% of |H| reordered by theta", frac=frac, lam=lam)

    # -- 190: the regulatory prior, per-(target, gene) ---------------------------------
    import nets190
    gk = genes_full[keep]
    Hk = H0[:, keep]
    N, info = nets190.ct_prop(list(targets), list(gk), Hk, topk=25)
    R = C.agreement_R(N, Hk)
    for lam, nm in ((-0.25, "k190_km25"), (+1.00, "k190_p100")):
        mf = np.ones(H0.shape)
        mf[:, keep] = np.exp(np.clip(lam * R, -50, 50))
        add(nm, mf, kind="190 ct-prop-k25", lam=lam, nonzero=int((N != 0).sum()))
    for j, sd in enumerate((101, 202, 303, 404, 505, 606)):
        rng = np.random.default_rng(sd)
        Ns, _ = nets190.ct_prop(list(targets), list(gk), Hk, topk=25, rng=rng,
                                shuffle="signs")
        Rs = C.agreement_R(Ns, Hk)
        mf = np.ones(H0.shape)
        mf[:, keep] = np.exp(np.clip(-0.25 * Rs, -50, 50))
        add(f"shuf190_{j+1}", mf, kind="190 SHUFFLED signs", lam=-0.25, seed=sd)

    return keys, meta


def arms_for(keys):
    """(name, key_id, seed_index).  Three seeds for the three load-bearing arms."""
    three = ("absH", "k189_l050", "k190_km25")
    A = []
    for k in three:
        for s in range(len(SEEDS)):
            A.append((k if s == 0 else f"{k}_s{s+1}", k, s))
    for k in keys:
        if k in three:
            continue
        A.append((k, k, 0))
    return A


def main() -> None:
    t0 = time.time()
    ctl = H.load_controls()
    cs, tie, keep, cpm, ref_lin = (ctl["cs"], ctl["tie"], ctl["keep"], ctl["cpm"],
                                   ctl["ref_lin"])
    genes_full = ctl["genes"]
    G_full, G = genes_full.size, cs.shape[0]
    control = RC.RcControlBlock(keep=keep, cs=cs, tie=tie, ref_lin=ref_lin,
                                n_controls=cs.shape[1], min_cpm=H.MIN_CPM)
    call_set = RC.CallSet(N_HEAD, N_CPM, Z_CAP)

    D = recon196.load_h1()
    assert np.array_equal(D["genes"], genes_full)
    targets, sums = list(D["targets"]), D["sums"]
    TOT = H.load_totals()
    tm = TOT["targets"] != "non-targeting"
    assert list(TOT["targets"][tm]) == targets
    assert np.array_equal(TOT["totals"][tm].astype(np.int64), sums.astype(np.int64))
    Hd = H.load_H()
    hpos = {t: i for i, t in enumerate(Hd["targets"])}
    H0 = Hd["H"]
    assert list(map(str, Hd["genes"])) == list(map(str, genes_full))

    KEYS, KMETA = build_keys(H0, genes_full, keep, [str(t) for t in Hd["targets"]])
    ARMS = arms_for(KEYS)
    if os.environ.get("ARMS"):
        want = set(os.environ["ARMS"].split(","))
        ARMS = [a for a in ARMS if a[0] in want]
    names = [a[0] for a in ARMS]

    # -- the call-set confound, measured before anything is emitted --------------------
    head0 = np.zeros((len(targets), G), bool)
    base_abs = np.abs(KEYS["absH"][:, keep])
    for i in range(len(targets)):
        head0[i, np.argsort(-base_abs[i], kind="stable")[:N_HEAD]] = True
    headinfo = {}
    print("HEAD-SET MOVEMENT (rc_callset selects its call-set head as the top 400 by |key|)")
    print(f"{'key':<16}{'head jaccard vs |H|':>22}{'median rank corr in pool':>27}")
    for k in KEYS:
        a = np.abs(KEYS[k][:, keep])
        jj, sp = [], []
        for i in range(len(targets)):
            hd = np.zeros(G, bool)
            hd[np.argsort(-a[i], kind="stable")[:N_HEAD]] = True
            jj.append((hd & head0[i]).sum() / max((hd | head0[i]).sum(), 1))
        headinfo[k] = dict(head_jaccard=float(np.mean(jj)))
        print(f"{k:<16}{np.mean(jj):>22.5f}{'':>27}")

    if MAXT:
        targets, sums = targets[:MAXT], sums[:MAXT]
    T_n = len(targets)
    ref = np.load(REF191 / "de191.a15.npz", allow_pickle=True)
    assert list(map(str, ref["perts"]))[:T_n] == targets
    assert np.array_equal(ref["genes"], genes_full[keep])

    out = {a: {k: np.full((T_n, G), np.nan) for k in ("p", "q", "l", "z")} for a in names}
    npred = {a: [] for a in names}
    digest = {a: [] for a in names}
    diag = {a: [] for a in names}
    print(f"\nh1-test: {cs.shape[1]} controls, {G}/{G_full} gated, {T_n} perturbations; "
          f"{len(ARMS)} arms; call set top {N_HEAD} |key| U top {N_CPM} cpm = "
          f"{call_set.size}, z_cap {Z_CAP}; mass_fix None [{time.time()-t0:.0f}s]",
          flush=True)

    for i, t in enumerate(targets):
        Rec = recon196.Recon(sums[i], D["r"], D["gate"], n=n, seed=1000 + i)
        assert np.array_equal(Rec.totals(1.0), sums[i].astype(np.int64)), f"{t}: m=1 broken"
        Lp = Rec.emit_L(1.0)
        libs = np.full(n, Lp)
        k_art = hpos[t]
        x = Rec.x
        Tw = np.floor(x).astype(np.int64)
        prof = (Tw + 0.5) / (n * Lp)

        M1 = sums[i].astype(np.int64)
        pcpm = 1e6 * M1 / max(int(M1.sum()), 1)
        effr = np.abs(np.log2((pcpm + 1e-9) / (cpm + 1e-9)))
        rank = np.lexsort((np.arange(G_full), -np.where(keep, effr, -np.inf)))
        sv = np.full(G_full, 1.0)
        sv[rank[:min(KEEP_TOP, int(keep.sum()))]] = KEEP_SPREAD

        for name, kid, ks in ARMS:
            rng = np.random.default_rng(SEEDS[ks] + i)
            B, dg = RC.rc_hybrid_block(
                prof, libs, DISPERSION, rng, pass2_rng=rng, spread=sv,
                overdispersion=PHI, control=control, key=KEYS[kid][k_art],
                jitter_seed=1000 + i + 100_000 * ks,
                emit_counts=emit_counts, call_set=call_set)
            delivered = B.sum(0).astype(np.int64)
            assert np.array_equal(delivered, Tw), (name, t)
            digest[name].append(hashlib.sha256(delivered.tobytes()).hexdigest()[:16])
            zz, pv, qv, lfc = H.de_from_block(B, keep, cs, tie, ref_lin)
            out[name]["z"][i], out[name]["p"][i] = zz, pv
            out[name]["q"][i], out[name]["l"][i] = qv, lfc
            npred[name].append(int((qv < 0.05).sum()))
            diag[name].append(dg)
            del B
        if (i + 1) % 10 == 0 or i == 0:
            print(f"  {i+1}/{T_n} [{time.time()-t0:.0f}s]  " +
                  " ".join(f"{a}={npred[a][-1]}" for a in names[:3]), flush=True)

    # ---- known-answer row -----------------------------------------------------------
    ok = True
    if "absH" in names and KA_TABLE.exists():
        print("\nKNOWN-ANSWER ROW (rule 10): arm `absH` IS `rc_callset`.")
        z = np.load(KA_TABLE, allow_pickle=True)
        for mine, theirs in (("p", "pred_p"), ("q", "pred_padj"), ("l", "pred_lfc"),
                             ("z", "pred_z")):
            same = bool(np.array_equal(out["absH"][mine], z[theirs][:T_n]))
            worst = float(np.nanmax(np.abs(out["absH"][mine] - z[theirs][:T_n])))
            print(f"  {theirs:<10} bit-identical on {T_n} x {G}: "
                  f"{'YES' if same else 'NO'}   max |diff| {worst:.3g}")
            ok &= same
        sd = sum(a == b for a, b in zip(digest["absH"], list(map(str, z["digest"][:T_n]))))
        print(f"  per-gene panel-total digest identical on {sd}/{T_n}")
        ok &= sd == T_n
        print(f"  KNOWN-ANSWER ROW: {'PASS' if ok else 'FAIL'}")
    if not ok and not os.environ.get("NOABORT"):
        raise SystemExit("ABORT (rule 10): `absH` is not rc_callset")

    print("\nPER-GENE PANEL TOTAL DIGESTS vs `absH` (rule 13's class test)")
    base = digest.get("absH")
    if base:
        for a in names:
            same = sum(x == y for x, y in zip(digest[a], base))
            flag = "PROFILE-PRESERVING" if same == T_n else "PROFILE-CHANGING"
            print(f"  {a:<18} identical to absH on {same}/{T_n} blocks   {flag}")

    for a in names:
        np.savez_compressed(
            WORK / f"{TAG}.{a}.npz", genes=genes_full[keep], perts=np.array(targets),
            real_p=ref["real_p"][:T_n], real_padj=ref["real_padj"][:T_n],
            real_lfc=ref["real_lfc"][:T_n], pred_p=out[a]["p"], pred_padj=out[a]["q"],
            pred_lfc=out[a]["l"], pred_z=out[a]["z"], digest=np.array(digest[a]))
        print(f"wrote {TAG}.{a}.npz  median n_pred {np.median(npred[a]):.0f}", flush=True)
    json.dump(dict(npred=npred, arms=names, digest=digest, head=headinfo,
                   kmeta={k: {kk: str(vv) for kk, vv in v.items()}
                          for k, v in KMETA.items()},
                   diag={a: diag[a] for a in names}),
              open(WORK / f"{TAG}_diag.json", "w"), default=float)
    print(f"done [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
