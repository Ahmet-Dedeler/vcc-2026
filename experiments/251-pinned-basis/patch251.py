"""251 -- the builder change (and the def_postbuild change) as exact-anchor replacements.

Applies to git HEAD `89cb22a`'s `scripts/make_submission_pm_cd4.py` (builder sha256 f7a12a70..., =
248's `make_submission_pm_cd4.248.py`) and `scripts/def_postbuild.py` (97fcacf8...).  Every anchor
must occur EXACTLY once or the script refuses.  `--check <pre.py> <pre_postbuild.py>` applies the
patch to byte copies in memory and compares with the working tree (the reproducibility row).

    uv run python experiments/251-pinned-basis/patch251.py            # apply to the working tree
    uv run python experiments/251-pinned-basis/patch251.py --check experiments/251-pinned-basis/work/make_submission_pm_cd4.PRE.py experiments/251-pinned-basis/work/def_postbuild.PRE.py
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILDER = ROOT / "scripts/make_submission_pm_cd4.py"
POSTBUILD = ROOT / "scripts/def_postbuild.py"

ARTIFACT_SHA = "ba7ade3fc4f59de5c2de5d8bf270f243ab8e8ddfd0728be9e2c8c71c9b18122b"
V_SHA = "dabf88f7f95ffc4dac109de640ad29b943773ec4e45c5ddb2a4231b0003cce61"
U_SHA = "f06e0262cec2affb6c0594d2cebf8285888fd36ed73b518d514c7f0129b74e4c"

# ------------------------------------------------------------------------------ builder
B1_ANCHOR = "    CD4_DEFLATE_K = None\n# --- THE SOURCE REGISTRY (experiment 229)."
B1_NEW = '''    CD4_DEFLATE_K = None
# --- `VCC_CD4_DEFLATE_BASIS` -- deflate the CD4 residual with a PINNED gene-space basis (250/251).
#
# WHY (250 sections 4.3-4.5).  `VCC_CD4_DEFLATE_K` above FITS its projector on the rows it is
# applied to.  On A/B/C that is the object 249 built (`cd4svd5`).  On a disjoint 300-target
# panel 250 measured that a day-of re-fit removes panel-specific signal as well as the shared
# T-cell hub: its aggregate-cosine cost is 4-14x A/B/C's (-2.1 / -6.6 / -5.2 / -3.5% on 199's
# and 238's four draws against -0.47%) and its six-member price is <= 0 on both low-coverage
# draws.  The SAME projector expressed in GENE space -- `V_5 = R_c^T U_5 / s_5` from A/B/C's
# own production residual, the object whose `U_5` sha is `f06e0262...` in 249's shipped audit
# -- applied after centring on the new panel's own carried rows beats the re-fit on 4/4
# disjoint panels: `pds` +0.0153 / +0.0037 / +0.0073 / +0.0107 at only -0.1 ... -2.0% cosine,
# price +0.0061 / +0.0009 / +0.0014 / +0.0030, never negative.  A basis fitted on OTHER rows
# removes the hub and nothing else.  So October pins; it never re-fits.
#
# THE CONSTRUCTION (250 `lib250.deflate_gene_basis`, reproduced verbatim in
# `deflate_cd4_residual_pinned` below).  On the carried rows (`target_condition_count > 0`):
# subtract the per-gene mean over those rows, re-zero outside the gene support,
# `R_k = R_c - (R_c V) V^T`, re-zero outside the support.  No fit, no eigh, no seed: the only
# panel-dependent quantity is the centring.  On A/B/C's own rows this is 249's projector to
# floating-point order (251 measured max |d| 7.8e-15 on entries up to 2.8 -- the two matrix
# products are associated differently), so the shipped re-fit path is NOT rewritten in terms
# of it: `cd4svd5`'s bytes stay `cd4svd5`'s.
#
# THE FLAG.  Unset (the default) is BYTE-INERT: `apply_frozen_cd4_residual` takes the branch
# below only when the name is set, and every prior build reproduces bit-for-bit (251's
# known-answer rows: the staged effects `abeb0749...` and the package `8df4dadd...`).  Set to a
# REGISTRY NAME and the CD4 residual is replaced by its pinned-deflated form BEFORE
# `solve_cd4_dose` sizes it (247's carrying-row rule off A/B/C, the all-rows rule on A/B/C), so
# the sizing rule sizes the object actually added -- 249's "resolved" construction.  Rules, at
# import: `VCC_CD4_DEFLATE_K` and `VCC_CD4_DEFLATE_BASIS` are TWO variables (rule 10) and may
# not both be set; the name must be in `CD4_DEFLATION_BASIS_REGISTRY`; and at load the file must
# hash to its pin, its `V` must hash to its pin, its provenance `U_k` sha must be the pin's, its
# gene axis must be the bundle's official axis, and its columns must be orthonormal.
# `effects.cd4.deflation` (mode `pinned`, the name, path and shas, `energy_removed_frac`, and
# `captured_frac_of_centred` -- the fraction of the day's centred CD4 energy the pinned basis
# captures, 250's census diagnostic: 0.10-0.14 on the four disjoint panels, 0.198 on A/B/C) and
# `config.cd4_deflate_basis` are written only when on (the `h1_union_g` convention).
#
# THE REGISTRY.  Rule 11: a path AND a sha, never a description.  `path` is relative to the
# repository root and every value is a plain string so `scripts/def_postbuild.py` can read
# this table as a literal out of the source, as it reads `WHITENER_REGISTRY`.
CD4_DEFLATION_BASIS_REGISTRY: dict[str, dict[str, str]] = {
    # 251: the gene-space top-5 basis of 249's staged A/B/C CD4 residual (the production
    # complete-mean3 object, 300 rows, 291 carried; `experiments/249-cd4-deflate/work/stage/
    # cd4-off.npz`), == 250's `basis250-production.npz[V5_abc]` bit for bit; its row-space
    # `U_5` hashes to the `basis_u_k_sha256` in `submissions/249-cd4svd5-f034.h5ad.meta.json`.
    # Written by `experiments/251-pinned-basis/basis251.py`; provenance in the `.json` sidecar.
    "V5_abc": {
        "path": "data/donors/cd4_deflation_basis/V5_abc.npz",
        "sha256": "@ARTIFACT_SHA@",
        "v_sha256": "@V_SHA@",
        "u_k_sha256": "@U_SHA@",
        "k": "5",
    },
}
CD4_DEFLATE_BASIS = os.environ.get("VCC_CD4_DEFLATE_BASIS", "").strip() or None
if CD4_DEFLATE_BASIS is not None:
    if CD4_DEFLATE_K is not None:
        raise ValueError(
            "VCC_CD4_DEFLATE_K and VCC_CD4_DEFLATE_BASIS are both set: a day-of re-fit and a "
            "pinned basis are two variables (rule 10); set exactly one"
        )
    if CD4_DEFLATE_BASIS not in CD4_DEFLATION_BASIS_REGISTRY:
        raise ValueError(
            f"VCC_CD4_DEFLATE_BASIS names {CD4_DEFLATE_BASIS!r}, which is not in "
            f"CD4_DEFLATION_BASIS_REGISTRY ({sorted(CD4_DEFLATION_BASIS_REGISTRY)}).  Add it "
            "there, with its path and shas, first."
        )
# --- THE SOURCE REGISTRY (experiment 229).'''.replace("@ARTIFACT_SHA@", ARTIFACT_SHA).replace("@V_SHA@", V_SHA).replace("@U_SHA@", U_SHA)

B2_ANCHOR = "\n\ndef single_source_residual(\n"
B2_NEW = '''

_CD4_DEFLATION_BASIS_CACHE: dict[str, tuple[np.ndarray, dict[str, object], tuple[str, ...]]] = {}


def load_cd4_deflation_basis(
    name: str, *, official_genes: Sequence[str] | None = None
) -> tuple[np.ndarray, dict[str, object]]:
    """Load a registered pinned basis (251), fail-closed on every pin it carries.

    Refuses, in order: an unregistered name; a missing file; a file whose sha256 is not the
    registry's; a `k` that is not the registry's; a `V` that is not (n_genes, k), finite and
    orthonormal; a `V` whose array sha256 is not the registry's; a provenance `U_k` sha that is
    not the registry's; and -- when `official_genes` is given (`run()` passes the bundle's axis
    before the three-donor model is loaded) -- a gene axis that is not this bundle's.  The
    array is cached per name; EVERY call that names an axis verifies it (a cache verified
    against one bundle says nothing about another).
    """
    if name not in CD4_DEFLATION_BASIS_REGISTRY:
        raise ValueError(
            f"{name!r} is not in CD4_DEFLATION_BASIS_REGISTRY "
            f"({sorted(CD4_DEFLATION_BASIS_REGISTRY)})"
        )
    entry = CD4_DEFLATION_BASIS_REGISTRY[name]
    cached = _CD4_DEFLATION_BASIS_CACHE.get(name)
    if cached is not None:
        basis, info, genes = cached
        if official_genes is not None:
            # every call that names an axis verifies it against the basis's own -- a cache
            # verified against one bundle says nothing about another
            _verify_cd4_deflation_basis_axis(name, genes, official_genes)
            info = {**info, "gene_axis_verified": True}
            _CD4_DEFLATION_BASIS_CACHE[name] = (basis, info, genes)
        return basis, info
    path = ROOT / str(entry["path"])
    if not path.is_file():
        raise FileNotFoundError(f"CD4 deflation basis {name!r} is missing: {path}")
    digest = file_sha256(path)
    if digest != entry["sha256"]:
        raise ValueError(
            f"CD4 deflation basis {name!r} at {path} does not match "
            f"CD4_DEFLATION_BASIS_REGISTRY: {digest} != {entry['sha256']}"
        )
    with np.load(path, allow_pickle=False) as archive:
        basis = np.ascontiguousarray(archive["V"], dtype=np.float64)
        genes = tuple(str(gene) for gene in archive["genes"])
        k = int(archive["k"])
        provenance = json.loads(str(archive["provenance"]))
    if k != int(entry["k"]):
        raise ValueError(f"CD4 deflation basis {name!r}: k {k} != the registry's {entry['k']}")
    if basis.ndim != 2 or basis.shape != (len(genes), k) or k < 1:
        raise ValueError(
            f"CD4 deflation basis {name!r}: V has shape {basis.shape}; expected "
            f"({len(genes)}, {k})"
        )
    if not np.all(np.isfinite(basis)):
        raise ValueError(f"CD4 deflation basis {name!r}: V has non-finite entries")
    gram_error = float(np.abs(basis.T @ basis - np.eye(k)).max())
    if gram_error > 1e-9:
        raise ValueError(
            f"CD4 deflation basis {name!r}: columns are not orthonormal (max |V^T V - I| "
            f"{gram_error:.3e})"
        )
    v_sha = array_sha256(basis)
    if v_sha != entry["v_sha256"]:
        raise ValueError(
            f"CD4 deflation basis {name!r}: V hashes to {v_sha}, not the registry's "
            f"{entry['v_sha256']}"
        )
    if str(provenance.get("basis_u_k_sha256")) != entry["u_k_sha256"]:
        raise ValueError(
            f"CD4 deflation basis {name!r}: provenance U_k sha "
            f"{provenance.get('basis_u_k_sha256')!r} != the registry's {entry['u_k_sha256']}"
        )
    axis_verified = False
    if official_genes is not None:
        _verify_cd4_deflation_basis_axis(name, genes, official_genes)
        axis_verified = True
    info: dict[str, object] = {
        "mode": "pinned",
        "basis": name,
        "basis_path": str(entry["path"]),
        "basis_sha256": digest,
        "basis_v_sha256": v_sha,
        "basis_u_k_sha256": str(entry["u_k_sha256"]),
        "k": k,
        "n_genes": len(genes),
        "orthonormality_err": gram_error,
        "gene_axis_verified": axis_verified,
        "provenance": {
            key: provenance[key]
            for key in ("experiment", "date", "fitted_on", "source")
            if key in provenance
        },
    }
    _CD4_DEFLATION_BASIS_CACHE[name] = (basis, info, genes)
    return basis, info


def _verify_cd4_deflation_basis_axis(name: str, genes: tuple[str, ...], official_genes) -> None:
    """RAISE unless the basis's gene axis IS the bundle's official axis, symbol for symbol."""
    official = tuple(map(str, official_genes))
    if official != genes:
        n_common = sum(1 for a, b in zip(official, genes) if a == b)
        raise ValueError(
            f"CD4 deflation basis {name!r} was fitted on a {len(genes)}-gene axis that is not "
            f"this bundle's official axis ({len(official)} genes, {n_common} positions agree): "
            "the basis is a function of the gene ORDER and cannot be applied here"
        )


def deflate_cd4_residual_pinned(cd4, basis: np.ndarray, *, basis_info: Mapping[str, object]):
    """250's `lib250.deflate_gene_basis`, verbatim: deflate with a gene-space basis fitted ELSEWHERE.

    On the CARRIED rows (`target_condition_count > 0`, exactly `deflate_cd4_residual`'s rows):

        R   = residual on the carried rows, masked to its gene support
        R_c = R - mean over rows(R)          per-gene centring on THIS panel's rows, re-zeroed
                                             outside the support
        P   = R_c V                          (n_carried x k)
        R_k = R_c - P V^T                    re-zeroed outside the support

    Nothing is fitted here.  `V` carries coordinates on genes that are panel targets today (they
    were ordinary genes on the panel it was fitted on); the final re-zeroing outside the day's
    support removes them, and the caller asserts the panel coordinates are still zero.

    The audit carries the census diagnostic 250 section 6.2 names: `captured_frac_of_centred`,
    the fraction of the day's centred CD4 energy inside the pinned subspace (0.10-0.14 on the
    four disjoint panels, 0.198 on A/B/C; a reading < 0.05 or > 0.19 is unmeasured territory),
    plus the day's OWN top-k share (`own_top_k_energy_frac_of_centred`, which 250 measured NOT
    to discriminate) for the record.  Returns a `CD4Prepared` with the SAME mode, votes, counts
    and condition details, as `deflate_cd4_residual` does.
    """
    V = np.asarray(basis, dtype=np.float64)
    eff = np.asarray(cd4.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(cd4.gene_votes)
    sup = votes > 0
    eff[~sup] = 0.0
    if V.ndim != 2 or V.shape[0] != eff.shape[1] or V.shape[1] < 1:
        raise ValueError(
            f"pinned basis must be (n_genes={eff.shape[1]}, k >= 1); got {V.shape}"
        )
    k = int(V.shape[1])
    rows = np.flatnonzero(np.asarray(cd4.target_condition_count) > 0)
    if rows.size == 0:
        raise ValueError(
            "VCC_CD4_DEFLATE_BASIS is set but the CD4 residual carries no rows; nothing to deflate"
        )
    R = eff[rows]
    e0 = float((R * R).sum())
    if e0 <= 0.0:
        raise ValueError(
            "VCC_CD4_DEFLATE_BASIS is set but the CD4 residual has zero energy; nothing to deflate"
        )
    mean = R.mean(axis=0, keepdims=True)
    Rc = R - mean
    Rc[~sup[rows]] = 0.0
    ec = float((Rc * Rc).sum())
    proj = Rc @ V
    captured = float((proj * proj).sum() / max(ec, 1e-300))
    Rk = Rc - proj @ V.T
    Rk[~sup[rows]] = 0.0
    e1 = float((Rk * Rk).sum())
    # The day's own spectrum, for the record only: 250 section 6.2 measured the top-k share at
    # 0.198-0.232 on every panel, transferring or not, so it is recorded and never gated on.
    w = np.linalg.eigvalsh(Rc @ Rc.T)
    own_top_k = float(np.sort(w)[::-1][:k].sum() / max(w.sum(), 1e-300))
    audit: dict[str, object] = {
        **{key: value for key, value in basis_info.items() if key != "provenance"},
        "n_carried": int(rows.size),
        "energy_removed_frac": float(1.0 - e1 / e0),
        "captured_frac_of_centred": captured,
        "own_top_k_energy_frac_of_centred": own_top_k,
        "centred_energy_frac_of_raw": float(ec / e0),
        **({} if "provenance" not in basis_info else {"provenance": basis_info["provenance"]}),
    }
    eff[rows] = Rk
    return (
        EXP062.CD4Prepared(
            mode=str(cd4.mode),
            mean_effect=eff,
            gene_votes=np.asarray(cd4.gene_votes).copy(),
            target_condition_count=np.asarray(cd4.target_condition_count).copy(),
            condition_details=cd4.condition_details,
        ),
        audit,
    )


def single_source_residual(
'''

B3_ANCHOR = '''            f"U_k sha256 {str(deflation_audit['basis_u_k_sha256'])[:16]}",
            flush=True,
        )
    fraction, frame = solve_cd4_dose(whitened, cd4)
'''
B3_NEW = '''            f"U_k sha256 {str(deflation_audit['basis_u_k_sha256'])[:16]}",
            flush=True,
        )
    elif CD4_DEFLATE_BASIS is not None:
        # 251: the PINNED gene-space basis (250 sections 4.3-4.5) -- nothing is fitted on this
        # panel's rows.  Mutually exclusive with `VCC_CD4_DEFLATE_K` at import, so exactly one
        # of the two branches can run; unset, neither does and `cd4` is untouched.
        basis, basis_info = load_cd4_deflation_basis(CD4_DEFLATE_BASIS)
        cd4, deflation_audit = deflate_cd4_residual_pinned(cd4, basis, basis_info=basis_info)
        if np.any(cd4.mean_effect[:, panel_positions] != 0):
            raise AssertionError(
                "pinned deflation put a VCC panel coordinate into the CD4 residual"
            )
        print(
            f"[cd4] deflation pinned {CD4_DEFLATE_BASIS} (file sha256 "
            f"{str(deflation_audit['basis_sha256'])[:16]}, k={deflation_audit['k']}): "
            f"{deflation_audit['n_carried']} carried rows, energy removed "
            f"{deflation_audit['energy_removed_frac']:.6f}, captured "
            f"{deflation_audit['captured_frac_of_centred']:.6f} of the day's centred CD4 "
            f"energy (250: 0.10-0.14 on disjoint panels, 0.198 on A/B/C; < 0.05 or > 0.19 is "
            f"unmeasured), own top-k share "
            f"{deflation_audit['own_top_k_energy_frac_of_centred']:.6f}",
            flush=True,
        )
    fraction, frame = solve_cd4_dose(whitened, cd4)
'''

B4_ANCHOR = '''            # Written only when the CD4 residual is deflated (249), the same convention.
            **({} if CD4_DEFLATE_K is None else {"cd4_deflate_k": CD4_DEFLATE_K}),
'''
B4_NEW = '''            # Written only when the CD4 residual is deflated (249), the same convention.
            **({} if CD4_DEFLATE_K is None else {"cd4_deflate_k": CD4_DEFLATE_K}),
            # Written only when the deflation uses a PINNED basis (251), the same convention.
            **({} if CD4_DEFLATE_BASIS is None else {"cd4_deflate_basis": CD4_DEFLATE_BASIS}),
'''

B5_ANCHOR = '''    donors = {name: Donor.load(path) for name, path in donor_paths.items()}
    model = MultiDonorModel(
'''
B5_NEW = '''    if CD4_DEFLATE_BASIS is not None:
        # 251: the pinned basis is hashed against its pins and its gene axis is checked against
        # THIS bundle's official axis BEFORE the three-donor model is loaded, so a wrong file
        # fails in seconds and not after the 17 GiB stage.
        _basis, basis_info = load_cd4_deflation_basis(
            CD4_DEFLATE_BASIS, official_genes=tuple(map(str, spec.genes))
        )
        print(
            f"[cd4-basis] {CD4_DEFLATE_BASIS}: {basis_info['basis_path']} sha256 "
            f"{basis_info['basis_sha256']} (V {str(basis_info['basis_v_sha256'])[:16]}, U_k "
            f"{str(basis_info['basis_u_k_sha256'])[:16]}), k={basis_info['k']}, gene axis "
            f"verified against the bundle ({basis_info['n_genes']} genes)",
            flush=True,
        )
        del _basis
    donors = {name: Donor.load(path) for name, path in donor_paths.items()}
    model = MultiDonorModel(
'''

BUILDER_EDITS = [(B1_ANCHOR, B1_NEW), (B2_ANCHOR, B2_NEW), (B3_ANCHOR, B3_NEW),
                 (B4_ANCHOR, B4_NEW), (B5_ANCHOR, B5_NEW)]

# ------------------------------------------------------------------------------ def_postbuild
P1_ANCHOR = '''def builder_literal(name: str, source: Path = BUILDER) -> dict:
'''
P1_NEW = '''def builder_literal_or_empty(name: str, source: Path = BUILDER) -> dict:
    """`builder_literal`, but `{}` when the builder predates the table (251: the pinned-basis
    registry exists from `V5_abc` on; a build carrying `cd4_deflate_basis` against an older
    builder source then fails the 'is registered' row instead of crashing the check)."""
    try:
        return builder_literal(name, source)
    except SystemExit:
        return {}


def builder_literal(name: str, source: Path = BUILDER) -> dict:
'''

P2_ANCHOR = '''    # --- 240: the arms the runbook says we ship, BY VALUE ---------------------------------
'''
P2_NEW = '''    # --- 251: a PINNED CD4 deflation basis, if (and only if) the build carries one ---------
    #
    # `VCC_CD4_DEFLATE_BASIS=<name>` writes `config.cd4_deflate_basis` and an
    # `effects.cd4.deflation` block with `mode == "pinned"`.  Rule 10 would refuse the new config
    # key as an unintended variable; it is accepted here ONLY when every pin it carries checks
    # out against the builder's `CD4_DEFLATION_BASIS_REGISTRY` (read as a literal, never
    # imported) and the file on disk, and the census diagnostic 250 section 6.2 names is echoed.
    # A build without the key sees none of these rows, so 222/249/244's metas are untouched.
    extra_allowed: set[str] = set()
    basis_name = cfg.get("cd4_deflate_basis")
    defl = meta["effects"]["cd4"].get("deflation") or {}
    if basis_name is not None:
        basis_registry = builder_literal_or_empty("CD4_DEFLATION_BASIS_REGISTRY", args.builder)
        entry = basis_registry.get(basis_name)
        chk("config.cd4_deflate_basis is in the builder's CD4_DEFLATION_BASIS_REGISTRY",
            entry is not None, f"{basis_name!r} (registered: {sorted(basis_registry)})")
        chk("config.cd4_deflate_k ABSENT while a pinned basis is on (one deflation variable, rule 10)",
            "cd4_deflate_k" not in cfg, f"cd4_deflate_k={cfg.get('cd4_deflate_k', '<absent>')!r}")
        chk("effects.cd4.deflation.mode == 'pinned' and .basis == config.cd4_deflate_basis",
            defl.get("mode") == "pinned" and defl.get("basis") == basis_name,
            f"mode={defl.get('mode')!r} basis={defl.get('basis')!r}")
        chk(f"effects.cd4.deflation.basis_sha256 == CD4_DEFLATION_BASIS_REGISTRY[{basis_name!r}]['sha256'] (the pin)",
            entry is not None and defl.get("basis_sha256") == entry.get("sha256"),
            f"{str(defl.get('basis_sha256'))[:16]}… vs pin {entry['sha256'][:16] + '…' if entry else None}")
        chk(f"effects.cd4.deflation.basis_v_sha256 == the registry's v_sha256 and .k == its k",
            entry is not None and defl.get("basis_v_sha256") == entry.get("v_sha256")
            and str(defl.get("k")) == str(entry.get("k")),
            f"V {str(defl.get('basis_v_sha256'))[:16]}… k={defl.get('k')!r} vs pin "
            f"{entry['v_sha256'][:16] + '…' if entry else None} k={entry.get('k') if entry else None!r}")
        bpath = (ROOT / str(entry["path"])) if entry else None
        if bpath is not None and bpath.is_file():
            on_disk_b = sha256_file(bpath)
            chk("sha256(the file at the basis pin's path) == the pin (read the pin, hash the file)",
                on_disk_b == entry["sha256"], f"{on_disk_b[:16]}… vs pin {entry['sha256'][:16]}… ({bpath.name})")
        else:
            chk("the pinned basis file exists at the registry's path", False, f"MISSING {bpath}")
        chk("effects.cd4.deflation.gene_axis_verified (the builder checked the basis's gene ORDER against this bundle)",
            defl.get("gene_axis_verified") is True, f"{defl.get('gene_axis_verified')!r}")
        captured = defl.get("captured_frac_of_centred")
        chk("effects.cd4.deflation.captured_frac_of_centred is a number in (0, 1) (250's census diagnostic, echoed below)",
            isinstance(captured, (int, float)) and 0.0 < float(captured) < 1.0, f"{captured!r}")
        extra_allowed.add("cd4_deflate_basis")
    elif defl.get("mode") == "pinned":
        chk("effects.cd4.deflation.mode == 'pinned' requires config.cd4_deflate_basis", False,
            "the audit says pinned but the config carries no basis name")

    # --- 240: the arms the runbook says we ship, BY VALUE ---------------------------------
'''

P3_ANCHOR = '''    if new_panel:
        allowed |= ALLOWED_CONFIG_DIFFS_NEW_PANEL
    moved = sorted(k for k in set(meta["config"]) | set(inc["config"])
'''
P3_NEW = '''    if new_panel:
        allowed |= ALLOWED_CONFIG_DIFFS_NEW_PANEL
    allowed |= extra_allowed          # 251: `cd4_deflate_basis`, only when its rows above ran
    moved = sorted(k for k in set(meta["config"]) | set(inc["config"])
'''

P4_ANCHOR = '''    else:
        print(f"cd4       all-rows rule f={cfg.get('cd4_target_f')} solved {rf.get('solved_fraction')} "
              f"({got_cd4} complete targets = {got_cd4 / manifest['n_constructs']:.1%} coverage)")
    print()
'''
P4_NEW = '''    else:
        print(f"cd4       all-rows rule f={cfg.get('cd4_target_f')} solved {rf.get('solved_fraction')} "
              f"({got_cd4} complete targets = {got_cd4 / manifest['n_constructs']:.1%} coverage)")
    if basis_name is not None:
        # 251: the census diagnostic, echoed off the meta (250 section 6.2: 0.10-0.14 on the four
        # disjoint panels, 0.198 on A/B/C; < 0.05 or > 0.19 is unmeasured territory; the day's own
        # top-k share 0.198-0.232 everywhere and NOT a gate).
        print(f"cd4       deflation PINNED {basis_name!r} sha256 {str(defl.get('basis_sha256'))[:16]}… k={defl.get('k')}: "
              f"{defl.get('n_carried')} carried rows, energy removed {defl.get('energy_removed_frac')!r}, "
              f"captured {defl.get('captured_frac_of_centred')!r} of the day's centred CD4 energy "
              f"(250: 0.10-0.14 disjoint, 0.198 A/B/C), own top-k share "
              f"{defl.get('own_top_k_energy_frac_of_centred')!r} (recorded, not gated)")
    print()
'''

POSTBUILD_EDITS = [(P1_ANCHOR, P1_NEW), (P2_ANCHOR, P2_NEW), (P3_ANCHOR, P3_NEW), (P4_ANCHOR, P4_NEW)]


def apply(src: str, edits, label: str) -> str:
    out = src
    for i, (anchor, new) in enumerate(edits, 1):
        n = out.count(anchor)
        if n != 1:
            raise SystemExit(f"{label}: anchor {i} found {n} times; refusing")
        out = out.replace(anchor, new)
    return out


def main() -> None:
    if len(sys.argv) >= 2 and sys.argv[1] == "--check":
        pre_b = Path(sys.argv[2]).read_text()
        pre_p = Path(sys.argv[3]).read_text()
        want_b = apply(pre_b, BUILDER_EDITS, "builder")
        want_p = apply(pre_p, POSTBUILD_EDITS, "def_postbuild")
        have_b = BUILDER.read_text()
        have_p = POSTBUILD.read_text()
        for name, want, have in (("builder", want_b, have_b), ("def_postbuild", want_p, have_p)):
            same = want == have
            print(f"{name}: patched byte copy {'==' if same else '!='} working tree  "
                  f"(sha256 {hashlib.sha256(have.encode()).hexdigest()[:16]}…)")
        raise SystemExit(0 if (want_b == have_b and want_p == have_p) else 1)
    for path, edits, label in ((BUILDER, BUILDER_EDITS, "builder"), (POSTBUILD, POSTBUILD_EDITS, "def_postbuild")):
        src = path.read_text()
        new = apply(src, edits, label)
        compile(new, str(path), "exec")
        path.write_text(new)
        print(f"{label}: applied {len(edits)} edits -> {path.relative_to(ROOT)} sha256 "
              f"{hashlib.sha256(new.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
