"""240 row 1 — did the DRIVER's defaults rebuild the shipped `h1union` package?

Compares the artifact `scripts/def_pipeline.py --stages build` produced with nothing exported
against the shipped `submissions/222-h1union-g028.h5ad.meta.json`: the file digests, the whole
config block field by field, and the quantities the mandate names.

    uv run python experiments/240-driver-fix/verify240.py
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BUILT = HERE / "work/ka1/240-ka-h1union.h5ad"
SHIPPED_META = ROOT / "submissions/222-h1union-g028.h5ad.meta.json"
# The shipped package, from 222 and reproduced by 229 row 2 and 235 row 2.
WANT = {"h5ad_sha256": "8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6",
        "h5ad_bytes": 12_648_807_800,
        "vcc_sha256": "f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e",
        "vcc_bytes": 2_824_785_920}
# The config-block quantities the mandate names, and where they live in the meta.
NAMED = {
    "effects.cd4.fraction": ("effects", "cd4", "fraction"),
    "effects.cd4.residual_frame.target_f_requested": ("effects", "cd4", "residual_frame",
                                                      "target_f_requested"),
    "effects.cd4.residual_frame.f_per_unit_fraction": ("effects", "cd4", "residual_frame",
                                                       "f_per_unit_fraction"),
    "effects.cd4.fully_supported_targets": ("effects", "cd4", "fully_supported_targets"),
    "operator.sha256": ("operator", "sha256"),
    "config.whiten_gamma": ("config", "whiten_gamma"),
    "config.whiten_lambda": ("config", "whiten_lambda"),
    "config.alpha": ("config", "alpha"),
    "config.cd4_target_f": ("config", "cd4_target_f"),
    "config.h1_union_g": ("config", "h1_union_g"),
    "config.h1_union_sources": ("config", "h1_union_sources"),
    "config.emitter": ("config", "emitter"),
    "effects.h1_union.carrying_targets": ("effects", "h1_union", "carrying_targets"),
    "effects.h1_union.residual_frame.applied_fraction": ("effects", "h1_union", "residual_frame",
                                                         "applied_fraction"),
}
# Config keys the driver moves ON PURPOSE relative to a hand-run build222.sh (both are recorded
# in the meta only; neither reaches the emitted counts).
INTENDED = {"expect_contexts"}


def dig(d, path):
    for key in path:
        if not isinstance(d, dict) or key not in d:
            return "<absent>"
        d = d[key]
    return d


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    built_meta = json.loads(Path(str(BUILT) + ".meta.json").read_text())
    ship = json.loads(SHIPPED_META.read_text())
    out: dict[str, object] = {}

    out["digests"] = {"h5ad_bytes": BUILT.stat().st_size, "h5ad_sha256": sha256(BUILT),
                      "want_h5ad_sha256": WANT["h5ad_sha256"], "want_h5ad_bytes": WANT["h5ad_bytes"]}
    vcc = BUILT.with_suffix(".vcc")
    if vcc.exists():
        out["digests"].update({"vcc_bytes": vcc.stat().st_size, "vcc_sha256": sha256(vcc),
                               "want_vcc_sha256": WANT["vcc_sha256"],
                               "want_vcc_bytes": WANT["vcc_bytes"]})

    moved = sorted(k for k in set(built_meta["config"]) | set(ship["config"])
                   if built_meta["config"].get(k, "<absent>") != ship["config"].get(k, "<absent>"))
    out["config_diff_vs_shipped"] = {k: [ship["config"].get(k, "<absent>"),
                                         built_meta["config"].get(k, "<absent>")] for k in moved}
    out["config_diff_is_only_intended"] = set(moved) <= INTENDED
    out["named"] = {k: {"built": dig(built_meta, p), "shipped": dig(ship, p)}
                    for k, p in NAMED.items()}

    print(f"built   {BUILT}")
    print(f"shipped {SHIPPED_META}\n")
    d = out["digests"]
    for kind in ("h5ad", "vcc"):
        if f"{kind}_sha256" in d:
            ok = d[f"{kind}_sha256"] == d[f"want_{kind}_sha256"] and d[f"{kind}_bytes"] == d[f"want_{kind}_bytes"]
            print(f"  {'PASS' if ok else 'FAIL'}  .{kind}  {d[f'{kind}_sha256']}  {d[f'{kind}_bytes']:,} bytes")
            if not ok:
                print(f"        want  {d[f'want_{kind}_sha256']}  {d[f'want_{kind}_bytes']:,} bytes")
    print(f"\n  config differences vs the shipped package: {moved or 'NONE'}"
          f"  -> {'only the intended ones' if out['config_diff_is_only_intended'] else 'UNEXPECTED'}")
    for k, v in out["config_diff_vs_shipped"].items():
        print(f"      {k}: shipped {v[0]!r} -> built {v[1]!r}")
    print("\n  the named quantities (built vs shipped):")
    for k, v in out["named"].items():
        mark = "OK " if v["built"] == v["shipped"] else "NO "
        print(f"    {mark} {k:<48} {v['built']!r}")
    (HERE / "work/verify240.json").write_text(json.dumps(out, indent=1, default=str) + "\n")
    print(f"\nwrote {HERE / 'work/verify240.json'}")


if __name__ == "__main__":
    main()
