"""180 gate 4 — the built knee package's meta against a23's, key by key.

Rule 10: diff the finished artifact's config against the incumbent's and require EXACTLY the
intended number of differences. Not the command you typed -- the meta the build actually
wrote. The intended set here is the alpha keys and nothing else; everything outside `config`
that moved is listed rather than hidden.

    uv run python metadiff180.py [built.meta.json] [incumbent.meta.json]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "work"

EXPECTED_CONFIG_KEYS = {"alpha_mode", "alpha_ref", "alpha_per_context"}


def flatten(node, prefix=""):
    out = {}
    if isinstance(node, dict):
        for key, value in node.items():
            out.update(flatten(value, f"{prefix}.{key}" if prefix else str(key)))
    elif isinstance(node, list):
        out[prefix] = json.dumps(node)
    else:
        out[prefix] = node
    return out


def main(built_path, incumbent_path) -> int:
    built = json.loads(Path(built_path).read_text())
    incumbent = json.loads(Path(incumbent_path).read_text())

    got, want = flatten(built["config"]), flatten(incumbent["config"])
    moved = [k for k in sorted(set(got) | set(want))
             if got.get(k, "<absent>") != want.get(k, "<absent>")]
    top = sorted({k.split(".")[0] for k in moved})
    print("[config] keys that moved: " + (", ".join(top) or "<none>"))
    for key in moved:
        print(f"[config]   {key}: built {got.get(key, '<absent>')!r}  "
              f"incumbent {want.get(key, '<absent>')!r}")
    config_ok = set(top) == EXPECTED_CONFIG_KEYS
    print(f"[config] expected exactly {sorted(EXPECTED_CONFIG_KEYS)} -> "
          f"{'PASS' if config_ok else 'FAIL'}")

    print("\n[effects] the effect audit, which must be a23's own (the knee is applied by "
          "REBUILDING per context; `effects` stays the reference-alpha build):")
    ge, we = flatten(built["effects"], "effects"), flatten(incumbent["effects"], "effects")
    ediff = [k for k in sorted(set(ge) | set(we))
             if ge.get(k, "<absent>") != we.get(k, "<absent>")]
    for key in ediff:
        print(f"[effects]   {key}: built {ge.get(key, '<absent>')!r}  "
              f"incumbent {we.get(key, '<absent>')!r}")
    print(f"[effects] fields {len(set(ge) | set(we))}, differences {len(ediff)}  "
          f"final_sha256 {built['effects']['final_sha256'][:24]} vs "
          f"{incumbent['effects']['final_sha256'][:24]}")

    print("\n[other] blocks outside `config`/`effects` that moved:")
    other = []
    for block in ("operator", "hashes", "output", "donor_coverage", "main_source_support",
                  "runtime", "frozen_bundle", "knee", "peak_rss_gib", "model", "status",
                  "submission"):
        gb = flatten(built.get(block, {}), block)
        wb = flatten(incumbent.get(block, {}), block)
        for key in sorted(set(gb) | set(wb)):
            if gb.get(key, "<absent>") != wb.get(key, "<absent>"):
                other.append(key)
    for key in other:
        print(f"[other]   {key}")

    if "knee" in built:
        print("\n[knee] per-context terms recorded in the meta:")
        for context, row in built["knee"]["per_context"].items():
            print(f"[knee]   {context}: alpha {row['alpha']:.6f}  "
                  f"A {row['A']:.6f}  B {row['B']:.6f}  claim {row['claim']:.6f}  "
                  f"C_p {row['C_p']:.4f} [{row['C_p_source']}]  "
                  f"rho {row['rho_unclamped']:.6f} -> {row['after']['rho']:.6f}  "
                  f"effect_sha {row['effect_sha256'][:16]}  "
                  f"homogeneity_residual {row['homogeneity_residual']:.2e}")

    (WORK / "metadiff180.json").write_text(json.dumps({
        "config_keys_moved": top, "config_ok": config_ok,
        "config_diffs": {k: [got.get(k, "<absent>"), want.get(k, "<absent>")] for k in moved},
        "effects_differences": ediff,
        "other_keys_moved": other,
    }, indent=1))
    print(f"\n[gate4] {'PASS' if config_ok and not ediff else 'FAIL'}")
    return 0 if config_ok and not ediff else 1


if __name__ == "__main__":
    a = sys.argv[1] if len(sys.argv) > 1 else str(WORK / "abc-knee.h5ad.meta.json")
    b = sys.argv[2] if len(sys.argv) > 2 else str(REPO / "submissions/a23-massfix.h5ad.meta.json")
    sys.exit(main(a, b))
