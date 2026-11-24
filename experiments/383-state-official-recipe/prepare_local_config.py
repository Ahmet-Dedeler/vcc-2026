#!/usr/bin/env python3
"""Create a local-path copy of the frozen official STATE starter config."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


EXPECTED_SOURCE_SHA256 = "e9b3d0dc3be7b78002bce40e49d3a0f1f4d1ed17cf5af591786ade6538cf85a9"
COLAB_DATA_ROOT = "/content/state/competition_support_set"


def render_local_config(source: Path, support_dir: Path) -> str:
    source_bytes = source.read_bytes()
    observed_sha256 = hashlib.sha256(source_bytes).hexdigest()
    if observed_sha256 != EXPECTED_SOURCE_SHA256:
        raise ValueError(
            f"starter.toml SHA-256 mismatch: observed {observed_sha256}, "
            f"expected {EXPECTED_SOURCE_SHA256}"
        )

    source_text = source_bytes.decode("utf-8")
    if source_text.count(COLAB_DATA_ROOT) != 1:
        raise ValueError(
            f"expected exactly one frozen Colab data root {COLAB_DATA_ROOT!r}"
        )

    return source_text.replace(COLAB_DATA_ROOT, str(support_dir.resolve()))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("support_dir", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    rendered = render_local_config(args.source, args.support_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(args.output.resolve())


if __name__ == "__main__":
    main()
