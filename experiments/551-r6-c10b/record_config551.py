"""551 (C10b) -- record the final bulk configuration F BEFORE the access (PREREGISTRATION section 2).

Writes work/c10b-config.json (the coordinator's message verbatim, when it arrived, who recorded it, the resolved
designs551.CONFIGS key) and appends one line to docs/CONFIRMATION-LEDGER.md's C10b section with the sha256 of that file
and of PREREGISTRATION.md. Refuses to overwrite an existing record, and refuses a configuration outside the six.

    uv run python experiments/551-r6-c10b/record_config551.py --config M6+S_P2 --received 2026-10-04T04:31Z \
        --message-file /path/to/message.txt
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import designs551 as D  # noqa: E402

CFG = HERE / "work" / "c10b-config.json"
LEDGER = ROOT / "docs" / "CONFIRMATION-LEDGER.md"
ANCHOR = "### Claim C10b — timing amendment"


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, choices=sorted(D.CONFIGS))
    ap.add_argument("--received", required=True, help="UTC time the coordinator's message arrived, e.g. 2026-10-04T04:31Z")
    ap.add_argument("--message-file", required=True, type=Path, help="the message text, verbatim")
    ap.add_argument("--uncovered", default="", help="bulk switches in the message that no 551 design expresses (-> R7)")
    a = ap.parse_args()
    assert not CFG.exists(), f"{CFG} exists: F is recorded once"
    msg = a.message_file.read_text()
    base, step = D.CONFIGS[a.config]
    rec = dict(config=a.config, base=base, step=step, recorded_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               received=a.received, by='session "Explore contribution areas"', message=msg, uncovered=a.uncovered,
               accesses_r6=step is not None, preregistration_sha256=sha(HERE / "PREREGISTRATION.md"))
    CFG.parent.mkdir(parents=True, exist_ok=True)
    CFG.write_text(json.dumps(rec, indent=1) + "\n")
    s = LEDGER.read_text()
    assert s.count(ANCHOR) == 1, "C10b section not found once in the ledger"
    line = (f"- **F recorded {rec['recorded_at']}** (coordinator's message received {a.received}): `{a.config}` "
            f"(base {base}, step {step or 'none: F is a C10a design, R6 is not accessed'}). "
            f"`experiments/551-r6-c10b/work/c10b-config.json` sha256 `{sha(CFG)}`; PREREGISTRATION.md sha256 "
            f"`{rec['preregistration_sha256']}`." + (f" Not covered by 551's designs (-> R7): {a.uncovered}." if a.uncovered else "") + "\n")
    # append at the end of the C10b section (before the next "## " heading after the anchor, or at the end of the file)
    i = s.index(ANCHOR)
    j = s.find("\n## ", i)
    s = s[:j + 1] + line + s[j + 1:] if j >= 0 else s.rstrip("\n") + "\n" + line
    LEDGER.write_text(s)
    print(line.strip())


if __name__ == "__main__":
    main()
