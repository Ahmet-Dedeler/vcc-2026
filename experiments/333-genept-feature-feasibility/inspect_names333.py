"""Key/schema diagnostic only; no embedding statistics or biological outcomes."""

from collections import Counter
import zipfile
import acquire333 as A
import convert_recovery333 as R


def main():
    R.verify()
    path = A.WORK / "name-audit.json"
    if path.exists():
        raise FileExistsError("name audit exists")
    with zipfile.ZipFile(A.ZIP) as archive:
        A.MEMBER = R.ACTUAL_MEMBER
        member = A.select_member(archive.infolist())
        with archive.open(member) as stream:
            mapping = A.RestrictedNumpyUnpickler(stream).load()
    assert isinstance(mapping, dict)
    keys = list(mapping)
    strings = [k for k in keys if isinstance(k, str)]
    mixed = [k for k in strings if k != k.upper()]
    spaces = [k for k in strings if k != k.strip()]
    counts = Counter(k.upper() for k in strings)
    collisions = [k for k, n in counts.items() if n > 1]
    A.write(path, {
        "archive_sha256": A.digest(A.ZIP), "n_keys": len(keys),
        "non_string_count": len(keys)-len(strings), "empty_count": strings.count(""),
        "mixed_case_count": len(mixed), "mixed_case_examples": sorted(mixed)[:30],
        "whitespace_count": len(spaces), "whitespace_examples": sorted(spaces)[:10],
        "uppercasing_collision_count": len(collisions), "collision_examples": sorted(collisions)[:20],
        "embedding_values_analyzed": False, "biological_outcomes_read": False})
    print(path.read_text())


if __name__ == "__main__":
    main()
