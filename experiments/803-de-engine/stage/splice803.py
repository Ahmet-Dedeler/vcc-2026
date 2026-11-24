#!/usr/bin/env python3
"""803 staging: ref_tv/reference.npz shares its two big members (cs.npy, emit_cs.npy, 3.17 GB, same CRC and size) byte for
byte with 323's reference.npz (the same h1-test controls). Upload only the rest and re-splice on EC2.
  cut    <donor.npz> <target.npz> <out.bin> <out.json>   (Mac) target minus the shared members' stored bytes
  splice <donor.npz> <rest.bin> <rest.json> <out.npz>     (EC2) rebuild; the caller checks sha256 against rest.json
"""
import hashlib, json, sys, zipfile

MEMBERS = ("cs.npy", "emit_cs.npy")


def data_ranges(path):
    out = {}
    with zipfile.ZipFile(path) as z, open(path, "rb") as f:
        for i in z.infolist():
            if i.filename in MEMBERS:
                assert i.compress_type == zipfile.ZIP_STORED
                f.seek(i.header_offset)
                h = f.read(30)
                n, m = int.from_bytes(h[26:28], "little"), int.from_bytes(h[28:30], "little")
                out[i.filename] = (i.header_offset + 30 + n + m, i.compress_size, i.CRC)
    return out


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


if sys.argv[1] == "cut":
    donor, tgt, ob, oj = sys.argv[2:6]
    rd, rt = data_ranges(donor), data_ranges(tgt)
    segs, pos = [], 0
    with open(tgt, "rb") as f, open(ob, "wb") as g:
        for m in sorted(MEMBERS, key=lambda k: rt[k][0]):
            assert rd[m][1:] == rt[m][1:], m
            off, size, _ = rt[m]
            f.seek(pos); g.write(f.read(off - pos))
            segs.append(dict(member=m, target_off=off, donor_off=rd[m][0], size=size))
            pos = off + size
        f.seek(pos); g.write(f.read())
    json.dump(dict(segments=segs, target_sha256=sha(tgt), donor_sha256=sha(donor)), open(oj, "w"), indent=1)
    print(json.dumps(segs))
else:
    donor, rb, rj, out = sys.argv[2:6]
    spec = json.load(open(rj))
    assert sha(donor) == spec["donor_sha256"], "donor sha differs"
    with open(rb, "rb") as r, open(donor, "rb") as d, open(out, "wb") as g:
        pos = 0
        for s in sorted(spec["segments"], key=lambda s: s["target_off"]):
            g.write(r.read(s["target_off"] - pos))
            d.seek(s["donor_off"]); left = s["size"]
            while left:
                b = d.read(min(left, 1 << 24)); g.write(b); left -= len(b)
            pos = s["target_off"] + s["size"]
        g.write(r.read())
    ok = sha(out) == spec["target_sha256"]
    print(f"splice {out}: sha {'OK' if ok else 'MISMATCH'}")
    sys.exit(0 if ok else 1)
