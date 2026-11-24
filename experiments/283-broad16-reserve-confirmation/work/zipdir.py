"""List the central directory of a remote ZIP via HTTP range requests (no full download)."""
import struct, sys, time, urllib.request, zlib
url = sys.argv[1]; total = int(sys.argv[2]); tail = int(sys.argv[3]) if len(sys.argv) > 3 else 2_000_000
def rng(a, b):
    req = urllib.request.Request(url, headers={"Range": f"bytes={a}-{b}"})
    with urllib.request.urlopen(req, timeout=60) as r:
        assert r.status == 206, r.status
        return r.read()
t0 = time.time(); data = rng(total - tail, total - 1); dt = time.time() - t0
print(f"tail {len(data)} bytes in {dt:.1f}s -> {len(data)/dt/1e6:.1f} MB/s")
# find EOCD
i = data.rfind(b"PK\x05\x06")
assert i >= 0
sig, disk, cdisk, nrec_disk, nrec, cdsize, cdoff, clen = struct.unpack("<IHHHHIIH", data[i:i+22])
if nrec == 0xFFFF or cdoff == 0xFFFFFFFF:
    j = data.rfind(b"PK\x06\x06")
    (sig64, sz, vmade, vneed, d1, d2, nrec_disk, nrec, cdsize, cdoff) = struct.unpack("<IQHHIIQQQQ", data[j:j+56])
print("entries", nrec, "cd size", cdsize, "cd offset", cdoff)
if cdoff < total - tail:
    cd = rng(cdoff, cdoff + cdsize - 1)
else:
    cd = data[cdoff - (total - tail): cdoff - (total - tail) + cdsize]
p = 0; rows = []
while p < len(cd) and cd[p:p+4] == b"PK\x01\x02":
    (sig, vm, vn, fl, meth, mt, md, crc, csz, usz, nlen, elen, clen2, dn, ia, ea, lho) = struct.unpack("<IHHHHHHIIIHHHHHII", cd[p:p+46])
    name = cd[p+46:p+46+nlen].decode("utf-8", "replace")
    extra = cd[p+46+nlen:p+46+nlen+elen]
    # zip64 extra
    if usz == 0xFFFFFFFF or csz == 0xFFFFFFFF:
        q = 0
        while q + 4 <= len(extra):
            hid, hsz = struct.unpack("<HH", extra[q:q+4]); body = extra[q+4:q+4+hsz]
            if hid == 1:
                vals = struct.unpack("<" + "Q" * (len(body)//8), body[:8*(len(body)//8)])
                k = 0
                if usz == 0xFFFFFFFF: usz = vals[k]; k += 1
                if csz == 0xFFFFFFFF: csz = vals[k]; k += 1
            q += 4 + hsz
    rows.append((name, usz, csz, meth))
    p += 46 + nlen + elen + clen2
print(f"{len(rows)} entries")
for name, usz, csz, meth in rows:
    print(f"{usz:>14,d} {csz:>14,d} m{meth} {name}")
