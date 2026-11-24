"""Fetch named members of a remote ZIP by HTTP range requests, inflate, verify CRC32/size. Writes into outdir.
usage: fetchzip.py <url> <total_bytes> <outdir> <member> [<member> ...]"""
import os, struct, sys, time, urllib.request, zlib
from pathlib import Path
url, total, outdir = sys.argv[1], int(sys.argv[2]), Path(sys.argv[3]); members = sys.argv[4:]
outdir.mkdir(parents=True, exist_ok=True)
def rng(a, b):
    for attempt in range(6):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={a}-{b}"})
            with urllib.request.urlopen(req, timeout=120) as r:
                assert r.status == 206, r.status
                return r.read()
        except Exception as e:  # noqa
            print("retry", attempt, e, flush=True); time.sleep(3 * (attempt + 1))
    raise RuntimeError("range fetch failed")
def rng_stream(a, b, sink):
    for attempt in range(6):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={a}-{b}"})
            with urllib.request.urlopen(req, timeout=180) as r:
                assert r.status == 206, r.status
                got = 0
                while True:
                    chunk = r.read(4 << 20)
                    if not chunk: break
                    sink(chunk); got += len(chunk)
                assert got == b - a + 1, (got, b - a + 1)
                return
        except Exception as e:  # noqa
            print("stream retry", attempt, e, flush=True); time.sleep(3 * (attempt + 1))
            raise
    raise RuntimeError("stream failed")
tail = rng(total - 2_000_000, total - 1)
i = tail.rfind(b"PK\x05\x06"); (_, _, _, _, nrec, cdsize, cdoff, _) = struct.unpack("<IHHHHIIH", tail[i:i+22])
cd = tail[cdoff - (total - 2_000_000): cdoff - (total - 2_000_000) + cdsize]
entries = {}; p = 0
while p < len(cd) and cd[p:p+4] == b"PK\x01\x02":
    (sig, vm, vn, fl, meth, mt, md, crc, csz, usz, nlen, elen, clen2, dn, ia, ea, lho) = struct.unpack("<IHHHHHHIIIHHHHHII", cd[p:p+46])
    name = cd[p+46:p+46+nlen].decode(); entries[name] = (meth, crc, csz, usz, lho); p += 46 + nlen + elen + clen2
for m in members:
    meth, crc, csz, usz, lho = entries[m]
    out = outdir / Path(m).name
    if out.exists() and out.stat().st_size == usz:
        print("skip (exists)", m, flush=True); continue
    lh = rng(lho, lho + 29); (sig, ver, flag, meth2, t, d, crc2, csz2, usz2, nlen, elen) = struct.unpack("<IHHHHHIIIHH", lh)
    assert sig == 0x04034b50, hex(sig)
    start = lho + 30 + nlen + elen
    t0 = time.time(); dec = zlib.decompressobj(-15) if meth == 8 else None; c = 0; n = 0
    tmp = out.with_suffix(out.suffix + ".part")
    with open(tmp, "wb") as f:
        def sink(chunk):
            global c, n
            data = dec.decompress(chunk) if dec else chunk
            c = zlib.crc32(data, c); n += len(data); f.write(data)
        rng_stream(start, start + csz - 1, sink)
        if dec:
            data = dec.flush(); c = zlib.crc32(data, c); n += len(data); f.write(data)
    dt = time.time() - t0
    ok = (c & 0xFFFFFFFF) == crc and n == usz
    print(f"{m}: {csz/1e6:.1f} MB compressed in {dt:.1f}s ({csz/dt/1e6:.2f} MB/s) -> {n:,} bytes crc_ok={ok}", flush=True)
    if not ok:
        tmp.unlink(); raise SystemExit(f"CRC/size mismatch for {m}")
    os.replace(tmp, out)
