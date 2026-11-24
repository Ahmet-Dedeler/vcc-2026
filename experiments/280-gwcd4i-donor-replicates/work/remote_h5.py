"""A minimal read-only, block-caching file object over HTTP range requests, for h5py.

Used to inspect (and selectively extract from) the GWCD4i S3 objects without downloading them.
No third-party dependency: urllib only.  Blocks are cached in memory (LRU) and, optionally, on disk.

    from remote_h5 import RemoteFile
    import h5py
    f = h5py.File(RemoteFile(url), "r")
"""
from __future__ import annotations

import io
import os
import urllib.request
from collections import OrderedDict

BLOCK = 1 << 22  # 4 MiB


class RemoteFile(io.RawIOBase):
    def __init__(self, url: str, block: int = BLOCK, max_blocks: int = 256, cache_dir: str | None = None):
        self.url = url
        self.block = block
        self.max_blocks = max_blocks
        self.cache_dir = cache_dir
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)
        self._cache: OrderedDict[int, bytes] = OrderedDict()
        self._pos = 0
        self.n_requests = 0
        self.n_bytes = 0
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=60) as r:
            self.size = int(r.headers["Content-Length"])

    # ------------------------------------------------------------------ io plumbing
    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def seek(self, offset: int, whence: int = 0) -> int:
        if whence == 0:
            self._pos = offset
        elif whence == 1:
            self._pos += offset
        else:
            self._pos = self.size + offset
        return self._pos

    def tell(self) -> int:
        return self._pos

    # ------------------------------------------------------------------ block fetch
    def _fetch(self, bi: int) -> bytes:
        b = self._cache.get(bi)
        if b is not None:
            self._cache.move_to_end(bi)
            return b
        path = None
        if self.cache_dir:
            path = os.path.join(self.cache_dir, f"{bi}.blk")
            if os.path.exists(path):
                with open(path, "rb") as fh:
                    b = fh.read()
        if b is None:
            lo = bi * self.block
            hi = min(lo + self.block, self.size) - 1
            req = urllib.request.Request(self.url, headers={"Range": f"bytes={lo}-{hi}"})
            for attempt in range(4):
                try:
                    with urllib.request.urlopen(req, timeout=120) as r:
                        b = r.read()
                    break
                except Exception:
                    if attempt == 3:
                        raise
            self.n_requests += 1
            self.n_bytes += len(b)
            if path:
                with open(path, "wb") as fh:
                    fh.write(b)
        self._cache[bi] = b
        self._cache.move_to_end(bi)
        while len(self._cache) > self.max_blocks:
            self._cache.popitem(last=False)
        return b

    def readinto(self, buf) -> int:  # type: ignore[override]
        n = len(buf)
        n = min(n, max(0, self.size - self._pos))
        if n == 0:
            return 0
        out = bytearray()
        pos = self._pos
        remaining = n
        while remaining > 0:
            bi = pos // self.block
            off = pos - bi * self.block
            blk = self._fetch(bi)
            take = min(remaining, len(blk) - off)
            out += blk[off:off + take]
            pos += take
            remaining -= take
        buf[:n] = out
        self._pos += n
        return n

    def read(self, size: int = -1) -> bytes:  # type: ignore[override]
        if size < 0:
            size = self.size - self._pos
        buf = bytearray(size)
        got = self.readinto(buf)
        return bytes(buf[:got])
