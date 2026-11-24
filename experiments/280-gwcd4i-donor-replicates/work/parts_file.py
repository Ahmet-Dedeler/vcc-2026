"""Read a byte-range-split download (work/parts/p0..p7) as one virtual file, so the 44.57 GB object
never has to exist twice on a disk that is 95 % full.  Passed straight to h5py.File(...)."""
from __future__ import annotations

import io
import os


class PartsFile(io.RawIOBase):
    def __init__(self, paths):
        self.fhs = [open(p, "rb") for p in paths]
        self.sizes = [os.path.getsize(p) for p in paths]
        self.offsets = []
        o = 0
        for s in self.sizes:
            self.offsets.append(o)
            o += s
        self.size = o
        self._pos = 0

    def readable(self):
        return True

    def seekable(self):
        return True

    def seek(self, offset, whence=0):
        if whence == 0:
            self._pos = offset
        elif whence == 1:
            self._pos += offset
        else:
            self._pos = self.size + offset
        return self._pos

    def tell(self):
        return self._pos

    def readinto(self, buf):
        n = min(len(buf), max(0, self.size - self._pos))
        if n == 0:
            return 0
        got = 0
        pos = self._pos
        while got < n:
            i = 0
            for j in range(len(self.offsets)):
                if pos >= self.offsets[j]:
                    i = j
                else:
                    break
            local = pos - self.offsets[i]
            take = min(n - got, self.sizes[i] - local)
            self.fhs[i].seek(local)
            chunk = self.fhs[i].read(take)
            if not chunk:
                break
            buf[got:got + len(chunk)] = chunk
            got += len(chunk)
            pos += len(chunk)
        self._pos = pos
        return got

    def close(self):
        for fh in self.fhs:
            fh.close()
        super().close()
