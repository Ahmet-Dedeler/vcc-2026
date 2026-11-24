"""424 — a streaming reader for R's XDR serialization (.rds, gzip), for objects too big to load.

`readRDS` on one VIPerturb-seq bin peaks at 16.7 GB RSS (binC, the smallest; binB would be ~22 GB). This parser walks the
serialized object in one pass over the gzip stream and never holds a large vector in memory: any integer / logical /
double vector longer than `big` elements is streamed to a native-endian file under `spill_dir` (doubles optionally cast
to float32, with an exactness audit) and returned as a `Spilled` handle. Everything else is returned as plain Python /
numpy objects:

    NULL -> None; SYMSXP -> Sym; CHARSXP -> str | None; STRSXP -> list[str|None]; INTSXP/LGLSXP -> int32 array;
    REALSXP -> float64 array; VECSXP -> RList (list with .attrs); pairlists -> PairList ([(tag, value)] with .attrs);
    S4SXP -> S4 (.attrs dict); ENVSXP -> Env; ALTREP compact_intseq / compact_realseq / deferred_string / wrap_* are
    expanded; any vector with attributes carries them in `.attrs` (numpy vectors are wrapped in `Vec`).

Format reference: R's src/main/serialize.c (ReadItem / InItem). Only what a Seurat v5 object needs is implemented;
anything else (byte code, class refs) raises, so a silent mis-parse is impossible.
"""
from __future__ import annotations

import gzip
import io
import struct
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

NILVALUE, GLOBALENV, UNBOUNDVALUE, MISSINGARG, BASENAMESPACE = 254, 253, 252, 251, 250
NAMESPACESXP, PACKAGESXP, PERSISTSXP, EMPTYENV, BASEENV, REFSXP = 249, 248, 247, 242, 241, 255
ATTRLANGSXP, ATTRLISTSXP, ALTREP_SXP = 240, 239, 238
SYMSXP, LISTSXP, CLOSXP, ENVSXP, PROMSXP, LANGSXP, SPECIALSXP, BUILTINSXP, CHARSXP = 1, 2, 3, 4, 5, 6, 7, 8, 9
LGLSXP, INTSXP, REALSXP, CPLXSXP, STRSXP, DOTSXP, VECSXP, EXPRSXP = 10, 13, 14, 15, 16, 17, 19, 20
BCODESXP, EXTPTRSXP, WEAKREFSXP, RAWSXP, S4SXP = 21, 22, 23, 24, 25
PAIRS = (LISTSXP, LANGSXP, CLOSXP, PROMSXP, DOTSXP, ATTRLANGSXP, ATTRLISTSXP)
NA_INT = -2147483648


@dataclass
class Sym:
    name: str


@dataclass
class Special:
    code: int


@dataclass
class Env:
    attrs: dict = field(default_factory=dict)
    frame: object = None
    enclos: object = None
    hashtab: object = None


@dataclass
class S4:
    attrs: dict


@dataclass
class Spilled:
    """A big numeric vector streamed to disk (native endian)."""
    path: Path
    dtype: str
    n: int
    sxp: int
    audit: dict
    attrs: dict = field(default_factory=dict)

    def memmap(self):
        return np.memmap(self.path, dtype=self.dtype, mode="r", shape=(self.n,))


class Vec(np.ndarray):
    """numpy vector carrying R attributes."""
    attrs: dict = {}


class RList(list):
    attrs: dict = {}


class PairList(list):
    attrs: dict = {}
    kind: int = LISTSXP


def attrs_dict(pl) -> dict:
    if pl is None:
        return {}
    out = {}
    for tag, val in pl:
        out[tag.name if isinstance(tag, Sym) else str(tag)] = val
    return out


class RDSStream:
    def __init__(self, path, spill_dir, big=5_000_000, real_as_f32=True, chunk=1 << 24, log=print):
        self.f = io.BufferedReader(gzip.open(path, "rb"), buffer_size=1 << 24)
        self.spill_dir = Path(spill_dir)
        self.spill_dir.mkdir(parents=True, exist_ok=True)
        self.big, self.real_as_f32, self.chunk, self.log = big, real_as_f32, chunk, log
        self.refs: list = []
        self.n_spilled = 0
        self.bytes_read = 0

    # ------------------------------------------------------------------ primitives
    def _read(self, n):
        b = self.f.read(n)
        if len(b) != n:
            raise EOFError(f"wanted {n} bytes, got {len(b)}")
        self.bytes_read += n
        return b

    def int(self):
        return struct.unpack(">i", self._read(4))[0]

    def length(self):
        n = self.int()
        if n == -1:
            hi, lo = self.int(), self.int()
            n = (hi << 32) + lo
        return n

    def header(self):
        fmt = self._read(2)
        if fmt != b"X\n":
            raise ValueError(f"not an XDR serialization: {fmt!r}")
        version, writer, minreader = self.int(), self.int(), self.int()
        enc = None
        if version == 3:
            enc = self._read(self.int()).decode()
        return dict(version=version, writer=writer, min_reader=minreader, encoding=enc)

    # ------------------------------------------------------------------ vectors
    def _numeric(self, sxp, n):
        width, src = (8, ">f8") if sxp == REALSXP else (4, ">i4")
        if n <= self.big:
            a = np.frombuffer(self._read(width * n), dtype=src).astype(src[1:].replace(">", ""), copy=True)
            return a.astype(np.float64 if sxp == REALSXP else np.int32).view(Vec)
        self.n_spilled += 1
        to32 = sxp == REALSXP and self.real_as_f32
        dtype = "float32" if to32 else ("float64" if sxp == REALSXP else "int32")
        path = self.spill_dir / f"spill{self.n_spilled:02d}_{dtype}_{n}.bin"
        audit = dict(nonint=0, max=-np.inf, min=np.inf, na=0, sum=0.0)
        per = self.chunk // width
        with open(path, "wb") as out:
            left = n
            while left:
                k = min(per, left)
                a = np.frombuffer(self._read(width * k), dtype=src)
                if sxp == REALSXP:
                    a = a.astype(np.float64)
                    audit["na"] += int(np.isnan(a).sum())
                    audit["nonint"] += int((a != np.round(a)).sum())
                    audit["max"] = max(audit["max"], float(np.nanmax(a)) if k else -np.inf)
                    audit["min"] = min(audit["min"], float(np.nanmin(a)) if k else np.inf)
                    audit["sum"] += float(np.nansum(a))
                    (a.astype(np.float32) if to32 else a).tofile(out)
                else:
                    a = a.astype(np.int32)
                    audit["na"] += int((a == NA_INT).sum())
                    audit["max"] = max(audit["max"], int(a.max()) if k else -np.inf)
                    audit["min"] = min(audit["min"], int(a.min()) if k else np.inf)
                    a.tofile(out)
                left -= k
        if to32:
            audit["f32_exact"] = bool(audit["nonint"] == 0 and audit["max"] < 2 ** 24)
        self.log(f"  spilled {['', 'int', 'real'][1 + (sxp == REALSXP)]} vector n={n:,} -> {path.name} {audit}")
        return Spilled(path=path, dtype=dtype, n=n, sxp=sxp, audit=audit)

    # ------------------------------------------------------------------ items
    def item(self, flags=None):
        if flags is None:
            flags = self.int()
        t = flags & 0xFF
        has_attr = bool(flags & (1 << 9))
        has_tag = bool(flags & (1 << 10))
        if t == NILVALUE:
            return None
        if t in (EMPTYENV, BASEENV, GLOBALENV, UNBOUNDVALUE, MISSINGARG, BASENAMESPACE):
            return Special(t)
        if t == REFSXP:
            idx = flags >> 8
            if idx == 0:
                idx = self.int()
            return self.refs[idx - 1]
        if t == PERSISTSXP or t == PACKAGESXP or t == NAMESPACESXP:
            if self.int() != 0:
                raise ValueError("string vec with nonzero lead")
            n = self.int()
            v = (t, [self.item() for _ in range(n)])
            self.refs.append(v)
            return v
        if t == SYMSXP:
            s = Sym(self.item())
            self.refs.append(s)
            return s
        if t == ENVSXP:
            self.int()                                  # locked
            env = Env()
            self.refs.append(env)
            env.enclos, env.frame, env.hashtab = self.item(), self.item(), self.item()
            env.attrs = attrs_dict(self.item())
            return env
        if t in PAIRS:
            out = PairList()
            out.kind = t
            first_attr = None
            while True:
                attr = self.item() if has_attr else None
                if first_attr is None and attr is not None:
                    first_attr = attr
                tag = self.item() if has_tag else None
                car = self.item()
                out.append((tag, car))
                nf = self.int()
                nt = nf & 0xFF
                if nt == t and t == LISTSXP:            # iterate along the CDR chain (no recursion)
                    flags, has_attr, has_tag = nf, bool(nf & (1 << 9)), bool(nf & (1 << 10))
                    continue
                cdr = self.item(nf)
                if cdr is not None:
                    out.append(("__cdr__", cdr))
                break
            out.attrs = attrs_dict(first_attr)
            return out
        if t == EXTPTRSXP:
            obj = Env()
            self.refs.append(obj)
            self.item(); self.item()                     # prot, tag
            if has_attr:
                obj.attrs = attrs_dict(self.item())
            return obj
        if t == WEAKREFSXP:
            obj = Env()
            self.refs.append(obj)
            if has_attr:
                obj.attrs = attrs_dict(self.item())
            return obj
        if t in (SPECIALSXP, BUILTINSXP):
            n = self.int()
            return Special(t), self._read(n).decode()
        if t == CHARSXP:
            n = self.int()
            if n == -1:
                return None
            return self._read(n).decode("utf-8", errors="replace")
        if t == ALTREP_SXP:
            info, state, attr = self.item(), self.item(), self.item()
            cls = info[0][1].name if isinstance(info, PairList) and isinstance(info[0][1], Sym) else str(info)
            val = self._altrep(cls, state)
            a = attrs_dict(attr)
            if a:
                val = self._with_attrs(val, a)
            return val
        if t in (LGLSXP, INTSXP, REALSXP):
            val = self._numeric(t, self.length())
        elif t == CPLXSXP:
            val = np.frombuffer(self._read(16 * self.length()), dtype=">f8").astype(np.float64).view(Vec)
        elif t == STRSXP:
            val = RList(self.item() for _ in range(self.length()))
        elif t in (VECSXP, EXPRSXP):
            val = RList(self.item() for _ in range(self.length()))
        elif t == RAWSXP:
            val = self._read(self.length())
        elif t == S4SXP:
            val = S4(attrs={})
        else:
            raise NotImplementedError(f"SEXP type {t} (flags {flags:#x}) is not supported")
        if has_attr:
            val = self._with_attrs(val, attrs_dict(self.item()))
        return val

    @staticmethod
    def _with_attrs(val, a):
        if isinstance(val, S4):
            val.attrs = a
        elif isinstance(val, (Spilled, RList, PairList, Vec)):
            val.attrs = a
        elif isinstance(val, np.ndarray):
            val = val.view(Vec)
            val.attrs = a
        elif isinstance(val, list):
            val = RList(val)
            val.attrs = a
        return val

    def _altrep(self, cls, state):
        if cls == "compact_intseq":
            n, start, inc = (float(x) for x in np.asarray(state))
            return (start + inc * np.arange(int(n))).astype(np.int32).view(Vec)
        if cls == "compact_realseq":
            n, start, inc = (float(x) for x in np.asarray(state))
            return (start + inc * np.arange(int(n))).astype(np.float64).view(Vec)
        if cls == "deferred_string":
            arg = state[0][1] if isinstance(state, PairList) else state
            arr = np.asarray(arg)
            return RList(str(int(x)) if float(x).is_integer() else repr(float(x)) for x in arr)
        if cls.startswith("wrap_"):
            return state[0] if isinstance(state, (list, RList)) else state
        raise NotImplementedError(f"ALTREP class {cls}")

    def read(self):
        hdr = self.header()
        obj = self.item()
        rest = self.f.read(1)
        return hdr, obj, len(rest) == 0
