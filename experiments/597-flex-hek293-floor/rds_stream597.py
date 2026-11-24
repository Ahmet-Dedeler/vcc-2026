"""597 — 424's streaming RDS reader plus R byte code (BCODESXP, SEXP type 21), which 424's reader rejects.

The VIPerturb pilot object stores a compiled R function (a closure whose body is byte code). R's serialize.c writes it as
  ReadBC:      nreps = InInteger;  ReadBC1
  ReadBC1:     code = ReadItem (an INTSXP);  ReadBCConsts
  ReadBCConsts: n = InInteger; n times: type = InInteger;
                 BCODESXP -> ReadBC1;  LANGSXP/LISTSXP/BCREPDEF/BCREPREF/ATTRLANGSXP/ATTRLISTSXP -> ReadBCLang;
                 otherwise -> ReadItem
  ReadBCLang:  BCREPREF -> InInteger (an index);  BCREPDEF -> pos = InInteger, type = InInteger, then as below;
               LANGSXP/LISTSXP/ATTR*: [attributes if ATTR*] tag = ReadItem; car = ReadBCLang(InInteger);
               cdr = ReadBCLang(InInteger);  otherwise -> ReadItem
and, as for every non-CHARSXP item, the item's attributes follow when its flags say so. The byte code is parsed and
discarded (returned as a placeholder); nothing a Seurat count matrix or meta.data needs lives in it. 424's file is
unchanged: this module subclasses its reader.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "424-viperturb-k562"))
from rds_stream import (ATTRLANGSXP, ATTRLISTSXP, BCODESXP, LANGSXP, LISTSXP, RDSStream,  # noqa: E402
                        Special)

BCREPREF, BCREPDEF = 243, 244
LANGTYPES = (LANGSXP, LISTSXP, BCREPDEF, BCREPREF, ATTRLANGSXP, ATTRLISTSXP)


class RDSStream597(RDSStream):
    n_bytecode = 0

    def item(self, flags=None):
        if flags is None:
            flags = self.int()
        if (flags & 0xFF) == BCODESXP:
            has_attr = bool(flags & (1 << 9))
            self.int()                      # number of reps (BCREPDEF slots)
            self._bc1()
            if has_attr:
                self.item()                 # attributes of the BCODESXP itself
            self.n_bytecode += 1
            return Special(BCODESXP)
        return super().item(flags)

    def _bc1(self):
        self.item()                         # the code vector
        for _ in range(self.int()):         # constants
            typ = self.int()
            if typ == BCODESXP:
                self._bc1()
            elif typ in LANGTYPES:
                self._bclang(typ)
            else:
                self.item()

    def _bclang(self, typ):
        if typ == BCREPREF:
            self.int()
            return
        if typ in LANGTYPES:
            if typ == BCREPDEF:
                self.int()                  # pos
                typ = self.int()            # the real type
            if typ in (ATTRLANGSXP, ATTRLISTSXP):
                self.item()                 # attributes
            self.item()                     # tag
            self._bclang(self.int())        # car
            self._bclang(self.int())        # cdr
            return
        self.item()                         # a full item follows its type integer
