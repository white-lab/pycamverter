"""Utility functions used in other modules."""

from __future__ import absolute_import, division

import copy
from collections import OrderedDict, Callable
import difflib
import itertools
import math

from . import regexes


def nCr(n, r):
    f = math.factorial
    return f(n) / f(r) / f(n - r)


class LenGen(object):
    def __init__(self, gen, len):
        self.gen = gen
        self.len = len

    def __call__(self):
        return itertools.islice(self.gen(), self.len)

    def __iter__(self):
        return self.gen

    def __len__(self):
        return self.len


class DefaultOrderedDict(OrderedDict):
    # Source: http://stackoverflow.com/a/6190500/562769
    def __init__(self, default_factory=None, *a, **kw):
        if (
            default_factory is not None and
            not isinstance(default_factory, Callable)
        ):
            raise TypeError("first argument must be callable")

        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)

        self[key] = value = self.default_factory()

        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,

        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(
            self.default_factory,
            self,
        )

    def __deepcopy__(self, memo):
        return type(self)(
            self.default_factory,
            copy.deepcopy(self.items()),
        )

    def __repr__(self):
        return "OrderedDefaultDict({}, {})".format(
            self.default_factory,
            OrderedDict.__repr__(self),
        )


class StrToBin(object):
    def __init__(self, f, encoding="utf-8"):
        self.f = f
        self.encoding = encoding

    def write(self, data):
        self.f.write(data.encode(self.encoding))


SUPERSCRIPT_UNICODE_START = ord(u"\u2070")
SUBSCRIPT_UNICODE_START = ord(u'\u2080')
SCRIPT_MAPPING = {
    str(i): i
    for i in range(10)
}
SCRIPT_MAPPING["+"] = 10
SCRIPT_MAPPING["-"] = 11
SCRIPT_MAPPING["("] = 12
SCRIPT_MAPPING[")"] = 13

try:
    unichr
except NameError:
    unichr = chr


def rewrite_ion_name(name):
    m = regexes.RE_B_Y_IONS.match(name)

    if m:
        name = "".join(m.group(1, 2))

    ret = ""
    sup, sub = False, False
    paren = 0

    for char in name:
        if char in "^_{}":
            if char == "^":
                sup, sub = True, False
            elif char == "_":
                sup, sub = False, True
            elif char == "}":
                sup, sub = False, False
                paren -= 1
            elif char == "{":
                paren += 1
            continue

        if sup:
            if char == "1":
                ret += u"\u00B9"
            elif char == "2":
                ret += u"\u00B2"
            elif char == "3":
                ret += u"\u00B3"
            else:
                ret += unichr(SUPERSCRIPT_UNICODE_START + SCRIPT_MAPPING[char])
        elif sub:
            ret += unichr(SUBSCRIPT_UNICODE_START + SCRIPT_MAPPING[char])
        else:
            ret += char

        if sup or sub:
            if not paren:
                sup, sub = False, False

    return ret


def fuzzy_find(needle, haystack):
    """
    Find the longest matching subsequence of needle within haystack.

    Returns the corresponding index from the beginning of needle.

    Parameters
    ----------
    needle : str
    haystack : str

    Returns
    -------
    int
    """
    s = difflib.SequenceMatcher(a=haystack, b=needle)
    best = s.find_longest_match(0, len(haystack), 0, len(needle))
    return best.a - len(needle) + best.size
