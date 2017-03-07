"""Utility functions used in other modules."""

from __future__ import absolute_import, division

import copy
from collections import OrderedDict, Callable
import itertools
import math


def nCr(n, r):
    f = math.factorial
    return f(n) / f(r) / f(n - r)


class LenGen(object):
    def __init__(self, gen, length):
        self.gen = gen
        self.length = length

    def __call__(self):
        return itertools.islice(self.gen(), self.length)

    def __iter__(self):
        return self.gen

    def __len__(self):
        return self.length


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
