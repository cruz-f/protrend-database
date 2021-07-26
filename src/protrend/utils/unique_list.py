from collections import UserList
from dataclasses import field
from functools import partial
from typing import Sequence


class UniqueList(UserList):

    def __init__(self, sequence: Sequence = None, output: str = 'take_all'):

        super().__init__()

        self.output = output

        if sequence is not None:
            sequence = self._unique_list(sequence)

            self.data = sequence

    def __call__(self, output: str = None):
        if output is None:
            output = self.output

        if output == 'take_first':
            return self.take_first()

        elif output == 'take_last':
            return self.take_last()

        return self.take_all()

    @staticmethod
    def _unique_list(sequence):
        lst = []

        for element in sequence:
            if element not in lst:
                lst.append(element)

        return lst

    def __setitem__(self, i, item):

        if item not in self.data:
            self.data[i] = item

    def __delitem__(self, i):
        del self.data[i]

    def __add__(self, other):

        other = self._unique_list(other)

        return self.__class__(self.data + other)

    def __radd__(self, other):

        other = self._unique_list(other)

        return self.__class__(self.data + other)

    def __iadd__(self, other):

        other = self._unique_list(other)
        self.data += other
        return self

    def add(self, item):
        return self.append(item)

    def append(self, item):
        if item not in self.data:
            self.data.append(item)

    def insert(self, i, item):
        if item not in self.data:
            self.data.insert(i, item)

    def copy(self):
        return self.__class__(self)

    def update(self, other):
        return self.extend(other)

    def extend(self, other):

        other = self._unique_list(other)

        self.data.extend(other)

    def take_all(self):
        return self.data

    def take_first(self):
        if self.data:
            return self.data[0]

    def take_last(self):
        if self.data:
            return self.data[-1]


def unique_field(output: str = None, init: bool = False):

    unique_list = partial(UniqueList, output=output)

    return field(default_factory=unique_list, init=init)