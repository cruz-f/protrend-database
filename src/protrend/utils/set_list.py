from collections import UserList
from dataclasses import field
from functools import partial
from typing import Sequence, Any


class SetList(UserList):

    def __init__(self,
                 sequence: Sequence = None,
                 output: str = 'take_all'):

        super().__init__()

        self.output = output

        if sequence is None:
            sequence = []

        self._add_elements(sequence)

    def __call__(self, output: str = None):
        if output is None:
            output = self.output

        if output == 'take_first':
            return self.take_first()

        elif output == 'take_last':
            return self.take_last()

        return self.take_all()

    def _add_element(self, element: Any):
        if element not in self.data:
            self.data.append(element)

    def _add_elements(self, sequence: Sequence):

        for element in sequence:
            self._add_element(element)

    def __setitem__(self, i, item):

        if item not in self.data:
            self.data[i] = item

    def __delitem__(self, i):
        del self.data[i]

    def __add__(self, other):

        other = self.__class__(other)

        return self.__class__(self.data + other)

    def __radd__(self, other):

        other = self.__class__(other)

        return self.__class__(self.data + other)

    def __iadd__(self, other):

        other = self.__class__(other)

        self.data += other

        return self

    def append(self, item):
        self._add_element(item)

    def insert(self, i, item):
        if item not in self.data:
            self.data.insert(i, item)

    def copy(self):
        return self.__class__(self)

    def extend(self, other):

        other = self.__class__(other)

        self.data.extend(other)

    def take_all(self):
        return self.data

    def take_first(self):
        if self.data:
            return self.data[0]

    def take_last(self):
        if self.data:
            return self.data[-1]


def set_list_field(output: str = None, init: bool = False):

    unique_list = partial(SetList, output=output)

    return field(default_factory=unique_list, init=init)