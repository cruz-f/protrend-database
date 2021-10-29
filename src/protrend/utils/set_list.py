from collections import UserList
from dataclasses import field
from functools import partial
from typing import Sequence, Any, Iterator, Union


class SetList(UserList):

    def __init__(self,
                 sequence: Union[Iterator, Sequence] = None,
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

    def __add__(self, other: Sequence):

        new_instance = self.copy()
        new_instance.extend(other)

        return new_instance

    def __radd__(self, other: Sequence):

        new_instance = self.copy()
        new_instance.extend(other)

        return new_instance

    def __iadd__(self, other: Sequence):

        self._add_elements(other)

        return self

    def append(self, item):
        self._add_element(item)

    def insert(self, i, item):
        if item not in self.data:
            self.data.insert(i, item)

    def copy(self) -> 'SetList':
        return self.__class__(self)

    def extend(self, other: Sequence):

        self._add_elements(other)

    def take_all(self):
        return self.data

    def take_first(self):
        if self.data:
            return self.data[0]

    def take_last(self):
        if self.data:
            return self.data[-1]


def set_list_field(output: str = None, init: bool = False):
    set_list = partial(SetList, output=output)

    return field(default_factory=set_list, init=init)
