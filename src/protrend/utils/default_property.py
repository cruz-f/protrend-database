from typing import Any


class DefaultProperty:

    def __init__(self, default: Any = None):
        """
        The default property is a python data-descriptor that can be used to define a default value
        for a custom pipeline component.
        If the pipeline initializer receives a new value the default one is ignored.

        :type default: Any

        :param default: The default value for that attribute
        """
        self.default = default

    def set_default(self, default: Any):

        if default is None:
            return

        self.default = default

    def __set_name__(self, owner, name):
        self.name = name

    def __get__(self, instance, owner):
        if instance is None:
            return self

        return instance.__dict__.get(self.name, self.default)

    def __set__(self, instance, value):
        if instance is None:
            return self

        if value is not None:
            instance.__dict__[self.name] = value
