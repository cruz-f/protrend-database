from typing import Any, Type


class DefaultProperty:

    """
    The default property is a python data-descriptor that can be used to define a default value
    for a custom pipeline component.
    If the pipeline initializer receives a new value the default one is ignored.
    """
    def get_default(self, cls):
        return getattr(cls, f'_{self.name}', None)

    def set_default(self, cls: Type, default: Any):

        if default is None:
            return

        setattr(cls, f'_{self.name}', default)

    def __set_name__(self, owner, name):
        self.name = name

    def __get__(self, instance, owner):
        if instance is None:
            return self

        default = getattr(instance, f'_{self.name}', None)
        return instance.__dict__.get(self.name, default)

    def __set__(self, instance, value):
        if instance is None:
            return self

        if value is not None:
            instance.__dict__[self.name] = value
