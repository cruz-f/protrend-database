import time
from functools import wraps
from typing import Union


def sleep(sec: Union[int, float] = 0.25):
    def func_wrapper(fn):
        @wraps(fn)
        def wrapper(*args, **kwargs):
            time.sleep(sec)
            return fn(*args, **kwargs)

        return wrapper

    return func_wrapper
