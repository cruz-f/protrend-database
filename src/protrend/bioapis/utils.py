import time
import unicodedata
import re
from functools import wraps
from typing import Union


def sleep(sec: Union[int, float] = 0.25):
    def func_wrapper(fn):
        @wraps(fn)
        def wrapper(*args, **kwargs):
            cache = kwargs.get('cache', False)

            result, cached_result = fn(*args, **kwargs)

            if cache and cached_result:
                return result, cached_result

            time.sleep(sec)
            return result, cached_result

        return wrapper

    return func_wrapper


def slugify(url):
    """
    Taken from https://github.com/django/django/blob/master/django/utils/text.py
    Convert to ASCII if 'allow_unicode' is False. Convert spaces or repeated
    dashes to single dashes. Remove characters that aren't alphanumerics,
    underscores, or hyphens. Convert to lowercase. Also strip leading and
    trailing whitespace, dashes, and underscores.
    """
    url = str(url)
    url = unicodedata.normalize('NFKD', url).encode('ascii', 'ignore').decode('ascii')
    url = re.sub(r'[^\w\s-]', '', url.lower())
    return re.sub(r'[-\s]+', '-', url).strip('-_')
