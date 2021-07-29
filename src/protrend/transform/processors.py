import functools
import re

from w3lib.html import remove_tags

more_pattern = re.compile(r'\s\smore\s\s\s')
white_space_pattern = re.compile(r'\s')
multiple_white_space_pattern = re.compile(r' +')
pubmed_pattern = re.compile(r'\[[0-9]*]|\[[0-9]*,|\s[0-9]*,|\s[0-9]*]')


def handle_nan(fn):
    @functools.wraps
    def wrapper(item):
        if item:
            return fn(item)
        return item

    return wrapper


@handle_nan
def lstrip(item: str) -> str:
    return item.lstrip()


@handle_nan
def rstrip(item: str) -> str:
    return item.rstrip()


@handle_nan
def remove_white_space(item: str) -> str:
    return re.sub(white_space_pattern, repl='', string=item)


@handle_nan
def remove_multiple_white_space(item: str) -> str:
    return re.sub(multiple_white_space_pattern, repl=' ', string=item)


@handle_nan
def remove_regprecise_more(item: str) -> str:
    return re.sub(more_pattern, repl='', string=item)


@handle_nan
def remove_pubmed(item: str) -> str:
    return re.sub(pubmed_pattern, repl='', string=item)


@handle_nan
def remove_html_tags(item: str) -> str:
    return remove_tags(item)
