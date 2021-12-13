import os
from dataclasses import dataclass, field, replace
from typing import Dict, List, Callable, Optional

from protrend.utils import Settings


@dataclass
class MultiStack:
    stack: Optional[List[str]] = field(default_factory=list)
    taxa: Optional[List[str]] = field(default_factory=list)
    source: Optional[List[str]] = field(default_factory=list)
    reader: Optional[List[Callable]] = field(default_factory=list)


@dataclass
class WriteStack:
    transformed: Optional[str] = field(default_factory=str)
    integrated: Optional[str] = field(default_factory=str)
    nodes: Optional[str] = field(default_factory=str)
    connected: Optional[str] = field(default_factory=str)


def build_stack(source: str,
                version: str,
                stack_to_load: Dict[str, str],
                sa: bool = True,
                dl: bool = True) -> Dict[str, str]:
    loaded_stack = {}

    for key, file in stack_to_load.items():

        if sa:
            file = os.path.join(Settings.staging_area, source, version, file)
            loaded_stack[key] = file

        if dl:
            file = os.path.join(Settings.data_lake, source, version, file)
            loaded_stack[key] = file

    return loaded_stack


def build_multi_stack(source: str,
                      version: str,
                      stack_to_load: Dict[str, MultiStack],
                      sa: bool = True,
                      dl: bool = True) -> Dict[str, MultiStack]:
    loaded_stack = {}

    for key, multi_stack in stack_to_load.items():

        new_multi_stack = replace(multi_stack)

        new_stack = []
        for file in multi_stack.stack:

            if dl:
                file = os.path.join(Settings.data_lake, source, version, file)
                new_stack.append(file)
                continue

            if sa:
                file = os.path.join(Settings.staging_area, source, version, file)
                new_stack.append(file)

        new_multi_stack.stack = new_stack
        loaded_stack[key] = new_multi_stack

    return loaded_stack


def build_load_stack(source: str,
                     version: str,
                     stack_to_load: List[str]) -> List[str]:
    loaded_stack = []

    for file in stack_to_load:
        file = os.path.join(Settings.data_lake, source, version, f'{file}.json')
        loaded_stack.append(file)

    return loaded_stack