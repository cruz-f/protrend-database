import os
from dataclasses import dataclass, field
from typing import Optional

from protrend.utils import Settings


@dataclass
class Stack:
    transformed: Optional[str] = field(default_factory=str)
    integrated: Optional[str] = field(default_factory=str)
    nodes: Optional[str] = field(default_factory=str)
    connected: Optional[str] = field(default_factory=str)


def build_file_path(source: str,
                    version: str,
                    file: str) -> str:
    return os.path.join(Settings.data_lake, source, version, file)
