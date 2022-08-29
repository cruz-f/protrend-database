import json
from pathlib import Path
from typing import Union


def read_report(file_path: Union[str, Path]) -> dict:
    """
    It reads the report
    """
    with open(file_path, 'r') as f:
        report = json.load(f)

    return report
