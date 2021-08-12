import json
import os

import pandas as pd


def read_json(file_path: str) -> dict:

    if not os.path.exists(file_path):
        raise OSError(f'Invalid file path {file_path}')

    with open(file_path, 'r') as handle:
        return json.load(handle)


def write_json(file_path: str, content: dict):

    with open(file_path, 'w') as handle:
        return json.dump(content, handle)


def read_json_lines(file_path: str) -> pd.DataFrame:

    if not os.path.exists(file_path):
        raise OSError(f'Invalid file path {file_path}')

    lines = []

    with open(file_path, 'r') as file:

        for line in file:
            line = line.rstrip()

            if line:
                line_dict = json.loads(line)
                lines.append(line_dict)

    return pd.DataFrame(lines)