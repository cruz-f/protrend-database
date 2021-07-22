import json
import os

import pandas as pd


def read_json_lines(file_path: str) -> pd.DataFrame:

    if not os.path.exists(file_path):
        raise OSError(f'Invalid file path {file_path}')

    lines = []

    with open(file_path, 'r') as file:
        for line in file:

            line_dict = json.loads(line)

            lines.append(line_dict)

    return pd.DataFrame(lines)