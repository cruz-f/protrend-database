import os

from protrend.utils import Settings


def build_file_path(source: str,
                    version: str,
                    file: str) -> str:
    return os.path.join(Settings.data_lake, source, version, file)
