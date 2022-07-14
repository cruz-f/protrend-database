import ftplib
import time
from pathlib import Path
from typing import List

from protrend.utils import Settings


def fetch_sequences(ftp_path: str, file_name_or_path: Path):
    """
    Fetch sequences from NCBI FTP server.
    :param ftp_path: Path to the FTP directory.
    :param file_name_or_path: Output file path.
    :return: None
    """
    time.sleep(Settings.request_sleep)
    ftp_path = ftp_path.replace('ftp://ftp.ncbi.nlm.nih.gov/', '')
    accession = ftp_path.split('/')[-1]

    ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
    ftp.login()
    ftp.cwd(ftp_path)

    gb_file = None
    files = ftp.nlst()
    for file in files:

        if file.endswith(f'{accession}_genomic.gbff.gz'):
            gb_file = file
            continue

    if gb_file is not None:

        with open(file_name_or_path, 'wb') as f:
            ftp.retrbinary(f'RETR {gb_file}', f.write)

    ftp.quit()

