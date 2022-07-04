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


def fetch_sequences_to_ncbi_ftp(protrend_ids: List[str], ftp_paths: List[str], output_folder: str):
    """
    Fetch sequences from NCBI FTP server.
    :param ftp_paths: List of FTP paths.
    :param output_folder: Output folder.
    :param protrend_ids: List of protrend IDs.
    :return: List of output files.
    """
    if not protrend_ids:
        protrend_ids = [None] * len(ftp_paths)

    output_folder = Path(output_folder)

    output_sequences = []
    for ftp_path, protrend_id in zip(ftp_paths, protrend_ids):
        file_path = output_folder.joinpath(f'{protrend_id}.gbff.gz')
        fetch_sequences(ftp_path=ftp_path, file_name_or_path=file_path)
        output_sequences.append(file_path)
        time.sleep(Settings.request_sleep)

    return output_sequences


if __name__ == '__main__':
    fetch_sequences_to_ncbi_ftp(
        ftp_paths=['ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/021/485/GCA_000021485.1_ASM2148v1'],
        output_folder=r'C:\Users\BiSBII\OneDrive - Universidade do Minho\PhD\Protrend\main\protrend-database\useless',
        protrend_ids=['PRT.ORG.0000001']
    )
