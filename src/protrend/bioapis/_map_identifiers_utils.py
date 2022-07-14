import re
import time
from typing import List

import requests
from requests.adapters import HTTPAdapter, Retry

from protrend.log import ProtrendLogger
from protrend.utils import request


NEXT_LINK_PATTERN = re.compile(r'<(.+)>; rel="next"')

POLLING_INTERVAL = 3

MAPPING_TERMS = ('UniProtKB_AC-ID',
                 'RefSeq_Protein',
                 'RefSeq_Nucleotide',
                 'GI_number',
                 'GeneID')
API_URL = "https://rest.uniprot.org"

retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def submit_id_mapping(from_db: str, to_db: str, ids: List[str]):

    if from_db not in MAPPING_TERMS:
        raise ValueError(f'Invalid from {from_db}')

    if to_db not in MAPPING_TERMS:
        raise ValueError(f'Invalid from {to_db}')

    url = f"{API_URL}/idmapping/run"
    data = {"from": from_db, "to": to_db, "ids": ",".join(ids)}
    response = request(method='post', url=url, data=data)
    try:
        return response.json().get("jobId")
    except Exception as e:
        ProtrendLogger.log.info(f"Failed to submit id mapping: {e}")
        return


def get_next_link(headers):
    if "Link" in headers:
        match = NEXT_LINK_PATTERN.match(headers["Link"])
        if match:
            return match.group(1)


def check_id_mapping_status(job_id):
    if not job_id:
        ProtrendLogger.log.info(f"Job {job_id} failed")
        return False

    while True:
        response = session.get(f"{API_URL}/idmapping/status/{job_id}")
        job = response.json()
        if "jobStatus" in job:
            if job["jobStatus"] == "RUNNING":
                ProtrendLogger.log.info(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                ProtrendLogger.log.info(f"Job {job_id} failed")
        else:
            return bool(job["results"] or job["failedIds"])


def get_batch(batch_response):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        yield batch_response.json()
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results):
    for key in ("results", "failedIds"):
        if key in batch_results and batch_results[key]:
            all_results[key] += batch_results[key]
    return all_results


def get_id_mapping_results_link(job_id):
    url = f"{API_URL}/idmapping/details/{job_id}"
    response = session.get(url)
    return response.json().get("redirectURL")


def get_id_mapping_results_search(url):
    response = session.get(url)
    results = response.json()
    for batch in get_batch(response):
        results = combine_batches(results, batch)
    return results
