import io
import time

import pandas as pd
import requests

from .settings import Settings


def _request(url: str,
             method: str = 'get',
             params: dict = None,
             **kwargs):

    requests_kwargs = dict(data=None,
                           headers=None,
                           cookies=None,
                           files=None,
                           auth=None,
                           timeout=Settings.request_timeout,
                           allow_redirects=True,
                           proxies=None,
                           hooks=None,
                           stream=None,
                           verify=None,
                           cert=None,
                           json=None)

    for key, val in kwargs.items():
        requests_kwargs[key] = val

    try:

        response = requests.request(method=method,
                                    url=url,
                                    params=params,
                                    **requests_kwargs)

        try:
            response.raise_for_status()

        except requests.exceptions.HTTPError:
            response = None

    except requests.exceptions.RequestException:
        response = None

    if response is None:

        response = requests.Response()
        response.status_code = 404

    return response


def request(url: str,
            method: str = 'get',
            params: dict = None,
            retries: int = Settings.request_retries,
            **kwargs) -> requests.Response:

    if retries > 0:

        for i in range(retries):

            response = _request(url=url,
                                method=method,
                                params=params,
                                **kwargs)

            if response.status_code == 404:
                time.sleep(2)

            else:
                return response

        response = requests.Response()
        response.status_code = 404

        return response

    else:

        return _request(url=url,
                        method=method,
                        params=params,
                        **kwargs)


def read_response(response: requests.Response, **kwargs) -> pd.DataFrame:
    data = response.content

    if data:

        buffer = io.StringIO(data.decode('utf-8'))
        return pd.read_csv(buffer, **kwargs)

    return pd.DataFrame()
