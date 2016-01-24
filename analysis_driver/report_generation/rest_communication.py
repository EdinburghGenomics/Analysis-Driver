from urllib.parse import urljoin
import requests
from pprint import pformat

from analysis_driver.app_logging import get_logger

app_logger = get_logger(__name__)


def _req(*args, **kwargs):
    r = requests.request(*args, **kwargs)
    if r.status_code != 200:
        app_logger.debug('%s %s %s %s' % (r.request.method, r.request.path_url, r.status_code, r.reason))
        json = r.json()
        if json:
            app_logger.debug(pformat(json))
    return r


def get_documents(url, **kwargs):
    param = []
    for key in kwargs:
        param.append('"%s":"%s"' % (key, kwargs.get(key)))
    if param:
        url += '?where={%s}' % ','.join(param)
    r = _req('GET', url)
    return r.json().get('data')


def get_document(url, **kwargs):
    documents = get_documents(url, **kwargs)
    if len(documents)>0:
        return documents[0]
    else:
        app_logger.error('No document found for ' + url + ' kwargs='+ str(kwargs))
        return None


def post_entry(url, payload):
    """Upload to the collection."""
    r = _req('POST', url, json=payload)
    if r.status_code != 200:
        return False
    return True


def put_entry(url, element_id, payload):
    """Upload Assuming we know the id of this entry"""
    url = urljoin(url, element_id)
    r = _req('PUT', url, json=payload)
    if r.status_code != 200:
        return False
    return True


def patch_entry(url, payload, update_lists=None, **kwargs):
    """Upload Assuming we can get the id of this entry from kwargs"""
    doc = get_document(url.rstrip('/'), **kwargs)
    if doc:
        url = urljoin(url, doc.get('_id'))
        headers = {'If-Match': doc.get('_etag')}
        if update_lists:
            for l in update_lists:
                payload[l] = list(set(payload.get(l, []) + doc.get(l, [])))
        r = _req('PATCH', url, headers=headers, json=payload)
        if r.status_code == 200:
            return True
    return False
