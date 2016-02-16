from urllib.parse import urljoin
import requests
from pprint import pformat
from analysis_driver.config import default as cfg
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


def get_documents(url, limit=10000, **kwargs):
    param = []
    for key in kwargs:
        param.append('"%s":"%s"' % (key, kwargs.get(key)))
    if param:
        url += '?where={%s}' % ','.join(param)
    url += '&max_results=%s' % limit
    r = _req('GET', url)
    return r.json().get('data')


def get_document(url, idx=0, **kwargs):
    documents = get_documents(url, **kwargs)
    if documents:
        return documents[idx]
    else:
        app_logger.error('No document found for ' + url + ' kwargs=' + str(kwargs))


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
                payload[l] = sorted(set(doc.get(l, []) + payload.get(l, [])))
        r = _req('PATCH', url, headers=headers, json=payload)
        if r.status_code == 200:
            return True
    return False

def patch_entries(url, payload, update_lists=None, **kwargs):
    """Apply the same upload to all the documents retrieved using  **kwargs"""
    docs = get_documents(url.rstrip('/'), **kwargs)
    if docs:
        result = True
        nb_docs = 0
        for doc in docs:
            url = urljoin(url, doc.get('_id'))
            headers = {'If-Match': doc.get('_etag')}
            if update_lists:
                for l in update_lists:
                    payload[l] = sorted(set(doc.get(l, []) + payload.get(l, [])))
            r = _req('PATCH', url, headers=headers, json=payload)
            if r.status_code != 200:
                result = False
            else:
                nb_docs += 1
        app_logger.info("Updated %s documents matching %s"%(nb_docs, kwargs))
        return result
    return False

def post_or_patch(endpoint, input_json, elem_key=None, update_lists=None):
    """
    :param str endpoint:
    :param list input_json:
    :param str elem_key:
    """
    url = '/'.join((cfg.query('rest_api', 'url').rstrip('/'), endpoint, ''))
    success = True
    for payload in input_json:
        if not post_entry(url, payload):
            elem_query = {}
            if elem_key:
                elem_query = {elem_key: payload.pop(elem_key)}
            success_tmp = patch_entry(url, payload, update_lists, **elem_query)
            success = success and success_tmp
    return success
