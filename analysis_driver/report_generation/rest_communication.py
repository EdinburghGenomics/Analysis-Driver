from urllib.parse import urljoin
import requests
from pprint import pformat
from analysis_driver.config import default as cfg
from analysis_driver.app_logging import get_logger

app_logger = get_logger(__name__)


def api_url(endpoint):
    return '{base_url}/{endpoint}/'.format(
        base_url=cfg.query('rest_api', 'url').rstrip('/'), endpoint=endpoint
    )


def _req(method, url, **kwargs):
    r = requests.request(method, url, **kwargs)
    if r.status_code != 200:
        app_logger.debug('%s %s %s %s' % (r.request.method, r.request.path_url, r.status_code, r.reason))
        json = r.json()
        if json:
            app_logger.debug(pformat(json))
    return r

    
def get_documents(endpoint, limit=10000, **kwargs):
    url = api_url(endpoint)
    param = []
    for key in kwargs:
        param.append('"%s":"%s"' % (key, kwargs.get(key)))
    if param:
        url += '?where={%s}' % ','.join(param)
    url += '&max_results=%s' % limit
    r = _req('GET', url)
    return r.json().get('data')


def get_document(endpoint, idx=0, **kwargs):
    documents = get_documents(endpoint, **kwargs)
    if documents:
        return documents[idx]
    else:
        app_logger.error('No document found for ' + endpoint + '. kwargs: ' + str(kwargs))


def post_entry(endpoint, payload):
    """Upload to the collection."""
    r = _req('POST', api_url(endpoint), json=payload)
    if r.status_code != 200:
        return False
    return True


def put_entry(endpoint, element_id, payload):
    """Upload Assuming we know the id of this entry"""
    r = _req('PUT', urljoin(api_url(endpoint), element_id), json=payload)
    if r.status_code != 200:
        return False
    return True


def patch_entry(endpoint, payload, update_lists=None, **kwargs):
    """Upload Assuming we can get the id of this entry from kwargs"""
    doc = get_document(endpoint, **kwargs)
    if doc:
        url = urljoin(api_url(endpoint), doc.get('_id'))
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
        app_logger.info('Updated %s documents matching %s' % (nb_docs, kwargs))
        return result
    return False


def post_or_patch(endpoint, input_json, elem_key=None, update_lists=None):
    """
    :param str endpoint:
    :param list input_json:
    :param str elem_key:
    """
    success = True
    for payload in input_json:
        if not post_entry(endpoint, payload):
            elem_query = {}
            if elem_key:
                elem_query = {elem_key: payload.pop(elem_key)}
            success = success and patch_entry(endpoint, payload, update_lists, **elem_query)
    return success
