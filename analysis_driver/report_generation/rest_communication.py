import logging
from urllib.parse import urljoin
import requests
from pprint import pformat

app_logger = logging.getLogger(__name__)

def _req(*args, **kwargs):
    r = requests.request(*args, **kwargs)
    if r.status_code != 200:
        app_logger.debug('%s %s %s'%(r.request, r.status_code, r.reason))
        json = r.json()
        if json:
            app_logger.debug(pformat(json))
    return r

def get_documents(url, **kwargs):
    param = []
    for key in kwargs:
        param.append('"%s":"%s"'%(key,kwargs.get(key)))
    if param: url = url + '?where={%s}'%','.join(param)
    r = _req('GET', url)
    return r.json().get('data')

def get_document(url, **kwargs):
    return get_documents(url, **kwargs)[0]


def post_entry(url, payload):
    """Upload to the collection."""
    r = _req('POST', url, json=payload)
    if r.status_code != 200:
        return False
    return True


def put_entry(url, id, payload):
    """Upload Assuming we know the id of this entry"""
    url = urljoin(url, id)
    r = _req('PUT', url, json=payload)
    if r.status_code != 200:
        return False
    return True

def patch_entry(url, payload, **kwargs):
    """Upload Assuming we can get the id of this entry from kwargs"""
    doc = get_document(url.rstrip('/'), **kwargs)
    url = urljoin(url, doc.get('_id'))
    headers={'If-Match':doc.get('_etag')}
    r = _req('PATCH', url, headers=headers, json=payload)
    if r.status_code != 200:
        return False
    return True
