from urllib.parse import urljoin
import requests
from pprint import pprint
from report_generation.model import ALL_PIECES, Info


def get_documents(url, **kwargs):
    param = []
    for key in kwargs:
        param.append('"%s":"%s"'%(key,kwargs.get(key)))
    if param: url = url + '?where={%s}'%','.join(param)
    r = requests.request('GET', url)
    return r.json().get('data')

def get_document(url, **kwargs):
    return get_documents(url, **kwargs)[0]


def post_entry(url, payload):
    """Upload to the collection."""
    r = requests.request('POST', url, json=payload)
    if r.status_code != 200:
        print('POST', r.status_code, r.reason, url)
        pprint(r.json())
        return False
    return True


def put_entry(url, id, payload):
    """Upload Assuming we know the id of this entry"""
    url = urljoin(url, id)
    r = requests.request('PUT', url, json=payload)
    if r.status_code != 200:
        print('PUT', r.status_code, r.reason, url)
        return False
    return True


def patch_entry(url, payload, **kwargs):
    """Upload Assuming we can get the id of this entry from kwargs"""
    doc = get_document(url.rstrip('/'), **kwargs)
    url = urljoin(url, doc.get('_id'))
    headers={'If-Match':doc.get('_etag')}
    r = requests.request('PATCH', url, headers=headers, json=payload)
    if r.status_code != 200:
        print('PATCH', r.status_code, r.reason, url)
        return False
    return True

def json_to_info(json):
    key_to_piece = {}
    for piece in ALL_PIECES:
        key_to_piece[piece.key] = piece
    info = Info()

    for key in json:
        if key in key_to_piece:
            info[key_to_piece.get(key)]=json.get(key)
    return info