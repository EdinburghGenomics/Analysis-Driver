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

    
def get_documents(endpoint, limit=10000, **query_args):
    url = api_url(endpoint) + '?max_results=%s' % limit
    q_string = '&'.join(['%s=%s' % (k, v) for k, v in query_args.items()]).replace(' ', '').replace('\'', '"')
    if q_string:
        url += '&' + q_string
    r = _req('GET', url)
    return r.json().get('data')


def get_document(endpoint, idx=0, **query_args):
    documents = get_documents(endpoint, **query_args)
    if documents:
        return documents[idx]
    else:
        app_logger.error('No document found for ' + endpoint + '. kwargs: ' + str(query_args))


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


def patch_entry(endpoint, payload, id_field, element_id, update_lists=None):
    """Upload Assuming we can get the id of this entry from kwargs"""
    doc = get_document(endpoint, where={id_field: element_id})
    if doc:
        url = urljoin(api_url(endpoint), doc.get('_id'))
        headers = {'If-Match': doc.get('_etag')}
        if update_lists:
            for l in update_lists:
                content = doc.get(l, [])
                new_content = [x for x in payload.get(l, []) if x not in content]
                payload[l] = content + new_content
        r = _req('PATCH', url, headers=headers, json=payload)
        if r.status_code == 200:
            return True
    return False


def patch_entries(endpoint, payload, update_lists=None, **kwargs):
    """Apply the same upload to all the documents retrieved using  **kwargs"""
    docs = get_documents(endpoint, **kwargs)
    if docs:
        result = True
        nb_docs = 0
        for doc in docs:
            url = urljoin(api_url(endpoint), doc.get('_id'))
            headers = {'If-Match': doc.get('_etag')}
            if update_lists:
                for l in update_lists:
                    content = doc.get(l, [])
                    new_content = [x for x in payload.get(l, []) if x not in content]
                    payload[l] = content + new_content
            r = _req('PATCH', url, headers=headers, json=payload)
            if r.status_code != 200:
                result = False
            else:
                nb_docs += 1
        app_logger.info('Updated %s documents matching %s' % (nb_docs, kwargs))
        return result
    return False


def post_or_patch(endpoint, input_json, id_field=None, update_lists=None):
    """
    :param str endpoint:
    :param list input_json:
    :param str id_field:
    """
    success = True
    for payload in input_json:
        if get_document(endpoint, where={id_field: payload[id_field]}):
            elem_key = payload.pop(id_field)
            s = patch_entry(endpoint, payload, id_field, elem_key, update_lists)
        else:
            s = post_entry(endpoint, payload)
        success = success and s
    return success
