from urllib.parse import urljoin
import requests
from analysis_driver.config import default as cfg
from analysis_driver.app_logging import get_logger

app_logger = get_logger(__name__)


def api_url(endpoint, **query_args):
    url = '{base_url}/{endpoint}/'.format(
        base_url=cfg.query('rest_api', 'url').rstrip('/'), endpoint=endpoint
    )
    if query_args:
        url += '?' + '&'.join(['%s=%s' % (k, v) for k, v in query_args.items()]).replace(' ', '').replace('\'', '"')

    return url


def _req(method, url, **kwargs):
    r = requests.request(method, url, **kwargs)
    # e.g: 'POST <url> ({"some": "args"}) -> {"some": "content"}. Status code 201. Reason: CREATED
    report = '%s %s (%s) -> %s. Status code %s. Reason: %s' % (
        r.request.method, r.request.path_url, kwargs, r.content.decode('utf-8'), r.status_code, r.reason
    )
    if r.status_code in (200, 201):
        app_logger.debug(report)
    else:
        app_logger.error(report)
    return r


def depaginate_documents(endpoint, **queries):
    elements = []
    page_size = queries.pop('max_results', 100)
    url = api_url(endpoint, max_results=page_size, **queries)
    content = _req('GET', url).json()
    elements.extend(content['data'])

    if 'next' in content['_links']:
        next_href, next_query = content['_links']['next']['href'].split('?')
        next_query = dict([x.split('=') for x in next_query.split('&')])
        queries.pop('page', None)
        next_query.update(queries)
        elements.extend(depaginate_documents(next_href, **next_query))
    return elements

    
def get_documents(endpoint, limit=10000, **query_args):
    if 'max_results' not in query_args:
        query_args['max_results'] = limit
    url = api_url(endpoint, **query_args)
    r = _req('GET', url)
    return r.json().get('data')


def get_document(endpoint, idx=0, **query_args):
    documents = get_documents(endpoint, **query_args)
    if documents:
        return documents[idx]
    else:
        app_logger.warning('No document found in endpoint %s for %s' % (endpoint, str(query_args)))


def post_entry(endpoint, payload):
    """Post a new entry to the endpoint."""
    r = _req('POST', api_url(endpoint), json=payload)
    if r.status_code != 200:
        return False
    return True


def put_entry(endpoint, element_id, payload):
    """Upload, assuming we know the id of the entry."""
    r = _req('PUT', urljoin(api_url(endpoint), element_id), json=payload)
    if r.status_code != 200:
        return False
    return True


def _patch_entry(endpoint, doc, payload, update_lists=None):
    """
    Patch a specific database item (specified by doc) with the given data payload.
    :param str endpoint:
    :param dict doc: The entry in the database to patch (contains the relevant _id and _etag)
    :param dict payload: Data with which to patch doc
    :param list update_lists: Doc items listed here will be appended rather than replaced by the patch
    """
    url = urljoin(api_url(endpoint), doc.get('_id'))
    _payload = dict(payload)
    headers = {'If-Match': doc.get('_etag')}
    if update_lists:
        for l in update_lists:
            content = doc.get(l, [])
            new_content = [x for x in _payload.get(l, []) if x not in content]
            _payload[l] = content + new_content
    r = _req('PATCH', url, headers=headers, json=_payload)
    if r.status_code == 200:
        return True
    return False


def patch_entry(endpoint, payload, id_field, element_id, update_lists=None):
    """
    Retrieve a document at the given endpoint with the given unique ID, and patch it with some data.
    :param str endpoint:
    :param str payload:
    :param str id_field: The name of the unique identifier (e.g. 'run_element_id', 'proc_id', etc.)
    :param str element_id: The value of id_field to retrieve (e.g. '160301_2_ATGCATGC')
    """
    doc = get_document(endpoint, where={id_field: element_id})
    if doc:
        return _patch_entry(endpoint, doc, payload, update_lists)
    return False


def patch_entries(endpoint, payload, update_lists=None, **query_args):
    """
    Retrieve many documents and patch them all with the same data.
    :param str endpoint:
    :param str payload:
    :param query_args: Database query args to pass to get_documents
    """
    docs = get_documents(endpoint, **query_args)
    if docs:
        success = True
        nb_docs = 0
        for doc in docs:
            if _patch_entry(endpoint, doc, payload, update_lists):
                nb_docs += 1
            else:
                success = False
        app_logger.info('Updated %s documents matching %s' % (nb_docs, query_args))
        return success
    return False


def post_or_patch(endpoint, input_json, id_field=None, update_lists=None):
    """
    For each document supplied, either post to the endpoint if the unique id doesn't yet exist there, or
    patch if it does.
    :param str endpoint:
    :param input_json: A single document or list of documents to post or patch to the endpoint.
    :param str id_field: The field to use as the unique ID for the endpoint.
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
