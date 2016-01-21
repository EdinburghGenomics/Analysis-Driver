__author__ = 'mwham'
from random import randint
from datetime import datetime
import json as j


DB = {}
endpoints = {
    'analysis_driver_procs': 'proc_id',
    'lanes': 'lane_id',
    'projects': 'project_id',
    'run_elements': 'run_element_id',
    'runs': 'run_id',
    'samples': 'sample_id',
    'unexpected_barcodes': 'run_element_id'
}


def fake_request(req_type, url, **kwargs):
    if req_type == 'GET':
        r = fake_get(url)
    elif req_type == 'POST':
        r = fake_post(url, kwargs['json'])
    elif req_type == 'PATCH':
        r = fake_patch(url, kwargs['headers'], kwargs['json'])
    else:
        r = None
    return r


def fake_patch(url, headers, json):
    endpoint = url.rstrip('/').split('/')[-1]
    data = DB[endpoint]['data']
    docs = [doc for doc in data if doc['_etag'] == headers['If-Match']]
    assert len(docs) == 1
    doc = docs[0]
    data.remove(doc)
    doc.update(json)
    data.append(doc)
    assert doc in DB[endpoint]['data']
    return FakeResponse({'_msg': 'patched entry'})


def fake_get(url):
    if '?' in url:
        url, query_string = url.split('?')
    else:
        query_string = None

    endpoint = url.split('/')[-1]
    json_content = dict(DB[endpoint])
    assert json_content == DB[endpoint] and json_content is not DB[endpoint]

    if query_string:
        queries = query_string.split('&')
        for q in queries:
            if 'where=' in q:
                for k, v in j.loads(q.replace('where=', '')).items():
                    json_content['data'] = [x for x in json_content['data'] if x[k] == v]
            if 'sort=' in q:
                sortcol = q.replace('sort=', '')
                reverse = sortcol.startswith('-')
                sortcol = sortcol.lstrip('-')

                def _itemget(e):
                    return e.get(sortcol, '0')

                json_content['data'] = sorted(
                    json_content['data'],
                    key=_itemget,
                    reverse=reverse
                )

    return FakeResponse(json_content)


def fake_post(url, json):
    endpoint = url.rstrip('/').split('/')[-1]
    identifier = endpoints[endpoint]
    data = DB[endpoint]['data']
    if any([doc for doc in data if doc[identifier] == json[identifier]]):
        return FakeResponse({'_error': {'code': 202}})
    else:
        json['_etag'] = str(randint(1000, 9000))
        json['_created'] = str(datetime.utcnow().strftime('%d_%m_%Y_%H:%M:%S'))
        data.append(json)
        assert json in DB[endpoint]['data']
        return FakeResponse({'_msg': 'posted entry'})


class FakeRequest:
    method = 'a method'
    path_url = 'a path url'


class FakeResponse:
    def __init__(self, json_content=None):
        if json_content is None:
            json_content = {}
        self.content = json_content

    def json(self):
        return self.content

    @property
    def status_code(self):
        error_code = self.content.get('_error', {}).get('code')
        if error_code:
            return error_code
        else:
            return 200

    reason = 'a reason'
    request = FakeRequest()
