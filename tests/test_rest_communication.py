__author__ = 'mwham'
from unittest.mock import patch, Mock
from analysis_driver.report_generation import rest_communication


class FakeRestReponse(Mock):
    def json(self):
        return self.content


def rest_url(endpoint):
    return 'http://localhost:4999/api/0.1/' + endpoint + '/'


test_endpoint = 'an_endpoint'
test_request_content = {'data': ['some', {'test': 'content'}]}


def patched_response(status_code=200):
    return patch(
        'requests.request',
        return_value=FakeRestReponse(status_code=status_code, content=test_request_content)
    )


@patched_response()
def test_req(mocked_instance):
    json = ['some', {'test': 'json'}]

    response = rest_communication._req('METHOD', rest_url(test_endpoint), json=json)
    assert response.status_code == 200
    assert response.content == test_request_content
    mocked_instance.assert_called_with('METHOD', rest_url(test_endpoint), json=json)


@patched_response()
def test_get_documents(mocked_instance):
    where = {'a_field': 'thing'}
    data = rest_communication.get_documents(test_endpoint, limit=100, **where)
    assert data == test_request_content['data']
    mocked_instance.assert_called_with('GET', rest_url(test_endpoint) + '?where={"a_field":"thing"}&max_results=100')


@patched_response()
def test_get_document(mocked_instance):
    where = {'a_field': 'thing'}

    assert rest_communication.get_document(test_endpoint, limit=100, **where) == test_request_content['data'][0]
    mocked_instance.assert_called_with('GET', rest_url(test_endpoint) + '?where={"a_field":"thing"}&max_results=100')


@patched_response()
def test_post_entry(mocked_instance):
    rest_communication.post_entry(test_endpoint, payload=test_request_content)
    mocked_instance.assert_called_with('POST', rest_url(test_endpoint), json=test_request_content)


@patched_response()
def test_put_entry(mocked_instance):
    rest_communication.put_entry(test_endpoint, 'an_element_id', payload=test_request_content)
    mocked_instance.assert_called_with('PUT', rest_url(test_endpoint) + 'an_element_id', json=test_request_content)


test_patch_document = {'_id': '1337', '_etag': 1234567, 'list_to_update': ['this', 'that', 'other']}


@patch('analysis_driver.report_generation.rest_communication.get_document', return_value=test_patch_document)
@patched_response()
def test_patch_entry(mocked_request, mocked_get_doc):
    patching_payload = {'_id': 1337, 'list_to_update': ['another']}
    rest_communication.patch_entry(test_endpoint, payload=patching_payload, update_lists=['list_to_update'])

    mocked_get_doc.assert_called_with(test_endpoint.rstrip('/'))
    mocked_request.assert_called_with(
        'PATCH',
        rest_url(test_endpoint) + '1337',
        headers={'If-Match': 1234567},
        json={'_id': 1337, 'list_to_update': ['another', 'other', 'that', 'this']}
    )


test_post_or_patch_payload = [
    {'_id': '1337', '_etag': 1234567, 'list_to_update': ['more'], 'another_field': 'this'}
]


def patched_post(success):
    return patch('analysis_driver.report_generation.rest_communication.post_entry', return_value=success)


def patched_patch(success):
    return patch('analysis_driver.report_generation.rest_communication.patch_entry', return_value=success)


def test_post_or_patch():
    with patched_post(True) as mocked_post, patched_patch(True) as mocked_patch:
        rest_communication.post_or_patch(
            'an_endpoint', test_post_or_patch_payload, elem_key='_id', update_lists=['list_to_update']
        )
        mocked_post.assert_called_with(test_endpoint, test_post_or_patch_payload[0])
        assert not mocked_patch.called

    with patched_post(False) as mocked_post, patched_patch(True) as mocked_patch:
        rest_communication.post_or_patch(
            'an_endpoint', test_post_or_patch_payload, elem_key='_id', update_lists=['list_to_update']
        )
        mocked_post.assert_called_with(test_endpoint, test_post_or_patch_payload[0])
        mocked_patch.assert_called_with(
            test_endpoint,
            {'_etag': 1234567, 'list_to_update': ['more'], 'another_field': 'this'},
            ['list_to_update'],
            _id='1337'
        )
