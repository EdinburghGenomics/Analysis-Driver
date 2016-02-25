__author__ = 'mwham'
from unittest.mock import patch, Mock
from analysis_driver import rest_communication


class FakeRestReponse(Mock):
    def json(self):
        return self.content


def rest_url(endpoint):
    return 'http://localhost:4999/api/0.1/' + endpoint + '/'


test_endpoint = 'an_endpoint'
test_request_content = {'data': ['some', {'test': 'content'}]}


patched_response = patch(
    'requests.request',
    return_value=FakeRestReponse(status_code=200, content=test_request_content)
)


@patched_response
def test_req(mocked_instance):
    json = ['some', {'test': 'json'}]

    response = rest_communication._req('METHOD', rest_url(test_endpoint), json=json)
    assert response.status_code == 200
    assert response.content == test_request_content
    mocked_instance.assert_called_with('METHOD', rest_url(test_endpoint), json=json)


@patched_response
def test_get_documents(mocked_instance):
    data = rest_communication.get_documents(test_endpoint, limit=100, where={'a_field': 'thing'})
    assert data == test_request_content['data']
    mocked_instance.assert_called_with('GET', rest_url(test_endpoint) + '?max_results=100&where={"a_field":"thing"}')


@patched_response
def test_get_document(mocked_instance):
    assert rest_communication.get_document(test_endpoint, limit=100, where={'a_field': 'thing'}) == test_request_content['data'][0]
    mocked_instance.assert_called_with('GET', rest_url(test_endpoint) + '?max_results=100&where={"a_field":"thing"}')


@patched_response
def test_post_entry(mocked_instance):
    rest_communication.post_entry(test_endpoint, payload=test_request_content)
    mocked_instance.assert_called_with('POST', rest_url(test_endpoint), json=test_request_content)


@patched_response
def test_put_entry(mocked_instance):
    rest_communication.put_entry(test_endpoint, 'an_element_id', payload=test_request_content)
    mocked_instance.assert_called_with('PUT', rest_url(test_endpoint) + 'an_element_id', json=test_request_content)


test_patch_document = {'_id': '1337', '_etag': 1234567, 'list_to_update': ['this', 'that', 'other']}


@patch('analysis_driver.rest_communication.get_document', return_value=test_patch_document)
@patched_response
def test_patch_entry(mocked_request, mocked_get_doc):
    patching_payload = {'list_to_update': ['another']}
    rest_communication.patch_entry(
        test_endpoint,
        payload=patching_payload,
        id_field='_id',
        element_id=1337,
        update_lists=['list_to_update']
    )

    mocked_get_doc.assert_called_with(test_endpoint, where={'_id': 1337})
    mocked_request.assert_called_with(
        'PATCH',
        rest_url(test_endpoint) + '1337',
        headers={'If-Match': 1234567},
        json={'list_to_update': ['this', 'that', 'other', 'another']}
    )


test_post_or_patch_payload = {
    '_id': '1337', '_etag': 1234567, 'list_to_update': ['more'], 'another_field': 'this'
}


def patched_post(success):
    return patch('analysis_driver.rest_communication.post_entry', return_value=success)


def patched_patch(success):
    return patch('analysis_driver.rest_communication.patch_entry', return_value=success)


def patched_get(payload):
    return patch('analysis_driver.rest_communication.get_document', return_value=payload)


def test_post_or_patch():
    with patched_get(test_post_or_patch_payload) as mget, patched_patch(True) as mpatch:
        success = rest_communication.post_or_patch(
            'an_endpoint', [dict(test_post_or_patch_payload)], id_field='_id', update_lists=['list_to_update']
        )
        mget.assert_called_with('an_endpoint', where={'_id': '1337'})
        mpatch.assert_called_with(
            'an_endpoint',
            {'_etag': 1234567, 'list_to_update': ['more'], 'another_field': 'this'},
            '_id',
            '1337',
            ['list_to_update']
        )
        assert success is True

    with patched_get(None) as mget, patched_post(True) as mpost:
        success = rest_communication.post_or_patch(
            'an_endpoint', [dict(test_post_or_patch_payload)], id_field='_id', update_lists=['list_to_update']
        )
        mget.assert_called_with('an_endpoint', where={'_id': '1337'})
        mpost.assert_called_with('an_endpoint', test_post_or_patch_payload)
        assert success is True
