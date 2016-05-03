from unittest.mock import patch, Mock
import json
from tests import test_analysisdriver
from analysis_driver import rest_communication


helper = test_analysisdriver.TestAnalysisDriver()


class FakeRestResponse(Mock):
    def __init__(self, *args, **kwargs):
        content = kwargs.pop('content', None)
        if type(content) in (list, dict):
            content = json.dumps(content)
        content = content.encode()
        super().__init__(*args, **kwargs)
        self.content = content
        self.request = Mock(method='a method', path_url='a url')
        self.status_code = 200
        self.reason = 'a reason'

    def json(self):
        return json.loads(self.content.decode('utf-8'))

    @property
    def text(self):
        return self.content.decode('utf-8')


def rest_url(endpoint):
    return 'http://localhost:4999/api/0.1/' + endpoint + '/'


test_endpoint = 'an_endpoint'
test_request_content = {'data': ['some', {'test': 'content'}]}


patched_response = patch(
    'requests.request',
    return_value=FakeRestResponse(status_code=200, content=test_request_content)
)


def test_api_url_query_strings():
    assert rest_communication.api_url('an_endpoint') == rest_url('an_endpoint')
    exp = '?where={"this":"that"}&embedded={"things":1}&aggregate=True&sort=-_created'
    obs = rest_communication.api_url(
        'an_endpoint',
        where={'this': 'that'},
        embedded={'things': 1},
        aggregate=True,
        sort='-_created'
    ).replace(rest_url('an_endpoint'), '')
    assert sorted(obs.lstrip('?').split('&')) == sorted(exp.lstrip('?').split('&'))


def test_parse_query_string():
    query_string = 'http://a_url?this=that&other={"another":"more"}'
    no_query_string = 'http://a_url'
    dodgy_query_string = 'http://a_url?this=that?other=another'

    assert rest_communication.parse_query_string(query_string) == {'this': 'that', 'other': '{"another":"more"}'}
    assert rest_communication.parse_query_string(no_query_string) is None

    try:
        rest_communication.parse_query_string(dodgy_query_string)
    except AssertionError as e:
        assert e.args[0] == 'Bad query string: ' + dodgy_query_string

    try:
        rest_communication.parse_query_string(query_string, requires=['things'])
    except AssertionError as e2:
        assert e2.args[0] == query_string + ' did not contain all required fields: ' + str(['things'])


@patched_response
def test_req(mocked_instance):
    json_content = ['some', {'test': 'json'}]

    response = rest_communication._req('METHOD', rest_url(test_endpoint), json=json_content)
    assert response.status_code == 200
    assert json.loads(response.content.decode('utf-8')) == response.json() == test_request_content
    mocked_instance.assert_called_with('METHOD', rest_url(test_endpoint), json=json_content)


def test_get_documents_depaginate():
    docs = (
        FakeRestResponse(content={'data': ['this', 'that'], '_links': {'next': {'href': 'an_endpoint?max_results=101&page=2'}}}),
        FakeRestResponse(content={'data': ['other', 'another'], '_links': {'next': {'href': 'an_endpoint?max_results=101&page=3'}}}),
        FakeRestResponse(content={'data': ['more', 'things'], '_links': {}})
    )
    patched_req = patch(
        'analysis_driver.rest_communication._req',
        side_effect=docs
    )
    with patched_req as mocked_req:
        assert rest_communication.get_documents('an_endpoint', depaginate=True, max_results=101) == [
            'this', 'that', 'other', 'another', 'more', 'things'
        ]
        assert all([a[0][1].startswith(rest_url('an_endpoint')) for a in mocked_req.call_args_list])
        assert [helper.query_args_from_url(a[0][1]) for a in mocked_req.call_args_list] == [
            {'page': '1', 'max_results': '101'},
            {'page': '2', 'max_results': '101'},
            {'page': '3', 'max_results': '101'}
        ]


@patched_response
def test_get_documents(mocked_instance):
    data = rest_communication.get_documents(test_endpoint, max_results=100, where={'a_field': 'thing'})
    assert data == test_request_content['data']
    assert mocked_instance.call_args[0][1].startswith(rest_url(test_endpoint))
    assert helper.query_args_from_url(mocked_instance.call_args[0][1]) == {
        'max_results': '100', 'where': {'a_field': 'thing'}, 'page': '1'
    }


@patched_response
def test_get_document(mocked_instance):
    assert rest_communication.get_document(test_endpoint, max_results=100, where={'a_field': 'thing'}) == test_request_content['data'][0]
    assert mocked_instance.call_args[0][1].startswith(rest_url(test_endpoint))
    assert helper.query_args_from_url(mocked_instance.call_args[0][1]) == {
        'max_results': '100', 'where': {'a_field': 'thing'}, 'page': '1'
    }


@patched_response
def test_post_entry(mocked_instance):
    rest_communication.post_entry(test_endpoint, payload=test_request_content)
    mocked_instance.assert_called_with('POST', rest_url(test_endpoint), json=test_request_content)


@patched_response
def test_put_entry(mocked_instance):
    rest_communication.put_entry(test_endpoint, 'an_element_id', payload=test_request_content)
    mocked_instance.assert_called_with('PUT', rest_url(test_endpoint) + 'an_element_id', json=test_request_content)


test_patch_document = {
    '_id': '1337', '_etag': 1234567, 'uid': 'a_unique_id', 'list_to_update': ['this', 'that', 'other']
}


@patch('analysis_driver.rest_communication.get_document', return_value=test_patch_document)
@patched_response
def test_patch_entry(mocked_request, mocked_get_doc):
    patching_payload = {'list_to_update': ['another']}
    rest_communication.patch_entry(
        test_endpoint,
        payload=patching_payload,
        id_field='uid',
        element_id='a_unique_id',
        update_lists=['list_to_update']
    )

    mocked_get_doc.assert_called_with(test_endpoint, where={'uid': 'a_unique_id'})
    mocked_request.assert_called_with(
        'PATCH',
        rest_url(test_endpoint) + '1337',
        headers={'If-Match': 1234567},
        json={'list_to_update': ['this', 'that', 'other', 'another']}
    )


test_post_or_patch_payload = {'uid': '1337', 'list_to_update': ['more'], 'another_field': 'that'}
test_post_or_patch_payload_no_uid = {'list_to_update': ['more'], 'another_field': 'that'}
test_post_or_patch_doc = {
    'uid': 'a_uid', '_id': '1337', '_etag': 1234567, 'list_to_update': ['things'], 'another_field': 'this'
}


def patched_post(success):
    return patch('analysis_driver.rest_communication.post_entry', return_value=success)


def patched_patch(success):
    return patch('analysis_driver.rest_communication._patch_entry', return_value=success)


def patched_get(payload):
    return patch('analysis_driver.rest_communication.get_document', return_value=payload)


def test_post_or_patch():
    with patched_get(test_post_or_patch_doc) as mget, patched_patch(True) as mpatch:
        success = rest_communication.post_or_patch(
            'an_endpoint',
            [test_post_or_patch_payload],
            id_field='uid',
            update_lists=['list_to_update']
        )
        mget.assert_called_with('an_endpoint', where={'uid': '1337'})
        mpatch.assert_called_with(
            'an_endpoint',
            test_post_or_patch_doc,
            test_post_or_patch_payload_no_uid,
            ['list_to_update']
        )
        assert success is True

    with patched_get(None) as mget, patched_post(True) as mpost:
        success = rest_communication.post_or_patch(
            'an_endpoint', [test_post_or_patch_payload], id_field='uid', update_lists=['list_to_update']
        )
        mget.assert_called_with('an_endpoint', where={'uid': '1337'})
        mpost.assert_called_with('an_endpoint', test_post_or_patch_payload)
        assert success is True
