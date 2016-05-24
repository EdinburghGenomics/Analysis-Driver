import json
from unittest.mock import Mock


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
