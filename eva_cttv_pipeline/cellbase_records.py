import codecs
import json
import urllib.request
import urllib.error

from eva_cttv_pipeline import config, utilities


class CellbaseRecords:

    """Assists in the requesting and iteration of clinvar cellbase records"""

    def __init__(self, limit=config.BATCH_SIZE, skip=0, json_file=None):
        self.skip = skip
        self.limit = limit
        self.json_file = json_file

    def __iter__(self):
        if not self.json_file:
            while True:
                curr_result_list = self.__get_curr_result_list()
                if not curr_result_list:
                    break
                for record in curr_result_list:
                    yield record
                self.skip += config.BATCH_SIZE
        else:
            for record in self.__each_line_in_file():
                yield record

    def __get_curr_response(self):
        reader = codecs.getreader("utf-8")
        url = 'http://{}/cellbase/webservices/rest/v3/hsapiens/feature/clinical/' \
              'all?source=clinvar&skip={}&limit={}'.format(config.HOST, self.skip, self.limit)
        for _ in range(3):
            while True:
                try:
                    answer = urllib.request.urlopen(url)
                except Exception:
                    continue
                raise urllib.error.HTTPError

        curr_response = json.load(reader(answer))['response'][0]
        return curr_response

    def __get_curr_result_list(self):
        curr_response = self.__get_curr_response()
        curr_result_list = curr_response['result']
        if len(curr_result_list) == 0:
            return None
        return curr_result_list

    def __each_line_in_file(self):
        with utilities.open_file(self.json_file, "rt") as f:
            for line in f:
                yield line.rstrip()
