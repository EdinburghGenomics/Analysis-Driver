__author__ = 'mwham'
import os
import shutil
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver import process_trigger
from analysis_driver.dataset_scanner import RunScanner
from analysis_driver.config import default as cfg


class TestProcessTrigger(TestAnalysisDriver):
    dataset = 'test_dataset'
    scanner = RunScanner(cfg)
    @property
    def from_dir(self):
        return os.path.join(self.data_transfer, 'from')

    @property
    def to_dir(self):
        return os.path.join(self.data_transfer, 'to')

    @property
    def triggerignore(self):
        return os.path.join(self.from_dir, '.triggerignore')

    def setUp(self):
        #self.old_input_dir = cfg['input_dir']
        #self.old_int_dir = cfg.get('intermediate_dir')
        #self.old_lock_dir = cfg.get('lock_file_dir')
        #cfg.content['input_dir'] = self.from_dir
        #cfg.content['intermediate_dir'] = self.to_dir
        #cfg.content['lock_file_dir'] = self.from_dir
        #if not os.path.isdir(self.to_dir):
        #    os.mkdir(self.to_dir)
        #if os.path.isdir(os.path.join(self.to_dir, self.dataset)):
        #    shutil.rmtree(os.path.join(self.to_dir, self.dataset))
        #assert not os.listdir(self.to_dir)
        self.scanner.reset(self.dataset)

    def tearDown(self):
        #cfg.content['input_dir'] = self.old_input_dir
        #if self.old_int_dir:
        #    cfg.content['intermediate_dir'] = self.old_int_dir
        #else:
        #    cfg.content.pop('intermediate_dir')
        #cfg.content['lock_file_dir'] = self.old_lock_dir

        with open(os.path.join(self.from_dir, '.triggerignore'), 'w') as f:
            for d in ['dir_t?_be_ign*d\n', 'test_dataset\n']:
                f.write(d)


class TestTransfer(TestProcessTrigger):
    def test_transfer(self):
        self.scanner.reset(self.dataset)
        print(os.getcwd())
        print(os.path.exists("tests/assets/data_transfer/from/"))
        process_trigger._transfer_to_int_dir(self.dataset, self.from_dir, self.to_dir, 1)

        new_dataset = os.path.join(self.to_dir, self.dataset)
        observed = os.listdir(new_dataset)
        expected = ['RTAComplete.txt', 'thang', 'thing']
        print(observed, '\n', expected)
        assert observed == expected


class TestDatasetScanner(TestProcessTrigger):
    @property
    def rta(self):
        return os.path.join(cfg['lock_file_dir'], 'this', 'RTAComplete.txt')

    def setUp(self):
        super().setUp()
        self.scanner._rm(self.rta)
        for d in ['this', 'other']:
            if not os.path.isdir(os.path.join(self.from_dir, d)):
                os.mkdir(os.path.join(self.from_dir, d))

    def tearDown(self):
        self.scanner.reset('this')
        self.scanner.reset('that')
        self.scanner.switch_status('other', 'transferring')
        self.scanner.switch_status('another', 'transferring')
        self.scanner.switch_status('more', 'active')
        super().tearDown()

    def test_scan_datasets(self):
        datasets = self.scanner.scan_datasets()
        print(datasets)

        for observed, expected in (
            (datasets['new'], ['this']),
            (datasets['new, rta complete'], ['test_dataset', 'that']),
            (datasets['transferring'], ['other']),
            (datasets['transferring, rta complete'], ['another']),
            (datasets['active'], ['more']),
        ):
            assert observed == expected

    def test_triggerignore(self):
        with open(self.triggerignore, 'r') as f:
            assert f.readlines() == ['dir_t?_be_ign*d\n', 'test_dataset\n']

        expected = ['other\n', 'that\n', 'this\n']
        with open(self.triggerignore, 'w') as f:
            for d in expected:
                f.write(d)
        with open(self.triggerignore, 'r') as f:
            assert f.readlines() == expected
            for d in expected:
                assert d not in self._flatten(self.scanner.scan_datasets())

    def test_switch_status(self):
        sw = self._switch_and_assert
        d = 'this'

        self.scanner._rm(self.rta)

        sw(d, 'transferring')
        sw(d, 'another_status')

        self.scanner._touch(self.rta)
        sw(d, 'active')
        self.scanner.reset(d)
        assert self.scanner.dataset_status(d) == 'new, rta complete'
        sw(d, 'transferring')
        sw(d, 'another_status')

        self.scanner._rm(self.rta)
        self.scanner.reset(d)
        assert self.scanner.dataset_status(d) == 'new'

    @staticmethod
    def _flatten(d):
        l = list()
        for k, v in d.items():
            if type(v) is list:
                l.extend(v)
            else:
                l.append(v)
        return l

    def _switch_and_assert(self, d, s):
        self.scanner.switch_status(d, s)
        assert self.scanner.dataset_status(d) == s or self.scanner.dataset_status(d) == s + ', rta complete'
