__author__ = 'mwham'
from unittest import TestCase
import os.path


class TestAnalysisDriver(TestCase):
    assets_path = os.path.join(os.path.dirname(__file__), 'assets')
    fastq_path = os.path.join(assets_path, 'fastqs')
