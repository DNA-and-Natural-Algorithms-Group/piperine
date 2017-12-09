import os
import subprocess
import unittest
import sys
from tempfile import mkstemp
from .test_data import fixed_file
from time import time

from .. import designer, tdm
from ..Srinivas2017 import energetics

# From a stackoverflow, 16571150
from io import StringIO

class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stdout = self._stdout

class Test_gen_th(unittest.TestCase):

    def setUp(self):
       self.ef = energetics.energyfuncs(targetdG=7.7)

    def tearDown(self):
        pass

    def runTest(self):
        pass

    def test_get_toeholds_timeout(self):
        starttime = time()
        timeout = 1e-4
        while timeflag:
            out = gen_th.new_toeholds(thold_l=3, timeout=1)
            if out == -1:
                timeflag = False
            self.assertGreater(timeout, starttime-time(), "Timeout does not work")

    def test_get_toeholds(self):
        out = gen_th.new_toeholds(self.ef, thold_l=5)




def suite():
    tests = ['test_tdm']
    return unittest.TestSuite(list(map(TestTDM, tests)))
