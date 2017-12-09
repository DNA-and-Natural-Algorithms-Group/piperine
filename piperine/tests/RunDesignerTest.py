import os
import subprocess
import unittest
import sys
from tempfile import mkstemp
from .test_data import fixed_file
from time import time

from .. import designer, tdm
from ..Srinivas2017 import translation, energetics

# From a stackoverflow, 16571150
if sys.version_info >= (3,0):
    from io import StringIO
else:
    from StringIO import StringIO

class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stdout = self._stdout

class Test_run_designer(unittest.TestCase):

    def setUp(self):
        # CRN used in this test
        self.crn_rxn = 'A + B -> A + D\n'
        # Establish random temporary filenames
        fid, self.basename = mkstemp(suffix='')
        os.close(fid)
        endings = ['.crn', '.fixed', '.pil', '.mfe']
        self.filenames = [ self.basename + suf for suf in endings ]
        self.crn, self.fixed, self.pil, self.mfe = self.filenames
        # Write CRN to basename.crn
        with open(self.crn, 'w') as f:
            f.write(self.crn_rxn)
        # Modules and module strings for import tests
        self.ef = energetics.energyfuncs(targetdG=7.7)
        self.trans = translation
        self.ef_str = 'energyfuncs_james'
        self.trans_str = 'DSDClasses'

    def tearDown(self):
        for f in self.filenames:
            if os.path.isfile(f):
                os.remove(f)

    def test_run_designer_smallcrn(self):
        with Capturing() as output:
            out = designer.run_designer(designer.small_crn[:-4], quick=True)

def suite():
    tests = [#'test_run_designer_accepts_string_modules',
             'test_run_designer_smallcrn']
    return unittest.TestSuite(list(map(Test_run_designer, tests)))
