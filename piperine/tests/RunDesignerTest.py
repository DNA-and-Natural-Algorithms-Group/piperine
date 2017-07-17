import os
import subprocess
import unittest
import sys
from tempfile import mkstemp
from .test_data import fixed_file
from time import time

from .. import designer, tdm, energyfuncs_james, gen_th, DSDClasses

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
        self.ef = energyfuncs_james.energyfuncs(targetdG=7.7)
        self.trans = DSDClasses
        self.ef_str = 'energyfuncs_james'
        self.trans_str = 'DSDClasses'
        
    def tearDown(self):
        for f in self.filenames:
            if os.path.isfile(f):
                os.remove(f)
    
    def runTest(self):
        self.test_run_designer_accepts_string_modules()
    
    def test_run_designer_noargs(self):
        with Capturing() as output:
            out = designer.run_designer(quick=True)
        
    def test_run_designer_accepts_string_modules(self):
        try:
            with Capturing() as output:
                out = designer.run_designer(basename=self.basename,
                                         e_module=self.ef_str,
                                         trans_module=self.trans_str,
                                         quick=True)
        except AttributeError as e:
            argsplit = e.args[0].split("'")
            mod = argsplit[1]
            att = argsplit[3]
            if mod == 'str':
                print(e)
                self.fail("Failed to import module from string")
            if self.trans_str in mod or self.ef_str in mod:
                print(e)
                self.fail("Piperine found module {} but it had no attribute {}".format(mod,att))
            raise(e)
        except ImportError as e:
            argsplit = e.args[0].split(" ")
            mod = argsplit[-1]
            self.fail("Found no module named {}".format(mod))
    
    def test_run_designer_alerts_unfound_modules(self):
        fakemod1 = "notAmodule"
        fakemod2 = "notAmodule"
        with self.assertRaises(ModuleNotFoundError):
            out = designer.run_designer(basename=self.basename,
                                     e_module=fakemod1,
                                     trans_module=fakemod2,
                                     quick=True)

def suite():
    tests = ['test_run_designer_accepts_string_modules', 'test_run_designer_noargs',
             'test_run_designer_alerts_unfound_modules']
    return unittest.TestSuite(list(map(Test_run_designer, tests)))
