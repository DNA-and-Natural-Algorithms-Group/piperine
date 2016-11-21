from __future__ import division
import os
import unittest
import sys
import os
from tempfile import mkstemp

#from .. import sloth, tdm
# From a stackoverflow, 16571150
from cStringIO import StringIO

class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stdout = self._stdout

class TestMakePepperCompilerInputs(unittest.TestCase):
    crn = 'A -> B + B\nB + B -> A\nD -> A\n'
    
    th_params = {"thold_l":7, "thold_e":7.7, "e_dev":1, \
                 "m_spurious":0.5, "e_module":'energyfuncs_james'}
    
    design_params = (7, 15, 2)
    
    # Set true values
    true_crn = ([{'products': ['B'],
                'reactants': ['A'],
                'rate':1,
                'stoich_p':2,
                'stoich_r':1,
                }],
                [{'products': ['A'],
                 'reactants': ['B'],
                 'rate':1,
                 'stoich_p':1,
                 'stoich_r':2,
                  }],
                [{'products': ['A'],
                 'reactants': ['D'],
                 'rate':1,
                 'stoich_p':1,
                 'stoich_r':1,
                    }])

    correct_sys_file = os.path.join(os.path.dirname(__file__), 'test_data', 'correct.sys')
    fid = open(correct_sys_file)
    correct_sys = fid.readlines()
    fid.close()
    
    def setUp(self):
        fid, self.crn_file = mkstemp(suffix='.crn')
        os.close(fid)
        fid, self.sys_file = mkstemp(suffix='.sys')
        os.close(fid)
        fid, self.fixed_file = mkstemp(suffix='.fixed')
        os.close(fid)
        fid, self.pil_file = mkstemp(suffix='.pil')
        os.close(fid)
        fid, self.save_file = mkstemp(suffix='.save')
        os.close(fid)
        fid, self.mfe_file = mkstemp(suffix='.mfe')
        os.close(fid)
        fid, self.seq_file = mkstemp(suffix='.seqs')
        os.close(fid)
        self.filelist = [self.crn_file,
                         self.sys_file,
                         self.fixed_file,
                         self.pil_file,
                         self.save_file,
                         self.mfe_file,
                         self.seq_file]

        self.basename = self.sys_file[0:-4]
        # Generate crn file
        f = open(self.crn_file, 'w')
        f.write(self.crn)
        f.close()

    def tearDown(self):
        for f in self.filelist:
            os.remove(f)

    def runTest(self):
        self.test_compilation()

    def test_crn_file(self):
        import logging
        from .. import sloth
        rxns, spcs = sloth.read_crn(self.crn_file)
        for i in range(len(rxns)):
            for key in rxns[i].keys():
                true_rxn = self.true_crn[i]
                test_rxn = rxns[i]
                self.assertEqual(test_rxn[key], test_rxn[key])
    
    def test_process_rxns(self):
        from .. import sloth
        from .. import DSDClasses as trans_mod
        design_params = trans_mod.default_params
        abstract_rxns, spcs = sloth.read_crn(self.crn_file)
        rxn_list, signals = trans_mod.process_rxns(abstract_rxns, spcs, design_params)

    def test_sys_file(self):
        import logging
        from .. import sloth
        from .. import DSDClasses as trans_mod
        rxns, spcs = sloth.read_crn(self.crn_file)
        sloth.write_sys_file('test', rxns, self.sys_file, trans_mod)
        fid = open(self.sys_file)
        lines_list = fid.readline()[:-1]
        fid.close()
        for i in range(len(lines_list)):
            self.assertEqual(lines_list[i], self.correct_sys[i])

class Test_readcrn(unittest.TestCase):
    crn = 'A -> B + B\nB + B -> A\nD -> A\n'

    def setUp(self):
        fid, self.crn_file = mkstemp(suffix='.crn')
        os.close(fid)
        self.filelist = [self.crn_file]
        self.basename = self.crn_file[0:-4]
        # Generate crn file
        f = open(self.crn_file, 'w')
        f.write(self.crn)
        f.close()

        # Set true values
        true_crn ([{'products': ['B'],
                    'reactants': ['A'],
                    'rate':1,
                    'stoich_p':2,
                    'stoich_r':1,
                    }],
                  [{'products': ['A'],
                     'reactants': ['B'],
                     'rate':1,
                     'stoich_p':1,
                     'stoich_r':2,
                      }],
                  [{'products': ['A'],
                     'reactants': ['D'],
                     'rate':1,
                     'stoich_p':1,
                     'stoich_r':1,
                        }])

    def tearDown(self):
        for f in self.filelist:
            os.remove(f)

    def runTest(self):
        self.test_compilation()

    def test_compilation(self):
        import logging
        from .. import sloth as slothcomp
        logging.captureWarnings(False)
        with Capturing() as output:
            crn = slothcomp.read_crn(self.crn_file)

def suite():
    tests = ['test_crn_file', 'test_sys_file']
    return unittest.TestSuite(map(TestMakePepperCompilerInputs, tests))
