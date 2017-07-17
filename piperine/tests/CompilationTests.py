import os
import unittest
import sys
import os
from tempfile import mkstemp
import pkg_resources

import filecmp

# From a stackoverflow, 16571150
from io import StringIO

from .. import designer
from .. import DSDClasses as trans_mod

# Grab package data
correct_sys_file = pkg_resources.resource_filename('piperine', 'tests/test_data/correct.sys')

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
            continue
            os.remove(f)

    def runTest(self):
        self.test_compilation()

    def test_crn_file(self):
        import logging
        rxns, spcs = designer.read_crn(self.crn_file)
        for i in range(len(rxns)):
            for key in list(rxns[i].keys()):
                true_rxn = self.true_crn[i]
                test_rxn = rxns[i]
                self.assertEqual(test_rxn[key], test_rxn[key])
    
    def test_sys_file(self):
        import logging
        rxns, spcs = designer.read_crn(self.crn_file)
        gates, strands = trans_mod.process_rxns(rxns, spcs, self.design_params)
        designer.write_sys_file(self.basename, gates, self.sys_file, trans_mod)
        with open(self.sys_file) as f:
            sys_lines = f.readlines()
        with open(correct_sys_file) as f:
            cor_lines = f.readlines()
        for i in range(1,len(sys_lines)):
            self.assertEqual(sys_lines[i], cor_lines[i])

def suite():
    tests = ['test_crn_file', 'test_sys_file']
    return unittest.TestSuite(list(map(TestMakePepperCompilerInputs, tests)))
