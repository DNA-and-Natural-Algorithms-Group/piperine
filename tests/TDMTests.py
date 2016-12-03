from __future__ import division
import os
import unittest
import sys
from tempfile import mkstemp
import os
from test_data import fixed_file

from .context import designer, tdm, DSDClasses

class TestTDM(unittest.TestCase):
    crn = 'A -> B + A\n'
    
    true_val = [0.01,
                0.01,
                4.54,
                9.74,
                1.00,
                5.00,
                0.96,
                0.99,
                0.96,
                0.99,
                0.05,
                "p0-Trans_int",
                0.03,
                1168.12,
                15156.88,
                0.00,
                10323.75,
                5.56,
                17.33]
    
    fmts = ['{:.2f}', '{:.2f}', '{:.2f}', '{:.2f}', '{:.2f}', '{:.2f}',\
            '{:.2f}', '{:.2f}', '{:.2f}', '{:.2f}', '{:.2f}', '{}', '{:.2f}',\
            '{:.2f}', '{:.2f}', '{:.2f}', '{:.2f}', '{:.2f}', '{:.2f}',\
            '{:.2f}', '{:.2f}']
    
    score_val = []
    
    def setUp(self):
        fid, self.crn_file = mkstemp(suffix='.crn')
        os.close(fid)
        fid, self.sys_file = mkstemp(suffix='.sys')
        os.close(fid)
        self.fixed_file = fixed_file
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
        filelist = os.listdir(os.getcwd())
        for f in filelist:
            if self.basename in f:
                os.remove(f)
    
    def runTest(self):
        self.test_tdm()
    
    def test_tdm(self):
        # Score test sequences
        self.scores, names = designer.score_fixed(self.fixed_file, 
                 basename=self.basename, 
                 crn_file=self.crn_file, 
                 sys_file=self.sys_file, 
                 pil_file=self.pil_file, 
                 save_file=self.save_file, 
                 mfe_file=self.mfe_file, 
                 seq_file=self.seq_file,
                 design_params=(7, 15, 2),
                 trans_module=DSDClasses)
        iterate_list = zip(self.scores, self.true_val, self.fmts, names)
        for score, tru, fmt, name in iterate_list:
            score_str = fmt.format(score)
            truth_str = fmt.format(tru)
            self.assertEqual(score_str, truth_str, 'Score {}'.format(name))

def suite():
    tests = ['test_tdm']
    return unittest.TestSuite(map(TestTDM, tests))
