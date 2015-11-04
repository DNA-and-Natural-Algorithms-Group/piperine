from __future__ import division
import os
import unittest
import sys
import os

from .. import sloth, tdm

class TestTDM(unittest.TestCase):
    cwd = os.getcwd()
    basename = 'crn_for_tdm_test'
    fixedfile = '{0}/constrained.fixed'.format(cwd)
    crn_file = 'A -> B + A\n'
    
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
        # Generate crn file
        crn_filename = '{0}/{1}.crn'.format(self.cwd, self.basename)
        f = open(crn_filename, 'w')
        f.write(self.crn_file)
        f.close()
        
        #outlist = sloth.generate_scheme(self.basename, fixedfile=self.fixedfile)
        #sloth.generate_seqs(self.basename)
        #
        #(gates, strands, toeholds, th_scores) = outlist
        #self.outlist = outlist
        #
        #self.namespace = tdm.get_seq_lists(self.basename, gates, strands)
    
    def tearDown(self):
        filelist = os.listdir(os.getcwd())
        for f in filelist:
            if self.basename in f:
                os.remove(f)
    
    def runTest(self):
        self.test_tdm()
    
    def test_tdm(self):
        # Score test sequences
        self.scores, names = sloth.score_fixed(self.fixedfile, 
                                                        basename=self.basename, 
                                                        mod_str="DSDClasses")
        #(gates, strands, toeholds, th_scores) = self.outlist
        #self.scores, names = tdm.EvalCurrent(self.basename, gates, strands,\
        #                                     th_scores)
        iterate_list = zip(self.scores, self.true_val, self.fmts, names)
        for score, tru, fmt, name in iterate_list:
            score_str = fmt.format(score)
            truth_str = fmt.format(tru)
            self.assertEqual(score_str, truth_str, 'Score {}'.format(name))

def suite():
    tests = ['test_tdm']
    return unittest.TestSuite(map(TestTDM, tests))
