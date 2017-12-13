import os
import unittest
import sys
import os

from .. import designer, tdm

class TestSPW(unittest.TestCase):
    cwd = os.getcwd()
    basename = 'crn_for_spw_test'
    fixedfile = '{0}/constrained.fixed'.format(cwd)
    crn_file = 'A -> B + A\n'
    spw_true = [684.58333333333337, 6215828.75, 4.166666666666667,\
                10164.166666666666, 14.666666666666666, 476.3333333333333]
    def setUp(self):
        # Generate crn file
        crn_filename = '{0}/{1}.crn'.format(self.cwd, self.basename)
        f = open(crn_filename, 'w')
        f.write(self.crn_file)
        f.close()

        designer.generate_scheme(self.basename, fixedfile=self.fixedfile)
        designer.generate_seqs(self.basename)

        self.seq_dict = tdm.Read_Finished('{0}.seq'.format(self.basename))
        self.domains_list = list(self.seq_dict.keys())

    def tearDown(self):
        filelist = os.listdir(os.getcwd())
        for f in filelist:
            if self.basename in f:
                os.remove(f)

    def runTest(self):
        test_spw()

    def test_spw(self):
        scores = tdm.Spurious_Weighted_Score(self.basename, self.domains_list,\
                                             self.seq_dict)
        self.assertEqual(self.spw_true, scores, 'SPW Comparison')

def suite():
    tests = ['test_spw']
    return unittest.TestSuite(list(map(TestSPW, tests)))
