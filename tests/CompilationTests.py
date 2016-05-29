from __future__ import division
import os
import unittest
import sys
import os
from tempfile import mkstemp

#from .. import sloth, tdm

class TestCompilation(unittest.TestCase):
    crn = 'A -> B + B\nB + B -> A\nD -> A\n'
    
    th_params = {"thold_l":7, "thold_e":7.7, "e_dev":1, \
                 "m_spurious":0.5, "e_module":'energyfuncs_james'}
    
    design_params = (7, 15, 2)
    
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
        
        #outlist = sloth.generate_scheme(self.basename, fixedfile=self.fixedfile)
        #sloth.generate_seqs(self.basename)
        #
        #(gates, strands, toeholds, th_scores) = outlist
        #self.outlist = outlist
        #
        #self.namespace = tdm.get_seq_lists(self.basename, gates, strands)
    
    def tearDown(self):
        for f in self.filelist:
            os.remove(f)
    
    def runTest(self):
        self.test_compilation()
    
    def test_compilation(self):
        from sloth import sloth as slothcomp
        from sloth import DSDClasses as mod
        gates, strands = slothcomp.generate_scheme(self.basename, 
                                           self.design_params, 
                                           crn_file=self.crn_file, 
                                           systemfile=self.sys_file)
        
        toeholds, th_scores = slothcomp.generate_seqs(self.basename, 
                                            gates, 
                                            strands, 
                                            self.design_params, 
                                            mod.n_th, 
                                            self.th_params, 
                                            systemfile=self.sys_file,
                                            pilfile=self.pil_file,
                                            mfefile=self.mfe_file,
                                            seqfile=self.seq_file,
                                            savefile=self.save_file,
                                            fixedfile=self.fixed_file,
                                            extra_pars='bored=10')
        

def suite():
    tests = ['test_compilation']
    return unittest.TestSuite(map(TestCompilation, tests))
