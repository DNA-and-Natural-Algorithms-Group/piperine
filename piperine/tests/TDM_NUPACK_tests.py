from pkg_resources import resource_filename
import os
import unittest
import sys
from tempfile import mkstemp, mkdtemp
from .test_data import fixed_file

import numpy as np

from .. import designer, tdm, DSDClasses


class TestTDMNUPACK(unittest.TestCase):
    basename = resource_filename('piperine', 'tests/test_data/sequences9mut.mfe')[:-4]
    
    def setUp(self):
        # get information from sequences9mut
        gates, strands = designer.process_crn(self.basename)
        self.h_inputs = tdm.get_heuristics_inputs(gates, strands)
        self.sequence_dicts = tdm.get_seq_dicts(self.basename, self.h_inputs)
        self.th_strs = [s.get_ths() for s in strands]
        self.toeholds = [self.sequence_dicts[0][i] for ths in self.th_strs for i in ths]
    
    def tearDown(self):
        pass 
    
    def test_NUPACKIntScore(self, tmpdir=None):
        seq_names  = ['hairpin1', 'hairpin2']
        seq_values = ["TTTTTTTTTTCCCAAAAAAAAAA", 
                      "AAAAAAAAAAGGGTTTTTTTTTT"]
        
        seq_dict = dict(zip(seq_names, seq_values))
        
        true_score = 100 * sum([7.289922e-13, 1.146297e-14, 1.156180e-12]) / 2e-6
        out = tdm.NUPACKIntScore(
                           seq_names[0], 
                           seq_names[1], 
                           seq_dict,
                           clean=False,
                           quiet=True,
                           tmpdir=tmpdir
                          )
        self.assertTrue(np.isclose(out, true_score, rtol=1e-6),
                         msg="true:{}\ntest:{}\n".format(true_score, out))
    
    def test_NUPACK_Cmpx_Conc(self, tmpdir=None):
        seq_names  = ['base', 'top1', 'top2']
        seq_values = ["TTTTTTTTTTCCCAAAAAAAAAA", 
                      "AAAAAAAAAA","TTTTTTTTTT"]
        
        seq_dict = dict(zip(seq_names, seq_values))
        
        true_score = 7.750862e-12
        out = tdm.NUPACK_Cmpx_Conc(seq_values, 
                                   params=[3, 25, 'dna', True, 'ted_calc'],
                                   clean=False,
                                   tmpdir=tmpdir)
        self.assertEqual(out, true_score, 
                         msg="true:{}\ntest:{}\n".format(true_score, out))
    

    def test_NUPACK_Cmpx_Conc_c4(self, tmpdir=None):
        seq_dict = {"r0-a_struct" : "TTCAGCCACAGAGTGACGCC",
                    "r0-b_struct" : "GGCGTCACTCGACCACGGCA",
                    "r0-c_struct" : "TGCCGTGGTCCGCAAAGGCG",
                    "r0-d_struct" : "CGCCTTTGCGATAGCGGTTA"}
        cmpx_ideal = "..........((((((((((+))))))))))((((((((((+))))))))))((((((((((+)))))))))).........."
        true_score = 7.547559e-07
        out = tdm.NUPACK_Cmpx_Conc(list(seq_dict.values()), 
                                   params=[4, 25, 'dna', True, 'ted_calc'],
                                   clean=False,
                                   tmpdir=tmpdir)
        self.assertEqual(out, true_score, 
                         msg="true:{}\ntest:{}\n".format(true_score, out))
    
    def test_NUPACK_Cmpx_Defect(self, tmpdir=None):
        # using gate r0-Gate from sequence set 9mut
        seqs = self.sequence_dicts[0]['r0-Gate'].split("+")
        struct = self.sequence_dicts[1]['r0-Gate']
        
        true_score = 2.329e+00
        out = tdm.NUPACK_Cmpx_Defect(seqs, 
                                     struct,
                                     params=[3, 25, 'dna', True, 'ted_calc'],
                                     clean=False,
                                     tmpdir=tmpdir)
        self.assertEqual(out, true_score, 
                         msg="true:{}\ntest:{}\n".format(true_score, out))
    
    def test_NUPACK_Cmpx_Defect_c4(self, tmpdir=None):
        seq_dict = {"r0-a_struct" : "TTCAGCCACAGAGTGACGCC",
                    "r0-b_struct" : "GGCGTCACTCGACCACGGCA",
                    "r0-c_struct" : "TGCCGTGGTCCGCAAAGGCG",
                    "r0-d_struct" : "CGCCTTTGCGATAGCGGTTA"}
        cmpx_ideal = "..........((((((((((+))))))))))((((((((((+))))))))))((((((((((+)))))))))).........."
        true_score = 4.265e+00
        out = tdm.NUPACK_Cmpx_Defect(seq_dict.values(), 
                                     cmpx_ideal,
                                     params=[4, 25, 'dna', True, 'ted_calc'],
                                     clean=False,
                                     tmpdir=tmpdir)
        self.assertEqual(out, true_score, 
                         msg="true:{}\ntest:{}\n".format(true_score, out))
    
    def test_Spurious_Weighted_Score(self, tmpdir=None):
        beta = 5
        w_lin = np.concatenate([np.zeros((beta, )), np.arange(12-beta)+1, 7*np.ones((2,))])
        n_strands = 30
        true_mis_intra_score = (1) 
        mis_inter_hits = np.array([6413,2457,793,227,47,5,5,5,5,5])
        true_mis_inter_score = sum(w_lin[4:14] * mis_inter_hits)
        spur_intra_hits = np.array([771, 211, 46, 9, 6, 6, 0, 0, 0, 0])
        true_intra_score = sum(w_lin[2:12] * spur_intra_hits)
        spur_inter_hits = np.array([25337, 7154, 1497, 142, 52, 36, 12, 11, 8, 8])
        true_inter_score = sum(w_lin[2:12] * spur_inter_hits)
        true_verboten = 685.00000 
        true_spc_weighted = 9310.66837 
        out = tdm.Spurious_Weighted_Score(
                self.basename,
                self.sequence_dicts[2],
                self.sequence_dicts[0],
                clean=False,
                tmpdir=tmpdir,
                bored=1)
        true_scores = [true_intra_score, true_inter_score, true_mis_intra_score, true_mis_inter_score,
                       true_verboten, true_spc_weighted]
        score_names =["spc_intra_score", "spc_inter_score",
                      "mis_intra_score", "mis_inter_score",
                      "verboten_score", "wsi_score"]
        for i in range(6):
            self.assertTrue(np.isclose(out[i], true_scores[i]), msg="Score {}".format(score_names[i]))
    
    def test_NUPACK_Eval_bad_nucleotide(self, tmpdir=None):
        import re
        mfe_seqs = self.sequence_dicts[0]
        ideal_structs = self.sequence_dicts[1]
        
        cmpx_names = self.h_inputs[1][::10]
        bps = np.array(
            [len(re.findall('[()]', ideal_structs[nam])) for nam in cmpx_names]
            )
        cmpx_concs = np.array([
                               1.000000e-06,
                               9.999991e-07,
                               1.000000e-06
                               ])
        
        cmpx_defex = np.array([
                               2.329e+00,
                               5.758e+00,
                               2.382e+00
                              ])
        seq_lens = np.array([
                            124, 
                            141, 
                            97])
        score_names = ['BN Max', 'BN Component', 'BN Mean']
        ted_vec = np.array(
                [cmpx_defex[i] * min(cmpx_concs[i], 1e-6) + \
                seq_lens[i]*max(1e-6 - cmpx_concs[i],0) \
                for i in range(len(cmpx_names))]
            )
        denom = (1e-6 * bps + ted_vec)
        bad_nuc = 100*ted_vec/denom
        argm = np.argmax(bad_nuc)
        trues = [bad_nuc.max(), cmpx_names[argm], bad_nuc.mean()]
        out = tdm.NUPACK_Eval_bad_nucleotide(mfe_seqs, 
                                          ideal_structs, 
                                          cmpx_names, 
                                          tmpdir=tmpdir,
                                          clean=False)
        for i in [0,2]:
            self.assertTrue(
                np.isclose(out[i],trues[i], rtol=1e-3),
                msg="Test:{} True:{} {}".format(out[i], trues[i], score_names[i])
            )
        self.assertEqual(trues[1], out[1], msg="True:{} Test:{}".format(trues[1], out[1]))
    
    def test_NUPACK_Eval_bad_nucleotide_c4(self, tmpdir=None):
        seq_dict = {"r0-a_struct" : "TTCAGCCACAGAGTGACGCC",
                    "r0-b_struct" : "GGCGTCACTCGACCACGGCA",
                    "r0-c_struct" : "TGCCGTGGTCCGCAAAGGCG",
                    "r0-d_struct" : "CGCCTTTGCGATAGCGGTTA",
                    "test_cmpx" : "TTCAGCCACAGAGTGACGCC+GGCGTCACTCGACCACGGCA+TGCCGTGGTCCGCAAAGGCG+CGCCTTTGCGATAGCGGTTA"}
        cmpx_ideal = "..........((((((((((+))))))))))((((((((((+))))))))))((((((((((+)))))))))).........."
        cmpx_name = "test_cmpx"
        ideal_structs = dict([(cmpx_name, cmpx_ideal)])
        import re
        cmpx_names = [cmpx_name]
        bps = np.array(
            [len(re.findall('[()]', ideal_structs[nam])) for nam in cmpx_names]
            )
        cmpx_concs = np.array([7.547559e-07])
        
        cmpx_defex = np.array([
                               4.265e+00,
                              ])
        seq_lens = np.array([80])
        ted_vec = np.array(
                [cmpx_defex[i] * min(cmpx_concs[i], 1e-6) + \
                seq_lens[i]*max(1e-6 - cmpx_concs[i],0) \
                for i in range(len(cmpx_names))]
            )
        denom = (1e-6 * bps + ted_vec)
        bad_nuc = 100*ted_vec/denom
        argm = np.argmax(bad_nuc)
        trues = [bad_nuc.max(), cmpx_names[argm], bad_nuc.mean()]
        out = tdm.NUPACK_Eval_bad_nucleotide(seq_dict, 
                                          ideal_structs, 
                                          cmpx_names, 
                                          tmpdir=tmpdir,
                                          clean=False)
        for i in [0,2]:
            self.assertTrue(
                np.isclose(out[i],trues[i], rtol=1e-3),
                msg="Test:{} True:{}".format(out[i], trues[i])
            )
        self.assertEqual(trues[1], out[1], msg="True:{} Test:{}".format(trues[1], out[1]))
    
    def test_NUPACKSSScore(self, tmpdir=None):
        strand = self.h_inputs[0][5]
        t_regi = self.h_inputs[3][strand]
        true_sstu = 0.9028 #0.9059
        true_sstu_sum = 9.5695 # 9.5847
        true_ssu = 0.8524 # 0.8579
        true_ssu_sum = 42.0609 # 42.1328
        trues = [true_ssu, true_ssu_sum,
                 true_sstu, true_sstu_sum, 44]
        score_names = ["min_Unpaired", "sum_Unpaired", "min_Unpaired_toe", \
                       "sum_Unpaired_toe", "seq_len"]
        out = tdm.NUPACKSSScore(strand, 
                                self.sequence_dicts[0], 
                                toe_region=t_regi,
                                clean=False,
                                tmpdir=tmpdir)
        for i in range(4):
            self.assertTrue(np.isclose(out[i], trues[i]), 
                            msg="Score {} Test:{} True:{}".format(score_names[i], out[i], trues[i]))
    
    def test_BM_Eval(self):
        sd = ['qwertyuiop', 'asdfghjkl', 'zxcvbnm']
        bm = range(3)
        out = tdm.BM_Eval(sd, bm, [])
        trues = [0, 0]
        names = ['No match BM score', 'No match largest match']
        for i in [0,1]:
            self.assertEqual(
                out[i], 
                trues[i], 
                "{} incorrect, test:{} true:{}".format(names[i], out[i], trues[i]))
        
        sd = ['qwertyuiop', 'asdfghjklyuiop', 'zxcvbnm']
        bm = range(3)
        out = tdm.BM_Eval(sd, bm, [])
        trues = [1, 5]
        names = ['5 match BM score', '5 match largest match']
        for i in [0,1]:
            self.assertEqual(
                out[i], 
                trues[i], 
                "{} incorrect, test:{} true:{}".format(names[i], out[i], trues[i]))
        
        sd = ['qwertyuiop', 'asdfghjkrtyuiop', 'zrtyuiopxcvbnm']
        bm = range(3)
        out = tdm.BM_Eval(sd, bm, [])
        trues = [12, 7]
        names = ['Thrice 7 match BM score', 'Thrice 7 match largest match']
        for i in [0,1]:
            self.assertEqual(
                out[i], 
                trues[i], 
                "{} incorrect, test:{} true:{}".format(names[i], out[i], trues[i]))
        
        sd = ['rtyuiopqwertyuiop', 'asdfghjkrtyuiop', 'zrtyuiopxcvbnm']
        bm = range(3)
        out = tdm.BM_Eval(sd, bm, [])
        trues = [12, 7]
        names = ['Thrice 7 match BM score, one internal repeat', 
                 'Thrice 7 match largest match, one internal repeat']
        for i in [0,1]:
            self.assertEqual(
                out[i], 
                trues[i], 
                "{} incorrect, test:{} true:{}".format(names[i], out[i], trues[i]))
        

def suite():
    tests = ['test_NUPACKIntScore', 'test_NUPACK_Cmpx_Conc', 'test_NUPACK_Cmpx_Defect',
             'test_NUPACK_Cmpx_Conc_c4', 'test_NUPACK_Cmpx_Defect_c4',
             'test_Spurious_Weighted_Score', 'test_NUPACK_Eval_bad_nucleotide', 
             'test_NUPACK_Eval_bad_nucleotide_c4', 
             'test_NUPACKSSScore', 'test_BM_Eval']
    return unittest.TestSuite(list(map(TestTDMNUPACK, tests)))
