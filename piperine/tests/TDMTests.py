import pkg_resources
import os
import unittest
import sys
from tempfile import mkstemp
import os
from .test_data import fixed_file

from .. import designer, tdm, DSDClasses

def to_sequences(name_list, seq_dict):
    return [seq_dict[name] for name in name_list]

class TestTDM(unittest.TestCase):
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
        self.crn_file = pkg_resources.resource_filename('piperine', 'tests/test_data/test_tdm.crn')
        self.sys_file = pkg_resources.resource_filename('piperine', 'tests/test_data/test_tdm.sys')
        self.fixed_file = pkg_resources.resource_filename('piperine', 'tests/test_data/test_tdm.fixed')
        self.seqs_file = pkg_resources.resource_filename('piperine', 'tests/test_data/test_tdm.seq')
        self.pil_file = pkg_resources.resource_filename('piperine', 'tests/test_data/test_tdm.pil')
        self.mfe_file = pkg_resources.resource_filename('piperine', 'tests/test_data/test_tdm.mfe')
        
        self.basename = self.sys_file[0:-4]
        self.gates, self.strands = designer.process_crn(crn_file = self.crn_file,
                                                        trans_module = DSDClasses,
                                                        design_params = (7, 15, 2))
        self.h_inputs = tdm.get_heuristics_inputs(self.gates, self.strands)
        self.seq_dict, self.cmplx_dict, self.domains_list = \
            tdm.get_seq_dicts(self.basename, self.h_inputs)
        
        # True top strand list
        self.true_tsl =  ['r1-c',
                          'r2-c',
                          'r2-d',
                          'r1-d',
                          'r0-c',
                          'r0-d',
                          'r2-b',
                          'r0-Out',
                          'r0-Backward',
                          'r0-cat_helper',
                          'r0-helper',
                          'r1-Out',
                          'r1-Backward',
                          'r1-cat_helper',
                          'r1-helper',
                          'r2-Out',
                          'r2-Backward',
                          'r2-cat_helper',
                          'r2-helper']
        
        # Base strand list contains all toehold complements
        self.true_bsl =  ['r0-toe-fa*',
                          #'r0-toe-sa*', # A does not appear as 2nd input
                          'r0-toe-fb*', 
                          'r0-toe-sb*', 
                          'r0-toe-fc*', 
                          #'r0-toe-sc*', # C does not appear as 2nd input
                          'r0-toe-fd*', 
                          'r0-toe-sd*', 
                          'r2-toe-fb*',
                          'r2-toe-sb*']
        
        # Top strand dict lead to these subsequences
        sd = self.seq_dict
        self.true_tsd = [sd['r1-cm'][-1] + sd['r1-ccm'] + sd['r1-toe-fc'], #'r1-c',
                         sd['r2-cm'][-1] + sd['r2-ccm'] + sd['r2-toe-fc'], #'r2-c',
                         sd['r2-dm'][-1] + sd['r2-cdm'] + sd['r2-toe-fd'], #'r2-d',
                         sd['r1-dm'][-1] + sd['r1-cdm'] + sd['r1-toe-fd'], #'r1-d',
                         sd['r0-cm'][-1] + sd['r0-ccm'] + sd['r0-toe-fc'], #'r0-c',
                         sd['r0-dm'][-1] + sd['r0-cdm'] + sd['r0-toe-fd'], #'r0-d',
                         sd['r2-bm'][-1] + sd['r2-cbm'] + sd['r2-toe-fb'], #'r2-b',
                         sd['r0-ch'][-1] + sd['r0-cch'] + sd['r0-toe-sb'], #'r0-Out',
                         sd['r0-toe-fb-suffix'] + sd['r0-toe-sa'][:3], #'r0-Backward',
                         sd['r0-dh'][-1]+sd['r0-cdh']+sd['r0-toe-fc']+sd['r0-ch'][:3], #'r0-cat_helper',
                         sd['r0-dh'][-1] + sd['r0-cdh'] + sd['r0-toe-fc'], #'r0-helper',
                         sd['r1-ch'][-1] + sd['r1-cch'] + sd['r1-toe-sb'], #'r1-Out',
                         sd['r1-toe-fb-suffix'] + sd['r1-toe-sa'][:3], #'r1-Backward',
                         sd['r1-dh'][-1]+sd['r1-cdh']+sd['r1-toe-fc']+sd['r1-ch'][:3], #'r1-cat_helper',
                         sd['r1-dh'][-1] + sd['r1-cdh'] + sd['r1-toe-fc'], #'r1-helper',
                         sd['r2-ch'][-1] + sd['r2-cch'] + sd['r2-toe-sb'], #'r2-Out',
                         sd['r2-toe-fb-suffix'] + sd['r2-toe-sa'][:3], #'r2-Backward',
                         sd['r2-dh'][-1]+sd['r2-cdh']+sd['r2-toe-fc']+sd['r2-ch'][:3], #'r2-cat_helper',
                         sd['r2-dh'][-1] + sd['r2-cdh'] + sd['r2-toe-fc']] #'r2-helper',
        
        # NotToInteract should return this dictionary
        self.true_nti = {'r0-toe-fa*': ['r0-Out',  
                         'r0-Backward',            
                         'r0-cat_helper',          
                         'r0-helper',              
                         'r1-Out',                 
                         'r1-Backward',            
                         'r1-toe-fdr1-dhr1-cdh',    # cat parts, bc A is 1st product
                         'r1-chr1-cch',             # cat parts, bc A is 1st product
                         'r1-toe-fdr1-dhr1-cdh',    # help parts, bc A is 1st product
                         'r2-Out',                 
                         'r2-Backward',            
                         'r2-dhr2-cdh',             # A is both products
                         'r2-chr2-cch',             # A is both products
                         'r2-dhr2-cdh',             # A is both products
                         'r1-chr1-cch',             # A is both products
                         'r0-toe-sar0-amr0-cam',   #
                         'r2-chr2-cch',            #
                         'r0-toe-sar0-amr0-cam',   #
                         'r2-dhr2-cdh',            #
                         'r0-toe-sar0-amr0-cam',   #
                         'r1-d',                   
                         'r0-c',                   
                         'r0-d',                   
                         'r2-b'],                  
                        'r0-toe-fb*': ['r0-Out',
                         'r0-toe-sar0-amr0-cam',   #
                         'r0-cat_helper',
                         'r0-helper',
                         'r1-Out',
                         'r1-Backward',
                         'r1-dhr1-cdhr1-toe-fcr1-chr1-cch', #
                         'r1-dhr1-cdhr1-toe-fc',            #
                         'r2-Out',
                         'r2-Backward',
                         'r2-cat_helper',
                         'r2-helper',
                         'r1-c',
                         'r2-c',
                         'r2-d',
                         'r1-dhr1-cdh',            # 
                         'r0-toe-sbr0-bmr0-cbm',   # 
                         'r0-c',
                         'r0-d',
                         'r2-b'],
                        'r0-toe-fc*': ['r0-Out',
                         'r0-Backward',
                         'r0-toe-fdr0-dhr0-cdh',   #
                         'r0-chr0-cch',            #
                         'r0-toe-fdr0-dhr0-cdh',   #
                         'r1-Out',
                         'r1-Backward',
                         'r1-cat_helper',
                         'r1-helper',
                         'r2-Out',
                         'r2-Backward',
                         'r2-cat_helper',
                         'r2-helper',
                         'r1-c',
                         'r2-c',
                         'r2-d',
                         'r1-d',
                         'r0-chr0-cch',            #
                         'r0-toe-scr0-cmr0-ccm',   #
                         'r0-d',
                         'r2-b'],
                        'r0-toe-fd*': ['r0-Out',
                         'r0-Backward',
                         'r0-dhr0-cdhr0-toe-fcr0-chr0-cch', #
                         'r0-dhr0-cdhr0-toe-fc',            #
                         'r1-Out',
                         'r1-toe-sar1-amr1-cam',  #
                         'r1-cat_helper',
                         'r1-helper',
                         'r2-Out',
                         'r2-Backward',
                         'r2-cat_helper',
                         'r2-helper',
                         'r1-c',
                         'r2-c',
                         'r2-d',
                         'r1-d',
                         'r0-c',
                         'r0-dhr0-cdh',          #
                         'r0-toe-sdr0-dmr0-cdm', #
                         'r2-b'],
                        'r0-toe-sb*': ['r0-chr0-cch', #
                         'r0-bmr0-cbm',               #
                         'r0-Backward',
                         'r0-cat_helper',
                         'r0-helper',
                         'r1-Out',
                         'r1-Backward',
                         'r1-cat_helper',
                         'r1-helper',
                         'r2-Out',
                         'r2-Backward',
                         'r2-cat_helper',
                         'r2-helper',
                         'r1-c',
                         'r2-c',
                         'r2-d',
                         'r0-bmr0-cbmr0-toe-fbr1-dhr1-cdh', #
                         'r0-c',
                         'r0-d',
                         'r2-b'],
                        'r0-toe-sd*': ['r0-Out',
                         'r0-Backward',
                         'r0-cat_helper',
                         'r0-helper',
                         'r1-chr1-cch', #
                         'r1-bmr1-cbm', #
                         'r1-Backward',
                         'r1-cat_helper',
                         'r1-helper',
                         'r2-Out',
                         'r2-Backward',
                         'r2-cat_helper',
                         'r2-helper',
                         'r1-c',
                         'r2-c',
                         'r2-d',
                         'r1-d',
                         'r0-c',
                         'r0-dmr0-cdmr0-toe-fdr0-dhr0-cdh', #
                         'r2-b'],
                        'r2-toe-fb*': ['r0-Out',
                         'r0-Backward',
                         'r0-cat_helper',
                         'r0-helper',
                         'r1-Out',
                         'r1-Backward',
                         'r1-cat_helper',
                         'r1-helper',
                         'r2-Out',
                         'r2-toe-sar2-amr2-cam', #
                         'r2-cat_helper',
                         'r2-helper',
                         'r1-c',
                         'r2-c',
                         'r2-d',
                         'r1-d',
                         'r0-c',
                         'r0-d',
                         'r2-toe-sbr2-bmr2-cbm'], # 
                        'r2-toe-sb*': ['r0-Out',
                         'r0-Backward',
                         'r0-cat_helper',
                         'r0-helper',
                         'r1-Out',
                         'r1-Backward',
                         'r1-cat_helper',
                         'r1-helper',
                         'r2-chr2-cch',
                         'r2-bmr2-cbm',
                         'r2-Backward',
                         'r2-cat_helper',
                         'r2-helper',
                         'r1-c',
                         'r2-c',
                         'r2-d',
                         'r1-d',
                         'r0-c',
                         'r0-d',
                         'r2-bmr2-cbmr2-toe-fb']} # 
        
        self.bmlist = ['r0-amr0-cam',
                       'r1-chr1-cch',
                       'r2-chr2-cch',
                       'r2-dhr2-cdh',
                       'r0-bmr0-cbm',
                       'r1-dhr1-cdh',
                       'r0-cmr0-ccm',
                       'r0-chr0-cch',
                       'r0-dmr0-cdm',
                       'r0-dhr0-cdh',
                       'r2-bmr2-cbm']

        self.cmplx_names = ['r0-Gate',
                            'r0-Gate_int',
                            'r0-Gate_waste',
                            'r0-Trans',
                            'r0-Trans_int',
                            'r0-Trans_waste',
                            'r0-Trans_cat_waste',
                            'r1-Gate',
                            'r1-Gate_int',
                            'r1-Gate_waste',
                            'r1-Trans',
                            'r1-Trans_int',
                            'r1-Trans_waste',
                            'r1-Trans_cat_waste',
                            'r2-Gate',
                            'r2-Gate_int',
                            'r2-Gate_waste',
                            'r2-Trans',
                            'r2-Trans_int',
                            'r2-Trans_waste',
                            'r2-Trans_cat_waste']
                                  
    def tearDown(self):
        pass 
    
    def runTest(self):
        self.test_tdm()
    
    def test_TopStrandlist(self):
        # all strands intended to be in solution as solitons with no secondary structure
        set_true = set(self.true_tsl)
        set_test = set(self.h_inputs[0])
        self.assertTrue(set_true <= set_test, msg="True set <= test set")
        self.assertTrue(set_test <= set_true, msg="Test set <= true set")
    
    def test_BaseStrandlist(self):
        # Base strand list contains all toehold complements
        set_true = set(to_sequences(self.true_bsl, self.seq_dict))
        set_test = set(to_sequences(self.h_inputs[2], self.seq_dict))
        self.assertTrue(set_true <= set_test, msg="True set <= test set {}".format(set_true - set_test))
        self.assertTrue(set_test <= set_true, msg="Test set <= true set {}".format(set_test - set_true))
    
    def test_TopStranddict(self):
        # TopStranddict maps each top strand to the critical toehold binding area
        def tsd_to_sequence(tsd, sd):
            outlist = []
            for key in tsd.keys():
                top_strand = sd[key]
                critical_region = ''.join([top_strand[i-1] for i in tsd[key]])
                outlist.append(critical_region)
            return outlist
        
        set_true = set(self.true_tsd)
        set_test = set(tsd_to_sequence(self.h_inputs[3], self.seq_dict))
        self.assertTrue(set_true <= set_test, msg="True set <= test set {}".format(set_true - set_test))
        self.assertTrue(set_test <= set_true, msg="Test set <= true set {}".format(set_test - set_true))
    
    def test_NotToInteract(self):
        # NotToInteract 
        def nti_to_sequence(nti, sd):
            outdict = {}
            for key in nti.keys():
                outdict[sd[key]] = to_sequences(nti[key], sd)
            return outdict
        true_nti_seqs = nti_to_sequence(self.true_nti, self.seq_dict)
        test_nti_seqs = nti_to_sequence(self.h_inputs[5], self.seq_dict)
        l1, l2 = len(true_nti_seqs), len(test_nti_seqs)
        self.assertEqual(l1, l2, msg="Len true: {} Len test: {}".format(l1, l2))
        for key in true_nti_seqs.keys():
            set_true = set(true_nti_seqs[key])
            set_test = set(test_nti_seqs[key])
            self.assertTrue(set_true <= set_test, msg="True set <= test set {}".format(set_true - set_test))
            self.assertTrue(set_test <= set_true, msg="Test set <= true set {}".format(set_test - set_true))
    
    def test_BMlist(self):
        set_true = set(to_sequences(self.bmlist, self.seq_dict))
        set_test = set(to_sequences(self.h_inputs[4], self.seq_dict))
        self.assertTrue(set_true <= set_test, msg="True set <= test set {}".format(set_true - set_test))
        self.assertTrue(set_test <= set_true, msg="Test set <= true set {}".format(set_test - set_true))
    
    def test_complex_names(self):
        set_true = set(to_sequences(self.cmplx_names, self.seq_dict))
        set_test = set(to_sequences(self.h_inputs[1], self.seq_dict))
        self.assertTrue(set_true <= set_test, msg="True set <= test set {}".format(set_true - set_test))
        self.assertTrue(set_test <= set_true, msg="Test set <= true set {}".format(set_test - set_true))
    
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
        iterate_list = list(zip(self.scores, self.true_val, self.fmts, names))
        for score, tru, fmt, name in iterate_list:
            score_str = fmt.format(score)
            truth_str = fmt.format(tru)
            self.assertEqual(score_str, truth_str, 'Score {}'.format(name))

def suite():
    tests = ['test_TopStrandlist', 'test_BaseStrandlist', 'test_TopStranddict', 'test_NotToInteract',\
             'test_BMlist', 'test_complex_names']
    return unittest.TestSuite(list(map(TestTDM, tests)))
