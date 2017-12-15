import unittest
import sys
from tempfile import mkstemp
import os

from ..designer import read_crn, get_parameters_from_crn_file
from .. import Srinivas2017 as default_translation_scheme
default_energyfuncs = default_translation_scheme.energetics.energyfuncs()

class TestCRNImport(unittest.TestCase):

    def setUp(self):
        fid, testfile = mkstemp(suffix='.crn')
        os.close(fid)
        self.testfile = testfile
        test_crn = 'A + B -> C + D (1.1e5)\nB + B -> A (1.1e5)\n2B -> 3.2C + D (1/5)\n2B -> 2C (1)\n'
        terms = ["toehold_energy",
                 "toehold_deviation",
                 "toehold_spurious",
                 "toehold_length",
                 "spurious_design_parameters",
                 "translation_scheme",
                 "n",
                 "bm"]
        values = [7.5, 0.7, 0.5, 5, "imax=-1 bored=10", "Chen2013", 6, 14]
        params = '\n'.join(["{} = {}".format(term, val) for term, val in zip(terms, values)])
        with open(self.testfile, 'w') as f:
            f.write(test_crn + params)
            f.close()

        output = read_crn(self.testfile)
        self.reactions_in = output[0]
        self.species_in = output[1]
        self.reactions_T = [{'reactants':['A', 'B'],
                              'products':['C', 'D'],
                              'stoich_p':[1, 1],
                              'stoich_r':[1, 1],
                              'rate':1.1e5},
                             {'reactants':['B', 'B'],
                              'products':['A'],
                              'stoich_p':[4/3],
                              'stoich_r':[1, 1],
                              'rate':1.1e5},
                             {'reactants':['B'],
                              'products':['C', 'D'],
                              'stoich_p':[3.2, 1],
                              'stoich_r':[2],
                              'rate':0.2},
                             {'reactants':['B'],
                              'products':['C'],
                              'stoich_p':[2],
                              'stoich_r':[2],
                              'rate':1}]

        self.species_T = ['A', 'B', 'C', 'D']
        self.param_dict_T = dict(zip(terms, values))

    def tearDown(self):
        os.remove(self.testfile)

    def test_proper_characters(self):
        readin_b_b = self.reactions_in[1]
        true_b_b = self.reactions_T[1]
        self.assertEqual(readin_b_b['reactants'], true_b_b['reactants'],
            'Incorrect A+A-> reactants interpretation')
        self.assertEqual(readin_b_b['stoich_r'], true_b_b['stoich_r'],
            'Incorrect A+A-> coefficients interpretation')
        readin_2b = self.reactions_in[2]
        true_2b = self.reactions_T[2]
        self.assertEqual(readin_b_b['reactants'], true_b_b['reactants'],
            'Incorrect 2A-> reactants interpretation')
        self.assertEqual(readin_b_b['stoich_r'], true_b_b['stoich_r'],
            'Incorrect 2A-> coefficients interpretation')

    def test_integer_coefficients(self):
        readin_b_b = self.reactions_in[1]
        true_b_b = self.reactions_T[1]
        self.assertEqual(readin_b_b['reactants'], true_b_b['reactants'],
            'Incorrect A+A-> reactants interpretation')
        self.assertEqual(readin_b_b['stoich_r'], true_b_b['stoich_r'],
            'Incorrect A+A-> coefficients interpretation')
        readin_2b = self.reactions_in[2]
        true_2b = self.reactions_T[2]
        self.assertEqual(readin_b_b['reactants'], true_b_b['reactants'],
            'Incorrect 2A-> reactants interpretation')
        self.assertEqual(readin_b_b['stoich_r'], true_b_b['stoich_r'],
            'Incorrect 2A-> coefficients interpretation')

    def test_noninteger_stoichiometry(self):
        readin_frac = self.reactions_in[1]
        true_frac = self.reactions_T[1]
        self.assertEqual(readin_frac['stoich_r'], true_frac['stoich_r'],
            'Incorrect scientific-notation coefficients interpretation')
        readin_deci = self.reactions_in[2]
        true_deci = self.reactions_T[2]
        self.assertEqual(readin_frac['stoich_r'], true_frac['stoich_r'],
            'Incorrect fractional coefficients interpretation')
        readin_deci = self.reactions_in[3]
        true_deci = self.reactions_T[3]
        self.assertEqual(readin_frac['stoich_r'], true_frac['stoich_r'],
            'Incorrect integer coefficients interpretation')

    def test_parameter_reading(self):
        p_dict = get_parameters_from_crn_file(self.testfile)
        for term in self.param_dict_T:
            self.assertEqual(self.param_dict_T[term], p_dict[term], "Comparing {}".format(term))

    def test_reaction_rate(self):
        readin_scino = self.reactions_in[1]
        true_scino = self.reactions_T[1]
        readin_frac = self.reactions_in[2]
        true_frac = self.reactions_T[2]
        self.assertEqual(readin_scino['reactants'], true_scino['reactants'],
            'Incorrect decimal rate constant interpretation')
        self.assertEqual(readin_frac['rate'], true_frac['rate'],
            'Incorrect fractional rate constant interpretation')

    def runTest(self):
        self.test_integer_coefficients()
        self.test_noninteger_coefficients()
        self.test_reaction_rate()

def suite():
    tests = ['test_integer_coefficients', 'test_noninteger_stoichiometry', 'test_reaction_rate', 'test_parameter_reading']
    return unittest.TestSuite(list(map(TestCRNImport, tests)))
