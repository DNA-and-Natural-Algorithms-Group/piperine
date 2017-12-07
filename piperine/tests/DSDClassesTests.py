import unittest
import sys
from tempfile import mkstemp
import os

from .. import designer
from ..Srinivas2017 import translation

def format_list(templates, word):
    if type(templates) is list:
        return [format_list(temp, word) for temp in templates]
    else:
        return templates.format(word)

class TestTranslation(unittest.TestCase):
    reactions = [{'products': ['a'],
                           'rate': 1,
                           'reactants': [],
                           'stoich_p': [1],
                           'stoich_r': []},
                          {'products': [],
                           'rate': 1,
                           'reactants': ['b'],
                           'stoich_p': [],
                           'stoich_r': [1]},
                          {'products': ['a'],
                           'rate': 1,
                           'reactants': ['b'],
                           'stoich_p': [1],
                           'stoich_r': [1]},
                          {'products': ['b', 'b'],
                           'rate': 1,
                           'reactants': ['a'],
                           'stoich_p': [1, 1],
                           'stoich_r': [1]},
                          {'products': ['a'],
                           'rate': 1,
                           'reactants': ['b', 'c'],
                           'stoich_p': [1],
                           'stoich_r': [1, 1]},
                          {'products': ['a', 'a'],
                           'rate': 1,
                           'reactants': ['b', 'c'],
                           'stoich_p': [1, 1],
                           'stoich_r': [1, 1]}]
    species = ['a', 'b', 'c']
    nicknames = ['{0}-Out', '{0}-Backward', '{0}-cat_helper', '{0}-helper']
    specific_names = ['{0}-ch{0}-cch{0}-toe-sb{0}-bm{0}-cbm',
                      '{0}-toe-fb-suffix{0}-toe-sa{0}-am{0}-cam',
                      '{0}-toe-fd{0}-dh{0}-cdh{0}-toe-fc{0}-ch{0}-cch',
                      '{0}-toe-fd{0}-dh{0}-cdh{0}-toe-fc']

    F_all_identical_Af = ['{0}-Out',
                          '{0}-toe-sa{0}-am{0}-cam',
                          '{0}-dh{0}-cdh',
                          '{0}-ch{0}-cch',
                          '{0}-dh{0}-cdh']

    F_all_identical_As = ['{0}-ch{0}-cch', '{0}-bm{0}-cbm',
                          '{0}-toe-fb-suffix', '{0}-am{0}-cam',
                          '{0}-cat_helper',
                          '{0}-helper']

    F_all_identical_Bf = F_all_identical_Af[:]

    F_all_identical_Bs = F_all_identical_As[:]

    F_all_identical_Cf = F_all_identical_Af[:]

    F_all_identical_Cs = F_all_identical_As[:]

    F_all_identical_Df = F_all_identical_Af[:]

    F_all_identical_Ds = F_all_identical_As[:]

    F_all_different_Af = [
                          '{0}-Out',
                          '{0}-Backward',
                          '{0}-cat_helper',
                          '{0}-helper'
                          ]

    F_all_different_As = [
                          '{0}-Out',
                          '{0}-toe-fb-suffix', '{0}-am{0}-cam',
                          '{0}-cat_helper',
                          '{0}-helper'
                          ]

    F_all_different_Bf = [
                          '{0}-Out',
                          '{0}-toe-sa{0}-am{0}-cam',
                          '{0}-cat_helper',
                          '{0}-helper'
                          ]

    F_all_different_Bs = [
                          '{0}-ch{0}-cch', '{0}-bm{0}-cbm',
                          '{0}-Backward',
                          '{0}-cat_helper',
                          '{0}-helper'
                         ]

    F_all_different_Cf = [
                          '{0}-Out',
                          '{0}-Backward',
                          '{0}-toe-fd{0}-dh{0}-cdh', '{0}-ch{0}-cch',
                          '{0}-toe-fd{0}-dh{0}-cdh'
                          ]

    F_all_different_Cs = nicknames

    F_all_different_Df = [
                          '{0}-Out',
                          '{0}-Backward',
                          '{0}-dh{0}-cdh{0}-toe-fc{0}-ch{0}-cch',
                          '{0}-dh{0}-cdh{0}-toe-fc',
                          ]

    F_all_different_Ds = nicknames

    def test_F_all_identical(self):
        strand = translation.SignalStrand('A')
        strand.set_identity_domains(0, "r0")
        in_strands = [strand, strand]
        out_strands = [strand, strand]
        toe_no_interact_map = translation.F(in_strands + out_strands, 'r0')
        self.assertEqual(format_list(self.F_all_identical_Af, 'r0'), toe_no_interact_map[strand.th(0)])
        self.assertEqual(format_list(self.F_all_identical_As, 'r0'), toe_no_interact_map[strand.th(1)])
        self.assertEqual(format_list(self.F_all_identical_Bf, 'r0'), toe_no_interact_map[strand.th(0)])
        self.assertEqual(format_list(self.F_all_identical_Bs, 'r0'), toe_no_interact_map[strand.th(1)])
        self.assertEqual(format_list(self.F_all_identical_Cf, 'r0'), toe_no_interact_map[strand.th(0)])
        self.assertEqual(format_list(self.F_all_identical_Cs, 'r0'), toe_no_interact_map[strand.th(1)])
        self.assertEqual(format_list(self.F_all_identical_Df, 'r0'), toe_no_interact_map[strand.th(0)])
        self.assertEqual(format_list(self.F_all_identical_Ds, 'r0'), toe_no_interact_map[strand.th(1)])

    def test_F_all_different(self):
        strandA = translation.SignalStrand('A')
        strandA.set_identity_domains(0, "r0")
        strandB = translation.SignalStrand('B')
        strandB.set_identity_domains(1, "r0")
        strandC = translation.SignalStrand('C')
        strandC.set_identity_domains(2, "r0")
        strandD = translation.SignalStrand('D')
        strandD.set_identity_domains(3, "r0")
        in_strands = [strandA, strandB]
        out_strands = [strandC, strandD]
        toe_no_interact_map = translation.F(in_strands + out_strands, 'r0')
        self.assertEqual(format_list(self.F_all_different_Af, 'r0'), toe_no_interact_map[strandA.th(0)])
        self.assertEqual(format_list(self.F_all_different_As, 'r0'), toe_no_interact_map[strandA.th(1)])
        self.assertEqual(format_list(self.F_all_different_Bf, 'r0'), toe_no_interact_map[strandB.th(0)])
        self.assertEqual(format_list(self.F_all_different_Bs, 'r0'), toe_no_interact_map[strandB.th(1)])
        self.assertEqual(format_list(self.F_all_different_Cf, 'r0'), toe_no_interact_map[strandC.th(0)])
        self.assertEqual(format_list(self.F_all_different_Cs, 'r0'), toe_no_interact_map[strandC.th(1)])
        self.assertEqual(format_list(self.F_all_different_Df, 'r0'), toe_no_interact_map[strandD.th(0)])
        self.assertEqual(format_list(self.F_all_different_Ds, 'r0'), toe_no_interact_map[strandD.th(1)])

    def test_bimrxn_reaction_line(self):
        true_rxnline = "component r0 = bimrxn(<t>, <bm>, <c>): A + B -> C + D\n"
        strandA = translation.SignalStrand('A')
        strandA.set_identity_domains(0, "r0")
        strandB = translation.SignalStrand('B')
        strandB.set_identity_domains(1, "r0")
        strandC = translation.SignalStrand('C')
        strandC.set_identity_domains(2, "r0")
        strandC.add_instance(2, "r0")
        strandD = translation.SignalStrand('D')
        strandD.set_identity_domains(3, "r0")
        strandD.add_instance(3, "r0")
        bim = translation.Bimrxn('r0', [strandA, strandB], [strandC, strandD], (7, 15, 2))
        self.assertEqual(true_rxnline, bim.get_reaction_line())

    def runTest(self):
        pass

def suite():
    tests = ['test_F_all_identical', 'test_F_all_different']
    return unittest.TestSuite(list(map(TestTranslation, tests)))
