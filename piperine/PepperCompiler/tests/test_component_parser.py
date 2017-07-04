#! /usr/bin/env python

import string

import unittest

import component_parser
from component_parser import sequence_flag, nucleotide_flag, domains_flag

## Helper functions
def format_signal(signal):
  (seq, reverse), struct = signal
  if reverse:
    seq += "*"
  if struct:
    return "%s(%s)" % (seq, struct)
  else:
    return seq

def format_constraint(constraint):
  constraint_type, value = constraint
  if constraint_type == sequence_flag:
    name, reverse = value
    if reverse:
      name += "*"
    return name
  elif constraint_type == domains_flag:
    name, reverse = value
    if reverse:
      name += "*"
    return "domains(%s)" % name
  else:
    assert constraint_type == nucleotide_flag, constraint_type
    nucleotides = [str(num) + letter for num, letter in value]
    nucleotides = string.join(nucleotides)
    return '"%s"' % nucleotides


class TestComponentParser(unittest.TestCase):
  
  ## Declare statements tests
  example_declare = [
    ["NAME", [], [], []],
    ["cOmP_f__32", [], [], []],
    ["NAME", [], [[["x", False], "X"]], [[["y", False], "Y"]]],
    ["NAME", [], [[["x", False], "X"]], []],
    ["NAME", [], [], [[["y", False], "Y"]]],
    ["NAME", [], [[["x", True], "X"]], [[["y", False], "Y"]]],
    ["NAME", [], [[["x", False], None]], [[["y", False], None]]],
    ["NAME", [], [[["x", False], "X"], [["cow", True], None], [["BullDoGG_31", False], "bull__dogg32Strand"]], 
                 [[["x", False], None], [["y", True], "X"]]],
    ["NAME", ["a", "toe", "BoB"], [[["x", False], "X"]], [[["x", False], "Y"]]],
  ]
  
  def test00_declare_noerror(self):
    """Test simple Component Declare statement is accepted"""
    statement = "declare component NAME: x(X) -> y(Y)"
    component_parser.parse_declare_statement(statement)
  
  def test01_declare_simple(self):
    """Test simple Component Declare statement is parsed correctly"""
    statement = "declare component NAME: x(X) -> y(Y)"
    result = component_parser.parse_declare_statement(statement)
    self.assertEqual(["NAME", [], [[["x", False], "X"]], [[["y", False], "Y"]]], result)
  
  def test02_declare_no_sig(self):
    """Test a Component Declare statement with no signals is parsed correctly"""
    statement = "declare component NAME: ->"
    result = component_parser.parse_declare_statement(statement)
    self.assertEqual(["NAME", [], [], []], result)
  
  def test03_declare_examples(self):
    """Test example Component Declare statements are parsed correctly"""
    for name, params, inputs, outputs in self.example_declare:
      # Build the statement
      params_str = string.join(params, ", ")
      inputs_str = map(format_signal, inputs)
      inputs_str = string.join(inputs_str, " + ")
      outputs_str = map(format_signal, outputs)
      outputs_str = string.join(outputs_str, " + ")
      statement = "declare component %s(%s): %s -> %s" % (name, params_str, inputs_str, outputs_str)
      # Test the statement
      result = component_parser.parse_declare_statement(statement)
      self.assertEqual([name, params, inputs, outputs], result)
  
  
  ## Sequence statement tests (includes super-sequences which can be specified this way)
  example_sequence = [ 
    ["NAME", [[nucleotide_flag, [[5, "N"]]]], None],
    ["cOmP_f__32", [[nucleotide_flag, [[5, "N"]]]], None],
    ["NAME", [[nucleotide_flag, [[2, "S"], ["?", "W"], [1, "A"], [1, "C"], [12, "Y"], [1, "B"]]]], None],
    ["NAME", [[nucleotide_flag, [[5, "N"]]]], 5],
    ["NAME", [[sequence_flag, ["seq1", False]]], None],
    ["NAME", [[sequence_flag, ["seq1", True]]], None],
    ["NAME", [[sequence_flag, ["seq1", False]], [nucleotide_flag, [[5, "N"]]], 
              [sequence_flag, ["seq2", True]]], None],
    ["NAME", [[domains_flag, ["seq1", False]]], None],
  ]
  
  def test10_sequence_noerror(self):
    """Test simple Component Sequence statement is accepted"""
    statement = 'sequence NAME = "5N"'
    component_parser.parse_general_sequence_statement(statement)
  
  def test11_sequence_simple(self):
    """Test simple Component Sequence statement is parsed correctly"""
    statement = 'sequence NAME = "5N"'
    result = component_parser.parse_general_sequence_statement(statement)
    # We don't test constraints, because it could be parsed as 5N or N + N + N + N + N
    self.assertEqual( ["NAME", [[nucleotide_flag, [[5, "N"]]]], None], result)
  
  def test12_sequence_super(self):
    """Test simple Component Super Sequence statement is parsed correctly"""
    statement = 'sequence NAME = seq1'
    result = component_parser.parse_general_sequence_statement(statement)
    # We don't test constraints, because it could be parsed as 5N or N + N + N + N + N
    self.assertEqual( ["NAME", [[sequence_flag, ["seq1", False]]], None], result)
  
  def test13_sequence_examples(self):
    """Test example Component Sequence statements are parsed correctly"""
    for name, constraints, length in self.example_sequence:
      # Build the statement
      constr_str = map(format_constraint, constraints)
      constr_str = string.join(constr_str, " ")
      length_str = (": %d" % length if length != None else "")
      statement = "sequence %s = %s %s" % (name, constr_str, length_str)
      # Test the statement
      result = component_parser.parse_general_sequence_statement(statement)
      self.assertEqual([name, constraints, length], result)
  
  
  ## Strand statement tests
  example_strand = [ 
    ["NAME", False, [[sequence_flag, ["seq", False]]], None],
    ["cOmP_f__32", False, [[sequence_flag, ["name", False]]], None],
    ["NAME", True, [[sequence_flag, ["seq", False]]], None],
    ["NAME", False, [[sequence_flag, ["seq", False]]], 13],
    ["NAME", False, [[sequence_flag, ["seq", False]], [sequence_flag, ["OtherSequence__42", False]]], None],
    ["NAME", False, [[sequence_flag, ["seq", True]]], None],
    ["NAME", False, [[nucleotide_flag, [[5, "N"]]]], 5],
    ["NAME", False, [[nucleotide_flag, [[2, "S"], ["?", "W"], [1, "A"], [1, "C"], [12, "Y"], [1, "B"]]]], None],
  ]
  
  def test20_strand_noerror(self):
    """Test simple Component Strand statement is accepted"""
    statement = 'strand NAME = seq'
    component_parser.parse_strand_statement(statement)
  
  def test21_strand_simple(self):
    """Test simple Component Strand statement is accepted"""
    statement = 'strand NAME = seq'
    result = component_parser.parse_strand_statement(statement)
    self.assertEqual( [False, "NAME", [[sequence_flag, ["seq", False]]], None], result )
  
  def test22_strand_examples(self):
    """Test example Component Strand statements are parsed correctly"""
    for name, dummy, constraints, length in self.example_strand:
      # Build the statement
      dummy_str = ("[dummy]" if dummy else "")
      constr_str = map(format_constraint, constraints)
      constr_str = string.join(constr_str, " ")
      length_str = (": %d" % length if length != None else "")
      statement = "strand %s %s = %s %s" % (dummy_str, name, constr_str, length_str)
      # Test the statement
      result = component_parser.parse_strand_statement(statement)
      self.assertEqual([dummy, name, constraints, length], result)
  
  
  ## Structure statement tests
  example_structure = [
    ["NAME", 1.0, ["strand"], [False, ".."]],
    ["NAME", 3.4, ["strand"], [False, ".."]],
    ["cOmP_f__32", 1.0, ["strand"], [False, ".."]],
    ["NAME", 1.0, ["strand1", "St4R__and42"], [False, ".."]],
    ["NAME", 1.0, ["strand"], [True, ".."]],
    ["NAME", 1.0, ["strand"], [False, "..(((+))..)"]],
#    ["NAME", 1.0, ["strand"], [False, "5. 3( 7( + 7) .. 3)"]], #TODO: This probably won't work
#    ["NAME", 1.0, ["strand"], [False, "U5 H3(H7(+) U2)"]], #TODO: This probably won't work
  ]
  
  def test30_structure_noerror(self):
    """Test simple Component Structure statement is accepted"""
    statement = 'structure NAME = strand : ..'
    component_parser.parse_structure_statement(statement)
  
  def test31_structure_simple(self):
    """Test simple Component Structure statement is parsed correctly"""
    statement = 'structure NAME = strand : ..'
    result = component_parser.parse_structure_statement(statement)
    self.assertEqual( [1.0, "NAME", ["strand"], [False, ".."]], result )
  
  def test32_structure_no_opt(self):
    """Test a no-opt Component Structure statement is parsed correctly"""
    statement = 'structure [no-opt] NAME = strand : ..'
    result = component_parser.parse_structure_statement(statement)
    self.assertEqual( [0.0, "NAME", ["strand"], [False, ".."]], result )
  
  def test33_structure_examples(self):
    """Test example Component Structure statements are parsed correctly"""
    for name, opt, strands, (domain, struct) in self.example_structure:
      # Build the statement
      strands_str = string.join(strands, " + ")
      domain_str = ("domain" if domain else "")
      statement = "structure [%snt] %s = %s : %s %s" % (opt, name, strands_str, domain_str, struct)
      # Test the statement
      result = component_parser.parse_structure_statement(statement)
      self.assertEqual([opt, name, strands, [domain, struct]], result)
  
  
  ## Kinetic statement tests
  # TODO: kinetic constraints
  example_kinetic = [
    [None, None, ["A"], ["B"]],
    [None, None, [], ["B"]],
    [None, None, ["A"], []],
    [None, None, ["A"], ["B", "A", "foo3Stru__ct99"]],
    [1.0e5, None, ["A"], ["B"]],
    [None, 1.0e10, ["A"], ["B"]],
    [1.0e5, 1.0e10, ["A"], ["B"]],
  ]
  
  def test40_kinetic_noerror(self):
    """Test simple Component Kinetic statement is accepted"""
    statement = 'kinetic A -> B'
    component_parser.parse_kinetic_statement(statement)
  
  def test41_kinetic_simple(self):
    """Test simple Component Kinetic statement is accepted"""
    statement = 'kinetic A -> B'
    result = component_parser.parse_kinetic_statement(statement)
    self.assertEqual( [None, None, ["A"], ["B"]], result )
  
  def test42_structure_examples(self):
    """Test example Component Kinetic statements are parsed correctly"""
    for low, high, inputs, outputs in self.example_kinetic:
      # Build the statement
      if low and high:
        params = "[%f /M/s < k < %f /M/s]" % (low, high)
      elif low and not high:
        params = "[k > %f /M/s]" % low
      elif high and not low:
        params = "[k < %f /M/s]" % high
      else:
        params = ""
        
      inputs_str  = string.join(inputs, " + ")
      outputs_str = string.join(outputs, " + ")
      statement = "kinetic %s %s -> %s" % (params, inputs_str, outputs_str)
      # Test the statement
      result = component_parser.parse_kinetic_statement(statement)
      self.assertEqual([low, high, inputs, outputs], result)

if __name__ == '__main__':
  unittest.main()

