#! /usr/bin/env python

import string

import unittest

import system_parser

## Helper functions
def format_seq(signal):
  seq, reverse = signal
  if reverse:
    seq += "*"
  return seq


class TestSystemParser(unittest.TestCase):
  
  ## Declare statements tests
  example_declare = [
    ["NAME", [], [], []],
    ["cOmP_f__32", [], [], []],
    ["NAME", [], [["x", False]], [["y", False]]],
    ["NAME", [], [["x", False]], []],
    ["NAME", [], [], [["y", False]]],
    ["NAME", [], [["x", True]], [["y", False]]],
    ["NAME", [], [["x", False]], [["y", False]]],
    ["NAME", [], [["x", False], ["cow", True], ["BullDoGG_31", False]], 
                 [["x", False], ["y", True]]],
    ["NAME", ["a", "toe", "BoB"], [["x", False]], [["x", False]]],
  ]
  
  def test00_declare_noerror(self):
    """Test simple System Declare statement is accepted"""
    statement = "declare system NAME: x -> y"
    system_parser.parse_declare_statement(statement)
  
  def test01_declare_simple(self):
    """Test simple Component Declare statement is parsed correctly"""
    statement = "declare system NAME: x -> y"
    result = system_parser.parse_declare_statement(statement)
    self.assertEqual(["NAME", [], [["x", False]], [["y", False]]], result)
  
  def test02_declare_examples(self):
    """Test example System Declare statements are parsed correctly"""
    for name, params, inputs, outputs in self.example_declare:
      # Build the statement
      params_str = string.join(params, ", ")
      inputs_str = map(format_seq, inputs)
      inputs_str = string.join(inputs_str, " + ")
      outputs_str = map(format_seq, outputs)
      outputs_str = string.join(outputs_str, " + ")
      statement = "declare system %s(%s): %s -> %s" % (name, params_str, inputs_str, outputs_str)
      # Test the statement
      result = system_parser.parse_declare_statement(statement)
      self.assertEqual([name, params, inputs, outputs], result)
  
  
  ## Import statement tests
  example_import = [
    [["COMP_NAME", None]],
    [["COMP_PATH", "COMP_NAME"]],
    [["COMP/DIR/PATH", None]],
    [["COMP/DIR/PATH", "COMP_NAME"]],
    [["COMP_NAME", None], ["COMP2", None], ["FOO/COMP__3", None]],
    [["COMP/DIR/PATH", "COMP_NAME"], ["../COMP_PATH/2", "COMP2"]],
  ]
  
  def test10_import_noerror(self):
    """Test simple System Import statement is accepted"""
    statement = "import NAME"
    system_parser.parse_import_statement(statement)
  
  def test11_import_simple(self):
    """Test simple System Import statement is parsed correctly"""
    statement = "import NAME"
    result = system_parser.parse_import_statement(statement)
    self.assertEqual([["NAME", None]], result)
  
  def test12_import_examples(self):
    """Test example System Import statements are parsed correctly"""
    for imports in self.example_import:
      # Build the statement
      parts = [("%s as %s" % (path, name) if name else path) for path, name in imports]
      statement = "import " + string.join(parts, ", ")
      # Test the statement
      result = system_parser.parse_import_statement(statement)
      self.assertEqual(imports, result)
  
  
  ## Component statement tests
  example_component = [
    ["TYPE", "NAME", [], [], []],
    ["TYPE", "NAME", [], [["x", False]], []],
    ["TYPE", "NAME", [], [], [["y", False]]],
    ["TYPE", "NAME", [], [["x", True]], []],
    ["TYPE", "NAME", [], [["x", False], ["y", True]], [["z", True]]],
    ["TYPE", "NAME", [5, 10], [], []],
    ["TYPE", "NAME", [5, 10, "ATCGGTCA"], [], []],
    ["TYPE", "NAME", [5, 10, "ATCGGTCA"], [["x", False], ["y", True]], [["z", True]]],
  ]
  
  def test20_component_noerror(self):
    """Test simple System Component statement is accepted"""
    statement = "component TYPE = NAME: ->"
    system_parser.parse_component_statement(statement)
  
  def test21_component_simple(self):
    """Test simple System Component statement is accepted"""
    statement = "component TYPE = NAME: ->"
    result = system_parser.parse_component_statement(statement)
    self.assertEqual(["TYPE", "NAME", [], [], []], result)
  
  def test22_component_examples(self):
    """Test example System Component statements are parsed correctly"""
    for type, name, params, ins, outs in self.example_component:
      # Build the statement
      params_str = map(repr, params)
      params_str = string.join(params_str, ", ")
      ins_str = map(format_seq, ins)
      ins_str = string.join(ins_str, " + ")
      outs_str = map(format_seq, outs)
      outs_str = string.join(outs_str, " + ")
      statement = "component %s = %s(%s): %s -> %s" % (type, name, params_str, ins_str, outs_str)
      # Test the statement
      result = system_parser.parse_component_statement(statement)
      self.assertEqual([type, name, params, ins, outs], result)


if __name__ == '__main__':
  unittest.main()

