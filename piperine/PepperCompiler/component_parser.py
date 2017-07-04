"""DNA component design grammar"""

import re

from component_parser_regex import (parse_declare_statement, 
  parse_general_sequence_statement, parse_super_sequence_statement,
  parse_strand_statement, parse_structure_statement, parse_kinetic_statement,
  sequence_flag, domains_flag, nucleotide_flag)

from var_substitute import process_list
from utils import error


def load_component(filename, args, prefix):
  """Load component file"""
  from component_class import Component  # Imported in the function because of import loop problems in python
  
  # Find first statement
  f = open(filename, "r")
  for line in f:
    line = re.sub(r"#.*\n", "", line)
    line = line.strip()
    if line:
      break
  # Parse it as a declare statement
  name, param_names, component_inputs, component_outputs = parse_declare_statement(line)
  
  # Create component object
  component = Component(name, prefix, param_names)
  
  # Open file and do parameter substitution
  params = {}
  if len(param_names) != len(args):
    error("Component %s takes %d parameters %r, %d given %r." % (filename, len(param_names), tuple(param_names), len(args), tuple(args)))
  for name, val in zip(param_names, args):
    params[name] = val
  doc = process_list(f, params)
  
  # Add all the sequences, strands, etc. to the Component
  for line in doc.split("\n"):
    # Strip comments and whitespace
    line = re.sub(r"#.*\n", "", line)
    line = line.strip()
    # Skip empty lines
    if not line:
      continue
    
    # Read the command name off (if there is a command name)
    command = line.split()[0]
    if command == "declare":
      error("Multiple declare statements in file %s\n%s" % (filename, line))
    
    elif command == "sequence":
      name, const, length = parse_general_sequence_statement(line)
      # Is it a base-sequence
      if len(const) == 1 and const[0][0] == nucleotide_flag:
        component.add_sequence(name, const[0][1], length)
      # Or a super-sequence
      else:
        component.add_super_sequence(name, const, length)
    
    elif command in ("super-sequence", "sup-sequence"):
      name, const, length = parse_sequence_statement(line)
      component.add_super_sequence(name, const, length)
    
    elif command == "strand":
      name, dummy, constraints, length = parse_strand_statement(line)
      component.add_strand(name, dummy, constraints, length)
    
    elif command == "structure":
      name, opt, strands, struct = parse_structure_statement(line)
      component.add_structure(name, opt, strands, struct)
    
    elif command == "kinetic":
      low, high, inputs, outputs = parse_kinetic_statement(line)
      component.add_kinetic(low, high, inputs, outputs)
    
    elif command == "equal":
      #TODO: implement
      error("equal command not implemented yet (in file %s)\n%s" (filename, line))
    
    else:
      error("Parse Error in file '%s': Command '%s' not valid.\n%s" % (filename, command, line)) # TODO: more info
  
  # Add IO to connect everything up.
  component.add_IO(component_inputs, component_outputs)
  return component

