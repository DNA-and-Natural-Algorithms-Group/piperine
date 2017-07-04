"""DNA component design grammar"""

import re

from system_parser_pyparsing import (parse_declare_statement, 
  parse_import_statement, parse_component_statement)

from system_class import System
from var_substitute import process_list
from utils import error


def load_system(filename, args, prefix, path, includes=None):
  """Load system file"""
  # Find first statement
  f = open(filename, "r")
  
  line = "" # Default initialization for line, so that empty files don't crash.
  for line in f:
    line = re.sub(r"#.*\n", "", line)
    line = line.strip()
    if line:
      break
  # Parse it as a declare statement
  name, param_names, system_inputs, system_outputs = parse_declare_statement(line)
  
  # Create system object
  system = System(path, name, prefix, param_names, includes=includes)
  
  # Open file and do parameter substitution
  params = {}
  if len(param_names) != len(args):
    error("System %s takes %d parameters %r, %d given %r." % (filename, len(param_names), tuple(param_names), len(args), tuple(args)))
  for name, val in zip(param_names, args):
    params[name] = val
  doc = process_list(f, params)
  
  # Add all the sequences, strands, etc. to the System
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
    
    elif command == "import":
      imports = parse_import_statement(line)
      system.add_import(imports)
    
    elif command == "component":
      template, name, params, ins, outs = parse_component_statement(line)
      system.add_component(template, name, params, ins, outs)
    
    else:
      error("Parse Error in file '%s': Command '%s' not valid.\n%s" % (filename, command, line)) # TODO: more info
  
  # Add IO to connect everything up.
  system.add_IO(system_inputs, system_outputs)
  return system

